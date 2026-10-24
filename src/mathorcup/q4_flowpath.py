"""Q4 route QA. A calculated corridor is a *provisional candidate*, never a verified avalanche flowpath.

Topographic Dijkstra: Penalise up-gradient steps; named source and target require external verification.
A verified route can be supplied as events/verified_route.csv with 'lon','lat' ordered along the valley.
"""
from __future__ import annotations
import heapq
import numpy as np,pandas as pd,rasterio,matplotlib.pyplot as plt
from pathlib import Path
from pyproj import Transformer
from rasterio.warp import reproject,Resampling
from rasterio.transform import from_origin,rowcol,xy
from scipy.ndimage import gaussian_filter1d
from .prepare_geo import glacier_grid
from .common import COLORS,savefig,write_json,record_warning,find_first

UTM_CRS='EPSG:32645'

def glacier_source_candidate(ctx):
    # Prefer Q2's contiguous upper-score ice patch based on historical speed,
    # slope and thickness rather than the whole-glacier centroid.
    q2=ctx.folder('q2')/'Q2_candidate_release_region.json'
    g=glacier_grid(ctx);mask=g['glacier_mask'];ii,jj=np.where(mask)
    if q2.exists():
        import json
        d=json.loads(q2.read_text(encoding='utf-8'))
        x=float(d['candidate_centroid_utm_x']);y=float(d['candidate_centroid_utm_y'])
        wgs=Transformer.from_crs(UTM_CRS,'EPSG:4326',always_xy=True)
        lon,lat=wgs.transform(x,y)
        return {'lon':float(lon),'lat':float(lat),'x':x,'y':y,
           'candidate_inventory_ice_volume_m3':float(d['candidate_inventory_ice_volume_m3']),
           'candidate_patch_area_km2':float(d['candidate_area_km2']),
           'identification':'Q2 top-15% physical ice component centroid; candidate not observed rupture'}

    if len(ii)==0:raise ValueError('No valid glacier mask')
    x=float(np.mean(g['x'][jj]));y=float(np.mean(g['y'][ii]))
    wgs=Transformer.from_crs(UTM_CRS,'EPSG:4326',always_xy=True)
    lon,lat=wgs.transform(x,y)
    return {'lon':float(lon),'lat':float(lat),'x':x,'y':y,
            'identification':'Fallback mask centroid if Q2 unrun, NOT verified initiation'}

def local_dem(ctx,xy_points,buffer_m=3400,resolution_m=60):
    xs=[p[0] for p in xy_points];ys=[p[1] for p in xy_points]
    left=min(xs)-buffer_m;right=max(xs)+buffer_m
    bottom=min(ys)-buffer_m;top=max(ys)+buffer_m
    nx=int(np.ceil((right-left)/resolution_m));ny=int(np.ceil((top-bottom)/resolution_m))
    if nx*ny>1_500_000:raise ValueError('Geographic ROI too large; verify endpoints and grid first')
    dempath=find_first(ctx.raw/'unpacked'/'7','DEM.tif')
    transform=from_origin(left,top,resolution_m,resolution_m)
    target=np.full((ny,nx),np.nan,dtype=np.float32)
    with rasterio.open(dempath) as src:
        reproject(source=rasterio.band(src,1),destination=target,
                src_transform=src.transform,src_crs=src.crs,src_nodata=src.nodata,
                dst_transform=transform,dst_crs=UTM_CRS,dst_nodata=np.nan,
                resampling=Resampling.bilinear)
    return target,transform,UTM_CRS

def _pixel_from_utm(transform,x,y,shape):
    r,c=rowcol(transform,x,y)
    r,c=int(r),int(c)
    if not (0<=r<shape[0] and 0<=c<shape[1]):raise ValueError('Point outside DEM ROI')
    return r,c

def directed_cost_path(z,source,target,cell_size,uphill_penalty=45):
    """8-neighbour directed Dijkstra. Cost is length*(1+uphill_penalty*positive_slope).
    This is a topographic CANDIDATE, not observed runout; it may cut across ridges.
    """
    nr,nc=z.shape;sr,sc=source;tr,tc=target
    if not np.isfinite(z[sr,sc]) or not np.isfinite(z[tr,tc]):raise ValueError('Source or target on DEM NoData')
    inf=float('inf');dist=np.full((nr,nc),inf,np.float64)
    previous=np.full((nr,nc),-1,np.int64)
    start=sr*nc+sc;goal=tr*nc+tc;dist[sr,sc]=0
    heap=[(0.,start)]
    steps=[(dr,dc,cell_size*np.hypot(dr,dc)) for dr in [-1,0,1] for dc in [-1,0,1] if dr or dc]
    nvisit=0
    while heap:
        cost,idx=heapq.heappop(heap)
        r,c=divmod(idx,nc)
        if cost>dist[r,c]+1e-9:continue
        if idx==goal:break
        zi=float(z[r,c]);nvisit+=1
        for dr,dc,dlen in steps:
            rr=r+dr;cc=c+dc
            if rr<0 or cc<0 or rr>=nr or cc>=nc or not np.isfinite(z[rr,cc]):continue
            rise=max(0.0,float(z[rr,cc])-zi)
            # Directed steepness penalty, positive edge cost so Dijkstra is exact for THIS assumed metric
            stepcost=dlen*(1+uphill_penalty*rise/dlen)
            new=cost+stepcost
            if new+1e-8<dist[rr,cc]:
                dist[rr,cc]=new;previous[rr,cc]=idx
                heapq.heappush(heap,(new,rr*nc+cc))
    if not np.isfinite(dist[tr,tc]):raise RuntimeError('No DEM-connected path between source and target')
    path=[];idx=goal
    while idx!=-1:
        r,c=divmod(int(idx),nc);path.append((r,c))
        if idx==start:break
        idx=previous[r,c]
    if path[-1]!=(sr,sc):raise RuntimeError('Disconnected backtrace')
    path.reverse()
    return np.array(path,dtype=int),float(dist[tr,tc]),nvisit

def _to_profile(utm_xy,z,smooth=0.8,spacing_m=60):
    coords=np.asarray(utm_xy,float)
    d=np.hypot(np.diff(coords[:,0]),np.diff(coords[:,1]));s=np.r_[0,np.cumsum(d)]
    keep=np.r_[True,np.diff(s)>1e-8];s=s[keep];coords=coords[keep];z=np.asarray(z)[keep]
    if s[-1]<spacing_m:raise ValueError('Path too short')
    su=np.arange(0,s[-1]+.01,spacing_m)
    if su[-1]<s[-1]:su=np.r_[su,s[-1]]
    xx=np.interp(su,s,coords[:,0]);yy=np.interp(su,s,coords[:,1])
    zz=np.interp(su,s,z)
    # smoothing only numerical DEM roughness; no forced monotonic heights
    if smooth>0:zz=gaussian_filter1d(zz,sigma=smooth)
    grad=np.gradient(zz,su)
    return pd.DataFrame({'distance_m':su,'utm_x':xx,'utm_y':yy,'elevation_m':zz,
                         'downward_grade':-grad})

def execute_q4_path(ctx,allow_provisional=False):
    out=ctx.folder('q4');c=ctx.cfg['q4'];candidate=glacier_source_candidate(ctx)
    wgs_to_utm=Transformer.from_crs('EPSG:4326',UTM_CRS,always_xy=True)
    target_lon,target_lat=c['port_lonlat_provisional'];target_xy=wgs_to_utm.transform(target_lon,target_lat)
    verified_file=ctx.root/c['route_waypoints_file'];verification=bool(c['route_verified']) and verified_file.exists()
    if not verification and allow_provisional:
        record_warning('Q4 runs a physically constrained competition scenario on an estimated DEM corridor; event-specific geolocation remains a source of error.')
    if not verification and not allow_provisional:
        record_warning('Q4 verified route unavailable. Generate geographic QA map & SYNTHETIC solver benchmark only.')
    z,transform,crs=local_dem(ctx,[(candidate['x'],candidate['y']),target_xy],resolution_m=c['grid_step_m'])
    start=_pixel_from_utm(transform,candidate['x'],candidate['y'],z.shape)
    dest=_pixel_from_utm(transform,*target_xy,z.shape)
    line=None;path_cost=None
    if verification:
        points=pd.read_csv(verified_file)
        if not {'lon','lat'}.issubset(points):raise ValueError('verified_route.csv requires lon,lat')
        if len(points)<3:raise ValueError('Verified route needs at least 3 ordered surveyed/validated points')
        coord=np.array([wgs_to_utm.transform(row.lon,row.lat) for row in points.itertuples()])
        if not np.isfinite(coord).all():raise ValueError('non finite verified path')
        from scipy.ndimage import map_coordinates
        pix=np.array([_pixel_from_utm(transform,x,y,z.shape) for x,y in coord]);bed=z[pix[:,0],pix[:,1]]
        profile=_to_profile(coord,bed,spacing_m=c['grid_step_m'])
        status='VERIFIED_ROUTE_USER_SUPPLIED'
        line=coord
    else:
        pixels,path_cost,visited=directed_cost_path(z,start,dest,c['grid_step_m'],uphill_penalty=45)
        line=np.array([xy(transform,int(r),int(col),offset='center') for r,col in pixels])
        bed=z[pixels[:,0],pixels[:,1]]
        profile=_to_profile(line,bed,spacing_m=c['grid_step_m'])
        status='PROVISIONAL_UNVERIFIED_TOPOGRAPHIC_COST_PATH'
    positive_climb=float(np.sum(np.maximum(np.diff(profile.elevation_m),0)))
    stats={'route_status':status,'q2_source_candidate_ice_region':candidate,
           'destination_lonlat_provisional':c['port_lonlat_provisional'],
           'routing_method':'directed topographic path with uphill penalty (UNVERIFIED)' if not verification else 'team provided confirmed waypoints',
           'path_length_km':float(profile.distance_m.iloc[-1]/1000),
           'path_start_elevation_m':float(profile.elevation_m.iloc[0]),
           'path_end_elevation_m':float(profile.elevation_m.iloc[-1]),
           'net_fall_m':float(profile.elevation_m.iloc[0]-profile.elevation_m.iloc[-1]),
           'cumulative_uphill_m':positive_climb,
           'license_to_claim_physical_reconstruction':bool(verification),
           'observed_2026_event_route_verified':bool(verification),
           'note':'Source point taken as official glacier mask centroid; candidate corridor calculated from official DEM with an uphill penalty. Appropriate as competition working hypothesis, subject to remote-sensing refinement.'}
    write_json(out/'Q4_route_quality_gate.json',stats)
    profile.to_csv(out/'Q4_route_profile.csv',index=False)
    # Plot georeferenced ROI DEM with tentative/verified line
    fig,axes=plt.subplots(1,2,figsize=(11.7,4.8),gridspec_kw={'width_ratios':[1.3,1]})
    extent=[transform.c,transform.c+z.shape[1]*transform.a,
            transform.f+z.shape[0]*transform.e,transform.f]
    im=axes[0].imshow(z,extent=extent,cmap='terrain',vmin=np.nanpercentile(z,4),vmax=np.nanpercentile(z,97))
    axes[0].plot(line[:,0],line[:,1],color='#C75041',lw=1.7,ls='--' if not verification else '-',label='未经验证的候选路径' if not verification else '团队核查路径')
    axes[0].scatter([candidate['x'],target_xy[0]],[candidate['y'],target_xy[1]],c=['#073B60','#F4B33C'],s=[65,65],edgecolors='white',zorder=10)
    axes[0].set(xlabel='UTM45N 东向坐标 (m)',ylabel='UTM45N 北向坐标 (m)',title='传播路径空间分布' if not verification else '传播路径空间分布')
    axes[0].legend(frameon=True,fontsize=8,loc='upper right');fig.colorbar(im,ax=axes[0],shrink=.7,label='高程 (m)')
    axes[1].plot(profile.distance_m/1000,profile.elevation_m,color=COLORS[0],lw=2)
    axes[1].fill_between(profile.distance_m/1000,profile.elevation_m,alpha=.15,color=COLORS[0]);axes[1].grid(True,alpha=.16)
    axes[1].set(xlabel='候选路径弧长 (km)',ylabel='高程 (m)',title='传播路径地形纵剖面')
    fig.suptitle('冰崩—泥石流传播路径及地形纵剖面',fontsize=12)
    savefig(fig,out/'Q4_01_route_geolocation_quality_gate.png')
    if positive_climb>300:record_warning(f'Estimated route accumulates {positive_climb:.0f} m local uphill; include route assumption in sensitivity.')
    if not verification and not allow_provisional:
        return None
    return profile
