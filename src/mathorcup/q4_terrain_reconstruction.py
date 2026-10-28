"""Event reconstruction geometry derived from the official DEM and glacier grid."""
from __future__ import annotations
import heapq
import json
import numpy as np
import pandas as pd
import yaml
from scipy.ndimage import gaussian_filter1d, map_coordinates, binary_dilation
from scipy.spatial import cKDTree
from pyproj import Transformer
from rasterio.transform import rowcol, xy, Affine
from .q4_flowpath import local_dem, UTM_CRS, _to_profile, directed_cost_path
from .prepare_geo import glacier_grid
from .common import write_json


def priority_drainage(z, dx):
    """Priority-flood depression filling for routing only; raw elevations stay intact.

    Each receiver has a lower filled elevation, or an earlier flood rank on a flat.
    This gives an acyclic drainage graph without carving an invented DEM channel.
    """
    z=np.asarray(z,float);ny,nx=z.shape;valid=np.isfinite(z)
    filled=z.copy();seen=~valid;rank=np.full(z.shape,-1,int)
    parent=np.full(z.size,-1,int);heap=[];order=[]
    boundary=valid&binary_dilation(~valid)
    boundary[[0,-1],:]|=valid[[0,-1],:];boundary[:,[0,-1]]|=valid[:,[0,-1]]
    for r,c in np.argwhere(boundary):
        seen[r,c]=True;heapq.heappush(heap,(float(z[r,c]),int(r*nx+c)))
    neighbours=[(dr,dc,float(dx*np.hypot(dr,dc))) for dr in (-1,0,1) for dc in (-1,0,1) if dr or dc]
    while heap:
        level,idx=heapq.heappop(heap);r,c=divmod(idx,nx)
        rank[r,c]=len(order);order.append(idx)
        for dr,dc,_ in neighbours:
            rr,cc=r+dr,c+dc
            if not (0<=rr<ny and 0<=cc<nx) or seen[rr,cc]:continue
            seen[rr,cc]=True;filled[rr,cc]=max(float(z[rr,cc]),level)
            parent[rr*nx+cc]=idx
            heapq.heappush(heap,(float(filled[rr,cc]),int(rr*nx+cc)))
    receiver=parent.copy()
    for idx in order:
        r,c=divmod(idx,nx);best=0.;dest=receiver[idx]
        for dr,dc,length in neighbours:
            rr,cc=r+dr,c+dc
            if not (0<=rr<ny and 0<=cc<nx) or not valid[rr,cc]:continue
            slope=(filled[r,c]-filled[rr,cc])/length
            if slope>best+1e-12:best=slope;dest=rr*nx+cc
        receiver[idx]=dest
    accumulation=np.where(valid,dx*dx,0.).ravel()
    # Sorting by filled height then flood rank respects every receiver, including flats.
    ids=np.array(order,int)
    topo=ids[np.lexsort((rank.ravel()[ids],filled.ravel()[ids]))]
    for idx in topo[::-1]:
        dest=receiver[idx]
        if dest>=0:accumulation[dest]+=accumulation[idx]
    return filled,receiver.reshape(z.shape),accumulation.reshape(z.shape)


def trace_drainage(receiver,start):
    nx=receiver.shape[1];idx=int(start[0]*nx+start[1]);result=[];visited=set()
    while idx>=0:
        if idx in visited:raise ValueError('Drainage receiver cycle')
        visited.add(idx);r,c=divmod(idx,nx);result.append((r,c));idx=int(receiver[r,c])
    return np.array(result,int)


def section_widths(z,tr,profile,relief_m=20.,sample_step_m=30.,half_span_m=900.):
    """Contiguous DEM valley width at a declared relative height, plus section samples."""
    coord=profile[['utm_x','utm_y']].to_numpy(float)
    tangent=np.gradient(coord,axis=0)
    tangent/=np.maximum(np.linalg.norm(tangent,axis=1)[:,None],1e-12)
    normals=np.column_stack([-tangent[:,1],tangent[:,0]])
    offsets=np.arange(-half_span_m,half_span_m+sample_step_m*.1,sample_step_m)
    points=coord[:,None,:]+normals[:,None,:]*offsets[None,:,None]
    rr=(points[:,:,1]-tr.f)/tr.e-.5;cc=(points[:,:,0]-tr.c)/tr.a-.5
    elevation=map_coordinates(z,[rr,cc],order=1,mode='constant',cval=np.nan)
    centre=len(offsets)//2;width=[];beds=[]
    for line in elevation:
        # Use the connected trough near the drainage line, without jumping across a ridge.
        lo=max(0,centre-2);hi=min(len(line),centre+3)
        k=lo+int(np.nanargmin(line[lo:hi]));bed=float(line[k]);beds.append(bed)
        wet=np.isfinite(line)&(line<=bed+relief_m);left=right=k
        while left>0 and wet[left-1]:left-=1
        while right<len(line)-1 and wet[right+1]:right+=1
        width.append(max(sample_step_m,(right-left+1)*sample_step_m))
    # Smooth pixel-scale jumps, retaining the underlying sampled cross sections for audit.
    widths=gaussian_filter1d(np.array(width,float),sigma=max(1,90/np.median(np.diff(profile.distance_m))))
    return widths,offsets,elevation,np.array(beds)


def prepare_reconstruction_geometry(ctx,evidence):
    settings=evidence['reconstruction'];dx=float(settings['routing_resolution_m'])
    project=Transformer.from_crs('EPSG:4326',UTM_CRS,always_xy=True)
    mapped_source=project.transform(*evidence['locations']['source_lonlat']['value'])
    # Select the same real ice-pixel footprint as before, then start the route
    # on its upper-elevation part. This covers the physical release footprint
    # instead of shifting its entire PCA length down an already lower valley.
    g=glacier_grid(ctx);ii,jj=np.where(g['glacier_mask']&np.isfinite(g['thickness'])&np.isfinite(g['dem']))
    all_pixels=np.column_stack([g['x'][jj],g['y'][ii]])
    pixel_area=abs(float(np.median(np.diff(g['x'])))*float(np.median(np.diff(g['y']))))
    n=int(round(evidence['published_values']['source_area_m2']['value']/pixel_area))
    nearest=np.argsort(np.linalg.norm(all_pixels-np.asarray(mapped_source),axis=1))[:n]
    pixels=all_pixels[nearest];pixel_heights=g['dem'][ii[nearest],jj[nearest]]
    upper=pixels[pixel_heights>=np.quantile(pixel_heights,.9)].mean(axis=0)
    origin=int(np.argmin(np.linalg.norm(pixels-upper,axis=1)))
    source=tuple(pixels[origin])
    site=project.transform(*evidence['locations']['port_area_lonlat']['value'])
    z,tr,_=local_dem(ctx,[source,site],buffer_m=10000,resolution_m=dx)
    # The main river bends north of the straight source-to-port rectangle. The
    # previous 3.5 km buffer cut that river off and encouraged ridge shortcuts.
    rr,cc=np.where(np.isfinite(z))
    r0,r1,c0,c1=int(rr.min()),int(rr.max()+1),int(cc.min()),int(cc.max()+1)
    z=z[r0:r1,c0:c1];tr=tr*Affine.translation(c0,r0)
    filled,receiver,acc=priority_drainage(z,dx)
    if ctx.cfg.get('strict_route',False):
        from .q4_flowpath import require_verified_route
        waypoints=require_verified_route(ctx)
        xx,yy=project.transform(waypoints.lon.to_numpy(),waypoints.lat.to_numpy())
        coord=np.column_stack([xx,yy]);source=tuple(coord[0]);site=tuple(coord[-1])
        rr=(coord[:,1]-tr.f)/tr.e-.5;cc=(coord[:,0]-tr.c)/tr.a-.5
        raw=map_coordinates(z,[rr,cc],order=1,mode='constant',cval=np.nan)
        if not np.isfinite(raw).all():raise ValueError('Verified cascade waypoints leave valid official DEM')
        method='verified ordered waypoints sampled on official DEM'
    else:
        start=tuple(map(int,rowcol(tr,*source)))
        chain=trace_drainage(receiver,start)
        coord=np.array([xy(tr,int(r),int(c),offset='center') for r,c in chain])
        distance=np.linalg.norm(coord-np.asarray(site),axis=1);closest=int(np.argmin(distance))
        if closest<10 or distance[closest]>500:
            dest=tuple(map(int,rowcol(tr,*site)))
            chain,_,_=directed_cost_path(z,start,dest,dx,uphill_penalty=250)
            coord=np.array([xy(tr,int(r),int(c),offset='center') for r,c in chain])
            method='high-uphill-penalty DEM corridor; drainage trace did not approach site'
            closest=len(chain)-1
        else:method='priority-flood D8 drainage trace, clipped at nearest port-area river cell'
        chain=chain[:closest+1];coord=coord[:closest+1]
        raw=z[chain[:,0],chain[:,1]]
    profile=_to_profile(coord,raw,smooth=.5,spacing_m=float(settings['solver_dx_m']))
    # Strictly uniform cells are required by the finite-volume update.
    profile=profile.iloc[:-1].copy() if abs(np.diff(profile.distance_m)[-1]-settings['solver_dx_m'])>1e-7 else profile
    r,c=rowcol(tr,profile.utm_x.to_numpy(),profile.utm_y.to_numpy());r=np.asarray(r);c=np.asarray(c)
    profile['drainage_area_m2']=acc[r,c]
    widths,offsets,sections,beds=section_widths(z,tr,profile,settings['cross_section_relief_m'],dx)
    profile['valley_width_m']=widths
    wider,_,_,_=section_widths(z,tr,profile,60.,dx)
    side=np.maximum((wider-widths)/80.,0.)
    side=np.minimum(side,np.maximum(widths-dx,0)/40.)
    bottom_width=np.maximum(dx,widths-40*side)
    profile['bank_side_slope']=side
    # The transition is picked among DEM drainage-area jumps near the published river distance.
    total=float(profile.distance_m.iloc[-1]);expected=total-evidence['published_values']['confluence_to_port_distance_m']['value']
    aa=profile.drainage_area_m2.to_numpy();jump=aa/np.maximum(np.r_[aa[0],aa[:-1]],dx*dx)
    window=(profile.distance_m.between(expected-1500,expected+1500)).to_numpy()
    ids=np.flatnonzero(window)
    if not len(ids):raise ValueError('Published confluence-distance window outside DEM route')
    score=np.log1p(jump[ids]-1)-np.abs(profile.distance_m.to_numpy()[ids]-expected)/1500
    junction=int(ids[np.argmax(score)])
    # Metric 3-D projection preserves both actual pixel locations and heights.
    # Moving a PCA coordinate to zero would lower the initial material by
    # hundreds of metres and change its gravitational energy before release.
    along,projection=project_release_pixels(pixels,pixel_heights,profile)
    if (along>=junction).any():raise ValueError('A mapped source pixel projects downstream of the transition')
    counts=np.bincount(along,minlength=len(profile)).astype(float)
    counts=gaussian_filter1d(counts,max(1,dx/settings['solver_dx_m']),mode='nearest')
    counts[junction:]=0
    source_weights=counts/counts.sum()
    source_width=counts*pixel_area/float(settings['solver_dx_m'])
    profile['effective_width_m']=np.maximum(bottom_width,source_width)
    profile['source_area_weight']=source_weights
    # Save exact arrays used by plots/solver; these generated tables remain local.
    out=ctx.folder('q4');profile.to_csv(out/'Q4_reconstruction_route.csv',index=False)
    projection.to_csv(out/'Q4_release_pixel_projection.csv',index=False)
    np.savez_compressed(out/'Q4_reconstruction_sections.npz',offsets_m=offsets,elevation_m=sections,
                        thalweg_elevation_m=beds,source_pixels_utm=pixels,
                        background_elevation_m=z,background_transform=np.array(tuple(tr)))
    quality={'method':method,'source_lonlat':evidence['locations']['source_lonlat'],
       'port_area_reference':evidence['locations']['port_area_lonlat'],
       'port_snap_distance_m':float(np.linalg.norm(profile[['utm_x','utm_y']].iloc[-1].to_numpy()-np.asarray(site))),
       'length_m':total,'drop_m':float(profile.elevation_m.iloc[0]-profile.elevation_m.iloc[-1]),
       'raw_uphill_sum_m':float(np.maximum(np.diff(profile.elevation_m),0).sum()),
       'drainage_fill_max_m':float(np.nanmax(filled-z)),'raw_DEM_used_for_dynamics':True,
       'confluence_index':junction,'confluence_distance_m':float(profile.distance_m.iloc[junction]),
       'downstream_reach_length_m':total-float(profile.distance_m.iloc[junction]),
       'source_pixels':len(pixels),'source_pixel_area_m2':pixel_area,
       'route_source_origin_utm_xy':list(source),
       'source_origin_rule':'actual source pixel nearest centroid of upper elevation decile',
       'release_projection_method':'nearest point on actual route in metric UTM-x/UTM-y/elevation space',
       'selected_source_mean_elevation_m':float(pixel_heights.mean()),
       'assigned_release_mean_bed_elevation_m':float(np.dot(source_weights,profile.elevation_m)),
       'release_projection_mean_abs_elevation_error_m':float(np.abs(projection.route_elevation_m-projection.source_elevation_m).mean()),
       'projected_source_area_m2':len(pixels)*pixel_area,
       'section_reference_relief_m':settings['cross_section_relief_m'],
       'cross_section_model':'trapezoid fitted to connected DEM widths at 20 m and 60 m relative height',
       'bank_side_slope_median':float(np.median(side)),
       'width_min_m':float(profile.effective_width_m.min()),'width_median_m':float(profile.effective_width_m.median()),
       'numerical_dx_m':settings['solver_dx_m'],'DEM_sampling_m':dx}
    write_json(out/'Q4_reconstruction_geometry.json',quality)
    return profile,junction,z,tr,quality


def project_release_pixels(pixels_utm,heights,profile):
    points=np.column_stack([np.asarray(pixels_utm,float),np.asarray(heights,float)])
    route=profile[['utm_x','utm_y','elevation_m']].to_numpy(float)
    if not np.isfinite(points).all() or not np.isfinite(route).all():
        raise ValueError('Release projection requires finite actual coordinates and elevations')
    distance,indices=cKDTree(route).query(points)
    matched=profile.iloc[indices]
    table=pd.DataFrame({'source_utm_x':points[:,0],'source_utm_y':points[:,1],
        'source_elevation_m':points[:,2],'route_index':indices,
        'route_distance_m':matched.distance_m.to_numpy(),
        'route_elevation_m':matched.elevation_m.to_numpy(),
        'projection_3d_distance_m':distance})
    return indices,table
