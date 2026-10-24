"""Paper-oriented, geographically referenced Q4 visualizations (2026 MathorCup B).

All maps derive from official TPDC DEM #7; no decorative satellite texture is passed
as geographical evidence. The external USGS 2026-08-27 public-domain disaster map
is linked as an independent qualitative reference in the source notes, rather than
being geometrically conflated with our model's footprint.

Design goals: equal-sized image panels, one fixed flow-depth scale, event-centred
snapshot timing, no in-panel text overlaps, legible geographic coordinates and
explicit observed/simulated provenance.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import xy
from matplotlib import pyplot as plt
from matplotlib import colors,ticker,patheffects
from matplotlib.colors import LightSource
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter, distance_transform_edt
from pyproj import Transformer
from .common import savefig,write_json,COLORS,record_warning

USGS_URL='https://www.usgs.gov/media/images/2026-nepal-debris-avalanche-and-flash-flood-map'

def extent_utm(tr,shape):
    return [float(tr.c),float(tr.c+tr.a*shape[1]),
            float(tr.f+tr.e*shape[0]),float(tr.f)]

def hillshade(z,dx):
    z=np.asarray(z,float)
    mask=~np.isfinite(z)
    if mask.any():
        _,ind=distance_transform_edt(mask,return_indices=True)
        z=z[tuple(ind)]
    smoothing=gaussian_filter(z,sigma=.6)
    ls=LightSource(azdeg=310,altdeg=42)
    hs=ls.hillshade(smoothing,vert_exag=.85,dx=dx,dy=dx)
    return .78+.22*hs

def decorate_geo(ax,xlabel=True,ylabel=True):
    ax.set_aspect('equal',adjustable='box')
    ax.tick_params(axis='both',labelsize=8,colors='#4e5661',length=2)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos:f'{x/1000:.1f}'))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos:f'{y/1000:.1f}'))
    if xlabel:ax.set_xlabel('UTM45N 东向坐标 (km)',fontsize=9,labelpad=5)
    if ylabel:ax.set_ylabel('UTM45N 北向坐标 (km)',fontsize=9,labelpad=5)
    for spine in ax.spines.values():spine.set_color('#D2D7DA');spine.set_linewidth(.8)

def _geo_north_and_scale(ax,extent,scale_m=1000):
    x0,x1,y0,y1=extent
    # north points toward larger UTM northing; preserve map geometry
    x=x0+(x1-x0)*.90;y=y0+(y1-y0)*.72
    ax.annotate('N',xy=(x,y+.09*(y1-y0)),xytext=(x,y),
        fontsize=9,weight='bold',ha='center',va='center',color='#172E3F',
        arrowprops={'arrowstyle':'-|>','color':'#172E3F','lw':1.2})
    sx=x0+.07*(x1-x0);sy=y0+.10*(y1-y0)
    size=min(scale_m,(x1-x0)*.22)
    ax.plot([sx,sx+size],[sy,sy],lw=3,color='#293845',solid_capstyle='butt',zorder=30)
    ax.text(sx+size*.5,sy+.02*(y1-y0),f'{size/1000:g} km',fontsize=8,ha='center',
            color='#293845',bbox=dict(fc='white',ec='none',alpha=.85,pad=1),zorder=30)

def _common_map_bg(ax,z,tr,scaling=None,alpha=1,palette='gray'):
    ext=extent_utm(tr,z.shape)
    if scaling is None:
        zz=z[np.isfinite(z)];scaling=(np.percentile(zz,2),np.percentile(zz,98))
    if palette=='hypsometric':
        ax.imshow(z,origin='upper',extent=ext,cmap='terrain',vmin=scaling[0],vmax=scaling[1],alpha=.92*alpha,interpolation='bilinear')
        ax.imshow(hillshade(z,abs(tr.a)),origin='upper',extent=ext,cmap='gray',vmin=.75,vmax=1,
                  alpha=.28*alpha,interpolation='bilinear')
    else:
        ax.imshow(z,origin='upper',extent=ext,cmap='Greys_r',vmin=scaling[0],vmax=scaling[1],alpha=.52*alpha,interpolation='bilinear')
        ax.imshow(hillshade(z,abs(tr.a)),origin='upper',extent=ext,cmap='gray',vmin=.75,vmax=1,
                  alpha=.43*alpha,interpolation='bilinear')
    ax.set_xlim(ext[0],ext[1]);ax.set_ylim(ext[2],ext[3])
    return ext

def _depth_norm(max_depth):
    vals=max_depth[np.isfinite(max_depth)&(max_depth>.2)]
    vmax=float(np.percentile(vals,99.5)) if len(vals) else 1.
    return colors.LogNorm(vmin=.2,vmax=max(vmax,.4))

def _mask_flow(depth):
    return np.ma.masked_less_equal(np.asarray(depth,float),.2)

def _make_axes_equal(axes, box_aspect=None):
    for ax in np.ravel(axes):
        if box_aspect is not None:ax.set_box_aspect(box_aspect)
        ax.tick_params(axis='both',labelsize=8)
        for sp in ax.spines.values():sp.set_linewidth(.65);sp.set_color('#B7C1C7')

def plot_local_2d(ctx,ans,tr,scenario_label='ONE_WAY_1D_TO_2D_SCENARIO'):
    """Publishable Q4_05 (two panels), Q4_06 (10 frames), and inlet timing diagnostic."""
    out=ctx.folder('q4')
    z=ans['z'];hmax=ans['max_depth'];arrival=ans['arrival'];frames=ans['captures']
    coupling_label='一维计算的通量驱动二维模型' if 'ONE_WAY' in scenario_label else '独立给定入流的二维数值试验'
    extent=extent_utm(tr,z.shape);norm=_depth_norm(hmax)
    cm=plt.get_cmap('inferno').copy();cm.set_bad((1,1,1,0))
    wet_area=np.array([np.sum(f['h']>.2)*ans['dx_m']**2/1e6 for f in frames])
    if np.max(wet_area)<=0:raise RuntimeError('No flow shown: refusing a 10-frame empty figure')
    # All frames informative after inlet entry, displayed from an identical geographic view.
    assert len(frames)==10, 'Paper layout expects exactly ten snapshots'
    if np.any(wet_area<=0):
        record_warning('Q4 frame contains no wet cells; this is a genuine solver or capture-time issue')

    mark=ans['inlet']['snapped_rc'];sx,sy=xy(tr,int(mark[0]),int(mark[1]),offset='center')
    old=ans['inlet']['nominal_rc'];oldx,oldy=xy(tr,int(old[0]),int(old[1]),offset='center')
    
    # Publication dual map: identical axes allocation, independent aligned horizontal colorbars.
    fig,axs=plt.subplots(1,2,figsize=(13.8,6.2))
    fig.subplots_adjust(left=.055,right=.965,top=.82,bottom=.23,wspace=.15)
    _make_axes_equal(axs,box_aspect=1.)
    for a in axs:
        _common_map_bg(a,z,tr)
        decorate_geo(a)
        a.plot(sx,sy,'^',color='#2F6477',markersize=7,markeredgecolor='white',zorder=15)
    im=axs[0].imshow(_mask_flow(hmax),extent=extent,origin='upper',cmap=cm,norm=norm,
                     interpolation='nearest',alpha=.91)
    arrvals=arrival[np.isfinite(arrival)]
    armin=float(arrvals.min()/60) if len(arrvals) else 0
    armax=float(np.percentile(arrvals,99)/60) if len(arrvals) else 1
    arrnorm=colors.Normalize(vmin=armin,vmax=max(armax,armin+.1))
    am=plt.get_cmap('cividis').copy();am.set_bad((1,1,1,0))
    amim=axs[1].imshow(np.ma.masked_invalid(arrival/60),extent=extent,origin='upper',
         norm=arrnorm,cmap=am,alpha=.87,interpolation='nearest')
    axs[0].set_title('a  最大流深',fontsize=12,pad=10,loc='left',weight='medium')
    axs[1].set_title('b  首次淹没时间',fontsize=12,pad=10,loc='left',weight='medium')
    # One precisely placed horizontal colorbar per panel: no panel-width distortion.
    cb1=fig.add_axes([.095,.145,.355,.017]);cb2=fig.add_axes([.565,.145,.355,.017])
    cbar=fig.colorbar(im,cax=cb1,orientation='horizontal')
    cbar.set_label('最大流深 (m，按对数色标)',fontsize=9,labelpad=4)
    cbar.ax.tick_params(labelsize=8)
    cbar2=fig.colorbar(amim,cax=cb2,orientation='horizontal')
    cbar2.set_label('相对源区释放的分钟数 (min)',fontsize=9,labelpad=4)
    cbar2.ax.tick_params(labelsize=8)
    fig.suptitle('泥石流二维最大流深及到达时间',fontsize=15,y=.955)
    fig.text(.5,.048,f'▲  横断面调整后的二维入流点   |   {coupling_label}   |   情景模拟',
             ha='center',fontsize=9,color='#4b5965')
    savefig(fig,out/'Q4_05_local_2d_depth_and_arrival.png')

    # 10 consistent panels: reserved title, panel labels, a common depth legend and caption outside map grid.
    fig,axs=plt.subplots(2,5,figsize=(17.4,7.7))
    fig.subplots_adjust(left=.048,right=.985,top=.885,bottom=.185,wspace=.12,hspace=.24)
    _make_axes_equal(axs,box_aspect=.76)
    for i,(ax,fr) in enumerate(zip(axs.flat,frames)):
        _common_map_bg(ax,z,tr)
        ax.imshow(_mask_flow(fr['h']),extent=extent,origin='upper',
                  cmap=cm,norm=norm,alpha=.96,interpolation='nearest')
        ax.scatter([sx],[sy],marker='^',s=22,c='#126B79',edgecolors='white',linewidths=.5,zorder=10)
        absolute_min=fr['time_s']/60
        relative_min=max(0,(fr['time_s']-ans['first_inlet_s'])/60)
        ax.set_title(f'{chr(97+i)}  t = {absolute_min:.1f} min',
                     fontsize=10.1,pad=8,color='#203E56')
        if i>=5:ax.set_xlabel('东向坐标 (km)',fontsize=8,labelpad=3)
        if i%5==0:ax.set_ylabel('北向坐标 (km)',fontsize=8,labelpad=3)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos:f'{x/1000:.1f}'))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos:f'{y/1000:.1f}'))
        ax.tick_params(labelsize=7,length=2)
        if i%5!=0:ax.tick_params(labelleft=False)
        if i<5:ax.tick_params(labelbottom=False)
    cax=fig.add_axes([.19,.102,.63,.022])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cm),cax=cax,orientation='horizontal')
    cb.set_label('等效流深 h (m，对数刻度)',fontsize=10,labelpad=5)
    cb.ax.tick_params(labelsize=8)
    fig.suptitle('泥石流二维流深演进过程',fontsize=15,y=.972)
    fig.text(.5,.030,
       f'标注说明：0 min 为源区一维释放启动；二维入流始于约 {ans["first_inlet_s"]/60:.2f} min；▲ 为二维注入位置。官方 DEM（数据集⑦），60 m 栅格。',
       ha='center',fontsize=9,color='#50616A')
    savefig(fig,out/'Q4_06_local_2d_10frames.png')

    # Publish intermediate arrays: the picture can be checked quantitatively / re-rendered.
    np.savez_compressed(out/'Q4_2d_snapshots_for_replot.npz',
        depth_m=np.stack([f['h'] for f in frames]),time_s=np.array([f['time_s'] for f in frames]),
        z_m=z,max_depth_m=hmax,arrival_s=arrival,wet_area_km2=wet_area,
        affine_coeff=np.array(tuple(tr)[:6]),first_inlet_s=ans['first_inlet_s'])
    ans['timeline'].to_csv(out/'Q4_2d_wet_area_time_diagnostics.csv',index=False)
    pd.DataFrame({'time_s':[f['time_s'] for f in frames],
                  'time_since_2d_onset_s':[f['time_s']-ans['first_inlet_s'] for f in frames],
                  'wet_area_km2':wet_area,
                  'max_depth_m':[float(np.max(f['h'])) for f in frames],
                  'cumulative_inlet_m3':[f['cum_injected_m3'] for f in frames]}).to_csv(out/'Q4_2d_snapshot_metrics.csv',index=False)

    hist=ans['timeline'];tt=hist.time_s/60
    fig,axs=plt.subplots(1,2,figsize=(12.2,4.4))
    fig.subplots_adjust(left=.09,right=.965,bottom=.25,top=.75,wspace=.25)
    for a in axs:a.grid(True,alpha=.13);a.axvline(ans['first_inlet_s']/60,color='#AE7D42',lw=1.3,ls='--')
    axs[0].plot(tt,hist.inlet_m3_s/1000,color='#25638B',lw=1.7)
    axs[0].fill_between(tt,hist.inlet_m3_s/1000,alpha=.13,color='#25638B')
    axs[0].set(xlabel='源区释放后时间 (min)',ylabel='耦合断面流量 (千 m³/s)',title='a  耦合断面流量过程')
    axs[1].plot(tt,hist.wet_area_km2,color='#BF613E',lw=1.8)
    axs[1].set(xlabel='源区释放后时间 (min)',ylabel='二维湿润面积 (km², h>0.2 m)',title='b  二维淹没面积变化')
    for f in frames:
        axs[1].axvline(f['time_s']/60,color='#727F86',lw=.55,alpha=.33)
    fig.suptitle('耦合断面入流与淹没面积变化',fontsize=13.5,y=.94)
    fig.text(.5,.07,'竖线：二维首次显著入流时刻及重新选择的 10 个关键帧；空间扩展随流量连续演化，旧图为时间采样过稀。',
             fontsize=9,ha='center',color='#52636E')
    savefig(fig,out/'Q4_06A_onset_and_wet_area_diagnostic.png')

def plot_2d_inlet_diagnostic(ctx,updated,legacy,tr):
    """Compare 2D branch artifact under equal forcing but different inlet placements."""
    out=ctx.folder('q4');extent=extent_utm(tr,updated['z'].shape)
    # Use one common *global* visual scale for a strict controlled comparison.
    peak=np.maximum(updated['max_depth'],legacy['max_depth'])
    norm=_depth_norm(peak);cm=plt.get_cmap('inferno').copy();cm.set_bad((1,1,1,0))
    # choose representative same *absolute* source time just after entry
    target_time=updated['first_inlet_s']+70
    newfr=min(updated['captures'],key=lambda f:abs(f['time_s']-target_time))
    oldfr=min(legacy['captures'],key=lambda f:abs(f['time_s']-target_time))
    fig,axs=plt.subplots(1,2,figsize=(13.2,6.4))
    fig.subplots_adjust(left=.065,right=.965,bottom=.28,top=.81,wspace=.16)
    _make_axes_equal(axs,box_aspect=1)
    for ax,ans,record,title in [
       (axs[0],legacy,oldfr,'a  侧坡均匀入流'),
       (axs[1],updated,newfr,'b  谷底定向入流')]:
        _common_map_bg(ax,ans['z'],tr)
        ax.imshow(_mask_flow(record['h']),cmap=cm,norm=norm,origin='upper',
                  extent=extent,interpolation='nearest',alpha=.92)
        r,c=ans['inlet']['snapped_rc']
        x,y=xy(tr,int(r),int(c),offset='center')
        ax.scatter([x],[y],marker='^',s=68,facecolor='#1B7485',edgecolors='white',linewidths=.8,zorder=12)
        decorate_geo(ax)
        ax.set_title(title,fontsize=12,pad=9,loc='left')
        _geo_north_and_scale(ax,extent,scale_m=500)
    cax=fig.add_axes([.24,.175,.52,.019]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cm),
           cax=cax,orientation='horizontal');cb.set_label('相同时间的流深（m；统一对数色标）',fontsize=9)
    fig.suptitle('不同二维入流条件下的传播范围',fontsize=14,y=.96)
    fig.text(.5,.045,f'比较时刻 t≈{target_time/60:.1f} min；▲ 为二维入流点。侧向分叉主要来自旧边界条件；真实支流仍以河网与遥感判读为准。',
             ha='center',fontsize=9,color='#4c5b66')
    savefig(fig,out/'Q4_06B_inlet_position_ablation.png')
    oldr,oldc=legacy['inlet']['nominal_rc'];newr,newc=updated['inlet']['snapped_rc']
    snap_dist=float(np.hypot(newr-oldr,newc-oldc)*updated['dx_m'])
    write_json(out/'Q4_2d_inlet_branch_diagnostic.json',{
      'inlet_original_rc':[oldr,oldc],'inlet_valley_rc':[newr,newc],
      'lateral_snap_distance_m':snap_dist,
      'legacy_peak_depth_m':legacy['mass']['peak_depth_m'],
      'updated_peak_depth_m':updated['mass']['peak_depth_m'],
      'legacy_flooded_area_km2':legacy['mass']['wet_area_h_gt_0p2_km2'],
      'updated_flooded_area_km2':updated['mass']['wet_area_h_gt_0p2_km2'],
      'interpretation':'Local branch may be numerical inlet artefact when 3x3 injection intersects side slope; no factual claim of geological bifurcation.'})

def plot_geographic_reference(ctx,profile):
    """Two independent, strictly georeferenced real-terrain maps.

    Different map scales intentionally use separate figures rather than one oddly
    sized multi-panel canvas. A full-valley map is landscape; local DEM is near-square.
    Both use official TPDC DEM #7, not procedural terrain.
    """
    from .q4_flowpath import local_dem
    from .prepare_geo import glacier_grid
    import geopandas as gpd
    out=ctx.folder('q4');beg=profile.iloc[0];end=profile.iloc[-1]
    full,tr,_=local_dem(ctx,[(beg.utm_x,beg.utm_y),(end.utm_x,end.utm_y)],
                        buffer_m=2200,resolution_m=90)
    loc=profile.iloc[max(0,len(profile)-42)]
    local,tr2,_=local_dem(ctx,[(loc.utm_x,loc.utm_y),(end.utm_x,end.utm_y)],
                         buffer_m=500,resolution_m=60)
    ext=extent_utm(tr,full.shape);ext2=extent_utm(tr2,local.shape)
    xx=profile.utm_x.to_numpy();yy=profile.utm_y.to_numpy()

    # 1: regional geographic map (equal physical aspect, no squeezed panel).
    fig,ax=plt.subplots(figsize=(14.7,6.4))
    fig.subplots_adjust(left=.075,right=.89,top=.83,bottom=.18)
    terrain_vals=full[np.isfinite(full)]
    terrain_limits=(np.percentile(terrain_vals,2),np.percentile(terrain_vals,98))
    _common_map_bg(ax,full,tr,scaling=terrain_limits,palette='hypsometric')
    ax.plot(xx,yy,color='#C87551',lw=2.0,ls='--',zorder=8)
    ax.scatter([beg.utm_x,end.utm_x],[beg.utm_y,end.utm_y],
         s=[100,90],facecolor=['#0C7588','#D29336'],edgecolor='white',linewidth=1.25,zorder=12)
    try:
        gg=glacier_grid(ctx);mask=np.array(gg['glacier_mask'],dtype=float)
        xgrid=np.asarray(gg['x']);ygrid=np.asarray(gg['y'])
        if len(xgrid)==mask.shape[1] and len(ygrid)==mask.shape[0]:
            ax.contour(xgrid,ygrid,mask,levels=[.5],colors='#13718D',linewidths=1.5,alpha=.95)
    except Exception as exc:record_warning('Glacier overlay skipped: '+str(exc))
    try:
        historic=gpd.read_file(ctx.processed/'historic_gullies.gpkg').to_crs('EPSG:32645')
        clip=historic.cx[ext[0]:ext[1],ext[2]:ext[3]]
        if len(clip):clip.boundary.plot(ax=ax,color='#9E8977',lw=.45,alpha=.65)
    except Exception as exc:record_warning('Historic polygons skipped: '+str(exc))
    ax.add_patch(Rectangle((ext2[0],ext2[2]),ext2[1]-ext2[0],ext2[3]-ext2[2],
        ec='#168B91',fc='none',lw=1.4,zorder=10))
    ax.annotate('官方冰川掩膜中心\n（候选源区）',xy=(beg.utm_x,beg.utm_y),xytext=(6,-30),
        textcoords='offset points',fontsize=9,color='#164459',
        bbox={'facecolor':'white','alpha':.84,'edgecolor':'none','pad':2})
    ax.annotate('候选口岸位置',xy=(end.utm_x,end.utm_y),xytext=(8,-27),
        textcoords='offset points',fontsize=9,color='#6B3E17',
        bbox={'facecolor':'white','alpha':.86,'edgecolor':'none','pad':2})
    decorate_geo(ax);_geo_north_and_scale(ax,ext,scale_m=2000)
    ax.set_title('地形高程、冰川分布与历史灾点',loc='left',fontsize=12,pad=9)
    cbx=fig.add_axes([.903,.30,.015,.36])
    sm=plt.cm.ScalarMappable(norm=colors.Normalize(*terrain_limits),cmap='terrain')
    cb=fig.colorbar(sm,cax=cbx);cb.set_label('高程 (m)',fontsize=9);cb.ax.tick_params(labelsize=8)
    fig.suptitle('吉隆冰崩—泥石流研究区地形图',fontsize=14.5,y=.97)
    legend=[Line2D([0],[0],color='#C87551',lw=2,ls='--',label='推算流路（非实测线）'),
            Line2D([0],[0],color='#13718D',lw=1.5,label='官方冰川边界'),
            Line2D([0],[0],color='#9E8977',lw=1,label='历史沟谷'),
            Line2D([0],[0],color='#168B91',lw=1.5,label='二维计算范围')]
    fig.legend(handles=legend,loc='lower center',ncol=4,bbox_to_anchor=(.5,.098),frameon=False,fontsize=9)
    fig.text(.5,.055,'数据：国家青藏高原科学数据中心 DEM ⑦、冰川④、历史灾害⑥；真实灾害范围对照：USGS（2026-08-27，网址见 README）。',
             ha='center',fontsize=8.5,color='#596674')
    savefig(fig,out/'Q4_00_actual_DEM_geographic_context.png')

    # 2: local inset *in another equal-sized map canvas*, to protect map geometry.
    fig,ax=plt.subplots(figsize=(10.2,7.1))
    fig.subplots_adjust(left=.1,right=.95,top=.84,bottom=.18)
    local_vals=local[np.isfinite(local)]
    local_limits=(np.percentile(local_vals,2),np.percentile(local_vals,98))
    _common_map_bg(ax,local,tr2,scaling=local_limits,palette='hypsometric')
    ax.plot(xx,yy,color='#C87551',lw=2,ls='--',zorder=7)
    ax.scatter([loc.utm_x,end.utm_x],[loc.utm_y,end.utm_y],
               s=[100,85],marker='^',facecolor=['#0C7588','#D29336'],
               edgecolor='white',linewidth=1.2,zorder=9)
    decorate_geo(ax);_geo_north_and_scale(ax,ext2,scale_m=500)
    ax.set_title('下游局部地形及入流断面',fontsize=12,loc='left',pad=9)
    fig.suptitle('下游二维计算区域地形图',fontsize=14,y=.968)
    fig.text(.5,.055,'灰度及明暗地形均由官方高程栅格直接计算；红虚线为地形算法推算路径，而非卫星识别的实际灾害边界。',
             ha='center',fontsize=9,color='#586672')
    savefig(fig,out/'Q4_00B_local_true_DEM_detail.png')
    write_json(out/'Q4_00_map_sources.json',{
        'actual_geodata':['TPDC DEM dataset #7','TPDC glacier dataset #4','TPDC historical gullies dataset #6'],
        'independent_reference':USGS_URL,
        'route_line_is_observed':False,
        'real_georeferenced_topography_not_artificial_satellite_rgb':True,
        'local_map_resolution_m':60})

def plot_inlet_factorial(ctx,variants,tr):
    """2x2 causal/numerical ablation: source placement and initial momentum separately."""
    out=ctx.folder('q4')
    ordered=['侧坡+静止入流','侧坡+沿沟动量','谷底+静止入流','谷底+沿沟动量']
    # Compare all four at the very same absolute time, approx 70 seconds after onset.
    base=variants[ordered[0]]
    t0=base['first_inlet_s']+70.
    frames=[min(variants[key]['captures'],key=lambda f:abs(f['time_s']-t0)) for key in ordered]
    top=np.maximum.reduce([variants[k]['max_depth'] for k in ordered])
    norm=_depth_norm(top);cm=plt.get_cmap('inferno').copy();cm.set_bad((1,1,1,0))
    z=base['z'];ext=extent_utm(tr,z.shape)
    fig,axes=plt.subplots(2,2,figsize=(11.7,9.4))
    fig.subplots_adjust(left=.09,right=.965,top=.88,bottom=.235,wspace=.18,hspace=.30)
    for k,ax,fr in zip(ordered,axes.flat,frames):
        _common_map_bg(ax,z,tr)
        ax.imshow(_mask_flow(fr['h']),cmap=cm,norm=norm,extent=ext,origin='upper',
                  interpolation='nearest',alpha=.94)
        src=variants[k]['inlet']['snapped_rc']
        x,y=xy(tr,int(src[0]),int(src[1]),offset='center')
        ax.scatter([x],[y],marker='^',s=44,facecolor='#1B6D7C',edgecolor='white',linewidths=.8)
        decorate_geo(ax)
        area=variants[k]['mass']['wet_area_h_gt_0p2_km2']
        ax.set_title(f'{chr(97+ordered.index(k))}  {k}  |  全程面积 {area:.3f} km²',
                     fontsize=10.7,loc='left',pad=8)
    cax=fig.add_axes([.23,.12,.54,.017]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cm),
                  orientation='horizontal',cax=cax)
    cb.set_label('同一时刻的流深 (m，固定对数色标)',fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.suptitle('二维入流位置与动量的组合试验',fontsize=14,y=.965)
    fig.text(.5,.030,'使用同一 DEM、流量 Q(t)、摩擦与网格。此图检验入流位置与动量对分叉的独立作用；真实河网另需遥感证据。',
             ha='center',fontsize=8.8,color='#526572')
    savefig(fig,out/'Q4_06C_2x2_inlet_boundary_factorial.png')
    np.savez_compressed(out/'Q4_2d_factorial_render_data.npz',
        depth_m=np.stack([f['h'] for f in frames]),
        max_depth_m=np.stack([variants[k]['max_depth'] for k in ordered]),
        z_m=z,time_s=np.array([f['time_s'] for f in frames]),
        snapped_rc=np.stack([variants[k]['inlet']['snapped_rc'] for k in ordered]),
        area_km2=np.array([variants[k]['mass']['wet_area_h_gt_0p2_km2'] for k in ordered]),
        affine_coeff=np.array(tuple(tr)[:6]))
