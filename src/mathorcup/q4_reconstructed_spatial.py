"""Nested native-DEM 2D calculation driven by exact 1D cumulative face flux."""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,LightSource
from rasterio.transform import rowcol
from .q4_solver import simulate_2d
from .q4_flowpath import local_dem
from .common import savefig,write_json,COLORS


def execute_reconstructed_spatial(ctx,lower,report,evidence):
    out=ctx.folder('q4');cfg=evidence['reconstruction'];p=lower['profile'];sensor=lower['sensor_idx']
    at=p.iloc[sensor];target=p.iloc[-1];ahead=p.iloc[min(sensor+3,len(p)-1)]
    vector=np.array([ahead.utm_x-at.utm_x,ahead.utm_y-at.utm_y],float);vector/=np.linalg.norm(vector)
    history=lower['history'];k=int(history.sensor_discharge_m3_s.to_numpy().argmax())
    height=float(history.sensor_depth_m.iloc[k]);b=float(p.effective_width_m.iloc[sensor]);m=float(p.bank_side_slope.iloc[sensor])
    inlet_width=b+2*m*height
    forcing=history[['time_s','sensor_discharge_m3_s','sensor_velocity_m_s','sensor_cumulative_m3']].copy()
    forcing.to_csv(out/'Q4_reconstruction_actual_2d_inflow.csv',index=False)
    points=[(float(at.utm_x),float(at.utm_y)),(float(target.utm_x),float(target.utm_y))]
    resolution=float(cfg['local_2d_resolution_m']);rows=[];main=None;main_transform=None
    for dx in [resolution,2*resolution]:
        z,tr,_=local_dem(ctx,points,buffer_m=550,resolution_m=dx)
        if not np.isfinite(z).all():raise ValueError('Nested 2D area contains DEM gaps')
        r,c=rowcol(tr,*points[0]);params=report['stage2_parameters']
        a=simulate_2d(z,dx,(int(r),int(c)),forcing,total_seconds=cfg['total_seconds'],
            mu=params['mu'],xi=params['xi'],cfl=.22,frames=10,flow_direction_xy=vector,
            density_kg_m3=1000.,inlet_width_m=inlet_width,inlet_depth_m=height)
        rr,cc=np.indices(z.shape);xs=tr.c+(cc+.5)*dx;ys=tr.f-(rr+.5)*dx
        radius=120.;near=(xs-points[1][0])**2+(ys-points[1][1])**2<=radius**2
        arrival=a['arrival'][near];valid=arrival[np.isfinite(arrival)]
        row={'resolution_m':dx,'mass_relative_error':a['mass']['relative_mass_error'],
           'input_m3':a['mass']['input_cumulative_m3'],'expected_input_from_1d_m3':float(forcing.sensor_cumulative_m3.iloc[-1]),
           'wet_area_km2':a['mass']['wet_area_h_gt_0p2_km2'],'max_depth_m':a['mass']['peak_depth_m'],
           'max_wet_speed_m_s':a['mass']['max_wet_cell_velocity_m_s'],
           'port_neighbourhood_arrival_s':float(valid.min()) if len(valid) else None,
           'port_evaluation_radius_m':radius,'inlet_cells':len(a['inlet']['cells']),
           'stored_final_m3':a['mass']['stored_final_m3'],'outflow_m3':a['mass']['signed_boundary_outflow_m3']}
        rows.append(row)
        if main is None:main=a;main_transform=tr
        print('[Q4 nested spatial]',row,flush=True)
    pd.DataFrame(rows).to_csv(out/'Q4_reconstruction_2d_grid_experiment.csv',index=False)
    write_json(out/'Q4_reconstruction_2d_summary.json',{
      'model_role':'native-DEM nested 2D event reconstruction','grid_experiments':rows,
      'input_source':'accepted 1D finite-volume interface flux and cumulative volume; no synthetic inlet',
      'inlet_width_from_calibrated_1d_section_m':inlet_width,'inlet_stage_from_1d_peak_flow_m':height,
      'reference_density_for_unit_conversion_kg_m3':1000,
      'displayed_intensity':'specific kinetic energy 0.5*u^2 in J/kg, not a measured mixture density',
      'stage2_friction_parameters':report['stage2_parameters'],
      'reverse_flux_at_1d_section_m3':float(history.sensor_reverse_cumulative_m3.iloc[-1])})
    tr=main_transform;a=main;z=a['z'];extent=[tr.c,tr.c+z.shape[1]*tr.a,tr.f+z.shape[0]*tr.e,tr.f]
    np.savez_compressed(out/'Q4_reconstruction_2d_arrays.npz',depth_m=a['max_depth'],arrival_s=a['arrival'],
        max_velocity_m_s=a['max_wet_velocity_m_s'],elevation_m=z,transform=np.array(tuple(tr)),
        time_s=[x['time_s'] for x in a['captures']],frames_m=[x['h'] for x in a['captures']])
    plot_reconstructed_fields(ctx,a,tr,p,points,resolution,rows,report)
    reconstruct_warning_range(ctx,lower,report,a,tr)
    return rows


def plot_reconstructed_fields(ctx,a,tr,p,points,resolution,rows,report):
    out=ctx.folder('q4');z=a['z']
    extent=[tr.c,tr.c+z.shape[1]*tr.a,tr.f+z.shape[0]*tr.e,tr.f]
    shade=LightSource(315,45).hillshade(z,dx=resolution,dy=resolution)
    fig,axs=plt.subplots(2,2,figsize=(11,9),layout='constrained')
    layers=[(a['max_depth'],'最大流深 (m)','turbo'),(a['arrival']/60,'首次到达 (min)','viridis'),
            (a['max_wet_velocity_m_s'],'最大流速 (m/s)','plasma'),
            (.5*a['max_wet_velocity_m_s']**2,'单位质量动能 (J/kg)','magma')]
    wet=a['max_depth']>.2
    for ax,(field,label,cmap),letter in zip(axs.flat,layers,'abcd'):
        ax.imshow(shade,extent=extent,cmap='Greys',alpha=.8)
        field=np.ma.masked_where(~wet,field)
        high=float(np.nanquantile(field.compressed(),.99))
        im=ax.imshow(field,extent=extent,cmap=cmap,vmin=0,vmax=max(high,.1),alpha=.85)
        ax.plot(p.utm_x,p.utm_y,color='#FFFFFF',lw=.8,alpha=.7)
        ax.scatter(*points[1],s=30,c='#B63126',marker='s',edgecolors='white',zorder=4)
        ax.set(xlabel='UTM东向 (m)',ylabel='UTM北向 (m)',title=f'{letter}  {label}')
        ax.set_xlim(extent[0],extent[1]);ax.set_ylim(extent[2],extent[3])
        fig.colorbar(im,ax=ax,shrink=.77,label=label,extend='max' if field.max()>high else 'neither')
    fig.suptitle('吉隆口岸河段的二维传播强度与到达时间')
    savefig(fig,out/'Q4_18_reconstructed_2d_depth_speed_arrival.png')
    fig,axs=plt.subplots(2,3,figsize=(12.2,7),sharex=True,sharey=True,layout='constrained')
    frames=a['captures'];indices=np.linspace(0,len(frames)-1,6).round().astype(int)
    vmax=max(1,float(np.quantile(a['max_depth'][wet],.99)))
    for ax,kframe in zip(axs.flat,indices):
        frame=frames[kframe];ax.imshow(shade,extent=extent,cmap='Greys')
        im=ax.imshow(np.ma.masked_where(frame['h']<=.2,frame['h']),extent=extent,cmap='turbo',norm=LogNorm(.2,vmax),alpha=.87)
        ax.scatter(*points[1],s=15,c='#B63126',marker='s',edgecolors='white')
        ax.set_title(f"源区释放后 {frame['time_s']/60:.2f} min",fontsize=10)
    fig.colorbar(im,ax=axs.ravel().tolist(),shrink=.8,label='等效流深 (m)')
    fig.suptitle('真实通量驱动的六个二维传播阶段')
    savefig(fig,out/'Q4_19_reconstructed_2d_six_frames.png')
    fig,axs=plt.subplots(1,3,figsize=(14,4.2),layout='constrained')
    vals=[report['upstream']['stored_m3'],report['deposit_partition_m3'],report['downstream']['stored_m3'],report['downstream']['outflow_m3']]
    axs[0].barh(['上游滞留','汇入口沉积','下游域内','下游出流'],np.array(vals)/1e6,color=COLORS[:4])
    axs[0].set(xlabel='体积 (百万m³)',title='a  全程一维体积收支')
    axs[1].barh(['二维累计输入','二维域内','二维边界出流'],np.array([rows[0]['input_m3'],rows[0]['stored_final_m3'],rows[0]['outflow_m3']])/1e6,color=COLORS[:3])
    axs[1].set(xlabel='体积 (百万m³)',title='b  局部二维体积收支')
    axs[2].plot([r['resolution_m'] for r in rows],[r['wet_area_km2'] for r in rows],'-o',color=COLORS[2])
    axs[2].set(xlabel='二维网格 (m)',ylabel='累计湿润面积 (km²)',title='c  原始DEM上的网格试验')
    fig.suptitle('级联体积守恒、实际入流传递与网格比较')
    savefig(fig,out/'Q4_20_cascade_volume_and_grid_experiments.png')


def reconstruct_warning_range(ctx,lower,report,two,transform):
    """Variable-width model corridor plus computed 2D wet cells, in one metric CRS."""
    import geopandas as gpd
    from shapely.geometry import LineString,shape,box
    from shapely.ops import unary_union
    from rasterio.features import shapes
    out=ctx.folder('q3');q4=ctx.folder('q4')
    route=pd.read_csv(q4/'Q4_reconstruction_route.csv')
    up=np.load(q4/'Q4_reconstruction_upper_arrays.npz');low=np.load(q4/'Q4_reconstruction_lower_arrays.npz')
    upper_max=up['max_depth_m'];lower_max=low['max_depth_m']
    j=len(upper_max)-1;depth=np.r_[upper_max,lower_max[1:]]
    base=route.effective_width_m.to_numpy(float,copy=True);side=route.bank_side_slope.to_numpy(float,copy=True)
    scale=report['stage2_parameters']['section_capacity_scale'];base[j:]*=scale;side[j:]*=scale
    width=base+2*side*depth;coords=route[['utm_x','utm_y']].to_numpy(float)
    pieces=[LineString(coords[i:i+2]).buffer(max(width[i],width[i+1])/2,cap_style=2)
            for i in range(len(coords)-1)]
    corridor=unary_union(pieces)
    wet=two['max_depth']>.2
    local=unary_union([shape(g) for g,v in shapes(wet.astype('uint8'),mask=wet,transform=transform) if v])
    source=np.load(q4/'Q4_reconstruction_sections.npz')['source_pixels_utm']
    release=unary_union([box(x-15,y-15,x+15,y+15) for x,y in source])
    union=unary_union([release,corridor,local])
    geometries=gpd.GeoDataFrame({'component':['source_grid','variable_width_corridor','nested_2d_wet','model_warning_envelope'],
      'geometry':[release,corridor,local,union]},crs='EPSG:32645')
    path=out/'Q3_reconstructed_warning_range.gpkg';path.unlink(missing_ok=True);geometries.to_file(path,driver='GPKG')
    result={'method':'model-based spatial range in retrospective warning replay',
      'source_to_port_corridor_km':float(route.distance_m.iloc[-1]/1000),
      'warning_envelope_km2':float(union.area/1e6),'release_grid_km2':float(release.area/1e6),
      'variable_width_corridor_km2':float(corridor.area/1e6),'nested_2d_wet_km2':float(local.area/1e6),
      'CRS':'EPSG:32645','width_basis':'computed maximum stages in DEM-derived, calibrated cross sections',
      'end_of_model_domain':'port river reach; length is a protected corridor distance, not final runout beyond the port'}
    write_json(out/'Q3_reconstructed_spatial_answer.json',result)
    fig,ax=plt.subplots(figsize=(11.5,6),layout='constrained')
    samples=coords[::max(1,len(coords)//60)].tolist()+[coords[-1].tolist()]
    dem,demtr,_=local_dem(ctx,samples,buffer_m=900,resolution_m=120)
    valid=np.isfinite(dem);surface=np.where(valid,dem,float(np.nanmedian(dem)))
    hill=LightSource(315,45).hillshade(surface,dx=120,dy=120)
    background=np.ma.masked_where(~valid,hill)
    background_extent=[demtr.c,demtr.c+dem.shape[1]*demtr.a,demtr.f+dem.shape[0]*demtr.e,demtr.f]
    ax.imshow(background,extent=background_extent,cmap='Greys',alpha=.75)
    geometries.iloc[[3]].plot(ax=ax,color='#D68B50',alpha=.4,edgecolor='#B57345',linewidth=.4)
    geometries.iloc[[0]].plot(ax=ax,color=COLORS[2],alpha=.8)
    geometries.iloc[[2]].plot(ax=ax,color=COLORS[0],alpha=.8)
    ax.plot(coords[:,0],coords[:,1],color='#3D5565',lw=1)
    ax.scatter(coords[[0,-1],0],coords[[0,-1],1],s=42,c=[COLORS[2],COLORS[1]],edgecolors='white')
    ax.set(xlabel='UTM东向 (m)',ylabel='UTM北向 (m)',title='源区、变断面廊道与口岸二维传播区')
    ax.text(.02,.98,f"源区至口岸 {result['source_to_port_corridor_km']:.2f} km\n模型警戒包络 {result['warning_envelope_km2']:.2f} km²",transform=ax.transAxes,
      va='top',bbox=dict(fc='white',ec='#cccccc',alpha=.9),fontsize=10)
    savefig(fig,out/'Q3_12_reconstructed_spatial_warning.png')
