"""Spatial warning envelopes from reachable propagation scenarios.
An uncoupled numerical 2D test is never included in the primary hazard footprint.
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import LineString,GeometryCollection,shape
from shapely.ops import unary_union
from rasterio.features import shapes
from rasterio.transform import Affine
from matplotlib.colors import LightSource
from scipy.ndimage import gaussian_filter
from .q4_flowpath import local_dem
from .common import savefig,write_json


def execute_q3_spatial(ctx):
    out=ctx.folder('q3');q4=ctx.folder('q4')
    p=pd.read_csv(q4/'Q4_route_profile.csv')
    main=json.loads((q4/'Q4_1d_settings_and_summary.json').read_text())['mass']
    status=json.loads((q4/'Q4_2d_status.json').read_text())
    xy=p[['utm_x','utm_y']].to_numpy(float);full=LineString(xy)
    runout=float(main['max_runout_m']);width=float(ctx.cfg['q4']['width_m'])
    active=p.loc[p.distance_m<=runout].copy()
    if len(active)<2:raise ValueError('No physically propagated 1D footprint')
    current_line=LineString(active[['utm_x','utm_y']].to_numpy(float));corridor=current_line.buffer(width/2,cap_style=2,join_style=2)
    local=GeometryCollection()
    coupled=bool(status.get('primary_1d_coupled',False)) and status.get('status')!='NUMERICAL_FAILURE'
    if coupled:
        arr=np.load(q4/'Q4_2d_snapshots_for_replot.npz');depth=arr['max_depth_m'];T=Affine(*arr['affine_coeff'])
        wet=np.isfinite(depth)&(depth>.2)
        local=unary_union([shape(geom) for geom,v in shapes(wet.astype('uint8'),mask=wet,transform=T) if v])
    union=unary_union([corridor,local])
    sweep=pd.read_csv(q4/'Q4_parameter_sensitivity.csv').drop_duplicates(['varied_parameter','value'])
    phase=pd.read_csv(q4/'Q4_volume_mu_phase_diagram.csv')
    extent=float(max(sweep.runout_distance_km.max(),phase.runout_km.max())*1000)
    potential_p=p.loc[p.distance_m<=extent+1e-6]
    potential=LineString(potential_p[['utm_x','utm_y']].to_numpy(float)).buffer(width/2,cap_style=2,join_style=2)
    geos=gpd.GeoDataFrame({'zone':['Primary_reachable_corridor','Primary_coupled_2d','Primary_union','Scenario_potential_envelope'],
       'meaning':['1D h>0.2m reach, constant assumed width','empty when no physical 1D inflow','only physically linked primary components','parameter-range scenario, no calibrated probability'],
       'geometry':[corridor,local,union,potential]},crs='EPSG:32645')
    name=out/'Q3_flow_connected_warning_footprint.gpkg';name.unlink(missing_ok=True);geos.to_file(name,driver='GPKG')
    report={'status':'REACHABILITY_GATED_SCENARIO_SPATIAL_RANGE','route_length_km':float(p.distance_m.max()/1000),
       'mapped_polyline_length_km':float(full.length/1000),'primary_reached_distance_km':runout/1000,
       'sensitivity_runout_range_km':[float(sweep.runout_distance_km.min()),float(sweep.runout_distance_km.max())],
       'local_2d_inundation_area_km2':float(local.area/1e6),'upstream_assumed_corridor_width_m':width,
       'combined_warning_envelope_area_km2':float(union.area/1e6),'parameter_envelope_area_km2':float(potential.area/1e6),
       'local_2d_east_west_extent_m':float(local.bounds[2]-local.bounds[0]) if not local.is_empty else None,
       'local_2d_north_south_extent_m':float(local.bounds[3]-local.bounds[1]) if not local.is_empty else None,
       '2d_included_only_if_primary_1d_inflow':coupled,'standalone_2d_excluded':not coupled,
       'impact_geometry':str(name),'verified_disaster_boundary':False,'pre_event_dem_provenance_verified':False,
       'geometry_interpretation':'retrospective DEM scenario; reachable 1D corridor + physically coupled 2D only; potential envelope separately identified'}
    write_json(out/'Q3_spatial_range_results.json',report)
    coords=xy[::max(1,len(xy)//40)];dem,tr,_=local_dem(ctx,coords,buffer_m=1800,resolution_m=120)
    dem=np.where(np.isfinite(dem),dem,np.nanmedian(dem));hill=LightSource(315,45).hillshade(gaussian_filter(dem,1),dx=120,dy=120)
    ex=[tr.c,tr.c+dem.shape[1]*tr.a,tr.f+dem.shape[0]*tr.e,tr.f]
    fig,ax=plt.subplots(figsize=(10,6.3),layout='constrained');ax.imshow(hill,cmap='Greys',extent=ex)
    gpd.GeoSeries([potential],crs=geos.crs).boundary.plot(ax=ax,color='#c39c65',linewidth=1.,label='参数情景包络')
    gpd.GeoSeries([union],crs=geos.crs).plot(ax=ax,color='#d96a49',alpha=.5)
    ax.plot(xy[:,0],xy[:,1],color='#888',lw=1,ls='--',label='完整候选通道')
    ac=active[['utm_x','utm_y']].to_numpy(float);ax.plot(ac[:,0],ac[:,1],color='#1c648b',lw=1.8,label='基准情景实际传播段')
    ax.scatter([xy[0,0],xy[-1,0]],[xy[0,1],xy[-1,1]],c=['#164c6a','#cb5b3e'],s=55,marker='s')
    ax.text(.02,.98,f'基准传播 {runout/1000:.2f} km\n参数范围 {report["sensitivity_runout_range_km"][0]:.2f}-{report["sensitivity_runout_range_km"][1]:.2f} km\n基准条件范围 {union.area/1e6:.2f} km²',transform=ax.transAxes,va='top',fontsize=9,bbox=dict(fc='white',ec='#ccc',alpha=.9))
    ax.set(xlabel='UTM45N 东向坐标 (m)',ylabel='UTM45N 北向坐标 (m)',title='按物理可达性约束的空间预警情景');ax.set_aspect('equal');ax.legend(loc='lower right',fontsize=8)
    savefig(fig,out/'Q3_R05_single_stage_spatial_reference.png')
    return report
