"""Q4 grid-resolution convergence / discretisation check with real DEM.
Uses same upstream 1D discharge curve, not additional data. Computes physically
separate 60m and 120m ROI tests; reports area/depth sensitivity honestly.
"""
from __future__ import annotations
import numpy as np,pandas as pd,matplotlib.pyplot as plt,json
from rasterio.transform import rowcol
from rasterio.warp import reproject,Resampling
from .q4_flowpath import local_dem
from .q4_solver import simulate_2d
from .common import savefig,write_json,COLORS

def execute_q4_grid_review(ctx,profile=None):
    out=ctx.folder('q4');cfg=ctx.cfg['q4'];fn=out/'Q4_1d_discharge_and_mass_timeseries.csv'
    if not fn.exists():return {'status':'NO_ONE_D_HYDROGRAPH'}
    h=pd.read_csv(fn)
    geo=json.loads((out/'Q4_2d_inlet_geometry.json').read_text())
    coupled=bool(geo['coupled_to_primary_1d'])
    actual=out/'Q4_actual_2d_inflow_hydrograph.csv'
    if actual.exists():h=pd.read_csv(actual)
    if 'sensor_discharge_m3_s' not in h or h.sensor_discharge_m3_s.max()<=.01:
        fallback=out/'Q4_STANDALONE_SYNTHETIC_2d_inflow.csv'
        if not fallback.exists():return {'status':'NO_COUPLING_INFLOW'}
        h=pd.read_csv(fallback);coupled=False
    if profile is None:profile=pd.read_csv(out/'Q4_route_profile.csv')
    target_dist=profile.distance_m.iloc[-1]-2500
    sens=profile.iloc[(profile.distance_m-target_dist).abs().argmin()]
    dest=profile.iloc[-1]
    # Use the same local tangent as the primary solver, not a different inlet direction.
    sensor_idx=int((profile.distance_m-target_dist).abs().argmin())
    ahead=profile.iloc[min(sensor_idx+3,len(profile)-1)]
    pdir=np.array([ahead.utm_x-sens.utm_x,ahead.utm_y-sens.utm_y],float)
    pdir/=max(np.linalg.norm(pdir),1e-12)
    geo=json.loads((out/'Q4_2d_inlet_geometry.json').read_text())
    sensor_xy=tuple(geo['sensor_utm_xy']);target_xy=tuple(geo['target_utm_xy']);pdir=np.asarray(geo['flow_direction_xy'])
    rows=[];fields=[]
    for resolution in [60,120]:
        z,tr,_=local_dem(ctx,[sensor_xy,target_xy],buffer_m=500,resolution_m=resolution)
        r,c=rowcol(tr,*sensor_xy)
        ans=simulate_2d(z,resolution,(int(r),int(c)),h,total_seconds=cfg['total_seconds'],
             mu=cfg['mu'],xi=cfg['xi_m_s2'],cfl=.22,frames=3,flow_direction_xy=pdir,
             density_kg_m3=cfg.get('density_kg_m3',1800.))
        m=ans['mass']
        rows.append({'resolution_m':resolution,'primary_1d_coupled':coupled,'wet_area_km2':m['wet_area_h_gt_0p2_km2'],
             'max_depth_m':m['peak_depth_m'],'mass_balance_rel_error':m['relative_mass_error'],
             'first_input_time_s':float(ans['first_inlet_s']),'n_steps':int(m['n_steps'])})
        fields.append((ans['max_depth'],tr))
    d=pd.DataFrame(rows);d.to_csv(out/'Q4_grid_resolution_reliability.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(9.5,4),layout='constrained')
    axes[0].plot(d.resolution_m,d.wet_area_km2,'-o',color=COLORS[0],lw=1.8)
    axes[0].set(xlabel='网格分辨率 (m)',ylabel='模拟淹没面积 (km²)',title='a  淹没范围')
    axes[1].plot(d.resolution_m,d.max_depth_m,'-s',color=COLORS[1],lw=1.8)
    axes[1].set(xlabel='网格分辨率 (m)',ylabel='最大流深 (m)',title='b  最大流深')
    for a in axes:a.grid(axis='y',alpha=.17)
    fig.suptitle('二维传播模型网格分辨率敏感性',fontsize=12)
    savefig(fig,out/'Q4_09_grid_resolution_reliability.png')
    status={'status':'COMPLETED_60_AND_120_M' if coupled else 'SECONDARY_SCENARIO_60_AND_120_M','results':rows,'primary_1d_coupled':coupled,
       'limitations':'Two-grid sensitivity does not prove asymptotic grid convergence; model boundary and DEM coarsening remain important'}
    write_json(out/'Q4_grid_resolution_status.json',status)
    print('[Q4-grid] 60m/120m',rows)
    return status
