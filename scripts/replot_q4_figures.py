"""Redraw Q4 maps and all four boundary-condition panels from cached numerical arrays.

Use after a purely visual edit; DO NOT rerun expensive Q4 simulation & parameter sweep.
Run from project root: python scripts/replot_q4_figures.py
"""
from pathlib import Path
import sys, json
import pandas as pd, numpy as np
from rasterio.transform import Affine
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context,configure_plotting
from mathorcup.q4_visuals import plot_local_2d,plot_inlet_factorial,plot_geographic_reference

ctx=Context(ROOT);configure_plotting();q=ROOT/'outputs'/'q4'
raw=np.load(q/'Q4_2d_snapshots_for_replot.npz')
status=json.loads((q/'Q4_2d_status.json').read_text())
tr=Affine(*raw['affine_coeff'].tolist())
frames=[{'time_s':t,'h':h,'cum_injected_m3':0} for t,h in zip(raw['time_s'],raw['depth_m'])]
ans={'z':raw['z_m'],'captures':frames,'max_depth':raw['max_depth_m'],
     'arrival':raw['arrival_s'],'timeline':pd.read_csv(q/'Q4_2d_wet_area_time_diagnostics.csv'),
     'first_inlet_s':float(raw['first_inlet_s']),'dx_m':ctx.cfg['q4']['grid_step_m'],
     'inlet':status['inlet'],'mass':status['mass']}
# Preserve the real snapshot cumulative inlet and wet area for the CSV after re-rendering.
fr_data=pd.read_csv(q/'Q4_2d_snapshot_metrics.csv')
for f,row in zip(frames,fr_data.itertuples()):f['cum_injected_m3']=float(row.cumulative_inlet_m3)
plot_local_2d(ctx,ans,tr,'ONE_WAY_1D_TO_2D_SCENARIO')
fac=q/'Q4_2d_factorial_render_data.npz'
if fac.exists():
 data=np.load(fac)
 tf=Affine(*data['affine_coeff'].tolist())
 keys=['侧坡+静止入流','侧坡+沿沟动量','谷底+静止入流','谷底+沿沟动量']
 variants={}
 for i,key in enumerate(keys):
  variants[key]={'z':data['z_m'],'captures':[{'time_s':float(data['time_s'][i]),'h':data['depth_m'][i]}],
    'max_depth':data['max_depth_m'][i],
    'inlet':{'snapped_rc':data['snapped_rc'][i].tolist()},
    'mass':{'wet_area_h_gt_0p2_km2':float(data['area_km2'][i])},
    'first_inlet_s':float(raw['first_inlet_s'])}
 plot_inlet_factorial(ctx,variants,tf)
else:
 print('Factorial cache missing: skip Q4_06C. Run once with q4.run_inlet_ablation=true to create it.')
# GIS map can be re-rendered without running any fluid mechanics
prof=pd.read_csv(q/'Q4_route_profile.csv')
plot_geographic_reference(ctx,prof)
print('Replotted Q4 figures from cached arrays, without dynamics/sensitivity recomputation')
