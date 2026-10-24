"""Recompute only the four Q4 inlet scenarios; skip the expensive parameter sweep."""
from pathlib import Path
import sys,numpy as np,pandas as pd
from rasterio.transform import rowcol
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context,configure_plotting
from mathorcup.q4_flowpath import local_dem
from mathorcup.q4_solver import simulate_2d
from mathorcup.q4_visuals import plot_inlet_factorial
ctx=Context(ROOT);configure_plotting()
pdpath=ctx.out/'q4'/'Q4_route_profile.csv';prof=pd.read_csv(pdpath)
import json
sensor_dist=json.loads((ctx.out/'q4'/'Q4_1d_settings_and_summary.json').read_text())['mass']['coupling_sensor_distance_m']
pos=int(np.argmin(np.abs(prof.distance_m.to_numpy()-sensor_dist)))
sens=prof.iloc[pos];dest=prof.iloc[-1]
zz,tr,_=local_dem(ctx,[(sens.utm_x,sens.utm_y),(dest.utm_x,dest.utm_y)],buffer_m=500,resolution_m=ctx.cfg['q4']['grid_step_m'])
rc=tuple(int(v) for v in rowcol(tr,float(sens.utm_x),float(sens.utm_y)))
nextp=prof.iloc[min(pos+3,len(prof)-1)];dr=np.array([nextp.utm_x-sens.utm_x,nextp.utm_y-sens.utm_y]);dr/=np.linalg.norm(dr)
hydro=pd.read_csv(ctx.out/'q4'/'Q4_1d_discharge_and_mass_timeseries.csv')
four={}
for label,mode in [('侧坡+静止入流','legacy_3x3'),('侧坡+沿沟动量','legacy_directional'),('谷底+静止入流','valley_zero_momentum'),('谷底+沿沟动量','valley_cross_section')]:
 print('Computing:',label,flush=True)
 four[label]=simulate_2d(zz,ctx.cfg['q4']['grid_step_m'],rc,hydro,total_seconds=ctx.cfg['q4']['total_seconds'],
                         mu=ctx.cfg['q4']['mu'],xi=ctx.cfg['q4']['xi_m_s2'],cfl=.22,frames=10,inlet_mode=mode,flow_direction_xy=dr)
plot_inlet_factorial(ctx,four,tr)
print('Saved four-condition visual cache and figure',flush=True)
