"""Replot the calibrated reconstruction from recorded arrays and experiment tables."""
from pathlib import Path
import sys,json,yaml
import numpy as np
import pandas as pd
from affine import Affine
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context,configure_plotting
from mathorcup.q4_terrain_reconstruction import prepare_reconstruction_geometry
from mathorcup.q4_cascade_reconstruction import plot_reconstruction,reconstruction_deposition_geometry
from mathorcup.q4_reconstructed_spatial import plot_reconstructed_fields,reconstruct_warning_range

def main():
    ctx=Context(ROOT);configure_plotting();out=ctx.folder('q4')
    evidence=yaml.safe_load((ROOT/'configs/case_evidence.yaml').read_text())
    report=json.loads((out/'Q4_cascade_reconstruction.json').read_text())
    geometry=json.loads((out/'Q4_reconstruction_geometry.json').read_text())
    terrain=np.load(out/'Q4_reconstruction_sections.npz')
    if 'background_elevation_m' not in terrain:
        prepare_reconstruction_geometry(ctx,evidence);terrain=np.load(out/'Q4_reconstruction_sections.npz')
    p=pd.read_csv(out/'Q4_reconstruction_route.csv',float_precision='round_trip');j=geometry['confluence_index']
    answers=[]
    for name,part in [('upper',p.iloc[:j+1]),('lower',p.iloc[j:])]:
        data=np.load(out/f'Q4_reconstruction_{name}_arrays.npz')
        answers.append({'s':data['distance_m'],'z':data['elevation_m'],'width':data['width_m'],
          'captures':[{'time_s':t,'h':h} for t,h in zip(data['time_s'],data['depth_m'])],
          'history':pd.read_csv(out/f'Q4_reconstruction_{name}_history.csv',float_precision='round_trip'),
          'mass':report['upstream' if name=='upper' else 'downstream'],'profile':part})
    trials=pd.read_csv(out/'Q4_reconstruction_calibration_trials.csv',float_precision='round_trip').loc[lambda d:d.stage.eq(2)]
    tr=Affine(*terrain['background_transform'][:6])
    plot_reconstruction(ctx,p,j,terrain['background_elevation_m'],tr,geometry,*answers,report,trials)
    ar=np.load(out/'Q4_reconstruction_2d_arrays.npz');two={'z':ar['elevation_m'],'max_depth':ar['depth_m'],
      'arrival':ar['arrival_s'],'max_wet_velocity_m_s':ar['max_velocity_m_s'],
      'captures':[{'time_s':t,'h':h} for t,h in zip(ar['time_s'],ar['frames_m'])]}
    lower=p.iloc[j:];i=int(np.argmin(abs(lower.distance_m-(lower.distance_m.iloc[-1]-2500))))
    points=[tuple(lower[['utm_x','utm_y']].iloc[i]),tuple(lower[['utm_x','utm_y']].iloc[-1])]
    tr=Affine(*ar['transform'][:6]);summary=json.loads((out/'Q4_reconstruction_2d_summary.json').read_text())
    plot_reconstructed_fields(ctx,two,tr,lower,points,evidence['reconstruction']['local_2d_resolution_m'],summary['grid_experiments'],report)
    reconstruct_warning_range(ctx,None,report,two,tr);reconstruction_deposition_geometry(ctx)
    print('Reconstruction figures refreshed from recorded numerical arrays')

if __name__=='__main__':main()
