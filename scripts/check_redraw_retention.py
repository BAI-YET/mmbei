"""Verify every delivered figure has machine-readable, lossless plotting inputs."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
base=ROOT/'outputs/redraw_data';metas=sorted(base.rglob('figure_metadata.json'));expected=[p for p in (ROOT/'outputs').rglob('*.png') if 'redraw_data' not in p.parts]
assert len(metas)==len(expected)>0
assert {json.loads(p.read_text())['figure'] for p in metas}=={str(p.relative_to(ROOT/'outputs')) for p in expected}
resolved=0
for p in metas:
    d=json.loads(p.read_text());assert d['datasets']
    for item in d['datasets']:
        path=(base/item['file']) if item['file'].startswith('shared_rasters/') else p.parent/item['file']
        assert path.exists(),str(path);resolved+=1
    for ax in d['axes']:
        for im in ax['rasters']:
            a=np.load(base/im['shared_data']/'array.npz');assert list(a['values'].shape)==im['shape']
            channels=1 if a['values'].ndim==2 else a['values'].shape[-1]
            for c in range(channels):
                suffix='' if channels==1 else f'_ch{c}'
                df=pd.read_csv(base/im['shared_data']/f'values{suffix}.csv',float_precision='round_trip')
                values=a['values'] if channels==1 else a['values'][...,c]
                assert np.array_equal(df.iloc[:,1:].to_numpy(float),values,equal_nan=True)
            for k in ['GeoTIFF','GeoTIFF_display','GeoTIFF_native','GeoTIFF_native_display']:
                if im.get(k):assert (p.parent/im[k]).exists()
# Individual source checks guard against row/order loss and colorbar clipping.
q1=json.loads((base/'q1/Q1_11_all_candidate_ranking/figure_metadata.json').read_text())
bars=[]
for a in q1['axes']:
    f=base/'q1/Q1_11_all_candidate_ranking'/f"{a['axis']}_bars.csv"
    bars.extend(pd.read_csv(f,float_precision='round_trip').width)
ref=pd.read_csv(ROOT/'outputs/q1/Q1_all_candidate_ranking.csv',float_precision='round_trip')
assert len(bars)==56 and np.allclose(sorted(bars),sorted(ref.rarity_anomaly_score),rtol=0,atol=1e-15)
p=json.loads((base/'q4/Q4_10_impact_pressure_distribution/figure_metadata.json').read_text())
layer=p['axes'][0]['rasters'][1];pressure=np.load(base/layer['shared_data']/'array.npz')['values']
peak=json.loads((ROOT/'outputs/q4/Q4_2d_status.json').read_text())['mass']['max_dynamic_pressure_kPa']
assert np.isclose(np.nanmax(pressure),peak,rtol=1e-14)
assert np.nanmax(pressure)>layer['vmax'],'The original truncated color scale must not truncate the exported values'
for f in (base/'q4/Q4_07_parameter_sensitivity').glob('*scatter.csv'):
    df=pd.read_csv(f)
    if 'marker_role' in df:assert df.loc[df.marker_role.eq('symbolic_nonarrival'),'actual_arrival_time_min'].isna().all()
result={'figure_count':len(expected),'resolved_plot_datasets':resolved,'all_numeric_raster_cells_roundtrip_exact':True,
        'Q1_all_56_factor_values_present':True,'full_unclipped_dynamic_pressure_retained':True,
        'nonarrival_not_converted_to_time':True,'missing_references':0}
(base/'retention_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(result)
