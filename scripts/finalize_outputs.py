"""Refresh derived publication assets, scenario IDs, joint results and data/source audit."""
from pathlib import Path
import sys,json,hashlib,subprocess
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context,configure_plotting,write_json
from mathorcup.publication_assets import make_publication_assets
from mathorcup.q34_joint_arrival import execute_joint_arrival
from mathorcup.report import make_report
from mathorcup.catalog import prepare_official_data
ctx=Context(ROOT);configure_plotting()
prepare_official_data(ctx)
sweep=ctx.out/'q4/Q4_parameter_sensitivity.csv';d=pd.read_csv(sweep);cfg=ctx.cfg['q4']
base={'mu':cfg['mu'],'xi':cfg['xi_m_s2'],'width':cfg['width_m'],'volume':cfg['release_volume_m3'],
      'duration':cfg['release_duration_s'],'total_seconds':cfg['total_seconds'],'cfl':cfg['cfl'],'frames':3,
      'initial_velocity_m_s':cfg['initial_velocity_m_s']}
ids=[]
for r in d.itertuples():
    kw=base.copy();kw[r.varied_parameter]=r.value
    key=json.dumps({k:round(float(v),7) for k,v in kw.items()},sort_keys=True)
    ids.append(hashlib.sha256(key.encode()).hexdigest()[:12])
d['scenario_id']=ids;d.to_csv(sweep,index=False)
from mathorcup.q3_spatial_footprint import execute_q3_spatial
execute_q3_spatial(ctx);make_publication_assets(ctx);execute_joint_arrival(ctx);make_report(ctx)
# Every file used for computation is identified; avoid recording absolute runtime paths.
records=[]
for group,base in [('source',ROOT/'src'),('scripts',ROOT/'scripts'),('tests',ROOT/'tests'),('config',ROOT/'configs'),('official_data',ROOT/'data/raw/unpacked'),('external_data',ROOT/'data/external')]:
    for p in sorted(base.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        sha=hashlib.sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):sha.update(chunk)
        records.append({'kind':group,'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha.hexdigest()})
write_json(ctx.folder('quality')/'reproducibility_manifest.json',{'version':'v2.0','files':records,
    'forecast_run_utc':'2026-08-25T18:00','validation_split':{'fit':[2021,2022,2023],'calibration':[2024],'test':[2025]},
    'source_trigger_truth_independently_observed':False})
print('Current figures, joint results, gallery and source/data hashes refreshed.')
