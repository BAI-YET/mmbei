"""Trace official inputs, archived forecasts, fitted models and every plotted figure."""
from pathlib import Path
import hashlib
import json
import yaml
from .common import write_json


def audit_project(ctx):
    records=[]
    for kind,base in [('official_input',ctx.raw/'unpacked'),('external_input',ctx.root/'data/external'),
                      ('processed_input',ctx.processed),
                      ('source',ctx.root/'src'),('script',ctx.root/'scripts'),('test',ctx.root/'tests'),
                      ('configuration',ctx.root/'configs')]:
        for p in sorted(base.rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts:continue
            with p.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
            records.append({'kind':kind,'path':str(p.relative_to(ctx.root)),
               'bytes':p.stat().st_size,'sha256':digest})
    for name,kind in [('run.py','source'),('requirements.txt','dependencies'),('README.md','documentation'),('.gitignore','configuration')]:
        p=ctx.root/name
        if p.exists():records.append({'kind':kind,'path':name,'bytes':p.stat().st_size,
            'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    forecast=ctx.root/'data/external/q2_forecast/provenance.json'
    raw=ctx.root/'data/external/q2_forecast/ecmwf_source_pre_event_single_run.json'
    forecast_verified=False
    if forecast.exists() and raw.exists():
        p=json.loads(forecast.read_text());forecast_verified=hashlib.sha256(raw.read_bytes()).hexdigest()==p['raw_sha256']
        if not forecast_verified:raise ValueError('Archived forecast differs from its recorded raw SHA256')
    evidence=yaml.safe_load((ctx.root/'configs/case_evidence.yaml').read_text());v=evidence['published_values']
    assert v['release_volume_m3']['value']==v['retained_deposit_volume_m3']['value']+v['mobile_debris_volume_m3']['value']
    figures=[p for name in ['q1','q2','q3','q4','joint'] for p in (ctx.out/name).glob('*.png')]
    missing=[]
    for p in figures:
        metadata=ctx.out/'redraw_data'/p.parent.name/p.stem/'figure_metadata.json'
        if not metadata.exists():missing.append(str(p.relative_to(ctx.out)))
    if missing:raise ValueError('Figures missing recorded plotting inputs: '+', '.join(missing))
    q3=json.loads((ctx.out/'q3/Q3_primary_competition_answer.json').read_text())
    cascade=json.loads((ctx.out/'q4/Q4_cascade_reconstruction.json').read_text())
    spatial=json.loads((ctx.out/'q4/Q4_reconstruction_2d_summary.json').read_text())
    assert q3['hazard_base_rate']['assumed_disaster_count'] is None
    assert cascade['cascade_volume_residual']<1e-8
    for r in spatial['grid_experiments']:
        assert abs(r['input_m3']-r['expected_input_from_1d_m3'])<=max(1.,r['input_m3']*1e-10)
    summary={'version':'v2.2','files':records,'figure_count':len(figures),
       'all_figures_have_numeric_plot_inputs':True,'archived_forecast_raw_sha256_verified':forecast_verified,
       'observed_data_roles':['official measurements','actual archived forecast'],
       'analysis_roles':['derived indicators','observed-day/year resampling','published event reconstruction constraints'],
       'simulation_parameter_roles':['fitted effective parameters','declared sensitivity experiments','numerical settings'],
       'production_synthetic_inlet_used':False,'invented_disaster_count_used':False,
       'case_constraint_sources':evidence['references'],'cascade_volume_residual':cascade['cascade_volume_residual']}
    write_json(ctx.folder('quality')/'reproducibility_manifest.json',summary)
    return summary
