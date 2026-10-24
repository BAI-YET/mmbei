"""Independent task-contract tests: target evidence vs derived proxies, not only no-crash tests."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
ROOT=Path(__file__).resolve().parents[1]

def test_q1_no_redundancy_in_matched_baseline():
    from mathorcup.common import Context
    from mathorcup.q1_event_signatures import _window_table,_get_ref_values
    from scipy.stats import spearmanr
    ctx=Context(ROOT)
    rank=ROOT/'outputs/q1/Q1_mechanism_representatives.csv'
    fn=ROOT/'data/processed/features_hourly.csv.gz'
    if not rank.exists() or not fn.exists():pytest.skip('Need official inputs and q1 output')
    table=pd.read_csv(rank)
    cols=table.feature.astype(str).tolist()
    assert len(cols)==5 and len(set(cols))==5
    # No duplicate raw channels or the same 2m/20m temperature variables.
    assert len([x for x in cols if x in ('temp_2m','temp_20m')])<=1
    for x in cols:
        assert x in pd.read_csv(fn,nrows=1).columns
    full=pd.read_csv(fn,parse_dates=['timestamp_raw']).set_index('timestamp_raw').sort_index()
    rolling=_window_table(full,24,ctx.cfg['min_valid_fraction'])
    historical=_get_ref_values(rolling,pd.Timestamp(ctx.cfg['snapshot_raw']),ctx.cfg)
    for i in range(len(cols)):
        for j in range(i+1,len(cols)):
            x=historical[cols[i]];y=historical[cols[j]]
            good=x.notna() & y.notna()
            if good.sum()>=25:
                rho=spearmanr(x[good],y[good]).statistic
                assert not np.isfinite(rho) or abs(rho) < ctx.cfg['ranking_redundancy_spearman']

def test_q2_zero_precipitation_never_triggers_rain():
    from mathorcup.q2_trigger_evaluation import dynamic_conditions
    import yaml
    fn=ROOT/'data/processed/features_hourly.csv.gz'
    if not fn.exists():pytest.skip('Need official processed MET')
    h=pd.read_csv(fn,parse_dates=['timestamp_raw']).set_index('timestamp_raw').sort_index()
    d,_=dynamic_conditions(h,source_z=5297,station_z=4140)
    good=(d.rain_24h<=.01) & (d.rain_24h.notna())
    assert int(good.sum())>100
    assert (d.loc[good,'rain_gate']<=1e-10).all()
    assert (~d.loc[good,'rain_condition_met']).all()

def test_q2_external_scope_and_weather_integrity(tmp_path):
    from mathorcup.common import Context
    from mathorcup.q2_external import load_q2_forecast,q2_forecast_url,q2_external_status
    from types import SimpleNamespace
    ctx=Context(ROOT);url=q2_forecast_url()
    assert 'run=2026-08-25T18%3A00' in url
    data,status=load_q2_forecast(ctx)
    assert status in ('NOT_DOWNLOADED','ACTUAL_DOWNLOADED_ARCHIVED_FORECAST',
                      'RUN_NOT_AVAILABLE_AT_ISSUANCE','FORECAST_VARIABLE_MISSING')
    if status=='NOT_DOWNLOADED':assert data is None
    # The source-only checkout must still record an honest missing-data status.
    sandbox=SimpleNamespace(root=tmp_path,folder=lambda name:tmp_path/'outputs'/name)
    report=q2_external_status(sandbox)
    assert report['weather_historical_run_status']=='NOT_DOWNLOADED'
    assert report['external_data_scope']=='QUESTION_2_ONLY'
    assert (tmp_path/'outputs/q2/Q2_external_data_provenance.json').exists()

def test_q2_thresholds_and_spatial_scale():
    fp=ROOT/'outputs/q2/Q2_primary_trigger_assessment.json'
    if not fp.exists():pytest.skip('q2 has not run')
    data=json.loads(fp.read_text())
    assert data['thresholds']['rain_24h_nonzero_q90_mm']>0
    assert data['thresholds']['source_pdd7_warmseason_q90_Cday']>0
    assert data['thresholds']['zero_rain_cannot_activate_rain_trigger'] is True
    df=pd.DataFrame(data['spatial_resolution_reliability'])
    assert {30,60,120}.issubset(set(df.grid_resolution_m))
    assert df.high20percent_IoU.between(0,1).all()

def test_q3_consumes_q2_physics_and_outputs_real_spatial_polygon():
    from mathorcup.q3_task_aligned import prepare_joint
    from mathorcup.common import Context
    ctx=Context(ROOT)
    if not (ROOT/'outputs/q3/Q3_primary_competition_answer.json').exists():pytest.skip('q3 not run')
    table=prepare_joint(ctx)
    assert set(('rain_gate','melt_gate','trigger_condition_index')).issubset(table.columns)
    report=json.loads((ROOT/'outputs/q3/Q3_primary_competition_answer.json').read_text())
    assert report['external_Q3_download'] is False
    assert report['empirical_disaster_time_PICP'] is None
    low,high=report['conditional_95pct_interval_h_from_issue']
    assert 0<=low<high<=ctx.cfg['q3']['primary_horizon_hours']
    gis=ROOT/'outputs/q3/Q3_flow_connected_warning_footprint.gpkg'
    assert gis.exists() and gis.stat().st_size>1000
    import geopandas as gpd
    g=gpd.read_file(gis)
    assert not g.empty and g.geometry.is_valid.all()
    assert g.crs.to_epsg()==32645

def test_q4_grid_and_mass_reliable():
    q4=ROOT/'outputs/q4'
    fn=q4/'Q4_grid_resolution_reliability.csv'
    if not fn.exists():pytest.skip('q4 has not run')
    df=pd.read_csv(fn)
    assert {60,120}.issubset(set(df.resolution_m))
    assert (df.mass_balance_rel_error < 1e-8).all()
    assert (df.wet_area_km2>0).all()
    assert (df.max_depth_m>0).all()

def test_q3_q4_joint_horizon_and_no_fake_disaster_coverage():
    fn=ROOT/'outputs/joint/Q34_status.json'
    if not fn.exists():pytest.skip('joint stage not run')
    s=json.loads(fn.read_text())
    assert s.get('empirical_joint_coverage') is None
    if 'conditional_on_event_within_configured_horizon' in s:
        assert s['conditional_on_event_within_configured_horizon']==6

def test_report_status_matches_v17_current_output():
    fn=ROOT/'outputs/SUMMARY_STATUS.json'
    if not fn.exists():pytest.skip('report not built')
    s=json.loads(fn.read_text())
    assert s['version']=='v2.0'
    assert s['extra_data_only_question_2'] is True
    assert s['Q3_disaster_event_PICP'] is None
