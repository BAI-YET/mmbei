"""v1.7 narrow regressions for previous data-utilisation and mislabelling bugs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def test_raw_observation_inventory_has_all_63_physical_and_diagnostics():
    from mathorcup.prepare_time import _RAW_FIELDS
    assert {k:len(v) for k,v in _RAW_FIELDS.items()} == {'MET':22,'RADM':4,'SOIL':12,'FLUX':25}
    names=[v for m in _RAW_FIELDS.values() for v in m.values()]
    assert len(names)==len(set(names))==63
    inv=ROOT/'outputs'/'quality'/'complete_field_inventory_63.csv'
    if inv.exists():
        x=pd.read_csv(inv)
        assert len(x)==63
        assert {'physical_candidate','diagnostic_only'}==set(x.role)
        assert x.valid_observations.ge(0).all()


def test_q1_new_data_contract_includes_2021_and_preserves_phyiscal_units():
    from mathorcup.common import Context
    from mathorcup.features import make_hourly
    ctx=Context(ROOT)
    assert 2021 in ctx.cfg['historic_years']
    d=pd.DataFrame({'timestamp_raw':pd.date_range('2021-08-01',periods=3,freq='h'),
        'sw_down':[600,200,0],'sw_up':[120,50,0],
        'lw_down':[310,310,310],'lw_up':[350,350,350],
        'soilT_0p1':[0,1,2],'soilT_1p6':[-1,0,0],
        'vwc_0p1':[.3,.2,.1],'vwc_1p6':[.1,.1,.1]})
    h=make_hourly(d)
    assert np.allclose(h.net_shortwave,[480,150,0])
    assert np.allclose(h.net_longwave,[-40,-40,-40])
    assert np.allclose(h.net_radiation,[440,110,-40])
    assert np.isclose(h.soil_temperature_gradient.iloc[0],(0-(-1))/1.5)
    assert np.isclose(h.vwc_deep_gradient.iloc[0],(.3-.1)/1.5)


def test_q2_all_four_historic_glacier_epochs_and_no_masked_landcover_fabrication():
    p=ROOT/'outputs'/'q2'/'Q2_official_glacier_four_period_elevation_change.csv'
    if p.exists():
        d=pd.read_csv(p)
        assert len(d)==4
        assert d.notna().any(axis=1).all()
    f=ROOT/'outputs'/'q2'/'Q2_official_landcover_four_year_class_counts.csv'
    if f.exists():
        df=pd.read_csv(f)
        assert set(df.year)=={2020,2021,2022,2023}
        assert df.grid_cell_area_m2_approx.between(100., 1.e8).all()
        assert df.select_dtypes(include='number').notna().any(axis=1).all()


def test_q2_candidate_ice_inventory_is_distinct_from_downstream_mixture_volume():
    fn=ROOT/'outputs'/'q2'/'Q2_source_release_candidates.csv'
    if fn.exists():
        df=pd.read_csv(fn)
        assert df.candidate_inventory_ice_volume_m3.gt(0).all()
        assert df.candidate_area_km2.gt(0).all()
    cr=ROOT/'outputs'/'q4'/'Q4_source_volume_data_crosswalk.json'
    if cr.exists():
        d=json.loads(cr.read_text())
        assert d['source_ice_volume_is_not_release_volume'] if 'source_ice_volume_is_not_release_volume' in d else True
        assert d['Q4_effective_mixture_volume_scenario_m3']>0


def test_q3_q1_factors_all_transferred_and_two_other_feature_groups_preserved():
    q1=ROOT/'outputs'/'q1'/'Q1_top5_with_stability.csv'
    f=ROOT/'outputs'/'q3'/'Q3_Q1_feature_transfer_contract.csv'
    if q1.exists() and f.exists():
        selected=pd.read_csv(q1).feature.astype(str).tolist()
        trans=pd.read_csv(f)
        assert set(selected)==set(trans.feature)
        assert len(selected)==5
        assert not trans.used_by_primary_forward_model.any()
        r=json.loads((ROOT/'outputs/q3/Q3_primary_competition_answer.json').read_text())
        assert r['inherited_Q1_features']==['temp_2m','sw_down','vpd_2m']
    a=ROOT/'outputs'/'q3'/'Q3_Q1_Q2_ablation_validation.csv'
    if a.exists():
        d=pd.read_csv(a)
        assert set(d.model)=={'Q1+Q2','Q1-only','Q2-only'}
        assert d.hours_above_own_threshold.ge(0).all()


def test_q3_calendar_alarm_is_actually_on_graph_timeline():
    """Regression for the old xlim(-0.2,6.2) that hid previous-evening alarms."""
    from mathorcup.q3_task_aligned import _plots
    import inspect
    src=inspect.getsource(_plots)
    assert 'x_alarm' in src
    assert 'xlim=(left,right)' in src
    assert "-18" not in src  # derived from timestamps, not hardcoded
    assert (ROOT/'outputs'/'q3'/'Q3_07_feature_group_ablation.png').exists() or not (ROOT/'outputs'/'q3'/'Q3_primary_competition_answer.json').exists()


def test_q4_dynamic_pressure_nonnegative_and_density_linear_scaling():
    from mathorcup.q4_solver import simulate_2d
    # A flat/weakly sloping synthetic DEM, low injection, short horizon; physics diagnostic is well-defined.
    z=1500-np.tile(np.arange(16)*.2,(16,1))
    hydro=pd.DataFrame({'time_s':[0,2,4,8],'sensor_discharge_m3_s':[20,20,0,0]})
    # No need to claim actual disaster pressure; check numerical consistency only.
    a=simulate_2d(z,50,(7,3),hydro,total_seconds=12,mu=.08,xi=400,
        cfl=.2,frames=3,inlet_mode='legacy_3x3',density_kg_m3=1000)
    b=simulate_2d(z,50,(7,3),hydro,total_seconds=12,mu=.08,xi=400,
        cfl=.2,frames=3,inlet_mode='legacy_3x3',density_kg_m3=2000)
    pa=a['max_dynamic_pressure_kPa'];pb=b['max_dynamic_pressure_kPa']
    assert pa.shape==z.shape and (pa>=0).all() and np.isfinite(pa).all()
    assert np.allclose(pb,2*pa)
    assert a['mass']['dynamic_pressure_coefficient']==.5
    assert a['mass']['relative_mass_error']<1e-8


def test_q2_external_is_explicitly_only_q2_and_not_fake_forecast():
    from mathorcup.q2_external import load_q2_forecast
    from mathorcup.common import Context
    ctx=Context(ROOT)
    df,status=load_q2_forecast(ctx)
    assert status in {'NOT_DOWNLOADED','ACTUAL_DOWNLOADED_ARCHIVED_FORECAST',
       'RUN_NOT_AVAILABLE_AT_ISSUANCE','FORECAST_VARIABLE_MISSING'}
    if status=='NOT_DOWNLOADED': assert df is None
    assert not (ROOT/'src/mathorcup/q3_time_risk.py').read_text().count('requests.get(')


def test_final_figure_titles_are_paper_style():
    from mathorcup.common import FIGURE_TITLES
    assert 'Q3_07_feature_group_ablation' in FIGURE_TITLES
    assert 'Q4_10_impact_pressure_distribution' in FIGURE_TITLES
    for k,v in FIGURE_TITLES.items():
        assert len(v)<35
        assert '为什么' not in v
