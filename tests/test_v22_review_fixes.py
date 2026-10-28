"""Regression checks for review findings; analytic fixtures are testing inputs only."""
from pathlib import Path
from types import SimpleNamespace
import importlib.util
import numpy as np
import pandas as pd
import pytest
from mathorcup.q4_high_order import advance_channel_muscl
from mathorcup.q4_cascade_reconstruction import section_depth,resample_channel_profile
from mathorcup.q4_flowpath import require_verified_route
from mathorcup.q3_task_aligned import empirical_first_passage
from mathorcup.q2_trigger_evaluation import dynamic_conditions


def test_strict_route_rejects_direct_cascade_before_geometry_or_output_changes(tmp_path,monkeypatch):
    import mathorcup.q4_cascade_reconstruction as cascade
    saved=tmp_path/'outputs/q4/existing.png';saved.parent.mkdir(parents=True);saved.write_bytes(b'old result')
    ctx=SimpleNamespace(root=tmp_path,cfg={'strict_route':True,'q4':{'route_verified':False}})
    monkeypatch.setattr(cascade,'prepare_reconstruction_geometry',lambda *_:pytest.fail('unguarded cascade'))
    with pytest.raises(ValueError,match='Strict route'):cascade.execute_cascade_reconstruction(ctx)
    assert saved.read_bytes()==b'old result'


def test_strict_route_requires_real_waypoints_and_checks_conflicting_cli_flags(tmp_path):
    entry=Path(__file__).resolve().parents[1]/'run.py'
    spec=importlib.util.spec_from_file_location('review_cli',entry);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    ctx=SimpleNamespace(root=tmp_path,cfg={'q4':{'route_verified':True}})
    with pytest.raises(ValueError,match='cannot be combined'):module.validate_route_policy(ctx,True,True)
    path=tmp_path/'events/verified_route.csv';path.parent.mkdir();path.write_text('lon,lat\n85,28\n')
    with pytest.raises(ValueError,match='two'):require_verified_route(ctx)
    path.write_text('lon,lat\n85,28\n85.1,28.1\n')
    assert len(require_verified_route(ctx))==2


def test_muscl_variable_section_lake_remains_at_rest():
    x=np.arange(51);bed=.25*np.sin(x/6);b=25+5*np.cos(x/8);side=1+.2*np.sin(x/7)
    h=4-bed;A=b*h+side*h*h;Q=np.zeros_like(A);initial=A.copy()
    for _ in range(80):A,Q,out,clip=advance_channel_muscl(A,Q,bed,b,10,.01,0,1e20,bank_side_slope=side)
    assert np.max(abs(Q))<1e-10
    assert np.allclose(A,initial,atol=1e-11) and abs(out)<1e-10 and clip==0


def test_muscl_keeps_real_free_surface_pressure_gradient():
    # A reconstruction bug that corrects both faces to centre pressure cancels
    # this physical acceleration even while passing a lake-at-rest check.
    n=51;b=np.full(n,20.);h=np.linspace(3.,2.,n);A=b*h;Q=np.zeros(n)
    area,discharge,out,_=advance_channel_muscl(A,Q,np.zeros(n),b,10,.02,0,1e20)
    assert np.all(discharge[5:-5]>0)
    assert area.sum()*10+out==pytest.approx(A.sum()*10,rel=1e-12)


def test_muscl_advancing_wet_front_is_positive_and_conservative():
    n=100;b=np.full(n,20.);bed=np.zeros(n);h=np.where(np.arange(n)<30,2.,0.)
    A=b*h;Q=np.zeros(n);initial=A.sum()*10.;exited=0.
    for _ in range(160):
        A,Q,out,clip=advance_channel_muscl(A,Q,bed,b,10,.05,0,1e20);exited+=out
        assert np.min(A)>=0 and clip<1e-8
    assert (A>1e-3).sum()>30
    assert A.sum()*10.+exited==pytest.approx(initial,rel=1e-11)


def test_mesh_comparisons_preserve_both_endpoints_and_uniform_spacing():
    p=pd.DataFrame({'distance_m':[100.,111.35,122.7],'elevation_m':[8.,7.,6.],
                    'effective_width_m':[20.,22.,25.]})
    for dx in [3.75,7.5,15.]:
        a=resample_channel_profile(p,dx)
        assert a.distance_m.iloc[0]==100. and a.distance_m.iloc[-1]==122.7
        assert np.ptp(np.diff(a.distance_m))<1e-12


def test_empty_passage_ensemble_yields_no_warning_probability_or_fake_interval():
    result=empirical_first_passage(np.zeros((10,12)),.5,.8)
    assert result['p_in_horizon']==0 and result['nonoccurrence_probability']==1
    assert result['interval_h']==(None,None) and result['median_h'] is None


def test_passage_median_uses_actual_discrete_time_support():
    result=empirical_first_passage(np.array([[0.,0.,1.],[0.,0.,1.]]),.5,.8)
    assert result['median_h']==1.5


def _energy_fixture():
    a=pd.date_range('2023-05-01',periods=24*100,freq='h')
    b=pd.date_range('2024-05-01',periods=24*20,freq='h')
    index=a.append(b);n=len(index);phase=np.arange(n)
    return pd.DataFrame({'temp_2m':5.+np.sin(phase/20),
      'precip_mm':.5+.2*np.sin(phase/7),'rain_24h':12.+np.sin(phase/7),
      'pdd_7d':35.,'energy_balance_residual':80.+5*np.sin(phase/31)},index=index)


def test_energy_channel_enters_primary_index_without_changing_earlier_fit():
    h=_energy_fixture();a,limits=dynamic_conditions(h,[2023],source_z=4140,station_z=4140)
    h.loc[h.index.year==2024,'energy_balance_residual']=400.
    b,changed=dynamic_conditions(h,[2023],source_z=4140,station_z=4140)
    assert limits==changed
    assert b.energy_gate.iloc[-1]>.99
    assert b.trigger_condition_index.iloc[-1]==b.energy_gate.iloc[-1]
    assert a.trigger_condition_index.iloc[-1]<b.trigger_condition_index.iloc[-1]
    assert np.allclose(a.loc[a.index.year==2023,'trigger_condition_index'],
                       b.loc[b.index.year==2023,'trigger_condition_index'],equal_nan=True)


def test_missing_energy_remains_unknown_in_a_planned_three_channel_index():
    h=_energy_fixture();h.loc[h.index[-200:],'energy_balance_residual']=np.nan
    a,_=dynamic_conditions(h,[2023],source_z=4140,station_z=4140)
    assert pd.isna(a.energy_gate.iloc[-1]) and pd.isna(a.trigger_condition_index.iloc[-1])


def test_warning_decision_checks_expiry_at_issue_without_future_event_truth():
    from mathorcup.warning_validation import decide_warning_at_issue
    d=decide_warning_at_issue('2026-08-26 09:00',.1,.25,True,
        '2026-08-25 18:00','2026-08-26 08:30')
    assert d['warning_time'] is None and d['valid_until'] is None
    d=decide_warning_at_issue('2026-08-26 09:00',.5,.25,True,
        '2026-08-25 18:00','2026-08-26 10:00')
    assert d['warning_time']==pd.Timestamp('2026-08-25 18:00')
    assert d['forecast_warning_time']==pd.Timestamp('2026-08-26 09:00')
    assert d['valid_until']==pd.Timestamp('2026-08-26 15:00')


def test_actual_source_projection_keeps_high_material_at_high_route_positions():
    from mathorcup.q4_terrain_reconstruction import project_release_pixels
    p=pd.DataFrame({'distance_m':[0.,100.,200.,1000.],
        'utm_x':[0.,100.,200.,1000.],'utm_y':[0.,0.,0.,0.],
        'elevation_m':[5200.,5180.,5160.,4000.]})
    indices,table=project_release_pixels(np.array([[0.,10.],[100.,10.],[200.,10.]]),
        np.array([5200.,5180.,5160.]),p)
    assert indices.tolist()==[0,1,2]
    assert np.allclose(table.source_elevation_m,table.route_elevation_m)


def test_delivered_release_retains_source_gravitational_elevation():
    import json
    root=Path(__file__).resolve().parents[1]
    path=root/'outputs/q4/Q4_reconstruction_geometry.json'
    if not path.exists():pytest.skip('Q4 has not run')
    g=json.loads(path.read_text())
    if 'selected_source_mean_elevation_m' not in g:pytest.skip('new geometry has not run')
    # Three steep-terrain DEM cells allow local projection error; a hundreds-
    # of-metres systematic downward shift is a different release experiment.
    assert abs(g['assigned_release_mean_bed_elevation_m']-g['selected_source_mean_elevation_m'])<100


def test_saved_spatial_route_matches_the_two_accepted_numerical_meshes():
    root=Path(__file__).resolve().parents[1]/'outputs/q4'
    if not (root/'Q4_reconstruction_lower_arrays.npz').exists():pytest.skip('Q4 has not run')
    route=pd.read_csv(root/'Q4_reconstruction_route.csv')
    up=np.load(root/'Q4_reconstruction_upper_arrays.npz');low=np.load(root/'Q4_reconstruction_lower_arrays.npz')
    assert len(route)==len(up['max_depth_m'])+len(low['max_depth_m'])-1
    assert (np.diff(route.distance_m)>0).all()
