"""Analytic and causal checks for the shared Q2/Q3 initiation calculation."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from mathorcup.q2_trigger_evaluation import (
    advance_water_storage,critical_storage,source_safety_factor,
    mechanical_source_geometry,mechanical_candidate_region,source_history_states)
from mathorcup.q3_task_aligned import (
    forecast_source_stability,first_instability_time,mechanical_weather_members,
    empirical_first_passage)


def parameters():
    return pd.DataFrame({'scenario_id':[0],'cohesion_kPa':[10.],
        'friction_angle_deg':[45.],'storage_yield':[.01],'drainage_hours':[10.],
        'degree_day_factor_mm_Cday':[6.],'infiltration_fraction':[1.], 'lapse_C_km':[6.5]})


def geometry():
    return {'normal_stress_kPa':100.,'driving_stress_kPa':60.,'source_elevation_m':4000.}


def test_water_reservoir_exact_solution_and_complete_mass_balance_with_spill():
    s,drain,spill=advance_water_storage(np.array([2.,2.]),np.array([3.,3.]),4,np.array([10.,10.]),np.array([1000.,5.]))
    exact=2*np.exp(-.4)+30*(1-np.exp(-.4))
    assert s[0]==pytest.approx(exact)
    assert s[1]==5 and spill[1]>0
    assert np.allclose(s+drain+spill,14.)


def test_vertical_thickness_force_balance_includes_basal_area():
    g={'glacier_mask':np.ones((4,4),bool),'thickness':np.full((4,4),20.),
        'slope':np.full((4,4),30.),'dem':np.full((4,4),5000.),'transform':(0,30,0,0,0,-30)}
    result=mechanical_source_geometry(g,g['glacier_mask'],900)
    theta=np.pi/6
    assert result['normal_stress_kPa']==pytest.approx(900*9.81*20*np.cos(theta)**2/1000)
    assert result['driving_stress_kPa']==pytest.approx(900*9.81*20*np.cos(theta)*np.sin(theta)/1000)
    assert result['effective_slope_deg']==pytest.approx(30)


def test_critical_pressure_solves_the_same_q3_safety_factor():
    p=parameters();g=geometry();pressure,storage=critical_storage(p,g)
    assert pressure[0]==pytest.approx(50.)
    assert source_safety_factor(storage,p,g)[0]==pytest.approx(1.)
    assert source_safety_factor(storage*.9,p,g)[0]>1
    assert source_safety_factor(storage*1.1,p,g)[0]<1


def test_cohesion_above_driving_stress_has_no_water_only_failure():
    p=parameters();p.cohesion_kPa=61.
    _,storage=critical_storage(p,geometry())
    assert np.isinf(storage[0])
    assert source_safety_factor(np.array([1e8]),p,geometry())[0]>1


def test_exactly_limiting_cohesion_fails_only_at_full_normal_pressure():
    p=parameters();p.cohesion_kPa=60.
    pressure,storage=critical_storage(p,geometry())
    assert pressure[0]==pytest.approx(100.)
    assert source_safety_factor(storage,p,geometry())[0]==pytest.approx(1.)


def test_forecast_first_loss_of_stability_matches_analytic_reservoir_crossing():
    p=parameters();g=geometry();_,limit=critical_storage(p,g)
    leads=np.arange(.25,12.01,.25)
    # T=0 stops melt while positive temperature is needed for liquid rain.
    p.degree_day_factor_mm_Cday=0
    fs,_=forecast_source_stability([0.],p,g,np.ones((1,len(leads))),np.full((1,len(leads)),10.),leads,4000.)
    exact=-10*np.log(1-limit[0]/100)
    numerical=first_instability_time(fs,leads)[0]
    assert exact<=numerical<exact+.25


def test_nonfailure_mass_is_retained_without_inventing_a_time_interval():
    result=empirical_first_passage(-np.full((20,12),2.),.5,-1)
    assert result['p_in_horizon']==0
    assert result['interval_h']==(None,None)
    assert result['nonoccurrence_probability']==1


def weather_history():
    idx=pd.date_range('2024-06-01',periods=24*50,freq='h');x=np.arange(len(idx))
    return pd.DataFrame({'temp_2m':5+2*np.sin(x/4),'precip_mm':np.maximum(0,np.sin(x/5))},index=idx)


def test_weather_members_are_invariant_to_any_future_observation_change():
    h=weather_history();issue=pd.Timestamp('2024-07-10 09:00');leads=np.arange(.5,6.1,.5)
    before=mechanical_weather_members(h,issue,leads,4000.)
    h.loc[h.index>=issue,['temp_2m','precip_mm']]=99999
    after=mechanical_weather_members(h,issue,leads,4000.)
    assert np.array_equal(before[0],after[0]) and np.array_equal(before[1],after[1])


def test_missing_spinup_forcing_is_unknown_rather_than_zero_rain():
    h=weather_history();h.loc['2024-07-01 00:00','precip_mm']=np.nan
    with pytest.raises(ValueError,match='Incomplete'):
        source_history_states(h,'2024-07-10 09:00',parameters(),geometry(),4000.)


def test_geometric_screen_does_not_choose_a_large_flat_inventory():
    m=np.ones((8,8),bool);slope=np.full((8,8),10.);slope[4:8,4:8]=40.
    g={'glacier_mask':m,'thickness':np.full((8,8),30.),'slope':slope}
    selected=mechanical_candidate_region(g,30,10)
    assert selected.sum()==16 and np.all(slope[selected]>=30)


def test_primary_release_distinguishes_model_failure_from_real_disaster_coverage():
    root=Path(__file__).resolve().parents[1]
    path=root/'outputs/q3/Q3_primary_competition_answer.json'
    if not path.exists():pytest.skip('Q3 not run')
    report=json.loads(path.read_text())
    assert report['model_family']=='Q2_conditioned_event_time'
    assert report['Q4_results_used_as_Q3_inputs'] is False
    assert report['empirical_disaster_time_PICP'] is None
    assert report['event_time_used_to_fit_parameters'] is False
    assert report['assessment_cross_probability_is_disaster_probability'] is False
    assert report['event_occurrence_probability'] is None
    assert report['case_reference_coverage']['n_references']==1
    q2=json.loads((root/'outputs/q2/Q2_sequential_model_contract.json').read_text())
    assert report['inherited_Q1_features']==q2['Q1_features']
    assert report['source_geometry_sha256']==q2['source_geometry_sha256']


def test_noninitiation_joint_output_retires_previous_time_samples(tmp_path):
    from types import SimpleNamespace
    from mathorcup.q34_joint_arrival import execute_joint_arrival
    root=tmp_path/'outputs'
    for q in ['q3','q4','joint']:(root/q).mkdir(parents=True)
    pd.DataFrame({'conditional_time_probability':[0.,0.]}).to_csv(root/'q3/Q3_source_initiation_time_distribution.csv',index=False)
    (root/'q4/Q4_reconstruction_geometry.json').write_text('{}')
    (root/'q3/Q3_primary_competition_answer.json').write_text('{}')
    old=root/'joint/Q34_computed_joint_samples.csv.gz';old.write_bytes(b'old forecast')
    ctx=SimpleNamespace(out=root,folder=lambda q:root/q)
    execute_joint_arrival(ctx)
    assert not old.exists()
    assert json.loads((root/'joint/Q34_status.json').read_text())['nominal_95_interval_hours_from_raw_issuance'] is None


def test_conditional_time_model_has_analytic_uniform_quantiles():
    from mathorcup.q3_task_aligned import conditional_trigger_time
    r=conditional_trigger_time(np.ones((7,12)),.5)
    assert r['conditional_bin_mass'].sum()==pytest.approx(1)
    assert r['interval_h']==pytest.approx((.15,5.85))
    assert r['median_h']==pytest.approx(3.) and r['event_occurrence_probability'] is None


def test_conditional_time_density_follows_positive_intensity_support():
    from mathorcup.q3_task_aligned import conditional_trigger_time
    r=conditional_trigger_time(np.array([[1]*6+[0]*6]),.5)
    assert r['interval_h']==pytest.approx((.075,2.925))
    assert r['conditional_bin_mass'][6:].sum()==0


def test_zero_intensity_uses_explicit_conditional_prior_without_inventing_event_probability():
    from mathorcup.q3_task_aligned import conditional_trigger_time
    r=conditional_trigger_time(np.zeros((2,12)),.5)
    assert r['n_uniform_prior_trajectories']==2
    assert r['event_occurrence_probability'] is None
    assert r['interval_h']==pytest.approx((.15,5.85))


def test_case_coverage_keeps_single_reference_denominator_and_uncertainty():
    from mathorcup.q3_task_aligned import single_case_coverage
    r=single_case_coverage((0,6),1.5)
    assert r['n_references']==1 and r['n_covered']==1
    assert r['binomial_95_interval']==pytest.approx([.025,1])
    assert single_case_coverage((2,6),1.5)['binomial_95_interval']==pytest.approx([0,.975])


def test_conditional_time_sampler_uses_only_the_declared_positive_bin():
    from mathorcup.q34_joint_arrival import sample_discrete_event_time
    t=sample_discrete_event_time([.5,1.,1.5],[1.,0.,0.],500,seed=2,uniform_within_bin=True)
    assert np.all((t>=0)&(t<=.5)) and len(np.unique(t))>100


def test_saved_time_distribution_recomputes_from_saved_Q2_intensity_ensemble():
    from mathorcup.q3_task_aligned import conditional_trigger_time
    root=Path(__file__).resolve().parents[1]/'outputs/q3'
    if not (root/'Q3_assessment_ensemble.npz').exists():pytest.skip('Q3 ensemble not generated')
    with np.load(root/'Q3_assessment_ensemble.npz') as a:
        r=conditional_trigger_time(a['assessment_score'],float(a['lead_h'][0]))
    saved=pd.read_csv(root/'Q3_source_initiation_time_distribution.csv')
    report=json.loads((root/'Q3_primary_competition_answer.json').read_text())
    assert np.allclose(saved.conditional_time_probability,r['conditional_bin_mass'],atol=1e-14)
    assert report['conditional_95pct_interval_h_from_issue']==pytest.approx(r['interval_h'])


def test_Q4_uses_the_exact_source_pixels_prepared_in_Q2():
    root=Path(__file__).resolve().parents[1]/'outputs'
    if not (root/'q2/Q2_reconstruction_sections.npz').exists():pytest.skip('Q2 geometry not generated')
    with np.load(root/'q2/Q2_reconstruction_sections.npz') as a,np.load(root/'q4/Q4_reconstruction_sections.npz') as b:
        assert np.array_equal(a['source_pixels_utm'],b['source_pixels_utm'])
    q2=json.loads((root/'q2/Q2_reconstruction_geometry.json').read_text())
    q4=json.loads((root/'q4/Q4_reconstruction_geometry.json').read_text())
    assert q2['geometry_contract_sha256']==q4['geometry_contract_sha256']
