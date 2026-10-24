"""Behavioral regressions for task-level errors found in the supplied release."""
from pathlib import Path
import numpy as np
import pandas as pd
from mathorcup.q4_solver import advance_2d,advance_1d
from mathorcup.q2_trigger_evaluation import dynamic_conditions,velocity_strain_proxy
from mathorcup.warning_validation import warning_state_machine,conformal_expansion,backtest_metrics


def test_geographic_north_inlet_moves_north_on_a_north_up_raster():
    z=np.zeros((21,21));U=np.zeros((3,21,21))
    U,_,_=advance_2d(U,z,30,.1,0,1e9,10,10,900.,velocity_xy=(0,20),
                     inlet_cells=np.array([[10,10]]),inlet_weights=np.array([1.]))
    assert U[2].sum()<0
    for _ in range(15):U,_,_=advance_2d(U,z,30,.1,0,1e9,10,10,0)
    yy=np.indices(z.shape)[0]
    assert (U[0]*yy).sum()/U[0].sum()<10


def test_still_lake_on_variable_bed_has_no_artificial_2d_flow():
    y,x=np.mgrid[:17,:19];z=.15*np.sin(x/3)+.2*np.cos(y/4)
    U=np.stack([3-z,np.zeros_like(z),np.zeros_like(z)]);expected=U.copy()
    for _ in range(40):U,flux,cut=advance_2d(U,z,10,.05,0,1e12,8,8,0)
    assert np.allclose(U,expected,atol=1e-12)
    assert abs(flux)<1e-12 and cut==0


def test_still_lake_on_variable_bed_has_no_artificial_1d_flow():
    z=np.sin(np.arange(35)/6)*.4;h=3-z;q=np.zeros_like(z);initial=h.copy()
    for _ in range(50):h,q,_,_=advance_1d(h,q,np.zeros_like(z),np.ones_like(z),10,.05,0,1e12,0,10,bed=z)
    assert np.allclose(h,initial,atol=1e-12)
    assert np.max(abs(q))<1e-11


def _weather():
    idx=pd.date_range('2022-05-01',periods=3600,freq='h');rng=np.random.default_rng(9)
    rain=rng.choice([0.,0.,1.,3.],len(idx))
    h=pd.DataFrame({'temp_2m':12+rng.normal(0,2,len(idx)),'precip_mm':rain,
       'rain_24h':pd.Series(rain,index=idx).rolling(24,min_periods=18).sum(),'pdd_7d':10.},index=idx)
    return h


def test_missing_weather_remains_unknown_not_zero_risk():
    h=_weather();h.iloc[-190:]=np.nan
    d,_=dynamic_conditions(h,baseline_years=[2022])
    assert pd.isna(d.trigger_condition_index.iloc[-1])
    assert pd.isna(d.melt_condition_met.iloc[-1]) and pd.isna(d.rain_condition_met.iloc[-1])


def test_velocity_derivatives_respect_decreasing_northing():
    x=np.arange(41)*30.;y=5000-np.arange(41)*30.;xx,yy=np.meshgrid(x,y)
    a,b,c,d=.01,.02,-.01,.03
    g={'x':x,'y':y,'vx':a*xx+b*yy,'vy':c*xx+d*yy,'glacier_mask':np.ones_like(xx,dtype=bool)}
    s=velocity_strain_proxy(g,sigma_cells=1.)
    expected=np.sqrt(a*a+d*d+2*((b+c)/2)**2)
    assert np.allclose(s[8:-8,8:-8],expected,rtol=1e-5)


def test_warning_retracts_and_does_not_reuse_previous_day_lead():
    idx=pd.date_range('2026-08-25 17:00',periods=17,freq='h')
    series=pd.Series([.9,.9,.2,.2,.2]+[.2]*12,index=idx)
    states,episodes=warning_state_machine(series,.7)
    assert states.active.iloc[1] and not states.active.iloc[4]
    assert states.index[1]==idx[1]+pd.Timedelta(hours=1)
    assert episodes.end_reason.iloc[0]=='three_low_hours'
    assert not states.active.iloc[-1]


def test_missing_observation_suspends_warning_and_requires_reentry():
    s=pd.Series([.9,.9,np.nan,.9,.9],index=pd.date_range('2026-08-25',periods=5,freq='h'))
    d,e=warning_state_machine(s,.7)
    assert d.state.iloc[2]=='UNKNOWN' and not d.active.iloc[3]
    assert d.active.iloc[4] and len(e)==2


def test_false_alarm_rate_and_false_alarm_ratio_have_distinct_denominators():
    df=pd.DataFrame({'forecast_cross_probability':[.8,.8,.1,.1,.1], 'true_crossed':[True,False,True,False,False],
      'issue_bjt':pd.date_range('2025-06-01',periods=5).astype(str),'already_above_at_issue':[False]*5,
      'first_true_cross_h':[1,np.nan,2,np.nan,np.nan],'lower_h':[0]*5,'upper_h':[6]*5})
    m=backtest_metrics(df)
    assert m['proxy_false_alarm']==1/3 and m['proxy_false_alarm_ratio']==.5


def test_too_few_calibration_events_produces_wide_interval():
    d=pd.DataFrame({'true_crossed':[True]*5,'first_true_cross_h':[1]*5,'lower_h':[0]*5,'upper_h':[6]*5})
    delta,n=conformal_expansion(d)
    assert n==5 and delta==6


def test_raw_q1_top5_matches_full_order_while_representatives_are_separate():
    q=Path(__file__).resolve().parents[1]/'outputs/q1'
    if not (q/'Q1_all_candidate_ranking.csv').exists():return
    full=pd.read_csv(q/'Q1_all_candidate_ranking.csv');top=pd.read_csv(q/'Q1_top5_with_stability.csv')
    assert top.feature.tolist()==full.head(5).feature.tolist()
    assert full.selected_top5.sum()==5
    assert (q/'Q1_mechanism_representatives.csv').exists()


def test_slope_acceleration_error_decreases_with_grid_refinement():
    # First-order hydrostatic reconstruction must converge on a wet slope.
    expected=9.81*.5;errors=[]
    for dx in [2.,1.,.5]:
        z=1000-np.arange(35)*dx*.5;h=np.full(35,10.);q=np.zeros(35)
        hn,qn,_,_=advance_1d(h,q,np.zeros(35),np.ones(35),dx,.001,0,1e15,0,10,bed=z)
        acceleration=qn[8:25]/hn[8:25]/.001
        errors.append(float(np.mean(abs(acceleration-expected))))
    assert errors[2]<errors[1]<errors[0]
    assert errors[2]/expected<.02


def test_production_forecast_is_unchanged_when_future_observations_are_modified():
    from mathorcup.common import Context
    from mathorcup.q3_task_aligned import prepare_joint,anomaly_score,forecast_analogue
    root=Path(__file__).resolve().parents[1]
    if not (root/'outputs/q2/Q2_primary_trigger_assessment.json').exists():return
    h=prepare_joint(Context(root));issue=pd.Timestamp('2025-07-22 09:00')
    a,info=anomaly_score(h,[2021,2022,2023]);changed=h.copy()
    changed.loc[changed.index>=issue,'temp_2m']=45.
    changed.loc[changed.index>=issue,'precip_mm']=90.
    b,info2=anomaly_score(changed,[2021,2022,2023])
    assert info['threshold']==info2['threshold']
    fa,_,_=forecast_analogue(a.combined_condition_score.loc[a.index<issue],issue,n_samples=50,seed=17)
    fb,_,_=forecast_analogue(b.combined_condition_score.loc[b.index<issue],issue,n_samples=50,seed=17)
    assert np.allclose(fa,fb)


def test_secondary_2d_is_not_added_to_primary_warning_polygon():
    import json
    root=Path(__file__).resolve().parents[1]
    p=root/'outputs/q3/Q3_spatial_range_results.json'
    if not p.exists():return
    d=json.loads(p.read_text())
    if not d['2d_included_only_if_primary_1d_inflow']:
        assert d['local_2d_inundation_area_km2']==0
        assert d['standalone_2d_excluded'] is True


def test_current_wet_flow_speed_does_not_exceed_a_generous_energy_head_bound():
    import json
    root=Path(__file__).resolve().parents[1]/'outputs/q4'
    if not (root/'Q4_1d_settings_and_summary.json').exists():return
    profile=pd.read_csv(root/'Q4_route_profile.csv')
    drop=profile.elevation_m.max()-profile.elevation_m.min()
    for name in ['Q4_1d_settings_and_summary.json','Q4_secondary_coupling_case.json']:
        if not (root/name).exists():continue
        case=json.loads((root/name).read_text());m=case['mass'];p=case['parameters']
        bound=np.sqrt(p['initial_velocity_m_s']**2+2*9.81*(drop+2*m['max_depth_over_time_m']))
        assert m['max_speed_over_time_m_s']<=bound*1.1
