"""Q3+Q4 physical-scenario coupling: T_arrival = T_trigger + T_travel.

Uses real Q4 solver ensemble travel times conditioned on downstream reachability;
reserves event-level coverage claims for independent observed disasters.
"""
from __future__ import annotations
import numpy as np,pandas as pd,matplotlib.pyplot as plt,json
from .common import savefig,write_json,COLORS,record_warning

def joint_arrival_samples(trigger_hours,travel_minutes):
    ts=np.asarray(trigger_hours,float);tt=np.asarray(travel_minutes,float)
    if ts.shape!=tt.shape or ts.ndim!=1:raise ValueError('Paired samples must have same one-dimensional shape')
    if np.any(tt<0) or not np.isfinite(ts+tt).all():raise ValueError('Invalid time inputs')
    return ts+tt/60

def sample_discrete_event_time(bin_ends_hours,conditional_mass,n,seed=7):
    """Sample a time within the chosen event bin; conditional_mass must sum to one."""
    bins=np.asarray(bin_ends_hours,float);p=np.asarray(conditional_mass,float)
    if not np.isclose(p.sum(),1,atol=1e-9):raise ValueError('Conditional discrete time mass does not sum to one')
    if len(bins)<2 or not np.all(np.diff(bins)>0):raise ValueError('invalid time bins')
    rng=np.random.default_rng(seed)
    ii=rng.choice(len(bins),size=n,p=p)
    edges=np.r_[0,bins[:-1]]
    return edges[ii]+rng.random(n)*(bins[ii]-edges[ii])

def execute_joint_arrival(ctx):
    out=ctx.folder('joint')
    q3=ctx.out/'q3'/'Q3_physical_condition_event_time_distribution.csv'
    q4=ctx.out/'q4'/'Q4_route_quality_gate.json'
    massfile=ctx.out/'q4'/'Q4_1d_settings_and_summary.json'
    if not q3.exists() or not q4.exists():
        write_json(out/'Q34_status.json',{'status':'NO_INPUTS','claim':'not executable before Q3,Q4'});return
    a=pd.read_csv(q3);route=json.loads(q4.read_text())
    if massfile.exists():
        m=json.loads(massfile.read_text())
        endpoint_time=m['mass'].get('first_arrival_to_downstream_end_s')
    else:endpoint_time=None
    n=4000;rng=np.random.default_rng(421)
    trigger=sample_discrete_event_time(a.bin_end_hours_from_issue,
        a.conditional_time_probability,n,seed=421)
    # Propagation uncertainty now comes from ACTUAL solver outputs under scenario parameters,
    # not invented 20..60 m/s average-speed samples.
    sweep_file=ctx.out/'q4'/'Q4_parameter_sensitivity.csv'
    if sweep_file.exists():
        sweep=pd.read_csv(sweep_file)
        if 'scenario_id' in sweep:sweep=sweep.drop_duplicates('scenario_id')
        valid=sweep['arrival_time_min'].dropna().to_numpy(float)
    else:valid=np.array([])
    if len(valid)>=2:
        travel=rng.choice(valid,size=n,replace=True)
        source='computed 1D parameter-sensitivity ensemble; conditional on downstream reachability'
        status='PHYSICS_SCENARIO_ENSEMBLE_CONDITIONAL_ON_REACH'
    elif endpoint_time is not None:
        travel=np.full(n,float(endpoint_time)/60.)
        source='one physical parameter scenario travel time; interval is conditional on selected parameters'
        status='PHYSICS_SCENARIO_SINGLE_TRAVEL_VALUE'
    else:
        write_json(out/'Q34_status.json',{'status':'NO_Q4_TRAVEL_RESULT',
           'reason':'No numerically computed downstream arrival; did not fabricate speed samples'})
        return None
    arrival=joint_arrival_samples(trigger,travel)
    q=np.quantile(arrival,[.025,.5,.975])
    write_json(out/'Q34_status.json',{
        'status':status,'source':source,'speed_band_m_s':None,
        'propagation_path_verified':route['observed_2026_event_route_verified'],
        'q4_ensemble_n_scenarios_with_reach':int(len(valid)),
        'trigger_time_empirically_calibrated':False,
        'n_unique_Q4_scenarios':int(len(sweep)) if sweep_file.exists() else 1,
        'n_nonreaching_Q4_scenarios':int(sweep.arrival_time_min.isna().sum()) if sweep_file.exists() else int(endpoint_time is None),
        'ensemble_weights':'equal per distinct tested scenario; sensitivity-design weights, not calibrated probabilities',
        'nominal_95_interval_hours_from_raw_issuance':[float(q[0]),float(q[2])],
        'conditional_on_event_within_configured_horizon':ctx.cfg['q3']['primary_horizon_hours'],
        'empirical_joint_coverage':None,
        'note':'Physical-scenario arrival uncertainty combines Q3 assumed risk time with calculated Q4 travel-time sensitivity; conditional on input assumptions.'})
    pd.DataFrame({'trigger_time_relative_hours':trigger,'travel_minutes_from_q4_scenarios':travel,
                 'arrival_relative_hours':arrival}).to_csv(out/'Q34_demo_joint_samples.csv.gz',index=False,compression='gzip')
    fig,axs=plt.subplots(1,3,figsize=(12.1,3.8))
    for ax,arr,col,title,xlab in zip(axs,[trigger,travel,arrival],COLORS[:3],
        ['源区触发时间情景','数值求解的传播时间集合','下游到达时间推演'],
        ['相对发布时刻 / h','传播 / min','相对发布时刻 / h']):
        ax.hist(arr,bins=30,color=col,alpha=.88);ax.set(xlabel=xlab,ylabel='模拟次数',title=title);ax.grid(axis='y',alpha=.14)
    fig.suptitle('泥石流触发与下游到达时间分布',fontsize=12)
    savefig(fig,out/'Q34_joint_uncertainty_MATHEMATICAL_DEMO.png')
    record_warning('Q34 joint interval is conditional on selected Q3 hazard and Q4 propagation parameters; not an independent coverage claim.')
    print('[Q34] model-conditional 95% arrival-time interval:',q[[0,2]],'hours; scenario ensemble status:',status)
