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

def sample_discrete_event_time(bin_ends_hours,conditional_mass,n,seed=7,uniform_within_bin=False):
    """Resample calculated first-passage times on their actual numerical time grid."""
    bins=np.asarray(bin_ends_hours,float);p=np.asarray(conditional_mass,float)
    if not np.isclose(p.sum(),1,atol=1e-9):raise ValueError('Conditional discrete time mass does not sum to one')
    if len(bins)<2 or bins[0]<=0 or not np.all(np.diff(bins)>0):raise ValueError('invalid time bins')
    if bins.shape!=p.shape or not np.isfinite(p).all() or np.any(p<0):raise ValueError('invalid time probabilities')
    rng=np.random.default_rng(seed)
    chosen=rng.choice(len(bins),size=n,p=p)
    if uniform_within_bin:
        starts=np.r_[0,bins[:-1]]
        return starts[chosen]+rng.random(n)*(bins[chosen]-starts[chosen])
    return bins[chosen]

def execute_joint_arrival(ctx):
    out=ctx.folder('joint')
    # Retire the former demonstration filenames when upgrading an existing run.
    for name in ['Q34_joint_uncertainty_MATHEMATICAL_DEMO.png','Q34_joint_uncertainty_MATHEMATICAL_DEMO.svg','Q34_demo_joint_samples.csv.gz']:
        (out/name).unlink(missing_ok=True)
    q3=ctx.out/'q3'/'Q3_source_initiation_time_distribution.csv'
    q4=ctx.out/'q4'/'Q4_reconstruction_geometry.json'
    massfile=ctx.out/'q4'/'Q4_1d_settings_and_summary.json'
    if not q3.exists() or not q4.exists():
        write_json(out/'Q34_status.json',{'status':'NO_INPUTS','claim':'not executable before Q3,Q4'});return
    a=pd.read_csv(q3);route=json.loads(q4.read_text())
    primary=json.loads((ctx.out/'q3/Q3_primary_competition_answer.json').read_text())
    if primary.get('source_geometry_sha256') and primary['source_geometry_sha256']!=route.get('geometry_contract_sha256'):
        raise ValueError('Q4 must inherit the source/terrain contract prepared in Q2')
    if not np.isclose(a.conditional_time_probability.sum(),1):
        # A new non-initiation result must not leave a previous arrival ensemble.
        (out/'Q34_computed_joint_samples.csv.gz').unlink(missing_ok=True)
        write_json(out/'Q34_status.json',{'status':'NO_MODEL_INITIATION_WITHIN_HORIZON',
            'source':'Q3 source-stability model produced no first-passage times',
            'nominal_95_interval_hours_from_raw_issuance':None,'empirical_joint_coverage':None,
            'task_impact_reference_covered':None,'forecast_warning_issued':False,
            'scientific_consequence':'this source-initiation hypothesis does not support a case warning; propagation alone cannot supply an onset forecast'})
        fig,ax=plt.subplots(figsize=(10,3.7),layout='constrained');ax.axis('off')
        ax.text(.5,.70,'本次源区起动模型未产生失稳时刻',ha='center',fontsize=17,color=COLORS[0])
        ax.text(.5,.43,'Q4的传播结果继续保留；联合到达时间区间不生成',ha='center',fontsize=12)
        ax.text(.5,.21,'预警模型需通过机制与独立资料检验，不能由事后传播时间反推为已验证的预报',ha='center',fontsize=10,color='#666')
        savefig(fig,out/'Q34_trigger_and_arrival_time.png')
        return
    if massfile.exists():
        m=json.loads(massfile.read_text())
        endpoint_time=m['mass'].get('first_arrival_to_downstream_end_s')
    else:endpoint_time=None
    reconstructed=ctx.out/'q4/Q4_cascade_reconstruction.json'
    if reconstructed.exists():endpoint_time=json.loads(reconstructed.read_text())['downstream']['arrival_s']
    n=4000;rng=np.random.default_rng(421)
    trigger=sample_discrete_event_time(a.bin_end_hours_from_issue,
        a.conditional_time_probability,n,seed=421,uniform_within_bin=primary.get('model_family')=='Q2_conditioned_event_time')
    # Propagation uncertainty now comes from ACTUAL solver outputs under scenario parameters,
    # not invented 20..60 m/s average-speed samples.
    sweep_file=ctx.out/'q4'/'Q4_parameter_sensitivity.csv'
    reconstructed_sweep=ctx.out/'q4/Q4_reconstruction_sensitivity.csv'
    if reconstructed_sweep.exists():sweep_file=reconstructed_sweep
    if sweep_file.exists():
        sweep=pd.read_csv(sweep_file)
        if 'arrival_s' in sweep:sweep['arrival_time_min']=sweep.arrival_s/60
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
        'propagation_path_verified':False,
        'source_mapping_assumption':'Q2 fixed source and terrain contract inherited by Q3 and Q4',
        'q4_ensemble_n_scenarios_with_reach':int(len(valid)),
        'trigger_time_empirically_calibrated':False,
        'trigger_time_source':primary['forecast_target'],
        'trigger_scenario_fraction':primary.get('model_failure_scenario_fraction'),
        'trigger_time_sampling':'piecewise-uniform density explicitly specified by the Q3 conditional time model' if primary.get('model_family')=='Q2_conditioned_event_time' else 'resample numerical grid support',
        'n_unique_Q4_scenarios':int(len(sweep)) if sweep_file.exists() else 1,
        'n_nonreaching_Q4_scenarios':int(sweep.arrival_time_min.isna().sum()) if sweep_file.exists() else int(endpoint_time is None),
        'ensemble_weights':'equal per distinct tested scenario; sensitivity-design weights, not calibrated probabilities',
        'nominal_95_interval_hours_from_raw_issuance':[float(q[0]),float(q[2])],
        'conditional_on_event_within_configured_horizon':ctx.cfg['q3']['primary_horizon_hours'],
        'empirical_joint_coverage':None,
        'note':'Arrival-time sensitivity combines the specified Q3 conditional event-time density with calculated Q4 scenario travel times; the event-clock simplification follows Q3.'})
    pd.DataFrame({'trigger_time_relative_hours':trigger,'travel_minutes_from_q4_scenarios':travel,
                 'arrival_relative_hours':arrival}).to_csv(out/'Q34_computed_joint_samples.csv.gz',index=False,compression='gzip')
    fig,axs=plt.subplots(1,3,figsize=(12.1,3.8))
    for ax,arr,col,title,xlab in zip(axs,[trigger,travel,arrival],COLORS[:3],
        ['源区触发时间情景','数值求解的传播时间集合','下游到达时间推演'],
        ['相对发布时刻 / h','传播 / min','相对发布时刻 / h']):
        ax.hist(arr,bins=30,color=col,alpha=.88);ax.set(xlabel=xlab,ylabel='模拟次数',title=title);ax.grid(axis='y',alpha=.14)
    fig.suptitle('泥石流触发与下游到达时间分布',fontsize=12)
    savefig(fig,out/'Q34_trigger_and_arrival_time.png')
    record_warning('Q34 joint interval is conditional on selected Q3 hazard and Q4 propagation parameters; not an independent coverage claim.')
    print('[Q34] model-conditional 95% arrival-time interval:',q[[0,2]],'hours; scenario ensemble status:',status)
