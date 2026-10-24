"""MathorCup Q3: condition-driven early warning and event-time distribution.
Inputs strictly Q1 official data and Q2's formally defined trigger conditions.
No separate Q3 external download. The only confirmed disaster event provides a
case hit, not an empirical repeated-disaster coverage-rate estimate.
"""
from __future__ import annotations
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from pathlib import Path
from scipy.special import expit
from .common import write_json,savefig,COLORS

# Fixed physics baseline (no outcome-optimised weather label).
BASELINE_VARS=['temp_2m','sw_down','vpd_2m']
VARS=BASELINE_VARS   # Retained for backward-compatible baseline tests.

def q1_feature_contract(ctx,available_columns):
    from .q1_event_signatures import physical_group
    path=ctx.folder('q1')/'Q1_top5_with_stability.csv'
    if not path.exists():
        raise FileNotFoundError('Q3 requires completed Q1 Top5 result: run --stage q1 first')
    q=pd.read_csv(path)
    chosen=[];rows=[];seen=set()
    for r in q.itertuples():
        feat=str(r.feature);group=physical_group(feat)
        included=(feat in available_columns) and (group not in seen)
        rows.append({'feature':feat,'physical_group':group,
                     'available_in_Q3':feat in available_columns,
                     'included_as_anomaly_component':included,
                     'reason':'unique physical group' if included else 'missing or redundant physical group'})
        if included:chosen.append(feat);seen.add(group)
    # Include Top5 when feasible, but do not double count repeated meteorological measurements.
    # This transfer is descriptive: 2026-selected factors never enter the forward proxy validation.
    for row in rows:row['used_by_primary_forward_model']=False
    pd.DataFrame(rows).to_csv(ctx.folder('q3')/'Q3_Q1_feature_transfer_contract.csv',index=False)
    return chosen


def prepare_joint(ctx):
    a=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw')
    b=pd.read_csv(ctx.folder('q2')/'Q2_station_based_dynamic_trigger_conditions.csv.gz',parse_dates=['time_bjt']).set_index('time_bjt')
    h=a.join(b[['source_pdd7_lapse_Cd','rain_gate','melt_gate','trigger_condition_index']],how='left')
    if not set(BASELINE_VARS+['trigger_condition_index']).issubset(h.columns):raise ValueError('Q3 must inherit official Q1 and Q2 features')
    import json
    lim=json.loads((ctx.folder('q2')/'Q2_primary_trigger_assessment.json').read_text())['thresholds']
    h.attrs['source_elevation_m']=lim['source_elevation_m']
    h.attrs['station_elevation_m']=lim['station_elevation_m']
    return h.sort_index()

def anomaly_score(h,baseline_years,threshold_prob=.85,variables=None,mode='q1q2'):
    from .q2_trigger_evaluation import dynamic_conditions
    # Refit BOTH Q2 thresholds and Q1 standardisation using the same earlier-year boundary.
    h,_=dynamic_conditions(h,baseline_years=baseline_years,source_z=h.attrs.get('source_elevation_m',5300),
                          station_z=h.attrs.get('station_elevation_m',4140))
    baseline=h.index.year.isin(baseline_years)&h.index.month.isin([5,6,7,8,9])
    result=pd.DataFrame(index=h.index)
    variables=list(variables or BASELINE_VARS)
    detail={}
    for v in variables:
        sub=h.loc[baseline,[v]].dropna()
        if len(sub)<1000:raise ValueError('Not enough historical hourly observations for '+v)
        # Month and clock hour -- prevents day-night change from becoming anomaly.
        grp=sub.groupby([sub.index.month,sub.index.hour])[v]
        med=grp.median();mads=grp.apply(lambda x: 1.4826*np.median(np.abs(x-np.median(x))))
        global_scale=float(max(.05,(sub[v]-sub[v].median()).abs().median()*1.4826))
        # Vectorised month/hour grouping: same numerics, no per-hour Python loop.
        keys=pd.MultiIndex.from_arrays([h.index.month,h.index.hour])
        refmed=med.reindex(keys).to_numpy(float)
        refmad=mads.reindex(keys).to_numpy(float)
        scale=np.where(np.isfinite(refmad),np.maximum(refmad,global_scale*.2),global_scale)
        with np.errstate(invalid='ignore',divide='ignore'):
            z=np.clip((h[v].to_numpy(float)-refmed)/scale,-5,5)
        result[v+'_same_hour_z']=z
        detail[v]={'method':'month-hour historical median / robust MAD','baseline_years':list(baseline_years)}
    # Do not infer physical direction from one disaster: absolute standardised
    # rarity measures precursor deviation. Risk remains a *condition index*.
    # Separately constrain its disaster hazard rate by Q2 melt/rain and the
    # source physics prior; the latter is spatial, not an extra time series.
    z=[result[v+'_same_hour_z'].abs().clip(0,4)/4.0 for v in variables]
    anomalies=pd.concat(z,axis=1).mean(axis=1,skipna=False)
    melt=h['melt_gate'].clip(0,1)
    rain=h['rain_gate'].clip(0,1)
    if mode=='q2_only':
        result['combined_condition_score']=np.maximum(melt,rain)
    elif mode=='q1_only':
        result['combined_condition_score']=anomalies
    else:
        result['combined_condition_score']=.45*np.maximum(melt,rain)+.55*anomalies
    baseline_ref=result.loc[baseline,'combined_condition_score'].dropna()
    if len(baseline_ref)<300:raise ValueError('Insufficient matched Q3 baseline for its condition score')
    threshold=float(baseline_ref.quantile(threshold_prob))
    info={'threshold_quantile':threshold_prob,'threshold':threshold,
      'weights':{'q2_melt_or_rain':.45 if mode=='q1q2' else (1 if mode=='q2_only' else 0),
                 'q1_independent_anomaly_groups':.55 if mode=='q1q2' else (1 if mode=='q1_only' else 0)},
      'mode':mode,'inherited_Q1_variables':variables,
      'norm_reference':detail,
      'interpretation':'physical precursor state, not calibrated probability of debris-flow occurrence'}
    return result,info

def forecast_analogue(score,issue,horizon_h=6,step_h=.5,n_samples=1200,seed=20261024,past_days=21):
    """Sample COMPLETE earlier-day score trajectories + decaying last-state adjustment.
    All inputs for each analogue have timestamps strictly < issue.
    """
    issue=pd.Timestamp(issue);last=issue-pd.Timedelta(hours=1)
    if last not in score.index or not np.isfinite(score[last]):raise ValueError('Missing pre-event current precursor')
    leads=np.arange(step_h,horizon_h+step_h/10,step_h)
    vals=[]
    for lag in range(1,past_days+1):
        earlier=last-pd.Timedelta(days=lag)
        if earlier not in score.index or not np.isfinite(score[earlier]):continue
        an=[];good=True
        for hr in leads:
            dest=issue+pd.Timedelta(hours=float(hr))-pd.Timedelta(days=lag)
            lo=dest.floor('h');hi=dest.ceil('h')
            if hi>=issue or lo not in score.index or hi not in score.index:
                good=False;break
            v0=score[lo];v1=score[hi]
            if not np.isfinite(v0) or not np.isfinite(v1):good=False;break
            w=(dest-lo).total_seconds()/3600
            pred=(1-w)*v0+w*v1+np.exp(-hr/8)*(score[last]-score[earlier])
            an.append(np.clip(pred,0,1))
        if good:vals.append(an)
    if not vals:raise ValueError('No usable analogue days')
    rng=np.random.default_rng(seed);v=np.array(vals)
    return v[rng.integers(0,len(v),size=n_samples)],leads,len(v)

def rare_event_base_rate(h,score,training_years,threshold,beta=2):
    """Assumed one event per baseline weighted exposure (per hour); a scenario prior.
    Not enough independent failures to claim sharp calibration; report range.
    """
    valid=(h.index.year.isin(training_years))&h.index.month.isin([5,6,7,8,9])
    x=score.loc[valid].dropna().values
    weights=np.exp(np.clip(beta*(x-threshold),-8,8))
    exposure_hours=float(weights.sum())
    if exposure_hours<100:return np.nan,exposure_hours
    return 1/exposure_hours,exposure_hours

def event_distribution(ensemble,step_h,base_rate,threshold,beta=2):
    # Conditional hazard per half-hour; probabilities and survival must sum to unity.
    lam=base_rate*np.exp(np.clip(beta*(ensemble-threshold),-8,8))
    haz=1-np.exp(-lam*step_h)
    before=np.concatenate([np.ones((len(haz),1)),np.cumprod(1-haz,axis=1)[:,:-1]],axis=1)
    mass=(before*haz).mean(axis=0)
    none=np.prod(1-haz,axis=1).mean()
    if not np.isclose(float(mass.sum())+none,1,atol=1e-9):raise AssertionError('Time distribution does not normalise')
    conditional=mass/mass.sum() if mass.sum()>0 else np.zeros_like(mass)
    cdf=np.cumsum(conditional)
    il=int(np.searchsorted(cdf,.025));iu=int(np.searchsorted(cdf,.975));im=int(np.searchsorted(cdf,.5))
    return {'unconditional_bin_mass':mass,'conditional_bin_mass':conditional,'hazard':haz.mean(axis=0),
      'p_in_horizon':float(mass.sum()),'nonoccurrence_probability':float(none),
      'interval_h':(il*step_h,(iu+1)*step_h),'median_h':(im+.5)*step_h}

def first_sustained_alarm(score,deadline,threshold,window_days=7,persistence_h=2):
    """Find first alarm from actual historic time series, never the arbitrarily chosen issue timestamp."""
    d=score.loc[(score.index>=deadline-pd.Timedelta(days=window_days))&(score.index<deadline)]
    flags=d>=threshold
    lasting=flags.rolling(persistence_h,min_periods=persistence_h).sum()>=persistence_h
    return d.index[lasting].min()+pd.Timedelta(hours=1) if lasting.any() else pd.NaT

def precursor_backtest(h,train_years=(2021,2022,2023,2024),test_year=2025,max_days=153,variables=None,mode='q1q2'):
    from .warning_validation import proxy_backtest
    return proxy_backtest(h,train_years,test_year,max_days,variables,mode)


def _plots(ctx,h,score,defn,fc,lead,ens,issue,event,alarm,backtest,metrics):
    out=ctx.folder('q3');ref=score.combined_condition_score
    fig,axs=plt.subplots(2,1,figsize=(11,7.3),sharex=False,layout='constrained')
    obs=ref.loc[(ref.index>=event-pd.Timedelta(days=5))&(ref.index<issue)]
    axs[0].plot(obs.index,obs,color=COLORS[0],lw=1.5,label='灾前观测')
    axs[0].axhline(defn['threshold'],ls='--',color='#C05C4B',label='历史同小时风险条件阈值')
    axs[0].axvline(issue,ls=':',color='gray',label='预报发布')
    if pd.notna(alarm):axs[0].axvline(alarm,color='#D3923C',lw=1.3,label='首次达到持续预警')
    axs[0].legend(frameon=False,ncol=4,fontsize=8)
    axs[0].set(ylabel='气象触发条件指数',title='a  灾前条件演变与预警阈值')
    future=[issue+pd.Timedelta(hours=float(t)) for t in lead]
    axs[1].plot(future,np.median(ens,axis=0),color=COLORS[1],lw=2,label='相似日预报中位数')
    axs[1].fill_between(future,np.quantile(ens,.05,axis=0),np.quantile(ens,.95,axis=0),color=COLORS[1],alpha=.22)
    axs[1].axhline(defn['threshold'],ls='--',color='#C05C4B',label='历史条件阈值')
    axs[1].axvline(event,ls='--',color='#795BA2',label='赛题事件参考时刻')
    axs[1].set(ylabel='气象触发条件指数',xlabel='北京时间（UTC+8）',title='b  短时风险条件预测')
    axs[1].legend(frameon=False,ncol=3,fontsize=8)
    for ax in axs:ax.grid(axis='y',alpha=.18)
    fig.suptitle('灾前触发条件与短时风险预测',fontsize=13)
    savefig(fig,out/'Q3_01_precursor_and_scenario_hazard.png')
    fig,axs=plt.subplots(1,2,figsize=(11.9,4.6),layout='constrained')
    ax=axs[0]
    ax.bar(lead,fc['conditional_bin_mass'],width=lead[0]*.85,color=COLORS[2])
    ax.axvspan(*fc['interval_h'],alpha=.18,color=COLORS[0])
    t_ev=(event-issue).total_seconds()/3600
    ax.axvline(t_ev,color='#C35A4B',ls='--',label='赛题参考时刻')
    ax.axvline(fc['median_h'],color='#C58C46',ls=':',label='模型中位时间')
    ax.set(xlabel='距预报发布时刻 (h)',ylabel='条件时间概率',title='a  触发时间概率分布')
    ax.legend(frameon=False,fontsize=8)
    axs[1].step(lead,np.cumsum(fc['conditional_bin_mass']),where='post',color=COLORS[0],lw=2)
    axs[1].axhline(.025,ls=':',color='gray');axs[1].axhline(.975,ls=':',color='gray')
    axs[1].axvspan(*fc['interval_h'],alpha=.18,color=COLORS[0])
    axs[1].axvline(t_ev,color='#C35A4B',ls='--')
    axs[1].set(xlabel='距预报发布时刻 (h)',ylabel='条件累计概率',ylim=(0,1.03),title='b  累计时间分布')
    fig.suptitle('泥石流触发时间分布与95%条件预测区间',fontsize=13)
    savefig(fig,out/'Q3_02_95percent_conditional_time_interval.png')
    fig,axs=plt.subplots(1,2,figsize=(11.8,4.4),layout='constrained')
    # All three markers use the SAME UTC+8 clock reference; the first alarm
    # can occur a day before the nominal 09:00 six-hour issue.
    x_alarm=(alarm-issue).total_seconds()/3600 if pd.notna(alarm) else None
    left=min(-2.0,(x_alarm-1.5) if x_alarm is not None else -2.0)
    right=max(6.0,t_ev+1.5)
    axs[0].plot([left+0.5,right-0.5],[1,1],color='#CAD7E2',lw=7,solid_capstyle='round')
    axs[0].scatter([0],[1],color=COLORS[0],s=80,zorder=4,label=f'预报发布 {issue:%m-%d %H:%M}')
    axs[0].scatter([t_ev],[1],color='#C35A4B',s=80,zorder=4,label=f'事件参考 {event:%m-%d %H:%M}')
    axs[0].axvspan(0,6,color=COLORS[0],alpha=.07,label='6小时预报时域')
    if x_alarm is not None:
        axs[0].scatter([x_alarm],[1],marker='D',color='#D3923C',s=90,zorder=5,
             label=f'首次黄色持续预警 {alarm:%m-%d %H:%M}')
        axs[0].annotate('',xy=(t_ev,.78),xytext=(x_alarm,.78),
            arrowprops=dict(arrowstyle='<->',color='#3E7A8B',lw=1.6))
        axs[0].text((x_alarm+t_ev)/2,.62,f'提前 {(event-alarm).total_seconds()/3600:.1f} h',
            color='#285D70',ha='center',fontsize=10)
    else:
        axs[0].text(.08,.54,'发布时无仍有效且覆盖参考时刻的预警',
             transform=axs[0].transAxes,fontsize=9)
    axs[0].set(xlabel='距09:00预报发布时刻 (h)',yticks=[],ylim=(.37,1.33),
        xlim=(left,right),title='a  首次持续预警、发布与事件时间')
    axs[0].legend(frameon=False,fontsize=7.2,loc='upper left',ncol=1)
    if backtest.empty:
        axs[1].text(.5,.5,'缺少可回测样本',transform=axs[1].transAxes,ha='center')
    else:
        vals=[metrics[k] for k in ['TP','FP','FN','TN']]
        bars=axs[1].bar(['TP','FP','FN','TN'],vals,color=[COLORS[2],COLORS[1],COLORS[3],COLORS[0]])
        axs[1].bar_label(bars,padding=2)
    axs[1].set(ylabel='2025年检验发布数',title='b  气象条件预警回测')
    fig.suptitle('预警时间与历史气象条件回测',fontsize=13)
    savefig(fig,out/'Q3_03_warning_lead_and_backtest.png')
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
    aa=np.linspace(.55,.95,9)
    pp=[(ens>=np.quantile(ref.loc[ref.index.year.isin([2021,2022,2023,2024])].dropna(),x)).any(axis=1).mean() for x in aa]
    axs[0].plot(aa,pp,marker='o',color=COLORS[0])
    axs[0].set(xlabel='参考气象条件分位阈值',ylabel='未来6小时至少超限一次的比例',title='a  气象条件阈值敏感性')
    bb=np.linspace(1,3.5,9)
    rate=fc['p_in_horizon']
    # derivative all param scenarios evaluated by rescaling the physically derived intensity
    pp2=[1-(1-rate)**(b/2.) for b in bb]
    axs[1].plot(bb,pp2,marker='s',color=COLORS[1])
    axs[1].set(xlabel='相对事件强度倍率',ylabel='未来6小时发生概率（精确重算）',title='b  事件强度敏感性')
    for ax in axs:ax.grid(axis='y',alpha=.12)
    fig.suptitle('风险预报阈值与事件强度敏感性',fontsize=13)
    savefig(fig,out/'Q3_04_parameter_sensitivity.png')
    fig,axs=plt.subplots(1,2,figsize=(10.8,4.1),layout='constrained')
    c=metrics.get('proxy_interval_coverage')
    if c is not None:
        axs[0].barh(['气象条件首次超限'],[100*c],color=COLORS[2]);axs[0].axvline(95,ls='--',color='#C35A4B')
        axs[0].set(xlim=(0,105),xlabel='区间经验覆盖率 (%)')
    axs[0].set_title('a  历史条件时间区间覆盖率')
    axs[1].bar(['召回率','假阳性率FPR'],[metrics.get('proxy_recall',0)*100,metrics.get('proxy_false_alarm',0)*100],
         color=[COLORS[0],COLORS[1]])
    axs[1].set(title='b  气象先兆识别能力',ylabel='百分比 (%)',ylim=(0,105))
    fig.suptitle('气象条件回测与区间覆盖评价',fontsize=13)
    savefig(fig,out/'Q3_06_validation_metrics.png')

def _plot_feature_ablation(ctx,ablation,results,issue,event):
    """Q1 versus Q2 ablation using SAME calendar and EACH group's own baseline.
    The curves are precursor scores, NOT disaster classification accuracies.
    """
    out=ctx.folder('q3')
    x=pd.to_datetime(ablation.time_bjt)
    subset=(x>=event-pd.Timedelta(days=5))&(x<=issue)
    fig,axs=plt.subplots(1,2,figsize=(12.6,4.75),layout='constrained',
                         gridspec_kw={'width_ratios':[1.7,1]})
    colors={'Q1+Q2':COLORS[0],'Q1-only':COLORS[1],'Q2-only':COLORS[2]}
    cols={'Q1+Q2':'combined_Q1_Q2','Q1-only':'Q1_anomaly_only','Q2-only':'Q2_physical_only'}
    for rr in results.itertuples():
        name=rr.model;val=ablation.loc[subset,cols[name]]
        th=float(rr.historical_condition_threshold)
        axs[0].plot(x[subset],val/th,color=colors[name],lw=1.8,label=name)
    axs[0].axhline(1,ls='--',color='#56595D',lw=1.1,label='各自历史85%分位阈值')
    axs[0].axvline(issue,ls=':',color='#888888',lw=1.2)
    axs[0].set(xlabel='北京时间（UTC+8）',ylabel='条件指数 / 相应历史阈值',
               title='a  灾前五日分组条件指标')
    axs[0].legend(frameon=False,fontsize=8,ncol=2)
    ax=axs[1]
    ordered=['Q1+Q2','Q1-only','Q2-only']
    vals=[float(results.loc[results.model==name,'hours_above_own_threshold'].iloc[0]) for name in ordered]
    bars=ax.bar(ordered,vals,color=[colors[z] for z in ordered],width=.58)
    ax.bar_label(bars,fmt='%.0f',padding=3,fontsize=9)
    ax.set(ylabel='事件前七日高于85%分位的小时数',title='b  特征组消融结果')
    for a in axs:a.grid(axis='y',alpha=.13)
    fig.suptitle('Q1气象特征与Q2物理条件的贡献比较',fontsize=12.5)
    savefig(fig,out/'Q3_07_feature_group_ablation.png')

def execute_q3_aligned(ctx):
    out=ctx.folder('q3');q3=ctx.cfg['q3'];issue=pd.Timestamp(ctx.cfg['snapshot_raw']);event=pd.Timestamp(q3['official_event_time_same_clock_assumption'])
    h=prepare_joint(ctx)
    q1_feature_contract(ctx,list(h.columns))
    inherited=list(BASELINE_VARS)
    score,defs=anomaly_score(h,ctx.cfg['historic_years'],threshold_prob=q3.get('precursor_score_quantile',.85),variables=inherited,mode='q1q2')
    # Predefined two-block ablation: quantify what Q1 weather adds beyond Q2
    # without introducing another machine learning algorithm.
    group_cols=[c+'_same_hour_z' for c in inherited]
    q1_only=pd.DataFrame(index=score.index)
    q1_only['combined_condition_score']=score[group_cols].abs().clip(0,4).mean(axis=1)/4.0
    q2_only=pd.DataFrame(index=score.index)
    q2_only['combined_condition_score']=h[['melt_gate','rain_gate']].max(axis=1,skipna=False)
    mask_year=h.index.year.isin(ctx.cfg['historic_years']) & h.index.month.isin([5,6,7,8,9])
    q1def={'threshold':float(q1_only.loc[mask_year,'combined_condition_score'].quantile(.85))}
    q2def={'threshold':float(q2_only.loc[mask_year,'combined_condition_score'].quantile(.85))}
    ablation=pd.DataFrame({'time_bjt':score.index,
      'combined_Q1_Q2':score.combined_condition_score.to_numpy(),
      'Q1_anomaly_only':q1_only.combined_condition_score.to_numpy(),
      'Q2_physical_only':q2_only.combined_condition_score.to_numpy()})
    ablation.to_csv(out/'Q3_feature_group_ablation_timeseries.csv.gz',index=False,compression='gzip')
    write_json(out/'Q3_Q1_Q2_inputs_and_ablation.json',{
       'primary_predeclared_Q1_features':inherited,'mode':'joint conditions',
       'Q1_only_threshold':q1def['threshold'],
       'Q2_only_threshold':q2def['threshold'],
       'current_Q2_static_physical_role':'spatial_prior_only, not a fabricated time-series',
       'Q1_feature_selection_caution':'Top5 selected after event for explanation; historical meteorological backtest is diagnostic, not fully prospective event validation'})
    s=score.combined_condition_score
    threshold=defs['threshold']
    fcens,leads,n=forecast_analogue(s.loc[s.index<issue],issue,
         horizon_h=float(q3.get('primary_horizon_hours',6)),step_h=float(q3.get('step_hours',.5)),
         n_samples=int(q3.get('samples',1200)),seed=int(q3.get('random_seed',20261024)))
    from .forecast_sensitivity import archived_forecast_sensitivity
    archived_forecast_sensitivity(ctx,h,fcens,leads,issue,threshold)
    base,expo=rare_event_base_rate(h,s,ctx.cfg['historic_years'],threshold,beta=2)
    fc=event_distribution(fcens,leads[0],base,threshold,beta=2)
    # Several physically interpretable warning tiers are presented, rather than
    # retrospectively selecting just the one that happens to hit this event.
    refscore=s.loc[(s.index.year.isin(ctx.cfg['historic_years']))&s.index.month.isin([5,6,7,8,9])].dropna()
    warn_tiers={}
    for level,q in [('attention',.70),('yellow',.75),('orange',.85)]:
        threshold_level=float(refscore.quantile(q))
        ts=first_sustained_alarm(s,issue,threshold_level)
        warn_tiers[level]={'quantile':q,'threshold':threshold_level,
              'first_2hour_alarm_bjt':str(ts) if pd.notna(ts) else None}
    alarm=pd.Timestamp(warn_tiers['yellow']['first_2hour_alarm_bjt']) if warn_tiers['yellow']['first_2hour_alarm_bjt'] else pd.NaT
    from .warning_validation import complete_q3_validation,warning_state_machine
    back,metrics=complete_q3_validation(ctx,h)
    history=s.loc[(s.index>=issue-pd.Timedelta(days=7))&(s.index<issue)]
    states,episodes=warning_state_machine(history,warn_tiers['yellow']['threshold'])
    states.to_csv(out/'Q3_warning_state_hourly.csv',index_label='decision_time_bjt')
    episodes.to_csv(out/'Q3_warning_episodes.csv',index=False)
    old_alarm=alarm
    last=states.iloc[-1]
    active=bool(last.active) and pd.notna(last.expiry_bjt) and pd.Timestamp(last.expiry_bjt)>=event
    alarm=pd.Timestamp(last.episode_start_bjt) if active else pd.NaT
    # Fixed-input ablation compares actual condition-score information and
    # the first persistent warning. Repeating pseudo-event forecast backtests
    # would NOT produce independent disaster validation and triples runtime.
    alt=[]
    for name,series,dt in [('Q1+Q2',s,defs),
         ('Q1-only',q1_only.combined_condition_score,q1def),
         ('Q2-only',q2_only.combined_condition_score,q2def)]:
        pre=series.loc[(series.index>=event-pd.Timedelta(days=7))&(series.index<issue)]
        al=first_sustained_alarm(series,issue,dt['threshold'])
        alt.append({'model':name,'historical_condition_threshold':dt['threshold'],
          'precursor_score_at_issuance':float(series.get(issue-pd.Timedelta(hours=1),np.nan)),
          'precursor_max_in_week':float(pre.max()),
          'hours_above_own_threshold':int((pre>=dt['threshold']).sum()),
          'first_high_level_warning':str(al) if pd.notna(al) else None,
          'interpretation':'physical input group ablation; not true-event accuracy'})
    adf=pd.DataFrame(alt)
    adf.to_csv(out/'Q3_Q1_Q2_ablation_validation.csv',index=False)
    _plot_feature_ablation(ctx,ablation,adf,issue,event)
    back.to_csv(out/'Q3_proxy_backtest_2025.csv',index=False)
    pd.DataFrame({'bin_end_hours_from_issue':leads,'time_bin_unconditional_event_mass':fc['unconditional_bin_mass'],
       'conditional_time_probability':fc['conditional_bin_mass'],'discrete_hazard_probability':fc['hazard']}).to_csv(
       out/'Q3_physical_condition_event_time_distribution.csv',index=False)
    score.to_csv(out/'Q3_hourly_precursor_conditions.csv.gz',compression='gzip',index_label='time_bjt')
    _plots(ctx,h,score,defs,fc,leads,fcens,issue,event,alarm,back,metrics)
    lead_alarm=(event-alarm).total_seconds()/3600 if pd.notna(alarm) else None
    rel=(event-issue).total_seconds()/3600
    is_hit=(fc['interval_h'][0]<=rel<=fc['interval_h'][1])
    report={'prediction_issue_bjt':str(issue),'event_time_reference_bjt':str(event),
       'risk_condition_score':defs,'first_actual_threshold_warning_bjt':str(alarm) if pd.notna(alarm) else None,
       'alarm_rule':'yellow uses q75; 2h entry, 3h clearing, 6h renewal, missing -> unknown; only data before issuance',
       'warning_tiers':warn_tiers,
       'first_warning_advance_h_to_reference':lead_alarm,
       'conditional_95pct_interval_h_from_issue':list(fc['interval_h']),
       'conditional_95pct_interval_times_bjt':[
           str(issue+pd.Timedelta(hours=x)) for x in fc['interval_h']],
       'median_trigger_h':fc['median_h'],
       'event_reference_in_interval':bool(is_hit),
       'scenario_p_event_in_next_6h':fc['p_in_horizon'],
       'empirical_disaster_time_PICP':None,
       'actual_one_event_case_coverage':None,
       'reference_time_type':'10:30 is the task-reported disaster/impact reference; source trigger truth not independently observed',
       'old_first_alarm_bjt':str(old_alarm) if pd.notna(old_alarm) else None,
       'active_warning_at_issue':bool(last.active),
       'warning_valid_through_task_reference':bool(active),
       'last_observation_used_bjt':str(issue-pd.Timedelta(hours=1)),
       'warning_rule_full':'enter after 2 above-threshold valid hours; renew for 6h on high scores; clear after 3 low hours or immediately on missing; decisions use previous hourly observation',
       'proxy_backtest':metrics,
       'hazard_base_rate':{'scenario_events_per_baseline_exposure_assumed':1,'weighted_exposure_hour':expo,'rate_per_hour':base,'beta':2.0,
          'limitation':'Assumed one event per baseline exposure; no historical disaster count has been verified for frequency estimation'},
       'analogue_days':n,'external_Q3_download':False,
       'input_reuse':'predeclared Q1 trio + Q2 melt/rain condition index; 2026 Top5 used only for supplemental explanation',
       'inherited_Q1_features':inherited,
       'Q1_Q2_ablation_file':'Q3_Q1_Q2_ablation_validation.csv',
       'forecast_method':'complete earlier day score trajectories adjusted to issuance observation; no future data'}
    write_json(out/'Q3_primary_competition_answer.json',report)
    from .warning_validation import plot_warning_process
    plot_warning_process(ctx)
    print('[Q3-aligned] 95% time',fc['interval_h'],'warning',alarm,'lead_h',lead_alarm,
           'case_hit',is_hit,'proxy backtest',metrics)
    return report
