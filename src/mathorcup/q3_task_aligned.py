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

def q1_feature_contract(ctx,available_columns,primary_features=None):
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
    for row in rows:row['used_by_primary_forward_model']=row['feature'] in set(primary_features or [])
    pd.DataFrame(rows).to_csv(ctx.folder('q3')/'Q3_Q1_feature_transfer_contract.csv',index=False)
    return chosen


def prepare_joint(ctx):
    a=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw')
    b=pd.read_csv(ctx.folder('q2')/'Q2_station_based_dynamic_trigger_conditions.csv.gz',parse_dates=['time_bjt']).set_index('time_bjt')
    columns=['source_pdd7_lapse_Cd','rain_gate','melt_gate','trigger_condition_index']
    columns += [c for c in ['energy_gate','station_positive_energy_7d_mmwe'] if c in b]
    h=a.join(b[columns],how='left')
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
    thermal_water=h['trigger_condition_index'].clip(0,1)
    if mode=='q2_only':
        result['combined_condition_score']=thermal_water
    elif mode=='q1_only':
        result['combined_condition_score']=anomalies
    else:
        result['combined_condition_score']=.45*thermal_water+.55*anomalies
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

def empirical_first_passage(ensemble,step_h,threshold):
    """First threshold passage in resampled observed trajectories.

    Event frequency is measured from the input trajectories. No disaster count,
    baseline hazard rate, or event time is invented to normalise the result.
    """
    a=np.asarray(ensemble,float)
    if a.ndim!=2 or not np.isfinite(a).all():raise ValueError('Finite trajectory matrix required')
    flags=a>=threshold;cross=flags.any(axis=1);mass=np.zeros(a.shape[1])
    if cross.any():mass=np.bincount(np.argmax(flags[cross],axis=1),minlength=a.shape[1])/len(a)
    probability=float(cross.mean());none=1-probability
    conditional=mass/probability if probability else np.zeros_like(mass)
    survival_before=1-np.r_[0,np.cumsum(mass)[:-1]]
    hazard=np.divide(mass,survival_before,out=np.zeros_like(mass),where=survival_before>0)
    if not probability:
        return {'unconditional_bin_mass':mass,'conditional_bin_mass':conditional,'hazard':hazard,
          'p_in_horizon':0.,'nonoccurrence_probability':1.,'interval_h':(None,None),
          'median_h':None,'time_target':'condition_threshold_first_passage'}
    cdf=np.cumsum(conditional);il=int(np.searchsorted(cdf,.025));iu=int(np.searchsorted(cdf,.975));im=int(np.searchsorted(cdf,.5))
    return {'unconditional_bin_mass':mass,'conditional_bin_mass':conditional,'hazard':hazard,
      'p_in_horizon':probability,'nonoccurrence_probability':none,
      'interval_h':(il*step_h,(iu+1)*step_h),'median_h':(im+1)*step_h,
      'time_target':'condition_threshold_first_passage'}

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
    savefig(fig,out/'Q3_R01_weather_precursor_reference.png')
    fig,axs=plt.subplots(1,2,figsize=(11.9,4.6),layout='constrained')
    ax=axs[0]
    ax.bar(lead,fc['conditional_bin_mass'],width=lead[0]*.85,color=COLORS[2])
    if fc['interval_h'][0] is not None:ax.axvspan(*fc['interval_h'],alpha=.18,color=COLORS[0])
    t_ev=(event-issue).total_seconds()/3600
    ax.axvline(t_ev,color='#C35A4B',ls='--',label='赛题参考时刻')
    if fc['median_h'] is not None:ax.axvline(fc['median_h'],color='#C58C46',ls=':',label='模型中位时间')
    ax.set(xlabel='距预报发布时刻 (h)',ylabel='条件时间概率',title='a  条件首次超限时间分布')
    ax.legend(frameon=False,fontsize=8)
    axs[1].step(lead,np.cumsum(fc['conditional_bin_mass']),where='post',color=COLORS[0],lw=2)
    axs[1].axhline(.025,ls=':',color='gray');axs[1].axhline(.975,ls=':',color='gray')
    if fc['interval_h'][0] is not None:axs[1].axvspan(*fc['interval_h'],alpha=.18,color=COLORS[0])
    axs[1].axvline(t_ev,color='#C35A4B',ls='--')
    axs[1].set(xlabel='距预报发布时刻 (h)',ylabel='条件累计概率',ylim=(0,1.03),title='b  累计时间分布')
    fig.suptitle('触发条件首次超限时间分布与95%预测区间',fontsize=13)
    savefig(fig,out/'Q3_R02_weather_passage_reference.png')
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
             label=f'模型预警发布 {alarm:%m-%d %H:%M}')
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
    savefig(fig,out/'Q3_R03_weather_warning_reference.png')
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
    aa=np.linspace(.55,.95,9)
    baseline=ref.loc[ref.index.year.isin(ctx.cfg['historic_years']) & ref.index.month.isin([5,6,7,8,9])].dropna()
    pp=[float((ens>=np.quantile(baseline,x)).any(axis=1).mean()) for x in aa]
    axs[0].plot(aa,pp,marker='o',color=COLORS[0])
    axs[0].set(xlabel='参考气象条件分位阈值',ylabel='未来6小时至少超限一次的比例',title='a  气象条件阈值敏感性')
    days=[7,14,21,28];pp2=[];counts=[]
    for window in days:
        try:
            trials,_,count=forecast_analogue(ref.loc[ref.index<issue],issue,
                horizon_h=ctx.cfg['q3']['primary_horizon_hours'],step_h=lead[0],
                n_samples=ctx.cfg['q3']['samples'],seed=ctx.cfg['q3']['random_seed'],past_days=window)
            pp2.append(float((trials>=defn['threshold']).any(axis=1).mean()));counts.append(count)
        except ValueError:pp2.append(np.nan);counts.append(0)
    axs[1].plot(days,pp2,marker='s',color=COLORS[1])
    axs[1].set(xlabel='相似日回看窗口 (天)',ylabel='基准模型6h条件超限频率',title='b  实际轨迹重算的窗口敏感性')
    rows=[{'experiment':'reference_quantile','setting':x,'cross_frequency':p,'n_analogue_days':None} for x,p in zip(aa,pp)]
    rows += [{'experiment':'analogue_lookback_days','setting':x,'cross_frequency':p,'n_analogue_days':n} for x,p,n in zip(days,pp2,counts)]
    pd.DataFrame(rows).to_csv(out/'Q3_sensitivity_recomputed_experiments.csv',index=False)
    for ax in axs:ax.grid(axis='y',alpha=.12)
    fig.suptitle('条件预报阈值与历史轨迹窗口敏感性',fontsize=13)
    savefig(fig,out/'Q3_R04_weather_sensitivity_reference.png')
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
    savefig(fig,out/'Q3_R06_weather_validation_reference.png')

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

def execute_q3_weather_reference(ctx):
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
    q2_only['combined_condition_score']=h.trigger_condition_index
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
    weather=archived_forecast_sensitivity(ctx,h,fcens,leads,issue,threshold)
    forecast_source='Q1 observed-analogue trajectories + Q2 actual archived weather'
    if weather is not None:fcens=weather['ensemble']
    else:forecast_source='complete earlier observed-condition trajectories'
    fc=empirical_first_passage(fcens,leads[0],threshold)
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
    from .warning_validation import complete_q3_validation,warning_state_machine,decide_warning_at_issue
    back,metrics=complete_q3_validation(ctx,h)
    history=s.loc[(s.index>=issue-pd.Timedelta(days=7))&(s.index<issue)]
    states,episodes=warning_state_machine(history,warn_tiers['yellow']['threshold'])
    states.to_csv(out/'Q3_warning_state_hourly.csv',index_label='decision_time_bjt')
    episodes.to_csv(out/'Q3_warning_episodes.csv',index=False)
    old_alarm=alarm
    last=states.iloc[-1]
    active=bool(last.active) and pd.notna(last.expiry_bjt) and pd.Timestamp(last.expiry_bjt)>=issue
    alarm=pd.Timestamp(last.episode_start_bjt) if active else pd.NaT
    # A forward warning is issued on predicted threshold passage, not only on
    # an observation alarm that happened to remain active from yesterday.
    probability_threshold=float(q3.get('warning_cross_probability_threshold',.25))
    decision=decide_warning_at_issue(issue,fc['p_in_horizon'],probability_threshold,
        active,alarm,last.expiry_bjt,float(q3.get('primary_horizon_hours',6)))
    forecast_alert=decision['forecast_issued']
    operational_alarm=decision['warning_time'] if decision['warning_time'] is not None else pd.NaT
    effective_reference=decision['valid_until'] is not None and decision['valid_until']>=event
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
       out/'Q3_weather_reference_time_distribution.csv',index=False)
    score.to_csv(out/'Q3_hourly_precursor_conditions.csv.gz',compression='gzip',index_label='time_bjt')
    _plots(ctx,h,score,defs,fc,leads,fcens,issue,event,operational_alarm,back,metrics)
    lead_alarm=(event-alarm).total_seconds()/3600 if pd.notna(alarm) else None
    rel=(event-issue).total_seconds()/3600
    is_hit=(fc['interval_h'][0]<=rel<=fc['interval_h'][1]) if fc['interval_h'][0] is not None else False
    report={'prediction_issue_bjt':str(issue),'event_time_reference_bjt':str(event),
       'risk_condition_score':defs,'first_actual_threshold_warning_bjt':str(alarm) if pd.notna(alarm) else None,
       'alarm_rule':'yellow uses q75; 2h entry, 3h clearing, 6h renewal, missing -> unknown; only data before issuance',
       'warning_tiers':warn_tiers,
       'first_warning_advance_h_to_reference':float((event-operational_alarm).total_seconds()/3600) if pd.notna(operational_alarm) and effective_reference else None,
       'first_monitoring_warning_advance_h_to_reference':lead_alarm,
       'conditional_95pct_interval_h_from_issue':list(fc['interval_h']),
       'conditional_95pct_interval_times_bjt':[
           str(issue+pd.Timedelta(hours=x)) if x is not None else None for x in fc['interval_h']],
       'median_trigger_h':fc['median_h'],
       'event_reference_in_interval':bool(is_hit),
       'forecast_target':'first passage of the multi-channel meteorological trigger-condition threshold',
       'predicted_condition_cross_probability':fc['p_in_horizon'],
       'predicted_no_condition_cross_probability':fc['nonoccurrence_probability'],
       'forecast_distribution_method':'empirical first passage of data-derived condition trajectories',
       'forecast_input_source':forecast_source,
       'interval_estimand':'first passage time conditional on meteorological condition crossing within the horizon',
       'current_forecast_validation_scope':'archived-weather hybrid case replay; historical proxy metrics belong to the analogue-reference model',
       'forecast_warning_probability_threshold':probability_threshold,
       'forecast_warning_issued':bool(forecast_alert),
       'forecast_warning_time_bjt':str(issue) if forecast_alert else None,
       'forecast_warning_lead_h_to_task_reference':float((event-issue).total_seconds()/3600) if forecast_alert else None,
       'monitoring_episode_start_bjt':str(alarm) if pd.notna(alarm) else None,
       'operational_warning_valid_until_bjt':str(decision['valid_until']) if decision['valid_until'] is not None else None,
       'warning_decision_uses_future_event_time':False,
       'operational_warning_time_bjt':str(operational_alarm) if pd.notna(operational_alarm) else None,
       'operational_warning_lead_h':float((event-operational_alarm).total_seconds()/3600) if pd.notna(operational_alarm) and effective_reference else None,
       'warning_decision_basis':'forecast P(condition passage within 6h) >= predeclared 0.25, or valid monitoring alarm',
       'empirical_disaster_time_PICP':None,
       'single_task_reference_hit':{'n_evaluable_references':1,'n_interval_hits':int(is_hit),
           'reference_kind':'task-reported impact clock, compared with the condition time window'},
       'reference_time_type':'10:30 is the task-reported disaster/impact reference; source trigger truth not independently observed',
       'old_first_alarm_bjt':str(old_alarm) if pd.notna(old_alarm) else None,
       'active_warning_at_issue':bool(last.active),
       'warning_valid_through_task_reference':bool(effective_reference),
       'last_observation_used_bjt':str(issue-pd.Timedelta(hours=1)),
       'warning_rule_full':'enter after 2 above-threshold valid hours; renew for 6h on high scores; clear after 3 low hours or immediately on missing; decisions use previous hourly observation',
       'proxy_backtest':metrics,
       'hazard_base_rate':{'method':'empirical trajectory passage frequencies','assumed_disaster_count':None},
       'analogue_days':n,'external_Q3_download':False,
       'input_reuse':'predeclared Q1 trio + Q2 rain/melt/station-energy condition index; 2026 Top5 used only for supplemental explanation',
       'inherited_Q1_features':inherited,
       'Q1_Q2_ablation_file':'Q3_Q1_Q2_ablation_validation.csv',
       'forecast_method':'complete earlier day score trajectories adjusted to issuance observation; no future data'}
    write_json(out/'Q3_weather_reference_assessment.json',report)
    from .question_completion import complete_q3
    complete_q3(ctx,report)
    from .warning_validation import plot_warning_process
    plot_warning_process(ctx)
    print('[Q3-aligned] 95% time',fc['interval_h'],'forecast warning',operational_alarm,
           'lead_h',report['operational_warning_lead_h'],
           'case_hit',is_hit,'proxy backtest',metrics)
    return report


def mechanical_weather_members(h,issue,leads,station_z,archived=None,model_z=None,past_days=21):
    """Complete earlier-day T/rain trajectories, never score trajectories.

    Archived weather, when present, supplies the centre. Earlier observed-day
    deviations supply a declared weather stress ensemble, NOT calibrated NWP
    errors. Without an archive use adjusted earlier-day trajectories directly.
    Hourly precipitation is piecewise constant over its ending hour; temperature
    is interpolated. A half-hour integration grid creates no extra observations.
    """
    issue=pd.Timestamp(issue);leads=np.asarray(leads,float);last=issue-pd.Timedelta(hours=1)
    if np.any(leads<=0) or np.any(np.diff(leads)<=0):raise ValueError('Invalid weather leads')
    if last not in h.index or not np.isfinite(h.loc[last,'temp_2m']):raise ValueError('Missing current temperature')
    temperatures=[];rainfalls=[];origins=[]
    for lag in range(1,past_days+1):
        anchor=last-pd.Timedelta(days=lag)
        if anchor not in h.index or not np.isfinite(h.loc[anchor,'temp_2m']):continue
        t=[];r=[];valid=True
        for lead in leads:
            dest=issue-pd.Timedelta(days=lag)+pd.Timedelta(hours=float(lead))
            lo=dest.floor('h');hi=dest.ceil('h')
            if hi>=issue or lo not in h.index or hi not in h.index:valid=False;break
            values=[h.loc[lo,'temp_2m'],h.loc[hi,'temp_2m'],h.loc[hi,'precip_mm']]
            if not np.isfinite(values).all():valid=False;break
            weight=(dest-lo).total_seconds()/3600
            t.append((1-weight)*values[0]+weight*values[1]+np.exp(-lead/8)*(h.loc[last,'temp_2m']-h.loc[anchor,'temp_2m']))
            r.append(max(0,values[2]))
        if valid:temperatures.append(t);rainfalls.append(r);origins.append(str(anchor))
    if not temperatures:raise ValueError('No complete earlier observed weather trajectories')
    temp=np.asarray(temperatures);rain=np.asarray(rainfalls);source='pre-issue observed-day weather trajectories'
    ref_z=float(station_z)
    if archived is not None:
        f=archived.copy();index=pd.to_datetime(f.time_bjt,utc=True).dt.tz_convert('Asia/Shanghai').dt.tz_localize(None)
        f.index=index;f=f.sort_index();hour=(f.index-issue).total_seconds()/3600
        good=np.isfinite(f.temperature_2m)&np.isfinite(f.precipitation)
        if good.sum()<2 or hour[good].min()>0 or hour[good].max()<leads.max():raise ValueError('Archived weather does not cover forecast horizon')
        centre_t=np.interp(leads,hour[good],f.loc[good,'temperature_2m'])
        end=issue+pd.to_timedelta(np.ceil(leads),unit='h')
        centre_r=f.precipitation.reindex(end).to_numpy(float)
        if not np.isfinite(centre_r).all():raise ValueError('Missing archived hourly precipitation')
        temp=centre_t[None,:]+temp-np.median(temp,axis=0)[None,:]
        rain=np.maximum(0,centre_r[None,:]+rain-np.median(rain,axis=0)[None,:])
        ref_z=float(model_z);source='Q2 archived ECMWF centre + complete pre-issue observed-day deviations (uncalibrated weather stress ensemble)'
    return temp,rain,ref_z,{'source':source,'n_weather_members':len(temp),'past_day_anchors':origins,
        'temperature_reference_elevation_m':ref_z,'forecast_error_distribution_calibrated':False,
        'precipitation_time_support':'piecewise constant ending-hour rates; T linearly interpolated',
        'future_station_observations_used':False}


def forecast_source_stability(initial_storage,parameters,geometry,temp,rain,leads,reference_z):
    """Apply exactly Q2's water balance and Fs to each parameter/weather pair."""
    from .q2_trigger_evaluation import advance_water_storage,source_safety_factor
    temp=np.asarray(temp,float);rain=np.asarray(rain,float);leads=np.asarray(leads,float)
    if temp.ndim!=2 or rain.shape!=temp.shape or temp.shape[1]!=len(leads):raise ValueError('Weather matrix shape mismatch')
    if not np.isfinite(temp).all() or not np.isfinite(rain).all() or np.any(rain<0):raise ValueError('Invalid forecast forcing')
    if not len(parameters):return np.empty((0,len(leads))),np.array([],int)
    nweather=len(temp);p=parameters.iloc[np.repeat(np.arange(len(parameters)),nweather)].reset_index(drop=True)
    s=np.repeat(np.asarray(initial_storage,float),nweather);result=[]
    capacity=geometry['normal_stress_kPa']/9.81*1000*p.storage_yield.to_numpy()
    if len(s)!=len(p):raise ValueError('Initial state does not match parameter design')
    steps=np.diff(np.r_[0,leads])
    if np.any(steps<=0):raise ValueError('Forecast leads must increase')
    for k,dt in enumerate(steps):
        source_t=np.tile(temp[:,k],len(parameters))-(geometry['source_elevation_m']-reference_z)/1000*p.lapse_C_km.to_numpy()
        liquid=(np.tile(rain[:,k],len(parameters))*(source_t>0)+
                p.degree_day_factor_mm_Cday.to_numpy()*np.maximum(source_t,0)/24)*p.infiltration_fraction.to_numpy()
        s,_,_=advance_water_storage(s,liquid,float(dt),p.drainage_hours.to_numpy(),capacity)
        result.append(source_safety_factor(s,p,geometry))
    return np.asarray(result).T,p.scenario_id.to_numpy()


def first_instability_time(fs,leads):
    a=np.asarray(fs,float);leads=np.asarray(leads,float)
    flag=a<=1;cross=flag.any(axis=1);times=np.full(len(a),np.nan)
    if cross.any():times[cross]=leads[np.argmax(flag[cross],axis=1)]
    return times


def mechanics_forcing_backtest(ctx,h,parameters,geometry):
    """Held-out forcing experiment, with Fs-model outcomes explicitly labelled.

    A reference outcome is the SAME physical model under later *observed*
    weather. Parameter x issue rows are not independent observed disasters.
    No resulting coverage value is exported as empirical disaster coverage.
    """
    from .q2_trigger_evaluation import source_history_states
    cfg=ctx.cfg['q2']['initiation'];q3=ctx.cfg['q3'];out=ctx.folder('q3')
    leads=np.arange(q3['step_hours'],q3['primary_horizon_hours']+1e-8,q3['step_hours'])
    rows=[];outcomes=[];skipped=[]
    for year in q3['mechanics_backtest_years']:
        candidates=pd.date_range(f'{year}-06-01 09:00',f'{year}-09-29 09:00',freq='D')
        chosen=candidates[np.unique(np.round(np.linspace(0,len(candidates)-1,int(q3['mechanics_backtest_issues_per_year']))).astype(int))]
        for issue in chosen:
            try:
                s,minimum,_=source_history_states(h,issue,parameters,geometry,ctx.cfg['station_elevation_m'],cfg['spinup_days'],cfg['consistency_days'])
                accepted=minimum>1;p=parameters.loc[accepted].copy();s=s[accepted]
                if not len(p):raise ValueError('No survival-consistent scenarios')
                t,r,z,_=mechanical_weather_members(h.loc[h.index<issue],issue,leads,ctx.cfg['station_elevation_m'])
                future=issue+pd.to_timedelta(leads,unit='h');lo=future.floor('h');hi=future.ceil('h')
                truth=h.reindex(lo)[['temp_2m']].to_numpy(float).ravel()
                upper=h.reindex(hi)[['temp_2m','precip_mm']].to_numpy(float)
                w=(future-lo).total_seconds().to_numpy()/3600
                true_t=(1-w)*truth+w*upper[:,0];true_r=upper[:,1]
                if not np.isfinite(true_t+true_r).all():raise ValueError('Missing held-out weather')
                fs,_=forecast_source_stability(s,p,geometry,t,r,leads,z)
                actual,_=forecast_source_stability(s,p,geometry,true_t[None,:],true_r[None,:],leads,ctx.cfg['station_elevation_m'])
                forecast=empirical_first_passage(-fs,float(leads[0]),-1)
                times=first_instability_time(actual,leads);cross=np.isfinite(times)
                low,high=forecast['interval_h'];covered=cross&(times>=low)&(times<=high) if low is not None else np.zeros(len(times),bool)
                rows.append({'issue_bjt':str(issue),'year':year,'n_surviving_parameter_scenarios':len(p),
                    'n_model_reference_crossings':int(cross.sum()),'n_model_reference_times_covered':int(covered.sum()),
                    'forecast_model_cross_fraction':forecast['p_in_horizon'],
                    'realized_weather_model_cross_fraction':float(cross.mean()),'interval_lower_h':low,'interval_upper_h':high,
                    'interval_width_h':high-low if low is not None else np.nan,
                    'target':'Fs-model failure under held-out observed weather; not observed disaster'})
                outcomes.extend({'issue_bjt':str(issue),'scenario_id':int(i),'realized_model_first_failure_h':float(v) if np.isfinite(v) else np.nan,
                    'covered_by_conditional_interval':bool(c)} for i,v,c in zip(p.scenario_id,times,covered))
            except ValueError as exc:skipped.append({'issue_bjt':str(issue),'reason':str(exc)})
    frame=pd.DataFrame(rows);frame.to_csv(out/'Q3_mechanics_forcing_backtest.csv',index=False)
    pd.DataFrame(outcomes).to_csv(out/'Q3_mechanics_reference_outcomes.csv.gz',index=False,compression='gzip')
    total=int(frame.n_model_reference_crossings.sum()) if len(frame) else 0
    hits=int(frame.n_model_reference_times_covered.sum()) if len(frame) else 0
    report={'experiment':'held-out observed weather versus pre-issue observed-day forecast, with a fixed declared parameter design',
        'n_issues':len(frame),'years':q3['mechanics_backtest_years'],'skipped_issues':skipped,
        'n_model_reference_crossings':total,'n_model_reference_times_covered':hits,
        'model_reference_conditional_PICP':hits/total if total else None,
        'mean_finite_interval_width_h':float(frame.interval_width_h.mean()) if len(frame) and frame.interval_width_h.notna().any() else None,
        'empirical_disaster_PICP':None,'n_independently_observed_source_onsets':0,
        'sampling_unit':'parameter-scenario x issue; correlated numerical outcomes, not independent disasters',
        'current_archived_hybrid_model_validated_by_this_experiment':False,
        'calibration':'none; nominal 95% quantiles retained even if forcing experiment undercovers'}
    write_json(out/'Q3_mechanics_validation.json',report)
    fig,axs=plt.subplots(1,2,figsize=(11.6,4.5),layout='constrained')
    if len(frame):
        dates=pd.to_datetime(frame.issue_bjt)
        axs[0].plot(dates,frame.forecast_model_cross_fraction,'o-',ms=3,label='相似日预报驱动')
        axs[0].plot(dates,frame.realized_weather_model_cross_fraction,'s',ms=3,label='后续实测天气驱动同一模型')
        axs[0].set(xlabel='历史发布日',ylabel='模型失稳情景比例',ylim=(0,1));axs[0].legend(frameon=False,fontsize=8)
        for year,color in zip(q3['mechanics_backtest_years'],COLORS):
            v=frame.loc[frame.year.eq(year)];den=v.n_model_reference_crossings.cumsum();num=v.n_model_reference_times_covered.cumsum()
            axs[1].plot(pd.to_datetime(v.issue_bjt),num/den.replace(0,np.nan),label=str(year),color=color)
    axs[1].axhline(.95,color='#777',ls='--');axs[1].set(xlabel='历史发布日',ylabel='模型参照时间的累计覆盖率',ylim=(0,1.05));axs[1].legend(frameon=False,fontsize=8)
    fig.suptitle('固定力学模型的天气驱动检验；分母为数值情景，非真实灾害')
    savefig(fig,out/'Q3_13_mechanical_forcing_validation.png')
    return report


def conditional_trigger_time(ensemble,step_h):
    """Specified conditional time model: density proportional to Q2 intensity.

    CONDITION: a single event occurs in the stated horizon and fixed source.
    Each observed-day scenario defines a piecewise-uniform time density. This
    predicts event time under an explicit modelling assumption; it does not
    identify P(event) or equate an observed weather crossing with a disaster.
    All-zero trajectories supply the uniform conditional prior, with count logged.
    """
    a=np.asarray(ensemble,float)
    if a.ndim!=2 or min(a.shape)<1 or not np.isfinite(a).all() or np.any(a<0) or step_h<=0:
        raise ValueError('Nonnegative finite Q2 intensity trajectories and positive step required')
    totals=a.sum(axis=1);weights=np.divide(a,totals[:,None],out=np.full_like(a,1/a.shape[1]),where=totals[:,None]>0)
    mass=weights.mean(axis=0);mass/=mass.sum();cdf=np.cumsum(mass);cdf[-1]=1.
    def quantile(q):
        i=min(int(np.searchsorted(cdf,q)),len(mass)-1);previous=0. if i==0 else cdf[i-1]
        return float((i+(q-previous)/mass[i])*step_h)
    return {'conditional_bin_mass':mass,'interval_h':(quantile(.025),quantile(.975)),
        'median_h':quantile(.5),'mean_h':float(np.dot((np.arange(len(mass))+.5)*step_h,mass)),
        'n_uniform_prior_trajectories':int((totals<=0).sum()),'event_occurrence_probability':None,
        'time_target':'event time conditional on occurrence within horizon and fixed source',
        'within_bin_model':'uniform density; numerical interpolation, not extra time observations'}


def single_case_coverage(interval_h,reference_h):
    from scipy.stats import beta
    hit=bool(interval_h[0]<=reference_h<=interval_h[1]);x=int(hit)
    return {'n_references':1,'n_covered':x,'coverage':float(x),
        'binomial_95_interval':[float(beta.ppf(.025,1,1)) if x else 0.,1. if x else float(beta.ppf(.975,1,1))],
        'scope':'one competition-reported reference under the adopted event clock, not multi-disaster generalization'}


def execute_q3_aligned(ctx):
    """Q1 -> Q2 -> conditional event time; Q4 results are not inputs here."""
    import json,hashlib
    out=ctx.folder('q3');q3=ctx.cfg['q3'];issue=pd.Timestamp(ctx.cfg['snapshot_raw'])
    for pattern in ['Q3_mechanics_*','Q3_mechanical_*','Q3_source_failure_*','Q3_forecast_source_stability.*','Q3_forecast_weather_members.*','Q3_13_mechanical_*']:
        for old in out.glob(pattern):
            if old.is_file():old.unlink()
    contract=json.loads((ctx.folder('q2')/'Q2_sequential_model_contract.json').read_text())
    if contract['Q1_answer_sha256']!=hashlib.sha256((ctx.folder('q1')/'Q1_answer_mechanism_top5.csv').read_bytes()).hexdigest():
        raise ValueError('Q1 features changed; rerun Q2 before Q3')
    if contract['assessment_issue_bjt']!=str(issue):raise ValueError('Issuance changed; rerun Q1/Q2 before Q3')
    if ctx.cfg.get('reference_experiments',False):execute_q3_weather_reference(ctx)
    series=pd.read_csv(ctx.folder('q2')/'Q2_sequential_assessment_timeseries.csv.gz',parse_dates=['time_bjt']).set_index('time_bjt')
    q1_feature_contract(ctx,pd.read_csv(ctx.processed/'features_hourly.csv.gz',nrows=0).columns,contract['Q1_features'])
    horizon=float(q3['primary_horizon_hours']);step=float(q3['step_hours']);n=int(q3['samples']);seed=int(q3['random_seed'])
    history=series.loc[series.index<issue];score=history.assessment_score
    ensemble,leads,days=forecast_analogue(score,issue,horizon,step,n,seed)
    forcing='complete pre-issue Q2 assessment trajectories'
    archive=ctx.folder('q2')/'Q2_archived_forecast_trigger_conditions.csv'
    if archive.exists():
        f=pd.read_csv(archive,parse_dates=['forecast_time_bjt']).set_index('forecast_time_bjt')
        hours=(f.index-issue).total_seconds()/3600;good=f.trigger_condition_forecast.notna()&(hours>=0)
        if good.sum()>=2 and hours[good].min()<=0 and hours[good].max()>=horizon:
            anom,_,_=forecast_analogue(history.Q1_selected_feature_anomaly,issue,horizon,step,n,seed)
            water=np.interp(leads,hours[good],f.loc[good,'trigger_condition_forecast'])
            floor=contract['q1_modulation_floor'];ensemble=water[None,:]*(floor+(1-floor)*anom)
            forcing='Q1-selected anomaly trajectories + Q2 actual archived weather water/thermal gate, with unchanged Q2 equation'
    forecast=conditional_trigger_time(ensemble,step)
    np.savez_compressed(out/'Q3_assessment_ensemble.npz',lead_h=leads,assessment_score=ensemble)
    crossing=float((ensemble>=contract['critical_score']).any(axis=1).mean())
    from .warning_validation import warning_state_machine,decide_warning_at_issue
    states,episodes=warning_state_machine(score.loc[score.index>=issue-pd.Timedelta(days=7)],contract['monitoring_score'])
    states.to_csv(out/'Q3_sequential_warning_states.csv',index_label='decision_time_bjt');episodes.to_csv(out/'Q3_sequential_warning_episodes.csv',index=False)
    last=states.iloc[-1];active=bool(last.active) and pd.notna(last.expiry_bjt) and pd.Timestamp(last.expiry_bjt)>=issue
    start=pd.Timestamp(last.episode_start_bjt) if active else None
    decision=decide_warning_at_issue(issue,crossing,float(q3['warning_cross_probability_threshold']),active,start,last.expiry_bjt,horizon)
    # Only evaluation below this point reads the given event reference.
    reference=pd.Timestamp(q3['official_event_time_same_clock_assumption'])
    relative=(reference-issue).total_seconds()/3600;coverage=single_case_coverage(forecast['interval_h'],relative)
    low,high=forecast['interval_h'];interval=[str((issue+pd.Timedelta(hours=low)).floor('min')),str((issue+pd.Timedelta(hours=high)).ceil('min'))]
    warning=decision['forecast_warning_time'] or decision['warning_time']
    valid=decision['valid_until'] is not None and decision['valid_until']>=reference
    warning_lead=(reference-warning).total_seconds()/3600 if warning is not None and valid else None
    pd.DataFrame({'bin_end_hours_from_issue':leads,'conditional_time_probability':forecast['conditional_bin_mass']}).to_csv(out/'Q3_source_initiation_time_distribution.csv',index=False)
    pd.DataFrame({'lead_h':leads,'assessment_q025':np.quantile(ensemble,.025,axis=0),
        'assessment_median':np.median(ensemble,axis=0),'assessment_q975':np.quantile(ensemble,.975,axis=0)}).to_csv(out/'Q3_forecast_assessment.csv',index=False)
    comparisons=[]
    for name,array in [('Q2_intensity_proportional',ensemble),('uniform_time_prior',np.ones_like(ensemble)),('squared_intensity_sensitivity',ensemble**2)]:
        dist=conditional_trigger_time(array,step)
        comparisons.append({'method':name,'lower_h':dist['interval_h'][0],'upper_h':dist['interval_h'][1],
            'width_h':dist['interval_h'][1]-dist['interval_h'][0],'median_h':dist['median_h'],
            'case_reference_covered':single_case_coverage(dist['interval_h'],relative)['n_covered']})
    pd.DataFrame(comparisons).to_csv(out/'Q3_time_model_comparison.csv',index=False)
    geometry=json.loads((ctx.folder('q2')/contract['source_geometry_file']).read_text())
    report={'model_family':'Q2_conditioned_event_time','prediction_issue_bjt':str(issue),
        'forecast_target':'mudslide event time conditional on the identified case source and occurrence within the forecast horizon',
        'criterion':'Q2 statistical trigger assessment; conditional time density proportional to assessment intensity',
        'conditional_95pct_interval_h_from_issue':[low,high],'conditional_95pct_interval_times_bjt':interval,
        'median_trigger_h':forecast['median_h'],'median_event_time_bjt':str(issue+pd.Timedelta(hours=forecast['median_h'])),
        'interval_estimand':forecast['time_target'],'within_bin_time_assumption':forecast['within_bin_model'],
        'event_occurrence_probability':None,'conditional_on_occurrence':True,
        'assessment_cross_probability':crossing,'assessment_cross_probability_is_disaster_probability':False,
        'forecast_warning_probability_threshold':float(q3['warning_cross_probability_threshold']),
        'forecast_warning_issued':bool(decision['forecast_issued']),
        'forecast_warning_time_bjt':str(decision['forecast_warning_time']) if decision['forecast_warning_time'] is not None else None,
        'operational_warning_time_bjt':str(warning) if warning is not None else None,
        'operational_warning_lead_h':warning_lead,'operational_warning_valid_until_bjt':str(decision['valid_until']) if decision['valid_until'] is not None else None,
        'monitoring_episode_start_bjt':str(start) if start is not None else None,'active_warning_at_issue':active,
        'warning_decision_uses_future_event_time':False,'last_observation_used_bjt':str(history.index[-1]),
        'event_time_used_to_fit_parameters':False,'case_reference_coverage':coverage,
        'case_reference_median_error_h':float(forecast['median_h']-relative),
        'case_reference_interval_width_h':float(high-low),
        'single_task_reference_hit':coverage,'empirical_disaster_time_PICP':None,
        'reference_time_type':'competition event clock adopted as trigger-time benchmark; minute-scale source/impact distinctions are simplified',
        'event_time_reference_bjt':str(reference),'inherited_Q1_features':contract['Q1_features'],
        'input_contract':'q2/Q2_sequential_model_contract.json','source_geometry_sha256':contract['source_geometry_sha256'],
        'source_to_port_km':contract['source_to_port_km'],'spatial_distance_range_km':[0.,contract['source_to_port_km']],
        'spatial_scope':'Q2 static source-to-port DEM corridor, determined before Q4 simulation; not full 2D inundation',
        'forecast_input_source':forcing,'analogue_days':days,'external_Q3_download':False,
        'n_uniform_prior_trajectories':forecast['n_uniform_prior_trajectories'],
        'hazard_base_rate':{'assumed_disaster_count':None,'method':'conditional intensity normalization; no invented baseline event count'},
        'model_assumptions':['fixed case source and available material','one event in the specified horizon when reporting conditional time','time density proportional to Q2 assessment intensity'],
        'independent_multi_event_validation':False,'Q4_results_used_as_Q3_inputs':False}
    write_json(out/'Q3_primary_competition_answer.json',report)
    _plot_sequential_answer(ctx,series,ensemble,leads,forecast,report,comparisons)
    print('[Q3 sequential]',interval,'warning',report['operational_warning_time_bjt'],'case coverage',coverage['n_covered'],'/1','range',report['spatial_distance_range_km'])
    return report


def _plot_sequential_answer(ctx,series,ensemble,leads,forecast,report,comparisons):
    import json
    out=ctx.folder('q3');issue=pd.Timestamp(report['prediction_issue_bjt']);reference=pd.Timestamp(report['event_time_reference_bjt'])
    contract=json.loads((ctx.folder('q2')/'Q2_sequential_model_contract.json').read_text())
    past=series.loc[(series.index>=issue-pd.Timedelta(days=7))&(series.index<issue)]
    fig,ax=plt.subplots(figsize=(10.8,4.3),layout='constrained')
    ax.plot(past.index,past.assessment_score,label='Q2综合触发评估',color=COLORS[0])
    dates=issue+pd.to_timedelta(leads,unit='h');ax.plot(dates,np.median(ensemble,axis=0),color=COLORS[1],label='预报评估中位数')
    ax.fill_between(dates,np.quantile(ensemble,.025,axis=0),np.quantile(ensemble,.975,axis=0),color=COLORS[1],alpha=.2)
    ax.axhline(contract['critical_score'],ls='--',color='#555',label='Q2统计临界水平');ax.axvline(issue,ls=':',color='#777')
    ax.set(xlabel='北京时间',ylabel='综合评估指数');ax.legend(frameon=False,fontsize=8)
    savefig(fig,out/'Q3_01_precursor_and_scenario_hazard.png')
    fig,axs=plt.subplots(1,2,figsize=(11.8,4.4),layout='constrained');low,high=forecast['interval_h']
    axs[0].bar(leads-leads[0]/2,forecast['conditional_bin_mass'],width=leads[0]*.9,color=COLORS[0])
    axs[0].axvspan(low,high,color=COLORS[1],alpha=.14);axs[0].set(xlabel='距发布时刻 (h)',ylabel='条件时间概率质量',title='a  发生条件下的时间分布')
    axs[1].plot(np.r_[0,leads],np.r_[0,np.cumsum(forecast['conditional_bin_mass'])],color=COLORS[0])
    axs[1].axvline((reference-issue).total_seconds()/3600,color=COLORS[1],ls='--',label='题面时刻；仅作评价')
    axs[1].axvspan(low,high,color=COLORS[0],alpha=.12);axs[1].set(xlabel='距发布时刻 (h)',ylabel='累计条件概率',ylim=(0,1));axs[1].legend(frameon=False,fontsize=8)
    savefig(fig,out/'Q3_02_95percent_conditional_time_interval.png')
    fig,ax=plt.subplots(figsize=(10.4,4.3),layout='constrained')
    labels=['Q2评估强度主模型','均匀时间先验','强度平方敏感性']
    for i,(row,label) in enumerate(zip(comparisons,labels)):
        ax.plot([row['lower_h'],row['upper_h']],[i,i],lw=9,color=COLORS[i],alpha=.6)
        ax.plot(row['median_h'],i,'o',color=COLORS[i]);ax.text(row['upper_h']+.06,i,f"宽度{row['width_h']:.2f}h",va='center',fontsize=8)
    ax.axvline((reference-issue).total_seconds()/3600,color='#555',ls='--',label='题面时刻')
    ax.set(yticks=range(3),yticklabels=labels,xlabel='距发布时刻 (h)',xlim=(0,leads[-1]+.9));ax.legend(frameon=False,fontsize=8)
    savefig(fig,out/'Q3_03_warning_lead_and_backtest.png')
    p=pd.read_csv(ctx.folder('q2')/'Q2_reconstruction_route.csv');xy=p[['utm_x','utm_y']].to_numpy()
    from .q4_flowpath import local_dem
    from matplotlib.colors import LightSource
    z,tr,_=local_dem(ctx,xy[::max(1,len(xy)//50)],buffer_m=900,resolution_m=120)
    z=np.where(np.isfinite(z),z,np.nanmedian(z));hill=LightSource(315,45).hillshade(z,dx=120,dy=120)
    fig,ax=plt.subplots(figsize=(10,5.5),layout='constrained');ax.imshow(hill,cmap='Greys',extent=[tr.c,tr.c+z.shape[1]*tr.a,tr.f+z.shape[0]*tr.e,tr.f])
    ax.plot(xy[:,0],xy[:,1],color=COLORS[0],lw=1.8,label='Q2源区至口岸静态通道');ax.scatter(xy[[0,-1],0],xy[[0,-1],1],c=COLORS[:2],s=45)
    ax.set(xlabel='UTM东向坐标 (m)',ylabel='UTM北向坐标 (m)');ax.legend(frameon=False,fontsize=8)
    ax.text(.02,.98,f"空间距离范围 0-{contract['source_to_port_km']:.2f} km\n由源区和DEM确定；Q4随后模拟传播",transform=ax.transAxes,va='top',bbox=dict(fc='white',alpha=.9,ec='#ccc'))
    savefig(fig,out/'Q3_05_spatial_warning_distance.png')
