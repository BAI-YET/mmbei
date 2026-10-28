"""Forward-only meteorological proxy validation and auditable warning episodes.
No disaster labels are manufactured from weather thresholds.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import brier_score_loss,roc_auc_score
from .common import write_json,savefig,COLORS


def decide_warning_at_issue(issue,probability,threshold,monitoring_active=False,
                            episode_start=None,expiry=None,horizon_h=6):
    """An issue-time decision has no observed future event-time argument."""
    issue=pd.Timestamp(issue);forecast=bool(probability>=threshold)
    monitored=bool(monitoring_active and pd.notna(episode_start) and pd.notna(expiry)
                   and pd.Timestamp(expiry)>=issue)
    start=pd.Timestamp(episode_start) if monitored else None
    ends=[]
    if monitored:ends.append(pd.Timestamp(expiry))
    if forecast:ends.append(issue+pd.Timedelta(hours=horizon_h))
    warning_time=min(issue,start) if forecast and monitored else (issue if forecast else start)
    return {'forecast_issued':forecast,'monitoring_active':monitored,
            'warning_time':warning_time,'valid_until':max(ends) if ends else None,
            'forecast_warning_time':issue if forecast else None}


def warning_state_machine(score,threshold,enter_hours=2,clear_hours=3,ttl_hours=6):
    """Observation stamped t can first inform decision at t+1h.
    Invalid or nonconsecutive observations reset the state to UNKNOWN.
    Expiry is a forecast validity limit, not an assertion that risk persists.
    """
    score=score.sort_index();rows=[];episodes=[]
    active=False;high=low=0;start=None;expiry=None;last_obs=None
    for t,v in score.items():
        decision=t+pd.Timedelta(hours=1)
        gap=last_obs is not None and t-last_obs!=pd.Timedelta(hours=1)
        last_obs=t
        reason=None
        if gap or not np.isfinite(v):
            reason='missing_or_nonconsecutive_observation';high=low=0;state='UNKNOWN'
        else:
            if v>=threshold:high+=1;low=0
            else:low+=1;high=0
            if active and ((expiry is not None and decision>expiry) or low>=clear_hours):reason='expiry' if decision>expiry else 'three_low_hours'
            state='NORMAL'
        if reason and active:
            episodes.append({'start_bjt':start,'end_bjt':decision,'expiry_bjt':expiry,'end_reason':reason})
            active=False;start=None;expiry=None
        if reason is None:
            if not active and high>=enter_hours:active=True;start=decision
            if active and v>=threshold:expiry=decision+pd.Timedelta(hours=ttl_hours)
            state='ALERT' if active else ('WATCH' if high else 'NORMAL')
        rows.append({'decision_time_bjt':decision,'observation_time_bjt':t,'condition_score':v,'state':state,
           'active':active,'episode_start_bjt':start,'expiry_bjt':expiry,'consecutive_high_hours':high,'consecutive_low_hours':low})
    if active:episodes.append({'start_bjt':start,'end_bjt':None,'expiry_bjt':expiry,'end_reason':'open_at_issue'})
    return pd.DataFrame(rows).set_index('decision_time_bjt'),pd.DataFrame(episodes,columns=['start_bjt','end_bjt','expiry_bjt','end_reason'])


def _block_coverage_ci(df,lower='lower_h',upper='upper_h',seed=7):
    event=df.loc[df.true_crossed].copy()
    if event.empty:return None
    dates=pd.to_datetime(event.issue_bjt)
    event['block']=(dates-dates.min()).dt.days//7
    groups=[g for _,g in event.groupby('block')];rng=np.random.default_rng(seed);vals=[]
    for _ in range(800):
        b=pd.concat([groups[i] for i in rng.integers(0,len(groups),size=len(groups))])
        vals.append(float(((b.first_true_cross_h>=b[lower])&(b.first_true_cross_h<=b[upper])).mean()))
    return np.quantile(vals,[.025,.975]).tolist()


def backtest_metrics(df,lower='lower_h',upper='upper_h',threshold=.25):
    if df.empty:return {'n_backtest':0}
    flag=df.forecast_cross_probability.ge(threshold);truth=df.true_crossed.astype(bool)
    TP=int((flag&truth).sum());FP=int((flag&~truth).sum());FN=int((~flag&truth).sum());TN=int((~flag&~truth).sum())
    targets=df.loc[truth]
    hit=(targets.first_true_cross_h>=targets[lower])&(targets.first_true_cross_h<=targets[upper])
    onset=targets.loc[~targets.already_above_at_issue]
    ohit=(onset.first_true_cross_h>=onset[lower])&(onset.first_true_cross_h<=onset[upper])
    def div(a,b):return a/b if b else None
    timely=df.loc[flag&truth&~df.already_above_at_issue.astype(bool),'first_true_cross_h'].dropna()
    lead_summary=({'n':int(len(timely)),'minimum_h':float(timely.min()),'q25_h':float(timely.quantile(.25)),
        'median_h':float(timely.median()),'q75_h':float(timely.quantile(.75)),
        'maximum_h':float(timely.max())} if len(timely) else {'n':0})
    return {'n_backtest':len(df),'n_proxy_crossings':int(truth.sum()),'TP':TP,'FP':FP,'FN':FN,'TN':TN,
       'new_onset_warning_lead_distribution':lead_summary,
       'proxy_recall':div(TP,TP+FN),'proxy_false_alarm':div(FP,FP+TN),
       'proxy_false_alarm_ratio':div(FP,TP+FP),'proxy_precision':div(TP,TP+FP),
       'proxy_interval_coverage':float(hit.mean()) if len(targets) else None,
       'n_proxy_intervals_evaluable':len(targets),'n_new_onsets':len(onset),
       'new_onset_interval_coverage':float(ohit.mean()) if len(onset) else None,
       'proxy_coverage_7day_block_bootstrap_95':_block_coverage_ci(df,lower,upper),
       'mean_interval_width_h':float((df[upper]-df[lower]).mean()),
       'full_horizon_interval_width_h':6.,'full_horizon_interval_coverage':1. if len(targets) else None,
       'mean_width_reduction_vs_full_horizon':float(1-(df[upper]-df[lower]).mean()/6.),
       'brier_score':float(brier_score_loss(truth,df.forecast_cross_probability)),
       'persistence_brier':float(brier_score_loss(truth,df.already_above_at_issue.astype(float))),
       'proxy_roc_auc':float(roc_auc_score(truth,df.forecast_cross_probability)) if truth.nunique()==2 else None,
       'meaning':'out-of-time weather-condition exceedance; dates have dependence; NOT debris-flow occurrence',
       'FPR_definition':'FP/(FP+TN)','false_alarm_ratio_definition':'FP/(TP+FP)'}


def proxy_backtest(h,train_years=(2021,2022,2023),test_year=2025,max_days=153,variables=None,mode='q1q2'):
    from .q3_task_aligned import anomaly_score,forecast_analogue,BASELINE_VARS
    variables=variables or BASELINE_VARS
    score,defs=anomaly_score(h,train_years,variables=variables,mode=mode)
    x=score.combined_condition_score;thr=defs['threshold'];rows=[];excluded=[]
    dates=pd.date_range(f'{test_year}-05-01',f'{test_year}-09-30',freq='D')[:max_days]
    for i,date in enumerate(dates):
        issue=date+pd.Timedelta(hours=9);last=issue-pd.Timedelta(hours=1)
        if not np.isfinite(x.get(last,np.nan)):
            excluded.append({'issue_bjt':str(issue),'reason':'missing_pre_issue_score'});continue
        try:ensemble,leads,n=forecast_analogue(x.loc[x.index<issue],issue,n_samples=300,seed=10000+i)
        except (ValueError,KeyError):
            excluded.append({'issue_bjt':str(issue),'reason':'no_complete_earlier_day_analogues'});continue
        # Interpolate only one half-hour inside TWO valid adjacent hourly observations.
        times=pd.DatetimeIndex([issue+pd.Timedelta(hours=float(t)) for t in leads])
        union=x.index.union(times)
        obs=x.reindex(union).interpolate('time',limit=1,limit_area='inside').reindex(times).to_numpy(float)
        exact=x.reindex(pd.date_range(issue,issue+pd.Timedelta(hours=6),freq='h'))
        if not np.isfinite(obs).all() or not exact.notna().all():
            excluded.append({'issue_bjt':str(issue),'reason':'incomplete_verification_horizon'});continue
        cross=(ensemble>=thr).any(axis=1);pcross=float(cross.mean());interval=(0.,6.)
        if cross.any():
            bins=np.argmax(ensemble[cross]>=thr,axis=1)
            freq=np.bincount(bins,minlength=len(leads))/len(bins);cdf=np.cumsum(freq)
            interval=(float(np.searchsorted(cdf,.025)*.5),float((np.searchsorted(cdf,.975)+1)*.5))
        truth=bool((obs>=thr).any());ttrue=float(leads[np.flatnonzero(obs>=thr)[0]]) if truth else np.nan
        rows.append({'issue_bjt':str(issue),'forecast_cross_probability':pcross,'true_crossed':truth,
           'already_above_at_issue':bool(x[last]>=thr),'first_true_cross_h':ttrue,'lower_h':interval[0],'upper_h':interval[1],
           'covered':bool(interval[0]<=ttrue<=interval[1]) if truth else np.nan,'n_complete_analogues':n,
           'train_years':','.join(map(str,train_years)),'Q2_threshold_max_training_year':max(train_years),
           'features':'|'.join(variables),'model':mode})
    df=pd.DataFrame(rows)
    metrics=backtest_metrics(df);metrics['excluded_issues']=excluded;metrics['n_scheduled_issues']=len(dates)
    metrics['train_years']=list(train_years);metrics['test_year']=test_year
    return df,metrics


def conformal_expansion(calibration,alpha=.05):
    d=calibration.loc[calibration.true_crossed]
    errors=np.maximum.reduce([d.lower_h.to_numpy()-d.first_true_cross_h.to_numpy(),
                            d.first_true_cross_h.to_numpy()-d.upper_h.to_numpy(),np.zeros(len(d))])
    n=len(errors);k=int(np.ceil((n+1)*(1-alpha)))
    if not n or k>n:return 6.,n
    return float(np.partition(errors,k-1)[k-1]),n


def complete_q3_validation(ctx,h):
    out=ctx.folder('q3');frames=[];comparison=[];selected=None;selected_metrics=None
    for mode in ['q1_only','q2_only','q1q2']:
        cal,cm=proxy_backtest(h,(2021,2022,2023),2024,mode=mode)
        test,tm=proxy_backtest(h,(2021,2022,2023),2025,mode=mode)
        delta,ncal=conformal_expansion(cal)
        test['calibrated_lower_h']=np.maximum(0,test.lower_h-delta)
        test['calibrated_upper_h']=np.minimum(6,test.upper_h+delta)
        test['calibration_year']=2024;test['interval_expansion_h']=delta
        met=backtest_metrics(test,'calibrated_lower_h','calibrated_upper_h')
        met.update({'train_years':[2021,2022,2023],'calibration_year':2024,'test_year':2025,
          'excluded_issues':tm['excluded_issues'],'n_scheduled_issues':tm['n_scheduled_issues'],
          'raw_proxy_interval_coverage':tm['proxy_interval_coverage'],'calibration_proxy_event_count':ncal,
          'interval_expansion_h':delta,'nominal_target':.95,
          'calibration_caution':'conformal quantile is diagnostic under temporal dependence; no guaranteed disaster coverage'})
        for split,frame in [('calibration',cal),('test',test)]:
            frame['validation_split']=split;frames.append(frame)
        comparison.append({'model':mode,**{k:v for k,v in met.items() if isinstance(v,(int,float,type(None)))}})
        if mode=='q1q2':selected=test;selected_metrics=met
    allrows=pd.concat(frames,ignore_index=True)
    allrows.to_csv(out/'Q3_all_models_forward_backtest.csv',index=False)
    pd.DataFrame(comparison).to_csv(out/'Q3_all_models_validation_comparison.csv',index=False)
    write_json(out/'Q3_forward_validation_protocol.json',{'feature_selection':'fixed before evaluation: temp_2m, sw_down, vpd_2m',
       'fit_years':[2021,2022,2023],'calibration_year':2024,'test_year':2025,
       'Q2_threshold_refitted_on_fit_years':True,'future_observations_only_for_scoring':True,
       '2026_Q1_top5_used_for_validation':False,'not_an_independent_disaster_validation':True,
       'selected_metrics':selected_metrics})
    fig,axs=plt.subplots(3,1,figsize=(11.8,8.5),sharex=True,layout='constrained')
    dates=pd.to_datetime(selected.issue_bjt)
    axs[0].plot(dates,selected.forecast_cross_probability,color=COLORS[0],label='6h内气象超限概率')
    axs[0].scatter(dates,selected.true_crossed.astype(float),s=14,color=COLORS[1],label='实际是否超限',alpha=.65)
    axs[0].axhline(.25,color='#777',ls='--',lw=1);axs[0].set(ylabel='超限概率 / 状态',ylim=(-.05,1.12));axs[0].legend(frameon=False,ncol=2,fontsize=8)
    axs[1].fill_between(dates,selected.calibrated_lower_h,selected.calibrated_upper_h,color=COLORS[2],alpha=.28,label='名义95%气象时间区间')
    ev=selected.true_crossed;axs[1].scatter(dates[ev],selected.loc[ev,'first_true_cross_h'],s=16,color=COLORS[1],label='实际首个超限时刻')
    axs[1].set(ylabel='距发布时刻 (h)',ylim=(0,6.3));axs[1].legend(frameon=False,ncol=2,fontsize=8)
    for mode,col in zip(['q1_only','q2_only','q1q2'],COLORS):
        d=allrows.loc[(allrows.model==mode)&(allrows.validation_split=='test')].copy()
        c=(d.first_true_cross_h>=d.calibrated_lower_h)&(d.first_true_cross_h<=d.calibrated_upper_h)
        mask=d.true_crossed
        axs[2].plot(pd.to_datetime(d.loc[mask,'issue_bjt']),c[mask].expanding().mean(),label=mode,color=col,lw=1.4)
    axs[2].axhline(.95,color='#777',ls='--',lw=1);axs[2].set(ylabel='累计气象区间覆盖率',ylim=(0,1.05),xlabel='2025年完整独立检验期（北京时间）');axs[2].legend(frameon=False,ncol=3,fontsize=8)
    fig.suptitle('完整历史检验期：预警、时间区间及三组模型比较')
    savefig(fig,out/'Q3_08_full_forward_validation.png')
    return selected,selected_metrics


def plot_warning_process(ctx):
    out=ctx.folder('q3');s=pd.read_csv(out/'Q3_warning_state_hourly.csv',parse_dates=['decision_time_bjt'])
    info=__import__('json').loads((out/'Q3_primary_competition_answer.json').read_text())
    fig,axs=plt.subplots(2,1,figsize=(11.6,6),sharex=True,layout='constrained')
    t=s.decision_time_bjt;th=info['warning_tiers']['yellow']['threshold']
    axs[0].plot(t,s.condition_score,color=COLORS[0]);axs[0].axhline(th,color=COLORS[1],ls='--',label='黄色条件阈值')
    axs[0].set(ylabel='先兆条件指数');axs[0].legend(frameon=False)
    axs[1].step(t,s.active.astype(int),where='post',color=COLORS[2],lw=1.8)
    axs[1].fill_between(t,0,s.active.astype(int),step='post',color=COLORS[2],alpha=.2)
    axs[1].set(yticks=[0,1],yticklabels=['已解除 / 未激活','有效预警'],xlabel='决策时刻（使用前一小时观测）',ylabel='预警状态',ylim=(-.15,1.2))
    for ax in axs:ax.axvline(pd.Timestamp(info['prediction_issue_bjt']),color='#555',ls=':',lw=1);ax.grid(axis='y',alpha=.12)
    fig.suptitle('预警全过程：两小时进入、三小时解除、六小时有效期')
    savefig(fig,out/'Q3_09_warning_episodes.png')
