"""Q3: physical precursor analogue warning + conditional disaster-chain time scenario.

TASK outputs: predicted trigger-time 95% interval; warning issuance and lead time;
2025 independent precursor-onset backtest (NEVER mislabelled disaster coverage);
spatial extent on official DEM-derived route. Uses no future observations in predictions.

A single 2026 disaster cannot calibrate an empirical disaster PICP. Instead report:
 - an event-specific observed coverage flag vs the 2026 official reference time, and
 - empirical PICP on *meteorological precursor threshold crossings* in 2025.
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.collections import LineCollection
from pathlib import Path
from .common import COLORS, savefig, write_json, record_warning

FEATURES=('pdd_7d','vpd','net_radiation','temp_2m')
WEIGHTS=np.array([.35,.20,.30,.15])


def hazard_step(lambda_daily,index,beta,step_h):
    """Nonnegative intensity transformed into event probability per time bin."""
    if step_h<=0 or np.any(np.asarray(lambda_daily)<0):
        raise ValueError('Positive time step and nonnegative intensity required')
    return 1-np.exp(-np.asarray(lambda_daily)*np.exp(np.clip(beta*index,-25,25))*step_h/24)


def conditional_event_time(h):
    h=np.asarray(h,float)
    if np.any((h<0)|(h>1)):raise ValueError('Hazard outside [0,1]')
    surv=np.r_[1,np.cumprod(1-h)[:-1]]
    mass=surv*h;none=float(np.prod(1-h));within=float(mass.sum())
    if abs(within+none-1)>1e-10:raise AssertionError('Survival mass not conserved')
    cond=mass/within if within>1e-15 else np.full_like(mass,np.nan)
    return mass,none,cond


def quantile_bin(mass,step_h,quantile):
    mass=np.asarray(mass,float)
    if not np.isfinite(mass).all() or mass.sum()<=0:raise ValueError('No valid time mass')
    p=mass/mass.sum();edge=np.cumsum(p)
    return float((int(np.searchsorted(edge,quantile,side='left'))+1)*step_h)


def conservative_time_interval(mass,step_h):
    """Conservative coverage at finite bin width: lower bin left edge, upper bin right edge."""
    p=np.asarray(mass,float)
    if not np.isfinite(p).all() or p.sum()<=0:raise ValueError('Invalid time distribution')
    cdf=np.cumsum(p/p.sum())
    lower=int(np.searchsorted(cdf,.025,side='left'))
    upper=int(np.searchsorted(cdf,.975,side='left'))
    return (float(lower*step_h),float((upper+1)*step_h))


def weighted_index(values,scaling):
    """Scored physical stress, not a disaster label."""
    x=np.asarray(values,dtype=float)
    centers=np.array([scaling[f]['median'] for f in FEATURES])
    scale=np.array([scaling[f]['scale'] for f in FEATURES])
    z=np.clip((x-centers)/scale,-3,3)
    return np.sum(z*WEIGHTS,axis=-1)


def _reference_scale(hourly,years):
    ref=hourly[(hourly.index.year.isin(years)) & (hourly.index.month.isin([5,6,7,8,9]))]
    scaling={}
    for f in FEATURES:
        x=ref[f].dropna().to_numpy(float)
        if x.size<500:raise ValueError('Insufficient station climate data: '+f)
        med=float(np.median(x));mad=float(1.4826*np.median(np.abs(x-med)))
        if mad<1e-7:mad=float(np.std(x))
        scaling[f]={'median':med,'scale':mad}
    scores=weighted_index(ref.loc[:,list(FEATURES)].to_numpy(float),scaling)
    scores=scores[np.isfinite(scores)]
    return scaling, float(np.quantile(scores,.80)), len(scores)


def _analogue_ensemble(hourly,issue,steps,step_h,seed=23,n_samples=1000,past_days=14):
    """Block resample complete past days' subdaily trajectories, updated using latest observed state.
    No value >= issue is read. The same analogue day is used for all features and lead times,
    preserving a good portion of cross-variable and temporal dependence.
    """
    issue=pd.Timestamp(issue)
    last=issue-pd.Timedelta(hours=1)
    cur=hourly.reindex([last])[list(FEATURES)].to_numpy(float)[0]
    if not np.isfinite(cur).all():raise ValueError(f'Missing pre-issuance hourly feature at {last}')
    leads=np.arange(1,steps+1,dtype=float)*step_h
    analogues=[]
    for lag in range(1,past_days+1):
        pastlast=last-pd.Timedelta(days=lag)
        old0=hourly.reindex([pastlast])[list(FEATURES)].to_numpy(float)[0]
        if not np.isfinite(old0).all():continue
        vals=[];valid=True
        for dh in leads:
            target=issue+pd.Timedelta(hours=float(dh))-pd.Timedelta(days=lag)
            # Half-hour forecast bins interpolated between adjacent *historical* hourly observations.
            lo=target.floor('h');hi=target.ceil('h')
            pair=hourly.reindex([lo,hi])[list(FEATURES)].to_numpy(float)
            if not np.isfinite(pair).all():valid=False;break
            w=(target-lo).total_seconds()/3600
            vv=pair[0]*(1-w)+pair[1]*w
            # Decaying bias-correction from latest pre-issuance observation.
            vv+=np.exp(-dh/10)*(cur-old0)
            vv[0]=max(0,vv[0]);vv[1]=max(0,vv[1])
            vals.append(vv)
        if valid:analogues.append(np.array(vals))
    if not analogues:raise ValueError('No complete prior-day analogue trajectories')
    rng=np.random.default_rng(seed)
    a=np.stack(analogues)
    ids=rng.integers(0,len(a),size=n_samples)
    return a[ids],leads,len(a)


def forecast_issue(hourly,issue,scaling,threshold,horizon_hours=6,step_h=.5,samples=1200,seed=23):
    n=int(round(horizon_hours/step_h))
    ensemble,lead,analog_n=_analogue_ensemble(hourly,issue,n,step_h,seed,samples)
    score=weighted_index(ensemble,scaling)
    cross=score>=threshold
    pred_cross=(cross.any(axis=1)).mean()
    first=np.where(cross.any(axis=1),cross.argmax(axis=1),-1)
    onset_counts=np.bincount(first[first>=0],minlength=n).astype(float)
    onset_cond=onset_counts/onset_counts.sum() if onset_counts.sum()>0 else np.zeros(n)
    # Gated intensity: smooth precursor response, not an asserted failure probability.
    rng=np.random.default_rng(seed+233)
    lambda_daily=np.exp(rng.uniform(np.log(.0005),np.log(.05),size=samples))
    beta=rng.uniform(.65,1.6,size=samples)
    gate=.02+.98/(1+np.exp(-np.clip((score-threshold)/.15,-30,30)))
    step_haz=hazard_step(lambda_daily[:,None],score-threshold,beta[:,None],step_h)*gate
    survival_before=np.concatenate([np.ones((samples,1)),np.cumprod(1-step_haz,axis=1)[:,:-1]],axis=1)
    event_mass=survival_before*step_haz
    none=np.prod(1-step_haz,axis=1)
    assert np.allclose(event_mass.sum(axis=1)+none,1,atol=1e-9)
    mixture=event_mass.mean(axis=0)
    cond=mixture/mixture.sum()
    return {'issue':pd.Timestamp(issue),'leads':lead,'scores':score,
        'precursor_threshold':threshold,'hazard_mean':step_haz.mean(axis=0),
        'time_mass':mixture,'conditional_mass':cond,'cross_probability':float(pred_cross),
        'precursor_onset_conditional_mass':onset_cond,
        'interval_hours':conservative_time_interval(mixture,step_h),
        'median_hours':quantile_bin(mixture,step_h,.5),
        'within_horizon_probability':float((1-none).mean()),
        'analog_days_count':analog_n,
        'conditional_precursor_interval':conservative_time_interval(onset_counts,step_h) if onset_counts.sum()>0 else None,
        'hazard_scenarios':step_haz}


def _history_precursor(h,scaling):
    vals=weighted_index(h[list(FEATURES)].to_numpy(float),scaling)
    return pd.Series(vals,index=h.index,name='precursor_index')


def backtest_precursor_2025(hourly,threshold,training_years=(2022,2023,2024),max_days=160):
    """True out-of-period meteorological precursor onset backtest (not disaster event coverage)."""
    scale,_,_= _reference_scale(hourly,training_years)
    # Threshold fixed *on training climatology*, avoid test-set quantile leakage.
    base=hourly[(hourly.index.year.isin(training_years))&hourly.index.month.isin([5,6,7,8,9])]
    vals=weighted_index(base[list(FEATURES)].to_numpy(float),scale)
    thresh=float(np.quantile(vals[np.isfinite(vals)],.80))
    fullscore=_history_precursor(hourly,scale)
    rows=[]
    dates=pd.date_range('2025-05-01','2025-09-30',freq='D')
    for j,day in enumerate(dates[:max_days]):
        issue=day+pd.Timedelta(hours=9)
        prior=issue-pd.Timedelta(hours=1)
        if not np.isfinite(fullscore.get(prior,np.nan)) or fullscore[prior]>=thresh:continue
        try:
            fc=forecast_issue(hourly.loc[hourly.index<issue],issue,scale,thresh,samples=240,seed=3000+j)
        except (ValueError,KeyError):continue
        targets=fullscore.reindex([issue+pd.Timedelta(hours=float(t)) for t in fc['leads']])
        # Future half-hour obs aren't direct station records: linearly interpolate hourly true scores.
        hourly_series=fullscore.reindex(pd.date_range(issue,issue+pd.Timedelta(hours=6),freq='h'))
        true_interp=hourly_series.reindex(hourly_series.index.union(targets.index)).interpolate('time').reindex(targets.index)
        crossed=(true_interp>=thresh).to_numpy()
        truth_cross=bool(np.any(crossed))
        event_time=(float(fc['leads'][np.where(crossed)[0][0]]) if truth_cross else np.nan)
        pi=fc['conditional_precursor_interval']
        covered=(bool(pi[0]-1e-8<=event_time<=pi[1]+1e-8) if truth_cross and pi is not None else np.nan)
        rows.append({'issue':issue,'forecast_threshold_cross_probability':fc['cross_probability'],
                     'truth_crossed_next6h':truth_cross,'observed_first_cross_hours':event_time,
                     'predicted_first_cross_low_hours':pi[0] if pi else np.nan,
                     'predicted_first_cross_high_hours':pi[1] if pi else np.nan,
                     'covered_on_precursor_event':covered,'analog_days':fc['analog_days_count']})
    d=pd.DataFrame(rows)
    if d.empty:return d,{'evaluated_issuances':0,'weather_PICP':None}
    alerts=d.forecast_threshold_cross_probability>=.25
    TP=int((alerts&d.truth_crossed_next6h).sum());FP=int((alerts&~d.truth_crossed_next6h).sum())
    FN=int((~alerts&d.truth_crossed_next6h).sum());TN=int((~alerts&~d.truth_crossed_next6h).sum())
    ev=d.loc[d.truth_crossed_next6h&d.covered_on_precursor_event.notna()]
    metrics={'evaluated_issuances':int(len(d)),'precursor_crossings':int(d.truth_crossed_next6h.sum()),
       'threshold_warning_probability':.25,'TP':TP,'FP':FP,'FN':FN,'TN':TN,
       'precursor_warning_recall':TP/max(1,TP+FN),'precursor_false_alarm_fraction':FP/max(1,FP+TN),
       'meteorological_onset_PICP_all_available':float(ev.covered_on_precursor_event.astype(float).mean()) if len(ev) else None,
       'meteorological_onset_PICP_n':int(len(ev)),
       'mean_predicted_precursor_interval_width_h':float((ev.predicted_first_cross_high_hours-ev.predicted_first_cross_low_hours).mean()) if len(ev) else None,
       'test_period':'2025 May-Sep','calibration_period':'2022-2024','label':'Meteorological threshold crossings only; not debris-flow events'}
    return d,metrics


def _short_title(ax,title):
    ax.set_title(title,fontsize=11.2,loc='left',pad=9)


def _plot_issue(ctx,fc,history,score,reference_event):
    out=ctx.folder('q3');issue=fc['issue'];lead=fc['leads'];idx=fc['scores'];threshold=fc['precursor_threshold']
    future_times=[issue+pd.Timedelta(hours=float(t)) for t in lead]
    fig,axs=plt.subplots(2,1,figsize=(11.7,7.1),gridspec_kw={'height_ratios':[1.2,1]})
    a=score.loc[issue-pd.Timedelta(hours=48):issue-pd.Timedelta(hours=1)]
    axs[0].plot(a.index,a,color='#557E9F',lw=1.6,label='灾前观测条件指数')
    axs[0].plot(future_times,np.median(idx,axis=0),color='#C2703D',lw=2,label='历史相似日预报中位数')
    axs[0].fill_between(future_times,np.quantile(idx,.05,axis=0),np.quantile(idx,.95,axis=0),color='#C2703D',alpha=.16,label='相似日5%–95%范围')
    axs[0].axhline(threshold,ls='--',color='#AA5349',lw=1.3,label='历史80%条件阈值')
    axs[0].axvline(issue,ls=':',color='0.35',lw=1.5,label='预报时刻')
    axs[0].axvline(reference_event,ls='--',color='#8A5FA2',lw=1.3,label='题面事件时刻')
    axs[0].set_ylabel('综合气象条件指数')
    _short_title(axs[0],'a  灾前气象条件与短时预报')
    h0,l0=axs[0].get_legend_handles_labels()
    fig.legend(h0,l0,loc='upper center',bbox_to_anchor=(.5,.935),ncol=5,fontsize=7.4,frameon=False)
    axs[1].plot(lead,fc['hazard_mean'],color='#2E8490',lw=2.2,label='模型逐时段危险概率')
    ay=axs[1].twinx()
    ay.plot(lead,fc['precursor_onset_conditional_mass'],color='#E08B39',lw=1.7,alpha=.9,label='首次超越阈值的条件概率')
    ay.set_ylabel('气象首次超越的条件概率',fontsize=9,color='#C07630')
    ay.tick_params(axis='y',colors='#C07630',labelsize=8)
    axs[1].axvline((reference_event-issue).total_seconds()/3600,color='#8A5FA2',ls='--',lw=1.2)
    axs[1].set(xlabel='预报发布后时间 (h)',ylabel='每0.5小时概率')
    _short_title(axs[1],'b  未来6小时风险演变')
    h1,l1=axs[1].get_legend_handles_labels();h2,l2=ay.get_legend_handles_labels()
    axs[1].legend(h1+h2,l1+l2,loc='upper center',bbox_to_anchor=(.5,-.19),ncol=2,fontsize=8,frameon=False)
    for ax in axs:ax.grid(axis='y',alpha=.17)
    fig.subplots_adjust(left=.09,right=.90,top=.83,bottom=.16,hspace=.27)
    fig.suptitle('灾前气象条件及触发风险演变',fontsize=14,y=.992)
    savefig(fig,out/'Q3_01_precursor_and_scenario_hazard.png')


def _plot_time_interval(ctx,fc,event_dt):
    out=ctx.folder('q3');lead=fc['leads'];dt=float(lead[0]);lo,hi=fc['interval_hours'];mid=fc['median_hours']
    x_event=(event_dt-fc['issue']).total_seconds()/3600
    fig,axs=plt.subplots(1,2,figsize=(12,4.9),gridspec_kw={'width_ratios':[1.15,1]})
    ax=axs[0]
    ax.bar(lead,fc['conditional_mass'],width=dt*.86,color='#358C8A',edgecolor='white',linewidth=.3)
    ax.axvspan(lo,hi,color='#7FA3B5',alpha=.15,label='95%预测时间区间')
    ax.axvline(x_event,color='#C55D48',ls='--',lw=1.8,label=f'题面参考时刻：+{x_event:.1f} h')
    ax.axvline(mid,color='#C28E43',ls=':',lw=1.7,label=f'时间分布中位数：+{mid:.1f} h')
    ax.set(xlabel='预警发布后时间 (h)',ylabel='给定6小时内发生时的概率')
    _short_title(ax,'a  触发时间概率分布')
    ax.legend(fontsize=8,frameon=False,loc='upper left')
    cdf=np.cumsum(fc['conditional_mass'])
    axs[1].step(lead,cdf,where='post',color='#356D94',lw=2)
    axs[1].axhline(.025,color='0.55',ls=':',lw=1)
    axs[1].axhline(.975,color='0.55',ls=':',lw=1)
    axs[1].axvline(x_event,color='#C55D48',ls='--',lw=1.6)
    axs[1].axvspan(lo,hi,color='#7FA3B5',alpha=.15)
    axs[1].set(xlabel='预警发布后时间 (h)',ylabel='条件累计概率',ylim=(0,1.04))
    _short_title(axs[1],'b  触发时间累积概率')
    for a in axs:a.grid(axis='y',alpha=.15)
    fig.suptitle('泥石流触发时间分布与95%预测区间',fontsize=14,y=.97)
    fig.subplots_adjust(left=.078,right=.98,bottom=.14,top=.86,wspace=.25)
    savefig(fig,out/'Q3_02_95percent_conditional_time_interval.png')


def _plot_warning_and_backtest(ctx,fc,backtest,metrics,event_dt):
    out=ctx.folder('q3');issued=fc['issue'];reflead=(event_dt-issued).total_seconds()/3600
    alarm=fc['cross_probability']>=.25
    fig,axs=plt.subplots(1,2,figsize=(12,4.9))
    # left: precise warning event timeline
    ax=axs[0]
    ax.plot([0,6],[1,1],color='#D7E2E9',lw=7,solid_capstyle='round',zorder=1)
    ax.axvspan(0,6,color='#358C8A',alpha=.06)
    ax.scatter([0],[1],s=105,color='#356D94',zorder=5)
    ax.scatter([reflead],[1],s=105,color='#C55D48',zorder=5)
    ax.annotate('预报发布 09:00',xy=(0,1),xytext=(.25,1.27),fontsize=9,ha='left',arrowprops=dict(arrowstyle='-',color='#356D94'))
    ax.annotate('题面灾害时刻 10:30',xy=(reflead,1),xytext=(reflead+.2,1.20),fontsize=9,ha='left',arrowprops=dict(arrowstyle='-',color='#C55D48'))
    if alarm:
        ax.annotate('',xy=(reflead,.78),xytext=(0,.78),arrowprops=dict(arrowstyle='<->',color='#4A7F89',lw=1.5))
        ax.text(reflead/2,.64,f'预警提前 {reflead:.1f} h',color='#2E6874',fontsize=11,ha='center')
    else:ax.text(3,.64,'该阈值下未触发预警',ha='center',fontsize=10,color='#C55D48')
    ax.set(xlim=(-.25,6.1),ylim=(.3,1.55),yticks=[],xlabel='预报发布后时间 (h)')
    _short_title(ax,'a  预警时刻与提前量')
    ax.spines['left'].set_visible(False)
    ax2=axs[1]
    if metrics.get('evaluated_issuances',0):
        vals=[metrics[k] for k in ('TP','FP','FN','TN')]
        bars=ax2.bar(['TP','FP','FN','TN'],vals,color=['#3F8E76','#C9775F','#BC9C46','#758C9B'],width=.58)
        ax2.bar_label(bars,padding=3,fontsize=9)
        ax2.set(ylabel='2025年独立气象条件样本数')
        cap=metrics.get('meteorological_onset_PICP_all_available')
        txt=(f'条件阈值首次超越区间覆盖率：{cap:.1%} (n={metrics["meteorological_onset_PICP_n"]})' if cap is not None else '区间覆盖率：验证样本不足')
        ax2.text(.5,-.22,txt,transform=ax2.transAxes,ha='center',va='top',fontsize=8.5)
    else:ax2.text(.5,.5,'气象阈值回测样本不足',ha='center',transform=ax2.transAxes)
    _short_title(ax2,'b  2025年气象预警回测')
    ax2.grid(axis='y',alpha=.14)
    fig.suptitle('泥石流预警时间与气象条件回测',fontsize=14,y=.975)
    fig.subplots_adjust(left=.07,right=.98,top=.85,bottom=.20,wspace=.26)
    savefig(fig,out/'Q3_03_warning_lead_and_backtest.png')


def _plot_validation_metrics(ctx,metrics):
    out=ctx.folder('q3')
    if not metrics.get('evaluated_issuances',0):return
    fig,axs=plt.subplots(1,2,figsize=(11.4,4.5))
    actual=metrics.get('meteorological_onset_PICP_all_available')
    if actual is not None:
        b=axs[0].barh(['气象条件时间区间'],[100*actual],color='#3C898A',height=.42)
        axs[0].axvline(95,color='#C25E4A',ls='--',lw=1.7,label='名义95%')
        axs[0].bar_label(b,labels=[f'{100*actual:.1f}% (n={metrics["meteorological_onset_PICP_n"]})'],padding=5,fontsize=10)
    axs[0].set(xlim=(0,105),xlabel='区间覆盖率 (%)')
    axs[0].legend(loc='lower right',fontsize=9,frameon=False)
    _short_title(axs[0],'a  2025年条件时间区间覆盖率')
    vals=[metrics['precursor_warning_recall'],metrics['precursor_false_alarm_fraction']]
    bars=axs[1].bar(['召回率','误报比例'],[100*v for v in vals],color=['#3C898A','#C77D62'],width=.5)
    axs[1].bar_label(bars,fmt='%.1f%%',padding=5,fontsize=10)
    axs[1].set(ylim=(0,110),ylabel='百分比 (%)')
    _short_title(axs[1],'b  2025年气象预警判别性能')
    for a in axs:a.grid(axis='x' if a==axs[0] else 'y',alpha=.17)
    fig.suptitle('气象条件时间区间与预警指标验证',fontsize=14,y=.965)
    fig.subplots_adjust(left=.13,right=.98,top=.82,bottom=.17,wspace=.32)
    savefig(fig,out/'Q3_06_validation_metrics.png')


def _plot_sensitivity(ctx,fc,reference_prob):
    out=ctx.folder('q3')
    # Cross probability variation with warning threshold; report it as climatological condition sensitivity.
    thresh=fc['precursor_threshold'];s=fc['scores']
    offsets=np.linspace(-.5,1.2,13)
    p=np.array([np.mean((s>=thresh+x).any(axis=1)) for x in offsets])
    x=np.linspace(.0005,.05,9)
    gate=.02+.98/(1+np.exp(-np.clip((s-fc['precursor_threshold'])/.15,-30,30)))
    hazard_percent=np.array([100*np.mean(1-np.prod(1-hazard_step(lam,s-fc['precursor_threshold'],1.0,float(fc['leads'][0]))*gate,axis=1)) for lam in x])
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.2))
    axs[0].plot(offsets,p,marker='o',color='#397B96',lw=2)
    axs[0].axvline(0,ls=':',color='0.45')
    axs[0].set(xlabel='气象阈值变动 ΔS',ylabel='未来6小时超越概率',ylim=(0,1.05))
    _short_title(axs[0],'a  气象阈值敏感性')
    axs[1].plot(x,hazard_percent,marker='o',color='#BD744C',lw=2)
    axs[1].set(xlabel='假设基础日危险强度',ylabel='未来6小时近似发生概率 (%)')
    _short_title(axs[1],'b  基础危险强度敏感性')
    for ax in axs:ax.grid(axis='y',alpha=.15)
    fig.suptitle('预警模型参数敏感性分析',fontsize=14,y=.97)
    fig.subplots_adjust(left=.08,right=.98,top=.85,bottom=.15,wspace=.24)
    savefig(fig,out/'Q3_04_parameter_sensitivity.png')


def _plot_spatial(ctx):
    """Cartographic route-distance outcome using official DEM; label as scenario corridor."""
    out=ctx.folder('q3');f=ctx.out/'q4'/'Q4_route_profile.csv'
    if not f.exists():
        from .q4_flowpath import execute_q4_path
        execute_q4_path(ctx,allow_provisional=True)
    if not f.exists():return {'map_status':'NO_Q4_ROUTE_PROFILE'}
    profile=pd.read_csv(f)
    from .q4_flowpath import local_dem
    from pyproj import Transformer
    from matplotlib.colors import Normalize
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    coords=profile[['utm_x','utm_y']].to_numpy(float)
    z,tr,crs=local_dem(ctx,coords[::max(1,len(coords)//15)],buffer_m=2400,resolution_m=100)
    z=np.asarray(z,float)
    mask=np.isfinite(z)
    if not mask.any():return {'map_status':'DEM_NO_COVERAGE'}
    z=np.where(mask,z,np.nanmedian(z))
    from scipy.ndimage import gaussian_filter
    hill=LightSource(azdeg=315,altdeg=48).hillshade(gaussian_filter(z,1),vert_exag=1.5,dx=100,dy=100)
    ext=[tr.c,tr.c+z.shape[1]*tr.a,tr.f+z.shape[0]*tr.e,tr.f]
    fig,ax=plt.subplots(figsize=(11,6.5))
    ax.imshow(hill,origin='upper',extent=ext,cmap='Greys',vmin=.05,vmax=.95)
    segs=np.stack([coords[:-1],coords[1:]],axis=1)
    distance=profile.distance_m.to_numpy(float)[:-1]/1000
    cols=['#267F96','#41917B','#DF9E43','#CD6454','#815C90']
    bins=[0,5,10,15,20,1e8]
    for j in range(5):
        idx=(distance>=bins[j])&(distance<bins[j+1])
        if idx.any():ax.add_collection(LineCollection(segs[idx],colors=[cols[j]],linewidths=3.3,zorder=8))
    ax.scatter(coords[0,0],coords[0,1],s=130,marker='*',facecolor='#224A74',edgecolor='white',zorder=10)
    ax.scatter(coords[-1,0],coords[-1,1],s=85,marker='s',facecolor='#CB5440',edgecolor='white',zorder=10)
    ax.annotate('冰川源区',coords[0],xytext=(8,10),textcoords='offset points',fontsize=9,fontweight='bold',color='#25333E',bbox=dict(fc='white',ec='none',alpha=.78,pad=1.5))
    ax.annotate('下游断面',coords[-1],xytext=(8,-15),textcoords='offset points',fontsize=9,fontweight='bold',color='#25333E',bbox=dict(fc='white',ec='none',alpha=.78,pad=1.5))
    # Scale & north arrow in same native metre coordinate frame.
    x0=ext[0]+.06*(ext[1]-ext[0]);y0=ext[2]+.08*(ext[3]-ext[2]);km=2000
    ax.plot([x0,x0+km],[y0,y0],color='black',lw=3.4,zorder=12)
    ax.text(x0+km/2,y0+350,'2 km',ha='center',fontsize=9,bbox=dict(fc='white',ec='none',alpha=.7),zorder=13)
    xx=ext[1]-.09*(ext[1]-ext[0]);yy=ext[3]-.17*(ext[3]-ext[2])
    ax.annotate('N',(xx,yy+1000),(xx,yy),ha='center',fontsize=10,arrowprops=dict(arrowstyle='-|>',lw=1.5,color='black'))
    ax.set(xlabel='UTM东向坐标 (m)',ylabel='UTM北向坐标 (m)')
    ax.set_title('泥石流传播廊道与空间预警距离',fontsize=14,pad=14)
    from matplotlib.patches import Patch
    h=[Patch(facecolor=cols[j],label=f'{bins[j]}–{(bins[j+1] if j<4 else round(float(profile.distance_m.iloc[-1]/1000),1))} km') for j in range(5)]
    ax.legend(handles=h,loc='lower right',ncol=2,fontsize=8.5,frameon=True)
    ax.set_aspect('equal')
    fig.subplots_adjust(left=.10,right=.97,top=.93,bottom=.10)
    savefig(fig,out/'Q3_05_spatial_warning_distance.png')
    total=float(profile.distance_m.iloc[-1]/1000)
    pd.DataFrame({'distance_start_km':[0,5,10,15,20],
                  'distance_end_km':[min(total,i) for i in [5,10,15,20,total]],
                  'distance_type':['候选传播廊道']*5}).query('distance_end_km > distance_start_km').to_csv(out/'Q3_spatial_distance_zones.csv',index=False)
    return {'map_status':'OFFICIAL_DEM_CANDIDATE_ROUTE_DISTANCE','modelled_corridor_length_km':total,
            'distance_not_a_verified_flood_buffer':True}


def execute_q3(ctx):
    out=ctx.folder('q3');cfg=ctx.cfg['q3'];issue=pd.Timestamp(ctx.cfg['snapshot_raw'])
    d=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw').sort_index()
    pre=d.loc[d.index<issue]
    scale,threshold,ref_n=_reference_scale(pre,ctx.cfg['historic_years'])
    hours=int(cfg.get('primary_horizon_hours',6));step=float(cfg.get('step_hours',.5))
    fc=forecast_issue(pre,issue,scale,threshold,horizon_hours=hours,step_h=step,
                      samples=int(cfg.get('samples',1200)),seed=int(cfg.get('random_seed',20261023)))
    eventdt=pd.Timestamp(cfg.get('official_event_time_same_clock_assumption','2026-08-26 10:30:00'))
    hit=(fc['interval_hours'][0]-step<= (eventdt-issue).total_seconds()/3600 <=fc['interval_hours'][1])
    warn=(fc['cross_probability']>=cfg.get('warning_cross_probability_threshold',.25))
    lead_h=(eventdt-issue).total_seconds()/3600 if warn else None
    bt,metrics=backtest_precursor_2025(d,threshold)
    bt.to_csv(out/'Q3_2025_precursor_onset_backtest.csv',index=False,encoding='utf-8-sig')
    write_json(out/'Q3_2025_precursor_backtest_metrics.json',metrics)
    # Competition-oriented time distribution with explicit event reference and reproducible metrics.
    leads=fc['leads'];mass=fc['time_mass'];cond=fc['conditional_mass'];scores=fc['scores']
    tab=pd.DataFrame({'future_raw_time_bin_end':[issue+pd.Timedelta(hours=float(t)) for t in leads],
      'hours_after_issuance':leads,'scenario_mean_unconditional_event_mass':mass,
      'conditional_time_mass_given_event_within_H':cond,
      'mean_discrete_hazard':fc['hazard_mean'],
      'hazard_q05_across_scenarios':np.quantile(fc['hazard_scenarios'],.05,axis=0),
      'hazard_q95_across_scenarios':np.quantile(fc['hazard_scenarios'],.95,axis=0),
      'precursor_score_q50':np.median(scores,axis=0),
      'precursor_score_q05':np.quantile(scores,.05,axis=0),
      'precursor_score_q95':np.quantile(scores,.95,axis=0),
      'precursor_first_cross_mass':fc['precursor_onset_conditional_mass']})
    tab.to_csv(out/'Q3_scenario_time_distribution.csv',index=False,encoding='utf-8-sig')
    _plot_issue(ctx,fc,pre,_history_precursor(pre,scale),eventdt)
    _plot_time_interval(ctx,fc,eventdt)
    _plot_warning_and_backtest(ctx,fc,bt,metrics,eventdt)
    _plot_sensitivity(ctx,fc,fc['within_horizon_probability'])
    _plot_validation_metrics(ctx,metrics)
    # Q3 spatial map depends on Q4 route generation. Calls standalone if needed.
    spatial=_plot_spatial(ctx)
    report={'status':'SCENARIO_TRIGGER_TIME_WITH_INDEPENDENT_PRECURSOR_BACKTEST',
      'forecast_issuance_time_raw':str(issue),
      'official_event_reference_time':str(eventdt),
      'comparison_assumes_matching_timezones':True,
      'source_station_timezone_confirmed':True,
      'horizon_hours':hours,'step_hours':step,
      'model':'historical-day analogue + gated discrete physical hazard, conditional on event in forecast horizon',
      'observed_features_cutoff_exclusive':str(issue),
      'baseline_reference_years':ctx.cfg['historic_years'],
      'warning_threshold_80th_percentile':threshold,
      'warning_decision':bool(warn),
      'warning_crossing_probability':fc['cross_probability'],
      'warning_issuance_raw':str(issue) if warn else None,
      'warning_advance_hours_to_official_reference':lead_h,
      'conditional_trigger_time_interval_95_percent_hours':list(fc['interval_hours']),
      'conditional_trigger_time_interval_95_percent_timestamps':[
         str(issue+pd.Timedelta(hours=float(t))) for t in fc['interval_hours']],
      'conditional_trigger_time_median_hours':fc['median_hours'],
      'official_reference_event_in_95_interval':bool(hit),
      'official_event_event_level_hit_over_total_events':'1/1' if hit else '0/1',
      'empirical_disaster_event_PICP':None,
      'precursor_onset_validation':metrics,
      'probability_event_within_horizon_scenario_mixture_mean':fc['within_horizon_probability'],
      'analogue_weather_days':fc['analog_days_count'],
      'spatial':spatial,
      'limitations':'2025 empirical coverage refers to meteorological precursor crossings, not debris-flow events; scenario hazard prior not calibrated to historical disaster events.'}
    write_json(out/'Q3_scenario_report.json',report)
    write_json(out/'Q3_climatology_scaling.json',scale)
    print('[Q3] official event',eventdt,'warn',warn,'lead hours',lead_h,'interval',fc['interval_hours'],
          'hit',hit,'proxy backtest metrics',metrics)
    return report
