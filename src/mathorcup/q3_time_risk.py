"""Shared time-distribution maths and observed-day analogue helpers.
The production forecast entry point is q3_task_aligned.py.
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from pathlib import Path

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
    # Frequencies come from complete observed analogue trajectories.
    from .q3_task_aligned import empirical_first_passage
    distribution=empirical_first_passage(score,step_h,threshold)
    mixture=distribution['unconditional_bin_mass']
    cond=distribution['conditional_bin_mass']
    step_haz=np.zeros_like(score)
    valid=first>=0
    step_haz[np.flatnonzero(valid),first[valid]]=1.
    none=(~valid).astype(float)

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

