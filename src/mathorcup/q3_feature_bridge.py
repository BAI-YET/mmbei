"""Q1 -> Q3 feature-use audit and fixed-weight sensitivity.

Q1 selects features for the 2026 EVENT case: consequently a 2025 backtest that
uses those post-hoc selected features must NOT be called genuinely out-of-time.
This supplemental module shows whether additional physical groups change Q3's
precursor score without turning proxy-condition targets into observed disasters.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from .common import savefig,write_json
from .q3_task_aligned import prepare_joint,anomaly_score

BASE=('temp_2m','sw_down','vpd_2m')

def execute_q3_feature_bridge(ctx):
    q3=ctx.folder('q3');p=ctx.folder('q1')/'Q1_top5_with_stability.csv'
    if not p.exists():
        raise FileNotFoundError('Q3 feature audit requires Q1 top5 output. Run q1 first.')
    chosen=pd.read_csv(p)
    h=prepare_joint(ctx)
    current=pd.Timestamp(ctx.cfg['snapshot_raw'])-pd.Timedelta(hours=1)
    baseline=h.index.year.isin(ctx.cfg['historic_years']) & h.index.month.isin([5,6,7,8,9])
    original,info=anomaly_score(h,ctx.cfg['historic_years'])
    local=original['combined_condition_score']
    # For each candidate compute its season/month/hour robust deviation and
    # partial incremental novelty against the baseline trio, not a disaster AUC.
    rows=[];indices={}
    for feat in chosen.feature:
        if feat not in h:rows.append({'feature':feat,'status':'not_in_processed_shared_features'});continue
        series=h[feat]
        reference=series.loc[baseline].dropna()
        if len(reference)<250 or current not in h.index:
            rows.append({'feature':feat,'status':'insufficient_reference'});continue
        data=reference.groupby([reference.index.month,reference.index.hour])
        med=data.median()
        mad=data.apply(lambda d:1.4826*np.median(np.abs(d-np.median(d))))
        med_idx=pd.MultiIndex.from_arrays([series.index.month,series.index.hour]);means=med.reindex(med_idx).to_numpy()
        scale=mad.reindex(med_idx).to_numpy()
        backup=max(.05,float(np.nanmedian(np.abs(reference-reference.median()))*1.4826*.2))
        scale=np.maximum(np.where(np.isfinite(scale),scale,backup),backup)
        z=pd.Series(np.clip((series.to_numpy(float)-means)/scale,-5,5),index=series.index)
        # magnitude rather than sign; a large negative anomaly may be equally important.
        rare=1/(1+np.exp(-(z.abs()-1)))
        indices[feat]=rare
        m=rare.notna() & local.notna() & baseline
        corr=float(spearmanr(rare.loc[m],local.loc[m]).statistic) if m.sum()>100 else np.nan
        rows.append({'feature':feat,'status':'evaluated','event_same_hour_z':float(z.loc[current]),
          'baseline_valid_hours':int(m.sum()),'spearman_abs_anomaly_vs_Q3_reference':corr,
          'incremental_information_interpretation':'descriptive only; not cross-event validation'})
    table=pd.DataFrame(rows);table.to_csv(q3/'Q3_Q1_feature_inheritance_audit.csv',index=False)
    if indices:
        mat=pd.DataFrame(indices,index=h.index)
        extra=mat.mean(axis=1,skipna=False)
        # Prespecified blending .75/.25 for exploration, never fit to one disaster.
        extended=.75*local+.25*extra
        valid=extended.notna()&local.notna()&baseline
        diff=extended.loc[valid]-local.loc[valid]
        result={'Q1_features_list':chosen.feature.tolist(),'baseline_Q3_features':list(BASE),
          'n_Q1_physical_features_calculable':int(len(indices)),
          'n_reference_hour_pairs':int(valid.sum()),
          'mean_difference_extended_minus_original':float(diff.mean()) if len(diff) else None,
          'event_ordinary_Q3_precursor':float(local.loc[current]),
          'event_extended_precursor':float(extended.loc[current]) if np.isfinite(extended.loc[current]) else None,
          'not_independent_disaster_validation':True,
          'warning':'Q1 feature selection used 2026 event; sensitivity is descriptive not 2025 blind model selection.'}
        write_json(q3/'Q3_feature_inheritance_sensitivity.json',result)
        fig,axes=plt.subplots(1,2,figsize=(10.8,4.2),layout='constrained')
        month=(h.index>=current-pd.Timedelta(days=10))&(h.index<=current)
        axes[0].plot(h.index[month],local[month],label='基础先兆指数',lw=1.6)
        axes[0].plot(h.index[month],extended[month],label='加入Q1特征的敏感性指数',lw=1.4)
        axes[0].set(title='a  灾前十天先兆指标变化',ylabel='标准化条件分值',xlabel='北京时间')
        axes[0].legend(frameon=False,fontsize=8);axes[0].tick_params(axis='x',rotation=15)
        pp=table[table.status.eq('evaluated')]
        axes[1].barh(pp.feature,pp.event_same_hour_z,color='#497EAB')
        axes[1].axvline(0,color='#777',lw=.8)
        axes[1].set(title='b  Q1入选因素同小时异常',xlabel='稳健Z值')
        savefig(fig,q3/'Q3_07_Q1_feature_inheritance.png')
    return table
