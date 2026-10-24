"""Task-level evidence: export complete comparisons before selecting headline results.
All figures read numerical arrays that are also exported as CSV/JSON.
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from .common import savefig,write_json,LABELS,COLORS


def q1_complete_evidence(ctx,hourly,ranking,ranks,refs,audits):
    out=ctx.folder('q1')
    names=ranking.feature.tolist();ev=ranking.event_value.to_numpy(float)
    rng=np.random.default_rng(20261024);years=refs.reference_year.unique()
    import warnings
    scores=[];zs=[];rankdraw=[]
    for _ in range(600):
        block=pd.concat([refs.loc[refs.reference_year.eq(y),names] for y in rng.choice(years,len(years),replace=True)])
        a=block.to_numpy(float);n=np.isfinite(a).sum(axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            med=np.nanmedian(a,axis=0);mad=1.4826*np.nanmedian(np.abs(a-med),axis=0)
        z=np.divide(ev-med,mad,out=np.full(len(names),np.nan),where=mad>1e-12)
        z[n<35]=np.nan
        left=(1+np.sum(a<=ev,axis=0))/(n+1);right=(1+np.sum(a>=ev,axis=0))/(n+1)
        tail=np.minimum(1.,2*np.minimum(left,right))
        sc=.45*np.minimum(np.abs(z)/3,1)+.55*(1-tail)
        order=np.argsort(-np.nan_to_num(sc,nan=-np.inf));rr=np.empty(len(names));rr[order]=np.arange(1,len(names)+1)
        rr[~np.isfinite(sc)]=np.nan
        scores.append(sc);zs.append(z);rankdraw.append(rr)
    zz=np.array(zs);rr=np.array(rankdraw);ss=np.array(scores)
    d=ranking.copy()
    d['event_minus_reference_median']=d.event_value-d.reference_median
    d['year_block_z_lower']=np.nanquantile(zz,.025,axis=0)
    d['year_block_z_upper']=np.nanquantile(zz,.975,axis=0)
    d['year_block_score_lower']=np.nanquantile(ss,.025,axis=0)
    d['year_block_score_upper']=np.nanquantile(ss,.975,axis=0)
    d['year_block_top5_frequency']=np.sum(rr<=5,axis=0)/np.maximum(1,np.isfinite(rr).sum(axis=0))
    d['year_block_valid_replicates']=np.isfinite(rr).sum(axis=0)
    d['year_block_rank_median']=np.nanmedian(rr,axis=0)
    d.to_csv(out/'Q1_all_factors_effect_and_stability.csv',index=False,encoding='utf-8-sig')
    fig,axs=plt.subplots(1,2,figsize=(13,10),layout='constrained')
    chunks=np.array_split(np.arange(len(d)),2)
    for ax,ids in zip(axs,chunks):
        t=d.iloc[ids].iloc[::-1];y=np.arange(len(t))
        colors=[COLORS[1] if i<=5 else COLORS[0] for i in t['rank']]
        ax.barh(y,t.rarity_anomaly_score,color=colors,alpha=.88)
        ax.set(yticks=y,yticklabels=[f'{int(r.rank):02d}  {r.label}' for r in t.itertuples()],xlim=(0,1.04),xlabel='匹配异常-稀有性得分')
        ax.tick_params(axis='y',labelsize=8);ax.grid(axis='x',alpha=.15)
        for yi,r in enumerate(t.itertuples()):ax.text(r.rarity_anomaly_score+.008,yi,f'{r.rarity_anomaly_score:.3f}',va='center',fontsize=7)
    fig.suptitle(f'全部{len(d)}项有效候选因素及原始前五项',fontsize=13)
    savefig(fig,out/'Q1_11_all_candidate_ranking.png')
    for part,begin in enumerate(range(0,len(d),28),1):
        t=d.iloc[begin:begin+28].iloc[::-1];y=np.arange(len(t))
        fig,ax=plt.subplots(figsize=(8.5,8.3),layout='constrained')
        ax.barh(y,t.rarity_anomaly_score,color=[COLORS[1] if r<=5 else COLORS[0] for r in t['rank']])
        ax.set(yticks=y,yticklabels=[f'{int(r.rank):02d} {r.label}' for r in t.itertuples()],xlim=(0,1.04),xlabel='历史匹配异常-稀有性得分',title=f'全部候选因素排序：第{begin+1}-{min(begin+28,len(d))}项')
        ax.tick_params(axis='y',labelsize=10);ax.grid(axis='x',alpha=.15)
        for yi,r in enumerate(t.itertuples()):ax.text(r.rarity_anomaly_score+.008,yi,f'{r.rarity_anomaly_score:.3f}',va='center',fontsize=8)
        savefig(fig,out/f'Q1_11_part{part}_complete_ranking.png')
    m=pd.DataFrame({f'{w}h':r.set_index('feature').z_mad for w,r in ranks.items()}).reindex(names)
    m.to_csv(out/'Q1_all_factors_multiscale_z.csv',encoding='utf-8-sig')
    fig,ax=plt.subplots(figsize=(8.9,max(8,len(m)*.18)),layout='constrained')
    b=max(1,float(np.nanquantile(np.abs(m),.95)))
    cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#e2e6ea')
    im=ax.imshow(m,cmap=cm,vmin=-b,vmax=b,aspect='auto')
    ax.set(yticks=np.arange(len(m)),yticklabels=[LABELS.get(x,x) for x in m.index],xticks=np.arange(len(m.columns)),xticklabels=m.columns,title='全部有效因素的多时间尺度有向异常')
    ax.tick_params(axis='y',labelsize=7.3)
    fig.colorbar(im,ax=ax,shrink=.55,label='MAD Z；灰色为不可估计')
    savefig(fig,out/'Q1_12_all_factor_multiscale.png')
    corr=refs[names].corr(method='spearman',min_periods=25)
    corr.to_csv(out/'Q1_all_factor_reference_spearman.csv',encoding='utf-8-sig')
    fig,ax=plt.subplots(figsize=(11.5,10),layout='constrained')
    im=ax.imshow(corr,cmap='RdBu_r',vmin=-1,vmax=1)
    ax.set(xticks=np.arange(len(names)),xticklabels=[str(x+1) for x in range(len(names))],yticks=np.arange(len(names)),yticklabels=[f'{i+1:02d} {LABELS.get(x,x)}' for i,x in enumerate(names)])
    ax.tick_params(axis='y',labelsize=5.8);ax.tick_params(axis='x',labelsize=5.5,labelrotation=90)
    ax.set_title('全部候选因素的历史同期相关结构（编号同排序图）')
    fig.colorbar(im,ax=ax,shrink=.55,label='Spearman ρ')
    savefig(fig,out/'Q1_13_all_factor_redundancy.png')
    full=pd.DataFrame(audits)
    write_json(out/'Q1_complete_factor_contract.json',{'predefined':len(full),'eligible':len(d),
       'excluded':full.loc[full.status.ne('ranked')].to_dict('records'),
       'raw_top5':d.head(5).feature.tolist(),
       'mechanism_representatives':d.loc[d.selected_mechanism_representative].feature.tolist(),
       'bootstrap_unit':'whole reference year; focal event held fixed; descriptive stability only',
       'bootstrap_repeats':600,'association_estimand':'focal pre-event window vs calendar/hour-matched historical windows',
       'multi_event_predictive_importance_identified':False})


def q2_complete_evidence(ctx):
    from .q2_spatial_mechanism import FEATURES,_spatial_folds
    from .q2_trigger_evaluation import dynamic_conditions
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    import geopandas as gpd
    out=ctx.folder('q2');summary=json.loads((out/'Q2_model_comparison.json').read_text())
    models=[]
    for key,val in summary.items():
        if key.endswith('_summary') and isinstance(val,dict):models.append({'model':key.removesuffix('_summary'),**val})
    md=pd.DataFrame(models);md.to_csv(out/'Q2_all_models_comparison.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.8),layout='constrained')
    x=np.arange(len(md))
    for ax,key,title in [(axs[0],'auc_oof','空间留出ROC-AUC'),(axs[1],'balanced_accuracy_oof','空间留出平衡准确率')]:
        bars=ax.bar(x,md[key],color=[COLORS[i%len(COLORS)] for i in x]);ax.bar_label(bars,fmt='%.3f',padding=3,fontsize=8)
        ax.axhline(.5,color='#888',ls='--',lw=1);ax.set(xticks=x,xticklabels=['地形四项','形态四项','单高差','单沟口高程','多数类基线'][:len(md)],ylim=(0,1.08),ylabel=title)
        ax.tick_params(axis='x',labelsize=8)
    fig.suptitle('全部候选沟谷成因判别模型（目标为成因类别）')
    savefig(fig,out/'Q2_15_all_model_comparison.png')
    g=gpd.read_file(ctx.processed/'historic_gullies.gpkg');yy=(g['类型']!='降雨型').astype(int).to_numpy()
    groups=_spatial_folds(g,ctx.cfg['q2']['cv_splits']);rows=[]
    xx=g[FEATURES].to_numpy(float)
    for k in np.unique(groups):
        for fields in [FEATURES[:4],[FEATURES[i] for i in [0,1,2,4]]]:
            ids=[FEATURES.index(f) for f in fields];pipe=make_pipeline(StandardScaler(),LogisticRegression(C=ctx.cfg['q2']['C'],max_iter=1000))
            pipe.fit(xx[groups!=k][:,ids],yy[groups!=k])
            for name,c in zip(fields,pipe[-1].coef_[0]):rows.append({'excluded_spatial_block':int(k),'model':'4geo' if fields==FEATURES[:4] else '4geomorph','feature':name,'coefficient_per_training_SD':float(c)})
    pd.DataFrame(rows).to_csv(out/'Q2_all_model_fold_coefficients.csv',index=False,encoding='utf-8-sig')
    norm=pd.DataFrame(xx,columns=FEATURES);norm=(norm-norm.median())/norm.std()
    fig,ax=plt.subplots(figsize=(9.5,4.5),layout='constrained')
    for j,label in enumerate(['降雨型','冰相关型']):
        dat=[norm.loc[yy==j,c].to_numpy() for c in FEATURES]
        b=ax.boxplot(dat,positions=np.arange(len(FEATURES))*3+j*.8,widths=.65,patch_artist=True,showfliers=False)
        for patch in b['boxes']:patch.set_facecolor(COLORS[j]);patch.set_alpha(.65)
        ax.plot([],[],color=COLORS[j],lw=6,label=label)
    ax.set(xticks=np.arange(len(FEATURES))*3+.4,xticklabels=FEATURES,ylabel='标准化变量（用于对比，不是失稳阈值）',title='全部五项历史沟谷判别变量的类别分布');ax.legend(frameon=False)
    savefig(fig,out/'Q2_16_all_predictor_distributions.png')
    h=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw')
    phys=json.loads((out/'Q2_primary_trigger_assessment.json').read_text());z=phys['thresholds']['source_elevation_m'];snap=pd.Timestamp(ctx.cfg['snapshot_raw'])
    rows=[];cadence=[]
    for test in [2023,2024,2025,2026]:
        train=[y for y in ctx.cfg['historic_years'] if y<test]
        d,lim=dynamic_conditions(h,baseline_years=train,source_z=z,station_z=ctx.cfg['station_elevation_m'])
        sub=d.loc[(d.index.year==test)&d.index.month.isin([5,6,7,8,9])&(d.index<snap),'trigger_condition_index']
        valid=sub.notna();flag=sub.ge(.5);starts=flag&~flag.shift(1,fill_value=False)
        rows.append({'test_year':test,'train_years':','.join(map(str,train)),'rain_q90_mm':lim['rain_24h_nonzero_q90_mm'],
          'source_pdd_q90_Cday':lim['source_pdd7_warmseason_q90_Cday'],'valid_hours':int(valid.sum()),'missing_hours':int((~valid).sum()),
          'condition_exceedance_fraction':float(flag.loc[valid].mean()),'condition_episodes':int(starts.sum()),
          'event_accuracy_evaluable':False})
        for step in [1,3,6,24]:
            pos=np.arange(len(sub));reported=pd.Series(sub.iloc[::step].to_numpy(),index=sub.index[::step]).reindex(sub.index,method='ffill')
            ok=valid&reported.notna();coarse=reported.ge(.5)
            both=ok&(flag|coarse)
            delay=[]
            for t in sub.index[starts&ok]:
                window=coarse.loc[t:t+pd.Timedelta(hours=step)]
                yes=window[window].index
                if len(yes):delay.append((yes[0]-t).total_seconds()/3600)
            cadence.append({'test_year':test,'update_interval_h':step,'n_valid_hour_pairs':int(ok.sum()),
                'agreement_with_hourly':float(flag.loc[ok].eq(coarse.loc[ok]).mean()),
                'exceedance_IoU':float((ok&flag&coarse).sum()/max(1,both.sum())),
                'missed_exceedance_hour_fraction':float((ok&flag&~coarse).sum()/max(1,(ok&flag).sum())),
                'median_detected_onset_delay_h':float(np.median(delay)) if delay else np.nan,
                'interpretation':'monitoring cadence vs hourly condition, not disaster prediction accuracy'})
    td=pd.DataFrame(rows);cd=pd.DataFrame(cadence)
    td.to_csv(out/'Q2_rolling_origin_threshold_reliability.csv',index=False)
    cd.to_csv(out/'Q2_update_cadence_reliability.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.4),layout='constrained')
    axs[0].plot(td.test_year,td.rain_q90_mm,'-o',color=COLORS[0]);axs[0].set(xlabel='独立评估年份',ylabel='仅更早年份计算的24h非零降雨阈值 (mm)')
    for y,d in cd.groupby('test_year'):axs[1].plot(d.update_interval_h,d.missed_exceedance_hour_fraction*100,'-o',label=str(y))
    axs[1].set(xlabel='监测更新间隔 (h)',ylabel='相对小时基准漏掉的超限小时 (%)');axs[1].legend(frameon=False)
    fig.suptitle('逐年向前验证与监测时间尺度误差')
    savefig(fig,out/'Q2_17_forward_time_reliability.png')
    sensitivity=[]
    for lapse in [5.,6.5,8.]:
        for q in [.85,.9,.95]:
            d,lim=dynamic_conditions(h,baseline_years=ctx.cfg['historic_years'],source_z=z,station_z=ctx.cfg['station_elevation_m'],lapse=lapse,quantile=q)
            now=d.loc[d.index<snap].iloc[-1]
            sensitivity.append({'lapse_C_per_km':lapse,'reference_quantile':q,'source_pdd_threshold_Cday':lim['source_pdd7_warmseason_q90_Cday'],
               'rain_threshold_mm':lim['rain_24h_nonzero_q90_mm'],'pre_issue_condition_index':now.trigger_condition_index})
    pd.DataFrame(sensitivity).to_csv(out/'Q2_joint_threshold_lapse_sensitivity.csv',index=False)
    write_json(out/'Q2_task_evidence_status.json',{'trigger_model_calibrated_on_disaster_occurrence':False,
        'thresholds_are':'conditional climate exceedance candidates, not measured rock/ice shear strength',
        'time_reliability':'earlier-year threshold fit + later-year evaluation + update-delay errors',
        'space_reliability':'gully-type held-out classification + raster scale agreement; neither validates event occurrence',
        'external_forecast_status':phys['external_forecast_status']})
