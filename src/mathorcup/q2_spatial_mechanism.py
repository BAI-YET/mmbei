"""Q2: real 84-label classification (type, NOT event occurrence), weather quantile warning conditions,
     and physically justified glacier driving stress. No synthetic-negative susceptibility AUC.
"""
from __future__ import annotations
import numpy as np,pandas as pd,geopandas as gpd,rasterio
import matplotlib.pyplot as plt
from rasterio.warp import reproject,Resampling
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,balanced_accuracy_score,precision_recall_fscore_support,brier_score_loss,roc_curve,confusion_matrix
from sklearn.base import clone
from .prepare_geo import glacier_grid
from .common import savefig,write_json,find_first,COLORS,record_warning

FEATURES=['高差','主沟长','沟口H','沟比降','DisFrmFlt']
BINARY_LAKE='冰湖溃决型'

def _spatial_folds(g,n=4):
    lat=pd.to_numeric(g["y_coord"],errors="coerce").to_numpy()
    # Use latitude quantile bins to prevent nearby samples from being scattered across train and test
    groups=pd.qcut(pd.Series(lat).rank(method='first'),q=n,labels=False).to_numpy()
    return groups

def _cross_validate(X,y,groups,C=1.0):
    options={
      'L2_Logistic_4geo':['高差','主沟长','沟口H','沟比降'],
      'L2_Logistic_4geomorph':['高差','主沟长','沟口H','DisFrmFlt'],
      'Single_elevation_drop':['高差'],
      'Single_mouth_elevation':['沟口H']
    }
    rows=[];allpred={}
    for modelname,fields in options.items():
        positions=[FEATURES.index(f) for f in fields]
        if max(positions)>=X.shape[1]:
            # An optional richer model is unavailable when a caller supplies
            # only the official four-feature historical baseline matrix.
            continue
        pred=np.full(len(y),np.nan)
        for k in np.unique(groups):
            train=groups!=k;test=~train
            if len(np.unique(y[train]))<2:continue
            pipeline=Pipeline([('scale',StandardScaler()),('logreg',LogisticRegression(C=C,max_iter=1000,solver='lbfgs'))])
            pipeline.fit(X[train][:,positions],y[train]);pred[test]=pipeline.predict_proba(X[test][:,positions])[:,1]
        keep=np.isfinite(pred)
        if np.unique(y[keep]).size<2:raise ValueError('OOF labels single class; choose different spatial folds')
        allpred[modelname]=pred
        for k in np.unique(groups):
            m=keep & (groups==k)
            if not m.any():continue
            met={'model':modelname,'spatial_block':int(k),'n_test':int(m.sum()),'n_positive':int(y[m].sum())}
            if np.unique(y[m]).size==2:met['auc']=float(roc_auc_score(y[m],pred[m]))
            else:met['auc']=np.nan
            met['balanced_accuracy']=float(balanced_accuracy_score(y[m],pred[m]>=.5)) if np.unique(y[m]).size==2 else np.nan
            rows.append(met)
        predclass=(pred[keep]>=.5).astype(int)
        # All probabilities are strictly held-out; do not compare to train-set AUC.
        allpred[modelname+'_summary']={
           'auc_oof':float(roc_auc_score(y[keep],pred[keep])),
           'balanced_accuracy_oof':float(balanced_accuracy_score(y[keep],predclass)),
           'macro_f1_oof':float(precision_recall_fscore_support(y[keep],predclass,average='macro',zero_division=0)[2]),
           'brier_oof':float(brier_score_loss(y[keep],pred[keep])),'n':int(keep.sum()),
           'n_spatial_blocks':int(len(np.unique(groups)))}
    allpred['majority_baseline_summary']={
      'auc_oof':0.5,'balanced_accuracy_oof':0.5,
      'note':'majority class prediction has no discrimination; constant-score ROC-AUC=0.5'}
    return pd.DataFrame(rows),allpred

def _plot_training_results(ctx,g,y,groups,pred,fold_metrics):
    out=ctx.folder('q2')
    fig,axes=plt.subplots(1,2,figsize=(10,4.5))
    colors=np.where(y==1,COLORS[0],COLORS[1])
    xx=pd.to_numeric(g["x_coord"],errors="coerce").values;yy=pd.to_numeric(g["y_coord"],errors="coerce").values
    axes[0].scatter(xx,yy,c=colors,s=22,alpha=.84,edgecolors='white',linewidths=.35)
    axes[0].set(xlabel='经度 (°E)',ylabel='纬度 (°N)',title='历史泥石流沟谷的成因类型分布')
    axes[0].scatter([],[],c=COLORS[0],label='冰相关型')
    axes[0].scatter([],[],c=COLORS[1],label='降雨型');axes[0].legend(loc='best',frameon=False)
    for k in np.unique(groups):
        yi=y[groups==k]
        if yi.size:
            axes[0].text(np.median(xx[groups==k]),np.median(yy[groups==k]),str(k+1),color='#202831',fontsize=12,weight='bold')
    axes[0].grid(True,alpha=.15)
    fpr,tpr,_=roc_curve(y,pred)
    auc=roc_auc_score(y,pred)
    axes[1].plot(fpr,tpr,color=COLORS[0],lw=2,label=f'L2 Logistic OOF AUC={auc:.3f}')
    axes[1].plot([0,1],[0,1],'--',color='#777',lw=1)
    axes[1].set(xlabel='假阳性率',ylabel='真阳性率',xlim=[0,1],ylim=[0,1],title='空间留出交叉验证 ROC 曲线')
    axes[1].legend(frameon=False,loc='lower right');axes[1].grid(True,alpha=.16)
    savefig(fig,out/'Q2_01_history_gullies_and_oof_roc.png')
    f=fold_metrics[fold_metrics.model=='L2_Logistic_4geo']
    fig,ax=plt.subplots(figsize=(7.3,4))
    ax.bar(f.spatial_block.astype(str),f.auc,color=COLORS[0],alpha=.9)
    ax.axhline(0.5,color='#777',ls='--',lw=1)
    ax.set(ylabel='ROC-AUC',xlabel='留出块编号（按纬度分块）',title='模型在不同空间区域的判别性能',ylim=(0,1))
    for i,row in enumerate(f.itertuples()):
        if np.isfinite(row.auc):ax.text(i,row.auc+.025,f'{row.auc:.2f}',ha='center',fontsize=9)
    savefig(fig,out/'Q2_02_spatial_holdout_auc.png')
    # Conditional gradients are not physical avalanche thresholds
    fig,ax=plt.subplots(figsize=(8.2,4.2))
    b=g.copy();b['binary_label']=y
    a=b.loc[b.binary_label==0,'高差'];c=b.loc[b.binary_label==1,'高差']
    ax.hist(a,bins=15,alpha=.58,label='降雨型沟谷',density=True,color=COLORS[1])
    ax.hist(c,bins=15,alpha=.58,label='冰相关型沟谷',density=True,color=COLORS[0])
    ax.set(xlabel='沟谷高差 (m)',ylabel='密度',title='不同成因沟谷的高差分布')
    ax.legend(frameon=False);savefig(fig,out/'Q2_03_elevation_drop_distributions.png')

def _weather_limits(ctx,hourly):
    # Climate exceedances are computed from physically relevant nonzero rain periods,
    # and thermal energy over full 24-h windows rather than a single 08:00 reading.
    cfg=ctx.cfg;out=ctx.folder('q2');snap=pd.Timestamp(cfg['snapshot_raw'])
    ref=hourly[(hourly.index.year.isin(cfg['historic_years'])) & hourly.index.month.isin(range(5,11))]
    thresholds={}
    for c in ['pdd_7d','rain_24h','rain_72h','net_radiation','vwc_0p1','precip_mm']:
        if c not in ref:continue
        # Hourly radiative variability is diurnal: use 24-hour average energy for thresholding.
        b=(ref[c].rolling('24h',min_periods=18).mean() if c=='net_radiation' else ref[c])
        vals=b.dropna()
        if c in ('rain_24h','rain_72h','precip_mm'):vals=vals[vals>.01]
        if len(vals)<50:continue
        thresholds[c]={'q75':float(vals.quantile(.75)),'q90':float(vals.quantile(.9)),
                       'q95':float(vals.quantile(.95)),'n':int(len(vals)),
                       'conditioning':'nonzero rainfall only' if c.startswith(('rain','precip')) else 'warm-season baseline'}
    current=hourly.loc[hourly.index<snap].iloc[-1]
    event_state={}
    for c,v in thresholds.items():
        x=(float(hourly[c].loc[hourly.index<snap].tail(24).mean()) if c=='net_radiation' else float(current[c]))
        physical=(x>.01 if c.startswith(('rain','precip')) else np.isfinite(x))
        event_state[c]={'pre_issue_value':x,'q95':v['q95'],
             'exceeds_q95':bool(physical and x>=v['q95']) if np.isfinite(x) else None}
    write_json(out/'Q2_weather_reference_quantiles.json',{'baseline_years':cfg['historic_years'],
       'thresholds':thresholds,'snapshot_bjt':str(snap),'event_state':event_state,
       'claim':'Climate candidate thresholds not mechanically measured slope-failure limits.'})
    df=hourly.loc[(hourly.index>=snap-pd.Timedelta(days=55))&(hourly.index<snap)].resample('D').agg(
       {'pdd_7d':'max','rain_24h':'max','net_radiation':'mean'})
    fig,axes=plt.subplots(3,1,figsize=(10.2,7.3),sharex=True,layout='constrained')
    for ax,(c,title),col in zip(axes,[('pdd_7d','7日正积温 (°C·d)'),('rain_24h','24小时累计降水 (mm)'),
       ('net_radiation','日平均净辐射 (W/m²)')],COLORS[:3]):
        ax.plot(df.index,df[c],lw=1.8,color=col,label='2026灾前观测')
        if c in thresholds:
            ax.axhline(thresholds[c]['q90'],ls='--',color='#B36F31',label='历史暖季90%分位')
            ax.axhline(thresholds[c]['q95'],ls=':',color='#B44A4A',label='历史暖季95%分位')
        ax.set_ylabel(title);ax.grid(axis='y',alpha=.15);ax.legend(loc='upper left',fontsize=8,frameon=False)
    axes[0].set_title('历史气候阈值与灾前水热条件')
    axes[-1].set_xlabel('北京时间（UTC+8）')
    savefig(fig,out/'Q2_04_climate_condition_timeseries.png')
    return thresholds

def _glacier_map(ctx):
    out=ctx.folder('q2');d=glacier_grid(ctx)
    m=d['glacier_mask'] & np.isfinite(d['thickness']) & np.isfinite(d['slope'])
    stress=900*9.81*d['thickness']*np.sin(np.deg2rad(d['slope']))/1000
    stress=np.where(m,stress,np.nan)
    x=d['x'];y=d['y'];ext=[float(x.min()),float(x.max()),float(y.min()),float(y.max())]
    fig,axes=plt.subplots(1,3,figsize=(12,4.7),sharex=True,sharey=True)
    for ax,(a,label,cmap) in zip(axes,[(d['thickness'],'冰厚 (m)','Blues'),(d['v'],'历史冰川速度 (m/a)','viridis'),(stress,'重力驱动剪应力代理 (kPa)','magma')]):
        a=np.where(m,a,np.nan);im=ax.imshow(a,origin='upper',extent=ext,cmap=cmap)
        ax.set_title(label,fontsize=10);ax.set_xlabel('UTM45N 东向坐标 (m)')
        fig.colorbar(im,ax=ax,fraction=.04,pad=.03)
    axes[0].set_ylabel('UTM45N 北向坐标 (m)')
    fig.suptitle('冰川冰厚、速度与驱动剪应力分布',fontsize=12)
    savefig(fig,out/'Q2_05_glacier_physical_indicators.png')
    write_json(out/'Q2_glacier_uncertainty_note.json',{
        'stress_formula':'rho_ice*g*H*sin(theta)',
        'ice_density_assumed_kg_m3':900,
        'glacier_velocity_time':'2014-2024 climatological velocity, NOT live 2026',
        'stress_is_failure_probability':False,
        'released_fraction_identified_from_data':False})

def execute_q2(ctx):
    out=ctx.folder('q2');g=gpd.read_file(ctx.processed/'historic_gullies.gpkg')
    if '类型' not in g:raise ValueError('No real type labels in official historical gullies')
    y=(g['类型'].astype(str)!='降雨型').astype(int).to_numpy()
    X=g[FEATURES].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    if not np.isfinite(X).all():raise ValueError('Missing/non-finite geophysical features: pre-specify imputation in TRAIN folds before running')
    groups=_spatial_folds(g,ctx.cfg['q2']['cv_splits'])
    foldmetrics,res=_cross_validate(X,y,groups,ctx.cfg['q2']['C'])
    outpred=pd.DataFrame({'gully_id':g.gully_id,'true_ice_related':y,'spatial_block':groups,
         'heldout_pred_prob_ice_related':res['L2_Logistic_4geo']})
    outpred.to_csv(out/'Q2_real_label_spatial_holdout_predictions.csv',index=False,encoding='utf-8-sig')
    foldmetrics.to_csv(out/'Q2_spatial_cv_by_block.csv',index=False)
    summary={k:v for k,v in res.items() if k.endswith('_summary')}
    summary['classification_task']='P(ice-related type | known historical debris-flow gully), NOT event occurrence hazard'
    summary['n_true_ice_related']=int(y.sum());summary['n_true_rainfall']=int((1-y).sum())
    # sensitivity: iceberg lakes are a distinct causal trigger and can be kept or dropped
    small=g['类型'].astype(str).ne(BINARY_LAKE).to_numpy()
    yy=y[small];xx=X[small];gs=_spatial_folds(g.loc[small].reset_index(drop=True),ctx.cfg['q2']['cv_splits'])
    fm2,rs2=_cross_validate(xx,yy,gs,ctx.cfg['q2']['C'])
    summary['excluding_lake_sensitivity']=rs2['L2_Logistic_4geo_summary']
    write_json(out/'Q2_model_comparison.json',summary)
    _plot_training_results(ctx,g,y,groups,res['L2_Logistic_4geo'],foldmetrics)
    # Genuine historical 'severity' labels and post-event fan geometry are used
    # for outcome description, NEVER as advance prediction covariates.
    descriptive=g.groupby('类型').agg(n_gullies=('gully_id','size'),
          median_gully_drop_m=('高差','median'),
          median_fault_distance_m=('DisFrmFlt','median'),
          median_deposit_area_raw=('Area','median'),
          median_fan_width_m=('扇宽','median')).reset_index()
    descriptive.to_csv(out/'Q2_historical_gully_mechanism_evidence.csv',index=False)
    pd.crosstab(g['类型'],[g['规模'],g['频率']]).to_csv(out/'Q2_historical_scale_frequency_table.csv')
    h=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw').sort_index()
    _weather_limits(ctx,h)
    _glacier_map(ctx)
    from .extra_official_sources import plot_landcover_and_county_population,optional_glacier_thinning_2015_2020,landcover_four_year_comparison
    plot_landcover_and_county_population(ctx)
    optional_glacier_thinning_2015_2020(ctx)
    landcover_four_year_comparison(ctx)
    print('[Q2] Real-label spatial OOF AUC=',round(summary['L2_Logistic_4geo_summary']['auc_oof'],3))
    return summary
