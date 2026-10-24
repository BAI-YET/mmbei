"""Q1: matched event studies, no TPI self-label. Exports reproducible candidate Top5."""
from __future__ import annotations
import numpy as np,pandas as pd
from datetime import timedelta
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from .common import savefig,write_json,LABELS,COLORS,record_warning,robust_z,quantile_tail

# All admissible physics-based candidates, not every instrument diagnostic.
# Repeated heights and related moisture/energy measures are considered at this
# stage, then de-duplicated on the SAME season/window as their event ranking.
CANDIDATES=['temp_1m', 'temp_2m', 'temp_4m', 'temp_10m', 'temp_20m', 'temp_gradient', 'wind_1m', 'wind_2m', 'wind_4m', 'wind_10m', 'wind_20m', 'wind_gradient', 'vector_shear_2_20m', 'bulk_richardson', 'rh_1m', 'rh_2m', 'rh_4m', 'rh_10m', 'rh_20m', 'humidity_gradient', 'pressure', 'precip_mm', 'sw_down', 'sw_up', 'lw_down', 'lw_up', 'net_shortwave', 'net_longwave', 'net_radiation', 'albedo', 'soilT_0p1', 'soilT_0p2', 'soilT_0p4', 'soilT_0p8', 'soilT_1p6', 'soil_temp_gradient', 'soil_freeze_depth_layers', 'vwc_0p1', 'vwc_0p2', 'vwc_0p4', 'vwc_0p8', 'vwc_1p6', 'vwc_gradient', 'soil_heat', 'soil_heat_0p2', 'soil_heat_gradient', 'sensible_heat', 'latent_heat', 'tke', 'ustar', 'tau', 'tstar', 'sonic_temp', 'flux_uz', 'vpd_2m', 'vpd', 'vapour_pressure', 'water_vapour_density', 'bowen_ratio', 'energy_balance_residual', 'pdd_7d', 'rain_24h']
# Instrument diagnosis excluded from physical importance ranking: FETCH, CO2 mass flux, probe fields.

SUM_FEATURES=['precip_mm']
ENDPOINT_FEATURES=['pdd_7d','rain_24h','rain_72h','rain_7d']
DIRECTIONAL_OR_BINARY={'temperature_inversion','soil_surface_frozen','soil_freeze_depth_layers'}
MIN_REFERENCE_YEARS=3

def _window_table(hourly,window_hours,min_fraction):
    need=max(1,int(np.ceil(window_hours*min_fraction)))
    core=[c for c in CANDIDATES if c in hourly]
    roll=hourly[core].rolling(f'{window_hours}h',min_periods=need).mean()
    for c in SUM_FEATURES:
        if c in hourly:roll[c]=hourly[c].rolling(f'{window_hours}h',min_periods=need).sum()
    for c in ENDPOINT_FEATURES:
        if c in hourly:roll[c]=hourly[c].copy()
    return roll

def _get_ref_values(rolled,anchor,cfg,year_exclude=None):
    # same hour, calendar-near controls; 2026 event excluded by construction
    records=[]
    for y in cfg['historic_years']:
        if y==year_exclude:continue
        for dd in range(-cfg['match_doy_radius'],cfg['match_doy_radius']+1):
            # Year replace can fail for leap day only; event is Aug 26
            target=anchor.replace(year=y)+pd.Timedelta(days=dd)-pd.Timedelta(hours=1)
            if target in rolled.index:
                r=rolled.loc[target].astype(object).copy()
                r['reference_year']=y;r['reference_anchor']=target
                records.append(r)
    return pd.DataFrame(records)

def _rank_window(rolled,anchor,cfg,exclude_year=None):
    case_end=anchor-pd.Timedelta(hours=1)
    if case_end not in rolled.index:raise ValueError(f'No exact hourly sample available at {case_end}')
    event=rolled.loc[case_end]
    refs=_get_ref_values(rolled,anchor,cfg,exclude_year)
    rows=[];excluded=[]
    for col in CANDIDATES:
        if col not in rolled:
            excluded.append({'feature':col,'reason':'feature_not_available'});continue
        value=float(event[col])
        r=refs[col].to_numpy(dtype=float)
        r=r[np.isfinite(r)]
        ref_years=int(refs.loc[refs[col].notna(),'reference_year'].nunique())
        if not np.isfinite(value):
            excluded.append({'feature':col,'reason':'event_window_insufficient_valid_hours'});continue
        if len(r)<35 or ref_years<MIN_REFERENCE_YEARS:
            excluded.append({'feature':col,'reason':'historical_valid_windows_or_years_insufficient','n':len(r),'years':ref_years});continue
        # Binary prevalence is reported in boundary-layer figures, not ranked by
        # a continuous MAD score; zero MAD must never be converted to p≈0.
        if col in DIRECTIONAL_OR_BINARY:
            excluded.append({'feature':col,'reason':'binary_or_discrete_diagnostic'});continue
        z=robust_z(value,r);tail=quantile_tail(value,r)
        if not np.isfinite(z) or not np.isfinite(tail):
            excluded.append({'feature':col,'reason':'undefined robust variance/score','n':len(r)})
            continue
        rows.append({'feature':col,'label':LABELS.get(col,col),
                     'rarity_anomaly_score':.45*min(abs(float(z))/3,1)+.55*(1-float(tail)),
                     'event_value':value,'reference_median':float(np.median(r)),
                     'z_mad':z,'abs_z':abs(z),'two_sided_empirical_tail':tail,
                     'ref_n_windows':len(r),'ref_years':ref_years})
    out=pd.DataFrame(rows).sort_values(['rarity_anomaly_score','abs_z'],ascending=[False,False])
    out.attrs['excluded_candidates']=excluded
    return out,refs

def physical_group(feature):
    """Do not select two sensors estimating the same physical process as
    independent top factors, even if their historical Spearman < 0.85.
    """
    if feature in ('vpd','vpd_2m','water_density','vapour_pressure','water_vapour_density','humidity_gradient') or feature.startswith('rh_'):
        return 'atmospheric_moisture'
    if feature.startswith('temp_') or feature in ('temp_gradient','sonic_temp','pdd_7d'):
        return 'air_temperature'
    if feature in ('sw_down','sw_up','lw_down','lw_up','net_radiation','net_shortwave','net_longwave','albedo'):
        return 'radiation_energy'
    if feature.startswith('vwc_'):
        return 'soil_moisture'
    if feature.startswith('soilT_') or feature.startswith('soil_temperature') or feature in ('soil_heat','soil_heat_0p2','soil_temp_gradient','soil_heat_gradient'):
        return 'soil_thermal'
    if feature in ('tke','ustar','tau','tstar','sensible_heat','latent_heat',
                   'turbulent_total_heat','bowen_ratio','z_over_L_proxy'):
        return 'turbulence_exchange'
    if feature.startswith('wind_') or feature in ('vector_wind_shear','vector_shear_2_20m'):
        return 'wind_structure'
    return feature

def _distinct_top5(rank,historic,threshold=.85):
    selected=[]
    for row in rank.itertuples():
        if len(selected)>=5:break
        redundant=any(physical_group(row.feature)==physical_group(past) for past in selected)
        for past in selected:
            a=historic[row.feature];b=historic[past]
            mask=a.notna() & b.notna()
            if mask.sum()<25:continue
            rho=spearmanr(a[mask],b[mask]).statistic
            if np.isfinite(rho) and abs(rho)>=threshold:redundant=True;break
        if not redundant:selected.append(row.feature)
    return selected

def _plot_qc(out,ctx):
    # Compute available fraction within each source's actual *operational* span;
    # 2021 pre-April and 2026 post-August are not monitoring losses.
    from .prepare_time import _RAW_FIELDS
    src=ctx.raw/'unpacked'/'1'
    matrix=pd.DataFrame(index=list(_RAW_FIELDS),columns=range(2021,2027),dtype=float)
    log=[]
    for source,cols in _RAW_FIELDS.items():
        for yr in range(2021,2027):
            f=src/f'{source}_Gyirong_{yr}.csv'
            raw=pd.read_csv(f,low_memory=False)
            dat=pd.to_datetime(raw['Timestamp'],errors='coerce')
            valid_dates=dat.dropna()
            if valid_dates.empty:continue
            # Operational window is delimited by availability of at least one physical field,
            # not by a full-year file padded with missing sentinel rows.
            phys=[c for c in cols if c in raw.columns]
            numeric=raw[phys].apply(pd.to_numeric,errors='coerce').replace([9999,9999.9,-9999,-9999.9],np.nan)
            good=(numeric.notna() & numeric.abs().lt(9999)).any(axis=1)
            available=dat.loc[good]
            if available.empty:continue
            start_t=available.min();end_t=available.max()
            span=((end_t-start_t).total_seconds()/3600)+1
            idx=dat.between(start_t,end_t)
            vals=numeric.loc[idx].notna().mean(axis=1)
            fraction=float(vals.sum()/span)
            matrix.loc[source,yr]=min(1.0,fraction)
            log.append({'source':source,'year':yr,'start':str(start_t),'end':str(end_t),
                        'expected_hours_within_operational_span':span,'mean_valid_field_fraction':fraction})
    pd.DataFrame(log).to_csv(out/'Q1_effective_coverage_by_year.csv',index=False)
    fig,ax=plt.subplots(figsize=(8.5,3.4),constrained_layout=True)
    im=ax.imshow(matrix.to_numpy(float),cmap='YlGnBu',vmin=0,vmax=1,aspect='auto')
    ax.set(yticks=np.arange(len(matrix)),yticklabels=matrix.index,
           xticks=np.arange(len(matrix.columns)),xticklabels=matrix.columns,
           title='逐小时观测有效数据比例')
    for i in range(len(matrix)):
        for j in range(len(matrix.columns)):
            v=matrix.iloc[i,j]
            if np.isfinite(v):ax.text(j,i,f'{v:.0%}',ha='center',va='center',fontsize=8,
                                      color='white' if v>.65 else '#243342')
    fig.colorbar(im,ax=ax,label='数据有效率',shrink=.8)
    savefig(fig,out/'Q1_00_data_coverage.png')

def _plot_timeseries(hourly,ctx):
    out=ctx.folder('q1'); anchor=pd.Timestamp(ctx.cfg['snapshot_raw'])
    hist=hourly.loc[(hourly.index>=anchor-pd.Timedelta(days=10)) & (hourly.index<anchor)].copy()
    cols=[('temp_2m','2 m 气温 (°C)'),('net_radiation','净辐射 (W/m²)'),('vwc_0p1','表层体积含水率'),('precip_mm','降水量 (mm/h)')]
    fig,axes=plt.subplots(4,1,figsize=(11,8.4),sharex=True,gridspec_kw={'hspace':0.26})
    for k,((c,lab),ax) in enumerate(zip(cols,axes)):
        if c not in hist:continue
        data=hist[c]
        if c=='precip_mm':
            ax.bar(data.index,data,width=.028,color=COLORS[k],alpha=.75)
        else:ax.plot(data.index,data,lw=1.3,color=COLORS[k])
        ax.set_ylabel(lab,fontsize=9);ax.grid(True,axis='y',alpha=.18)
        ax.axvline(anchor,color='#CB5564',ls='--',lw=1.5)
    axes[0].set_title('灾前十天气象要素变化过程')
    axes[-1].set_xlabel('北京时间（UTC+8）')
    fig.autofmt_xdate();savefig(fig,out/'Q1_01_event_hourly_timeseries.png')

def _plot_profiles(hourly,ctx):
    anchor=pd.Timestamp(ctx.cfg['snapshot_raw']);hh=hourly.loc[hourly.index<anchor]
    month=hh.loc[hh.index.month.isin([6,7,8])]
    hourly_means=month.groupby(month.index.hour).mean(numeric_only=True)
    hours=[1,2,4,10,20];cmap=plt.cm.viridis
    fig,axes=plt.subplots(1,2,figsize=(9,4.7))
    for hr,color in zip([2,8,14,20],[cmap(0.08),cmap(.35),cmap(.65),cmap(.94)]):
        for ax,key,xlabel in [(axes[0],'temp','温度 (°C)'),(axes[1],'wind','风速 (m/s)')]:
            vs=[hourly_means.loc[hr,f'{key}_{h}m'] if f'{key}_{h}m' in hourly_means else np.nan for h in hours]
            ax.plot(vs,hours,marker='o',ms=3,lw=1.5,color=color,label=f'{hr:02d}:00')
            ax.set(xlabel=xlabel,ylabel='高度 (m)');ax.grid(True,alpha=.17)
    axes[0].legend(title='站点小时',frameon=False,loc='best');fig.suptitle('夏季近地层气象要素垂直廓线')
    savefig(fig,ctx.folder('q1')/'Q1_02_boundary_layer_profiles.png')

def _plot_seasonal_structure_and_extremes(hourly,ctx):
    from scipy.stats import rankdata
    out=ctx.folder('q1')
    h=hourly.loc[hourly.index.year.isin(ctx.cfg['historic_years'])]
    h=h.assign(season=np.where(h.index.month.isin([6,7,8]),'暖季',
                      np.where(h.index.month.isin([12,1,2]),'冬季','其他')))
    summaries=[]
    for name,d in h.groupby('season'):
        if name=='其他':continue
        tgrad=d['temp_20m']-d['temp_2m']
        shear=d['wind_20m']-d['wind_2m']
        summaries.append({'season':name,'valid_inversion_hours':int(tgrad.notna().sum()),
            'inversion_frequency':float((tgrad>0).sum()/max(1,tgrad.notna().sum())),
            'median_T20_minus_T2_C':float(tgrad.median()),
            'median_shear_m_s':float(shear.median()),
            'mean_net_radiation_Wm2':float(d['net_radiation'].mean())})
    pd.DataFrame(summaries).to_csv(out/'Q1_boundary_layer_seasonal_metrics.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(10.6,4.1),layout='constrained')
    for season,color in [('冬季','#6586A8'),('暖季','#D67A51')]:
        d=h[h.season.eq(season)]
        inv=(d.temp_20m>d.temp_2m).astype(float).where((d.temp_20m-d.temp_2m).notna())
        vals=inv.groupby(inv.index.hour).mean().reindex(range(24))
        axes[0].plot(range(24),100*vals.to_numpy(),label=season,lw=2,color=color)
        sh=(d.wind_20m-d.wind_2m).groupby(d.index.hour).median().reindex(range(24))
        axes[1].plot(range(24),sh.to_numpy(),label=season,lw=2,color=color)
    axes[0].set(xlabel='北京时间 (h)',ylabel='逆温出现频率 (%)',title='a  冬夏逆温日变化',ylim=(0,100))
    axes[1].set(xlabel='北京时间 (h)',ylabel='20–2 m风速差 (m/s)',title='b  近地层风切变日变化')
    for ax in axes:ax.grid(axis='y',alpha=.16);ax.legend(frameon=False)
    fig.suptitle('大气边界层季节结构特征',fontsize=13)
    savefig(fig,out/'Q1_08_boundary_layer_seasonal_structure.png')
    # Same-calendar-period extremes: comparison uses actual independently observed dates.
    anchor=pd.Timestamp(ctx.cfg['snapshot_raw']);hist=h[h.index.month.isin([8])]
    hbase=hist.groupby(hist.index.hour)[['temp_2m','precip_mm','sw_down','wind_20m']].quantile(.95)
    d=hourly.loc[(hourly.index>=anchor-pd.Timedelta(days=7))&(hourly.index<anchor),list(hbase.columns)]
    ex=[]
    for col in hbase:
        high=d[col] >= d.index.hour.map(hbase[col])
        ex.append({'feature':col,'n_valid_hours':int(d[col].notna().sum()),
                   'n_above_historical_same_hour_q95':int(high.sum()),
                   'fraction_of_valid_hours':float(high.sum()/max(1,d[col].notna().sum()))})
    pd.DataFrame(ex).to_csv(out/'Q1_weather_extreme_hours.csv',index=False)

def execute_q1(ctx):
    out=ctx.folder('q1');hourly=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw').sort_index()
    anchor=pd.Timestamp(ctx.cfg['snapshot_raw'])
    ranks={};refs_primary=None
    for w in ctx.cfg['windows_hours']:
        roll=_window_table(hourly,w,ctx.cfg['min_valid_fraction'])
        r,ref=_rank_window(roll,anchor,ctx.cfg)
        r['window_hours']=w;ranks[w]=r
        if w==ctx.cfg['ranking_primary_window_hours']:refs_primary=ref
        r.to_csv(out/f'Q1_window_{w:03d}h_ranked_anomaly.csv',index=False,encoding='utf-8-sig')
    primary=ctx.cfg['ranking_primary_window_hours'];ranking=ranks[primary].copy()
    pd.DataFrame(ranks[primary].attrs.get('excluded_candidates',[])).to_csv(out/'Q1_excluded_candidates.csv',index=False)
    # Explicitly report all candidates not assessed (insufficient QA, years, or
    # binary/diagnostic property), rather than silently dropping them.
    qc_inventory=set(hourly.columns)
    included=set(ranking.feature)
    audits=[]
    omitted={d['feature']:d['reason'] for d in ranks[primary].attrs.get('excluded_candidates',[])}
    for name in CANDIDATES:
        if name in included:reason='ranked'
        elif name not in qc_inventory:reason='not_in_official_processed_features'
        elif name in DIRECTIONAL_OR_BINARY:reason='binary_prevalence_diagnostic'
        else:reason=omitted.get(name,'unclassified_exclusion')
        audits.append({'feature':name,'label':LABELS.get(name,name),'status':reason})
    pd.DataFrame(audits).to_csv(out/'Q1_candidate_inclusion_exclusion.csv',index=False,encoding='utf-8-sig')
    selected=_distinct_top5(ranking, refs_primary,ctx.cfg['ranking_redundancy_spearman'])
    ranking['selected_mechanism_representative']=ranking.feature.isin(selected)
    ranking['selected_top5']=np.arange(len(ranking))<5
    ranking['rank']=np.arange(1,len(ranking)+1)
    ranking.to_csv(out/'Q1_all_candidate_ranking.csv',index=False,encoding='utf-8-sig')
    representatives=ranking[ranking.selected_mechanism_representative].copy()
    top5=ranking.head(5).copy()
    for w,r in ranks.items():
        positions={f:i+1 for i,f in enumerate(r.feature)}
        top5[f'rank_{w}h']=top5.feature.map(positions)
    top5.to_csv(out/'Q1_top5_with_stability.csv',index=False,encoding='utf-8-sig')
    representatives.to_csv(out/'Q1_mechanism_representatives.csv',index=False,encoding='utf-8-sig')
    from .task_evidence import q1_complete_evidence
    q1_complete_evidence(ctx,hourly,ranking,ranks,refs_primary,audits)
    # leave-one-reference-year-out sensitivity on the actual event
    rlo=[]
    roll=_window_table(hourly,primary,ctx.cfg['min_valid_fraction'])
    for y in ctx.cfg['historic_years']:
        r,_=_rank_window(roll,anchor,ctx.cfg,exclude_year=y)
        for z in r.itertuples():rlo.append({'omitted_reference_year':y,'feature':z.feature,'rank':int(np.where(r.feature==z.feature)[0][0]+1),'abs_z':z.abs_z})
    pd.DataFrame(rlo).to_csv(out/'Q1_leave_one_year_out_rank.csv',index=False)
    _plot_qc(out,ctx);_plot_timeseries(hourly,ctx);_plot_profiles(hourly,ctx)
    _plot_seasonal_structure_and_extremes(hourly,ctx)
    _plot_full_energy_soil_structure(hourly,ctx)
    from .extra_official_sources import plot_county_temperature_climate
    plot_county_temperature_climate(ctx)
    # primary anomaly bar chart
    r=ranking.head(min(15,len(ranking))).iloc[::-1]
    fig,ax=plt.subplots(figsize=(9.5,6.5))
    colors=['#D17744' if c in selected else '#39789A' for c in r.feature]
    ax.barh(r.label,r.rarity_anomaly_score,color=colors,alpha=.92)
    ax.axvline(0,lw=.9,color='#555')
    ax.set_xlabel('历史匹配异常稀有性评分（0–1）')
    ax.set_title('灾前多源气象因素异常稀有性排序')
    ax.grid(axis='x',alpha=.17);savefig(fig,out/'Q1_03_top15_anomaly_ranking.png')
    heat_cols=[f for f in ranking.feature.head(15)]
    hm=np.full((len(heat_cols),len(ranks)),np.nan)
    for j,w in enumerate(ranks):
        m=ranks[w].set_index('feature').z_mad
        for i,feat in enumerate(heat_cols):
            if feat in m:hm[i,j]=m[feat]
    fig,ax=plt.subplots(figsize=(8.2,6.5));bound=np.nanpercentile(np.abs(hm),95);im=ax.imshow(hm,cmap='RdBu_r',vmin=-bound,vmax=bound,aspect='auto')
    ax.set(yticks=np.arange(len(heat_cols)),yticklabels=[LABELS.get(q,q) for q in heat_cols],xticks=np.arange(len(ranks)),xticklabels=[str(k)+' h' for k in ranks],title='不同时间窗口的气象特征异常分布')
    fig.colorbar(im,ax=ax,label='MAD Z',shrink=.8)
    savefig(fig,out/'Q1_04_multiscale_anomaly_heatmap.png')
    yrtable=pd.DataFrame(rlo).pivot_table(index='feature',columns='omitted_reference_year',values='rank')
    subset=yrtable.reindex(selected)
    fig,ax=plt.subplots(figsize=(8.2,3.5))
    im=ax.imshow(subset.values,cmap='YlGnBu_r',vmin=1,vmax=max(10,np.nanmax(subset.values)),aspect='auto')
    ax.set(yticks=range(len(selected)),yticklabels=[LABELS.get(f,f) for f in selected],xticks=range(len(ctx.cfg['historic_years'])),xticklabels=[f'不含{y}' for y in ctx.cfg['historic_years']],title='关键气象特征排序的年份敏感性')
    for i in range(len(selected)):
        for j in range(len(ctx.cfg['historic_years'])):
            if np.isfinite(subset.values[i,j]):ax.text(j,i,int(subset.values[i,j]),ha='center',va='center',fontsize=9)
    fig.colorbar(im,ax=ax,label='名次',shrink=.8)
    savefig(fig,out/'Q1_05_reference_year_sensitivity.png')
    # Pseudo-event controls: same calendar season, each control year compared with OTHER years.
    # As Top5 was selected using the real event, this is a descriptive selection-biased diagnostic,
    # NOT a valid hypothesis-testing p-value or supervised disaster forecast.
    case_map=ranking.set_index('feature').abs_z
    case_score=float(np.nanmean([case_map.get(f,np.nan) for f in selected]))
    pseudo=[]
    # Only pre-selected independent factors are needed for a descriptive control
    # score. Precompute each control-year baseline once; do not re-rank all 70+
    # predictors thousands of times (and this is NOT a holdout significance test).
    yearref={}
    for yy in ctx.cfg['historic_years']:
        yearref[yy]=rolled_ref=_get_ref_values(roll,anchor,ctx.cfg,year_exclude=yy)
    ctrl=(anchor-pd.Timedelta(hours=1))
    for yy in ctx.cfg['historic_years']:
        ref=yearref[yy]
        for dd in range(-ctx.cfg['match_doy_radius'],ctx.cfg['match_doy_radius']+1):
            pend=ctrl.replace(year=yy)+pd.Timedelta(days=dd)
            if pend not in roll.index:continue
            values=[]
            for feat in selected:
                if feat not in roll or feat not in ref:continue
                z=robust_z(float(roll.loc[pend,feat]),ref[feat].to_numpy(float))
                if np.isfinite(z):values.append(abs(z))
            if len(values)>=4:
                pseudo.append({'pseudo_anchor':str(pend+pd.Timedelta(hours=1)),
                               'control_year':yy,'mean_selected_abs_z':float(np.mean(values)),
                               'matched_feature_n':len(values)})
    if pseudo:
        pc=pd.DataFrame(pseudo)
        pc.to_csv(out/'Q1_pseudo_event_control_scores.csv',index=False)
        frac=float(np.mean(pc.mean_selected_abs_z>=case_score))
        fig,ax=plt.subplots(figsize=(8.4,4.1))
        ax.hist(pc.mean_selected_abs_z,bins=18,color='#8AA7B6',edgecolor='white',alpha=.85)
        ax.axvline(case_score,color='#C45C41',lw=2,ls='--',label=f'真实事件对照得分 {case_score:.2f}')
        ax.set(xlabel='入选五项的平均 |MAD Z|',ylabel='历史伪事件窗口数',
               title='历史同期对照窗口异常得分分布')
        ax.legend(frameon=False);savefig(fig,out/'Q1_06_pseudo_event_controls.png')
        write_json(out/'Q1_pseudo_event_caution.json',{'historical_control_windows':int(len(pc)),
            'event_selected_feature_score':case_score,
            'fraction_of_control_windows_at_least_as_high':frac,
            'valid_hypothesis_test_pvalue':False,
            'reason':'Features selected on focal event; pseudo windows overlap and are temporally dependent.'})
    write_json(out/'Q1_method_status.json',{
       'analysis':'historical matched magnitude+empirical rarity, NOT supervised disaster importance',
       'source_clock_verified':True,
       'snapshot_raw':str(anchor),'primary_window_h':primary,'baseline_years':ctx.cfg['historic_years'],
       'n_candidate_predefined':len(CANDIDATES),'n_features':len(ranking),'top5':[{'feature':r.feature,'name':r.label,'z_mad':r.z_mad} for r in top5.itertuples()],
       'mechanism_representatives':selected,'top5_policy':'literal first five of all-factor ranking; mechanism representatives are a separate explanatory result',
       'caveat':'Only one confirmed disaster; this is physical association and rarity, not causal identification.'})
    print('[Q1] Raw score Top5:',top5.feature.tolist(),'; mechanism representatives:',selected)
    return top5


def _plot_full_energy_soil_structure(hourly,ctx):
    """Paper figures: source-specific turbulent/surface energy and five soil depths.
    These represent the official station, not measured glacier melt energy.
    """
    out=ctx.folder('q1');snap=pd.Timestamp(ctx.cfg['snapshot_raw'])
    h=hourly.loc[(hourly.index>=snap-pd.Timedelta(days=7))&(hourly.index<snap)]
    keys=[('net_shortwave','净短波'),('net_longwave','净长波'),
          ('sensible_heat','感热 H'),('latent_heat','潜热 LE'),('soil_heat','土壤热通量 G')]
    fig,axs=plt.subplots(2,1,figsize=(10.8,6.4),sharex=True,layout='constrained')
    for idx,(key,name) in enumerate(keys):
        if key not in h:continue
        axs[0 if idx<2 else 1].plot(h.index,h[key],lw=1.2,color=COLORS[idx],label=name)
    axs[0].set(ylabel='辐射通量 (W/m²)',title='a  净短波与净长波辐射')
    axs[1].set(xlabel='北京时间（UTC+8）',ylabel='热通量 (W/m²)',title='b  湍流与土壤热量交换')
    for ax in axs:ax.legend(ncol=3,fontsize=8,frameon=False,loc='upper left');ax.grid(axis='y',alpha=.16)
    fig.suptitle('灾前地表辐射与热量通量变化')
    savefig(fig,out/'Q1_09_radiation_and_turbulent_energy.png')
    depths=['0p1','0p2','0p4','0p8','1p6'];z=[.1,.2,.4,.8,1.6]
    fig,axs=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    # use 24h averages, suppress unreliable depth cells
    t=h.tail(24)
    for ax,prefix,title,unit in [(axs[0],'soilT_','a  土壤温度','°C'),
                                 (axs[1],'vwc_','b  体积含水率','m³/m³')]:
        values=[t[f'{prefix}{k}'].mean() if f'{prefix}{k}' in t else np.nan for k in depths]
        ax.plot(values,z,'-o',color=COLORS[0 if prefix=='soilT_' else 2],lw=2)
        ax.set(xlabel=unit,ylabel='土壤深度 (m)',title=title)
        ax.set_ylim(1.7,0);ax.grid(alpha=.18)
    fig.suptitle('灾前土壤多深度温湿廓线')
    savefig(fig,out/'Q1_10_multidepth_soil_profiles.png')
    # Instrument completeness by family explicitly grounded in cleaned data
    groups={'MET':['temp_20m','wind_20m','rh_20m','vector_wind_shear'],
            'RADM':['sw_down','net_radiation'],
            'SOIL':['vwc_0p1','vwc_1p6','soilT_0p1','soilT_1p6'],
            'FLUX':['tau','ustar','tke','sensible_heat','latent_heat']}
    rows=[]
    for group,cols in groups.items():
        for col in cols:
            if col in hourly:
                d=hourly[col].loc[(hourly.index.year.isin(ctx.cfg['historic_years'])) & hourly.index.month.isin([8])]
                rows.append({'group':group,'feature':col,'n_valid_august_hours':int(d.notna().sum()),
                             'n_valid_historical_years':int(d.dropna().index.year.nunique())})
    pd.DataFrame(rows).to_csv(out/'Q1_physical_family_evidence_coverage.csv',index=False)
