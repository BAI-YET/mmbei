"""Additional task answers computed directly from the official hourly observations."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .common import savefig,write_json,LABELS,COLORS


def _longest_run(mask):
    best=current=0
    for value in np.asarray(mask,bool):
        current=current+1 if value else 0;best=max(best,current)
    return best


def complete_q1(ctx,hourly=None):
    out=ctx.folder('q1');issue=pd.Timestamp(ctx.cfg['snapshot_raw'])
    if hourly is None:
        hourly=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw')
    from .q1_event_signatures import physical_group
    reps=pd.read_csv(out/'Q1_mechanism_representatives.csv')
    reps['mechanism']=reps.feature.map(physical_group)
    reps['answer_rank']=np.arange(1,len(reps)+1)
    reps.to_csv(out/'Q1_answer_mechanism_top5.csv',index=False,encoding='utf-8-sig')
    event=hourly.loc[(hourly.index>=issue-pd.Timedelta(days=7))&(hourly.index<issue)]
    reference=hourly.loc[hourly.index.year.isin(ctx.cfg['historic_years'])&
       (abs(hourly.index.dayofyear-issue.dayofyear)<=ctx.cfg['match_doy_radius'])]
    fields=['temp_2m','temp_20m','vpd','sw_down','energy_balance_residual','rain_24h']
    rows=[];percentiles=pd.DataFrame(index=event.index)
    for name in fields:
        if name not in hourly:continue
        result=[];limit=[];record=[]
        for t,v in event[name].items():
            ref=reference.loc[reference.index.hour==t.hour,name].dropna().to_numpy(float)
            if not np.isfinite(v) or len(ref)<25:result.append(np.nan);limit.append(np.nan);record.append(False);continue
            result.append(float((np.sum(ref<v)+.5*np.sum(ref==v))/len(ref)))
            limit.append(float(np.quantile(ref,.95)));record.append(bool(v>ref.max()))
        p=np.array(result);percentiles[name]=p;valid=np.isfinite(p)
        rows.append({'feature':name,'label':LABELS.get(name,name),'valid_event_hours':int(valid.sum()),
          'hours_above_matched_q95':int(np.sum(p>=.95)),'hours_above_matched_q99':int(np.sum(p>=.99)),
          'new_matched_record_hours':int(np.sum(record)),
          'longest_consecutive_q95_hours':_longest_run(np.where(valid,p>=.95,False)),
          'event_max':float(event[name].max()),'event_max_time_bjt':str(event[name].idxmax()),
          'maximum_matched_percentile':float(np.nanmax(p)) if valid.any() else None})
    extremes=pd.DataFrame(rows);extremes.to_csv(out/'Q1_explicit_extreme_weather.csv',index=False,encoding='utf-8-sig')
    percentiles.to_csv(out/'Q1_extreme_hourly_percentiles.csv',index_label='time_bjt')
    daily=hourly.loc[hourly.index<issue,'temp_2m'].resample('D').agg(['mean','count'])
    daily=daily.loc[daily['count']>=18];recent=daily.loc[daily.index>=issue.normalize()-pd.Timedelta(days=30)]
    streak=0
    for t,v in recent['mean'].iloc[::-1].items():
        if not np.isfinite(v) or v<=0:break
        if t!=recent.index[-1]-pd.Timedelta(days=streak):break
        streak+=1
    report={'event_window_start':str(event.index.min()),'event_window_end':str(event.index.max()),
      'historical_match':'same calendar ±15 days and same clock hour; observed valid years only',
      'mechanism_top5':reps[['answer_rank','feature','label','mechanism','rarity_anomaly_score']].to_dict('records'),
      'raw_top5_file':'Q1_top5_with_stability.csv','extreme_summary':rows,
      'positive_mean_temperature_days_at_end_of_available_record':streak,
      'positive_temperature_streak_search_days':30,'daily_minimum_valid_hours':18}
    write_json(out/'Q1_direct_answer.json',report)
    fig,axs=plt.subplots(1,2,figsize=(12.5,5),layout='constrained')
    r=reps.iloc[::-1];axs[0].barh(r.label,r.rarity_anomaly_score,color=COLORS[0])
    axs[0].set(xlabel='历史匹配关联评分',title='a  合并同机制后的五项因素')
    im=axs[1].imshow(percentiles.T,aspect='auto',origin='lower',vmin=.5,vmax=1,cmap='YlOrRd')
    axs[1].set(yticks=np.arange(len(fields)),yticklabels=[LABELS.get(x,x) for x in fields],
      xticks=np.arange(0,len(event),24),xticklabels=[f'{event.index[i]:%m-%d}' for i in range(0,len(event),24)],
      xlabel='事件前七日（北京时间）',title='b  同季同小时历史百分位')
    fig.colorbar(im,ax=axs[1],label='历史百分位',shrink=.8)
    fig.suptitle('五项互补机制与明确的极端天气识别')
    savefig(fig,out/'Q1_14_mechanisms_and_extreme_weather.png')
    return report


def energy_water_channels(hourly,source_z,station_z,lapse=6.5):
    """Unit-consistent energy equivalents; station energy is not labelled glacier melt."""
    h=hourly.sort_index();d=pd.DataFrame(index=h.index)
    source_t=h.temp_2m-(source_z-station_z)*lapse/1000
    d['source_PDD7_Cday']=(source_t.clip(lower=0)/24).rolling('168h',min_periods=126).sum()
    d['source_PDD30_Cday']=(source_t.clip(lower=0)/24).rolling('720h',min_periods=540).sum()
    # W/m² * s / (J/kg) gives kg/m², numerically mm of water equivalent.
    for name,col in [('net_shortwave','station_net_shortwave_7d_mmwe'),
                     ('energy_balance_residual','station_positive_energy_residual_7d_mmwe')]:
        if name in h:d[col]=(h[name].clip(lower=0)*3600/334000).rolling('168h',min_periods=126).sum()
    paired=source_t.notna()&source_t.shift(1).notna()
    changes=((source_t>=0)!=(source_t.shift(1)>=0)).astype(float).where(paired)
    d['source_freeze_thaw_transitions_7d']=changes.rolling('168h',min_periods=126).sum()
    d['rain24_mm']=h.rain_24h
    return d


def complete_q2(ctx):
    out=ctx.folder('q2');issue=pd.Timestamp(ctx.cfg['snapshot_raw'])
    h=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw')
    limits=json.loads((out/'Q2_primary_trigger_assessment.json').read_text())['thresholds']
    channels=energy_water_channels(h,limits['source_elevation_m'],limits['station_elevation_m'],limits['lapse_rate_C_per_km'])
    baseline=channels.index.year.isin(ctx.cfg['historic_years'])&channels.index.month.isin([5,6,7,8,9])
    current=channels.loc[channels.index<issue].iloc[-1];rows=[];ranks=pd.DataFrame(index=channels.index)
    for name in channels:
        ref=channels.loc[baseline,name].dropna()
        if name in ('rain24_mm','source_freeze_thaw_transitions_7d'):ref=ref[ref>0]
        if len(ref)<50:continue
        threshold=float(ref.quantile(.90));v=float(current[name])
        percentile=float((np.sum(ref<v)+.5*np.sum(ref==v))/len(ref)) if np.isfinite(v) else None
        rows.append({'channel':name,'value_at_issue':v,'historical_q90':threshold,
          'historical_percentile':percentile,'exceeds_q90':bool(v>=threshold and v>0) if np.isfinite(v) else None,
          'reference_n':len(ref),'reference_years':','.join(map(str,sorted(ref.index.year.unique())))})
        # Historical empirical rank from training values only, with missing values retained.
        ordered=np.sort(ref.to_numpy());values=channels[name].to_numpy(float)
        ranks[name]=np.where(np.isfinite(values),np.searchsorted(ordered,values,side='right')/len(ordered),np.nan)
    channels.to_csv(out/'Q2_energy_water_channels.csv.gz',compression='gzip',index_label='time_bjt')
    table=pd.DataFrame(rows);table.to_csv(out/'Q2_expanded_channel_thresholds.csv',index=False)
    # A fixed, equally weighted combination is evaluated against its own historical reference.
    compound=ranks.mean(axis=1,skipna=False);ref=compound.loc[baseline].dropna()
    compound_limit=float(ref.quantile(.90));value=float(compound.loc[compound.index<issue].iloc[-1])
    report={'baseline_years':ctx.cfg['historic_years'],'reference_quantile':.9,'channels':rows,
      'channels_exceeding_reference':[x['channel'] for x in rows if x['exceeds_q90']],
      'compound_mean_rank':value,'compound_mean_rank_q90':compound_limit,
      'compound_reference_exceeded':bool(value>=compound_limit),
      'energy_conversion':'positive measured station flux W/m² ×3600 /334000 J/kg; 1kg/m²=1mm water equivalent',
      'source_temperature':'station T minus declared lapse rate times source-station elevation difference',
      'quantity_scope':'station energy potential and lapse-adjusted thermal proxies; no fabricated source melt or pore-pressure measurement',
      'threshold_selection':'fixed historical q90; not adjusted to make the focal event cross'}
    write_json(out/'Q2_expanded_mechanism_answer.json',report)
    names={'source_PDD7_Cday':'源区7日正积温','source_PDD30_Cday':'源区30日正积温',
      'station_net_shortwave_7d_mmwe':'站点净短波能量当量','station_positive_energy_residual_7d_mmwe':'站点正能量余量当量',
      'source_freeze_thaw_transitions_7d':'源区冻融转换代理','rain24_mm':'24小时降水'}
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.8),layout='constrained')
    vals=table.historical_percentile.to_numpy(float);axs[0].barh([names[x] for x in table.channel],vals,
      color=[COLORS[1] if x else COLORS[0] for x in table.exceeds_q90])
    axs[0].axvline(.9,color='#666',ls='--');axs[0].set(xlabel='同一历史基准中的百分位',xlim=(0,1.05),title='a  本次事件的各通道位置')
    pre=channels.index.to_series().between(issue-pd.Timedelta(days=30),issue-pd.Timedelta(hours=1))
    for key in ['source_PDD7_Cday','source_PDD30_Cday','station_positive_energy_residual_7d_mmwe']:
        if key in ranks:axs[1].plot(ranks.index[pre],ranks.loc[pre,key],label=names[key],lw=1.3)
    axs[1].axhline(.9,color='#666',ls='--');axs[1].set(xlabel='北京时间',ylabel='历史百分位',title='b  前期水热与能量累积');axs[1].legend(frameon=False,fontsize=8)
    fig.suptitle('降雨、持续热量、辐射与冻融通道的共同评估')
    savefig(fig,out/'Q2_18_energy_water_mechanism_channels.png')
    return report


def complete_q3(ctx,report=None):
    out=ctx.folder('q3')
    if report is None:report=json.loads((out/'Q3_primary_competition_answer.json').read_text())
    if report.get('model_family')=='Q2_conditioned_event_time':
        answer={k:report[k] for k in ['forecast_target','conditional_95pct_interval_times_bjt',
            'operational_warning_time_bjt','operational_warning_lead_h','case_reference_coverage',
            'spatial_distance_range_km','inherited_Q1_features']}
        write_json(out/'Q3_warning_lead_answer.json',answer)
        contract=json.loads((ctx.folder('q2')/'Q2_sequential_model_contract.json').read_text())
        components=pd.read_csv(ctx.folder('q2')/'Q2_Q1_feature_components.csv.gz',parse_dates=['time_bjt']).set_index('time_bjt')
        current=components.loc[components.index<pd.Timestamp(report['prediction_issue_bjt'])].iloc[-1]
        fig,axs=plt.subplots(1,2,figsize=(11.4,4.4),layout='constrained')
        features=contract['Q1_features'];weights=np.array([contract['Q1_weights'][f] for f in features])
        axs[0].barh([LABELS.get(f,f) for f in features],weights,color=COLORS[0])
        axs[0].set(xlabel='由Q1关联评分归一化的设计权重',title='a  前问结果进入后问')
        axs[1].barh([LABELS.get(f,f) for f in features],weights*current[features].to_numpy(float),color=COLORS[2])
        axs[1].set(xlabel='当前异常调制项的加权分量',title='b  发布前最新观测的贡献')
        savefig(fig,out/'Q3_11_warning_time_and_lead_distribution.png')
        return answer
    if report.get('model_family')=='mechanical_source_stability':
        fraction=report['model_failure_scenario_fraction'];validation=report['mechanics_forcing_validation']
        answer={'prediction_target':report['forecast_target'],'warning_time_bjt':report['operational_warning_time_bjt'],
            'lead_h_to_task_impact_reference':report['operational_warning_lead_h'],
            'trigger_interval_bjt':report['conditional_95pct_interval_times_bjt'],
            'model_failure_scenario_fraction':fraction,'empirical_disaster_PICP':None,
            'model_forcing_experiment':validation,
            'reference_clock':'task impact clock is not a measured source-onset time'}
        write_json(out/'Q3_warning_lead_answer.json',answer)
        fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
        axs[0].bar(['未来时域内失稳','未来时域内未失稳'],[fraction,1-fraction],color=COLORS[:2])
        axs[0].axhline(report['forecast_warning_probability_threshold'],color='#777',ls='--',label='预设情景比例阈值')
        axs[0].set(ylabel='声明参数与天气集合中的比例',ylim=(0,1),title='a  发布时刻的预警判断');axs[0].legend(frameon=False,fontsize=8)
        rows=pd.read_csv(out/'Q3_mechanics_forcing_backtest.csv')
        if len(rows):
            good=rows.n_model_reference_crossings>0
            coverage=rows.loc[good,'n_model_reference_times_covered']/rows.loc[good,'n_model_reference_crossings']
            axs[1].hist(coverage,bins=np.linspace(0,1,11),color=COLORS[2],edgecolor='white')
        axs[1].set(xlabel='各发布日的模型参照时间覆盖率',ylabel='发布日数',title='b  天气驱动检验；无真实起动时刻分母')
        fig.suptitle('预警判断与对应验证口径')
        savefig(fig,out/'Q3_11_warning_time_and_lead_distribution.png')
        return answer
    data=pd.read_csv(out/'Q3_proxy_backtest_2025.csv');threshold=report['forecast_warning_probability_threshold']
    hit=data.true_crossed.astype(bool)&data.forecast_cross_probability.ge(threshold)
    new=hit&~data.already_above_at_issue.astype(bool)
    leads=data.loc[new,'first_true_cross_h'].dropna()
    summary=report['proxy_backtest']['new_onset_warning_lead_distribution']
    write_json(out/'Q3_weather_reference_warning_lead_answer.json',{'historical_new_onset_lead':summary,
      'current_warning_time_bjt':report['operational_warning_time_bjt'],
      'current_reference_lead_h':report['operational_warning_lead_h'],
      'current_condition_probability':report['predicted_condition_cross_probability'],
      'current_no_cross_probability':report['predicted_no_condition_cross_probability'],
      'forecast_alert_threshold':threshold,'denominator':'correctly alerted new meteorological threshold passages in 2025',
      'monitoring_alarm_first_bjt':report['old_first_alarm_bjt'],
      'monitoring_alarm_active_at_issue':report['active_warning_at_issue']})
    fig,axs=plt.subplots(1,2,figsize=(10.7,4.2),layout='constrained')
    axs[0].hist(leads,bins=np.arange(0,6.51,.5),color=COLORS[0],edgecolor='white')
    if len(leads):axs[0].axvline(leads.median(),color=COLORS[1],ls='--',label=f'中位提前量 {leads.median():.2f} h')
    axs[0].set(xlabel='预报发布至新超限的提前量 (h)',ylabel='正确预警次数',title='a  独立年份预警提前量分布');axs[0].legend(frameon=False,fontsize=8)
    probability=report['predicted_condition_cross_probability']
    axs[1].bar(['未来6h条件超限','未来6h不超限'],[probability,1-probability],color=COLORS[:2])
    axs[1].axhline(threshold,color='#777',ls='--',label='预设预报预警阈值')
    axs[1].set(ylabel='轨迹频率',ylim=(0,1),title='b  本次发布的前瞻预警判断');axs[1].legend(frameon=False,fontsize=8)
    fig.suptitle('预警发布时间、提前量与检验分布')
    savefig(fig,out/'Q3_R11_weather_warning_lead_reference.png')
