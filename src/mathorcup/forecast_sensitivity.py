"""Compare actual Q2 archived-weather forcing against the validated analogue reference.
No claim that one available forecast run has itself been externally validated.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .common import savefig,write_json,COLORS


def archived_forecast_sensitivity(ctx,h,baseline_ensemble,leads,issue,threshold):
    from .q3_task_aligned import anomaly_score,forecast_analogue,BASELINE_VARS
    path=ctx.folder('q2')/'Q2_archived_forecast_trigger_conditions.csv'
    if not path.exists():
        write_json(ctx.folder('q3')/'Q3_archived_forecast_sensitivity.json',{'status':'UNAVAILABLE'});return
    f=pd.read_csv(path,parse_dates=['forecast_time_bjt']).set_index('forecast_time_bjt')
    q1,_=anomaly_score(h,ctx.cfg['historic_years'],variables=BASELINE_VARS,mode='q1_only')
    q1ensemble,_,n=forecast_analogue(q1.combined_condition_score.loc[q1.index<issue],issue,n_samples=len(baseline_ensemble),seed=ctx.cfg['q3']['random_seed'])
    hours=(f.index-issue).total_seconds()/3600
    good=f.trigger_condition_forecast.notna()&(hours>=0)&(hours<=6)
    if good.sum()<6:return
    q2=np.interp(leads,hours[good],f.loc[good,'trigger_condition_forecast'])
    joint=.55*q1ensemble+.45*q2[None,:]
    d=pd.DataFrame({'lead_h':leads,'analogue_reference_median':np.median(baseline_ensemble,axis=0),
        'Q2_actual_archived_forecast_condition':q2,'archived_weather_joint_median':np.median(joint,axis=0),
        'archived_weather_joint_p025':np.quantile(joint,.025,axis=0),'archived_weather_joint_p975':np.quantile(joint,.975,axis=0)})
    d.to_csv(ctx.folder('q3')/'Q3_archived_forecast_vs_analogue.csv',index=False)
    report={'status':'ACTUAL_Q2_FORECAST_FORCING_SENSITIVITY','baseline_method':'all-condition trajectory analogues',
       'alternative_method':'.55*Q1 anomaly analogues + .45*Q2 archived-weather condition',
       'forecast_threshold':threshold,'alternative_proxy_exceedance_probability':float((joint>=threshold).any(axis=1).mean()),
       'reference_proxy_exceedance_probability':float((baseline_ensemble>=threshold).any(axis=1).mean()),
       'not_disaster_probability':True,'archived_run_out_of_time_accuracy_evaluable':False,
       'limitation':'one run only; model/site bias and future weather uncertainty not calibrated; retain analogue as validation reference'}
    write_json(ctx.folder('q3')/'Q3_archived_forecast_sensitivity.json',report)
    fig,ax=plt.subplots(figsize=(10.2,4.4),layout='constrained')
    ax.plot(leads,d.analogue_reference_median,'-o',color=COLORS[0],label='历史相似日基准')
    ax.plot(leads,d.archived_weather_joint_median,'-s',color=COLORS[1],label='Q2真实历史天气预报驱动')
    ax.fill_between(leads,d.archived_weather_joint_p025,d.archived_weather_joint_p975,color=COLORS[1],alpha=.18)
    ax.axhline(threshold,color='#777',ls='--',label='气象先兆阈值')
    ax.set(xlabel='距09:00发布时刻 (h)',ylabel='气象触发条件指数',title='实际历史天气预报与相似日外推的结果比较');ax.legend(frameon=False,ncol=3,fontsize=8)
    savefig(fig,ctx.folder('q3')/'Q3_10_archived_forecast_sensitivity.png')
    return {'ensemble':joint,'report':report}
