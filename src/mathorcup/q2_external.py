"""Q2-only external observation / archived-forecast ingestion.

Nothing here is silently synthesised. Historical forecast runs have a UTC INIT time
and may become available hours later; when unavailable, explicit status is recorded.
Other problems access only processed Q2 outputs, never pull an extra dataset.
"""
from __future__ import annotations
import json, hashlib
from pathlib import Path
from datetime import datetime, timezone, timedelta
from urllib.parse import urlencode
import numpy as np, pandas as pd
from .common import write_json

ICIMOD_URL='https://www.icimod.org/press-releases/second-major-himalayan-flood-in-seven-weeks-raises-new-questions-about-mountain-risk/'
ECMWF_URL='https://open-meteo.com/en/docs/single-runs-api'
ZENODO_URL='https://zenodo.org/records/22842744'
# Real reported Sentinel-1 *summary*, not an invented 2D displacement raster.
SATELLITE_OBS={
  'product':'Sentinel-1 InSAR LOS deformation, ICIMOD published synthesis',
  'source_url':ICIMOD_URL,
  'measurement_period_start':'2026-07-07', 'measurement_period_end':'2026-08-24',
  'los_max_displacement_m':0.40, 'duration_days':48,
  'average_los_motion_proxy_mm_day':round(400/48,3),
  'quantity':'maximum line-of-sight deformation in affected northern slope, not vertical fall or a grid',
  'research_role':'retrospective external Q2 corroboration only; published after disaster',
  'publication_date':'2026-10-13','available_at_2026_issue':False,
  'cannot_create_pixel_displacement':True
}

FORECAST_VARS=('temperature_2m','precipitation','shortwave_radiation','freezing_level_height')
# 2026-08-25 18:00 UTC + ~6h distribution -> available no later than 08:00 Beijing on Aug 26.
# 00 UTC Aug 26 would NOT be available by 09:00 Beijing for ECMWF; avoid hindsight.
DEFAULT_RUN_UTC='2026-08-25T18:00'
DEFAULT_ISSUE_BJT='2026-08-26 09:00:00'

def q2_forecast_url(lat=28.2877,lon=85.52723,run_utc=DEFAULT_RUN_UTC):
    args={'latitude':lat,'longitude':lon,'hourly':','.join(FORECAST_VARS),'models':'ecmwf_ifs',
          'run':run_utc,'timezone':'UTC','forecast_hours':48,'temperature_unit':'celsius',
          'precipitation_unit':'mm'}
    return 'https://single-runs-api.open-meteo.com/v1/forecast?'+urlencode(args)

def download_q2_archived_forecast(ctx,run_utc=DEFAULT_RUN_UTC,timeout=40):
    """Explicitly invoked from CLI only. Do not call on ordinary Q1/Q3/Q4 runs."""
    import requests
    output=ctx.root/'data/external/q2_forecast';output.mkdir(parents=True,exist_ok=True)
    url=q2_forecast_url(run_utc=run_utc)
    resp=requests.get(url,timeout=timeout,headers={'User-Agent':'MathorCupB-scientific-reproducibility/1.6'})
    resp.raise_for_status()
    data=resp.json()
    if data.get('error') or not data.get('hourly'):raise ValueError(f'Forecast API error: {data}')
    hourly=data['hourly'];requested=set(FORECAST_VARS)
    if not requested.issubset(hourly):raise ValueError('Missing requested forecast columns: '+str(requested-set(hourly)))
    t=pd.to_datetime(hourly['time'],utc=True)
    if len(t)==0:raise ValueError('Empty forecast times')
    raw=output/'ecmwf_source_pre_event_single_run.json';raw.write_text(json.dumps(data,ensure_ascii=False),encoding='utf-8')
    cols={k:hourly[k] for k in FORECAST_VARS}
    df=pd.DataFrame(cols);df.insert(0,'time_utc',t);df['time_bjt']=df.time_utc.dt.tz_convert('Asia/Shanghai')
    df.to_csv(output/'source_forecast_hourly.csv',index=False)
    sha=hashlib.sha256(raw.read_bytes()).hexdigest()
    write_json(output/'provenance.json',{'provider':'Open-Meteo (archived ECMWF initialised run)',
          'docs':ECMWF_URL,'download_url':url,'run_initialised_utc':run_utc,
          'minimum_scheduled_latency_h':4,'conservative_publication_latency_h':6,
          'source_coordinates':[28.2877,85.52723], 'raw_sha256':sha,
          'retrieved_utc':datetime.now(timezone.utc).isoformat(),
          'do_not_use_as_observed_station_data':True})
    return df

def load_q2_forecast(ctx,issue_bjt=DEFAULT_ISSUE_BJT):
    root=ctx.root/'data/external/q2_forecast';path=root/'source_forecast_hourly.csv';meta=root/'provenance.json'
    if not path.exists() or not meta.exists():return None,'NOT_DOWNLOADED'
    m=json.loads(meta.read_text()); run=pd.Timestamp(m['run_initialised_utc'],tz='UTC')
    issue=pd.Timestamp(issue_bjt).tz_localize('Asia/Shanghai').tz_convert('UTC')
    # The issue-time check is conservative, model run available ~4-6h after initialisation.
    if run+pd.Timedelta(hours=m.get('conservative_publication_latency_h',6))>issue:
        return None,'RUN_NOT_AVAILABLE_AT_ISSUANCE'
    data=pd.read_csv(path,parse_dates=['time_utc','time_bjt'])
    if not set(FORECAST_VARS).issubset(data.columns):return None,'FORECAST_VARIABLE_MISSING'
    return data,'ACTUAL_DOWNLOADED_ARCHIVED_FORECAST'

def q2_external_status(ctx):
    base=ctx.root/'data/external';base.mkdir(parents=True,exist_ok=True)
    obs=base/'q2_satellite'/'ICIMOD_Sentinel1_2026_reported_measurement.json'
    obs.parent.mkdir(parents=True,exist_ok=True)
    write_json(obs,SATELLITE_OBS)
    forecast,code=load_q2_forecast(ctx)
    report={'external_data_scope':'QUESTION_2_ONLY','satellite_status':'REPORTED_MEASUREMENT_FROM_PUBLISHED_SENTINEL_1_ANALYSIS',
       'satellite_summary':SATELLITE_OBS,
       'weather_historical_run_status':code,
       'forecast_retrieval_command':'python run.py --fetch-q2-external',
       'forecast_url':q2_forecast_url(),
       'forecast_publication_timing':'2026-08-25 18:00 UTC init + 6h conservative latency before Beijing 09:00 issuance',
       'post_event_8m_sar_dem_dataset_url':ZENODO_URL,
       'post_event_dem_use':'Q2 / Q4 retrospective geometric checking only; never forecasting input',
       'forecast_data_rows':int(len(forecast)) if forecast is not None else 0,
       'interpolation_not_performed_if_missing':True}
    write_json(ctx.folder('q2')/'Q2_external_data_provenance.json',report)
    return report


def source_forecast_conditions(ctx,h,limits,forecast):
    """Propagate actual archived temperatures/rain through Q2's condition equations.
    Pre-issue station history + future source grid form a declared spatial proxy.
    Missing freezing-level output is retained, never replaced by zero.
    """
    from .q2_trigger_evaluation import physical_exceedance
    raw=ctx.root/'data/external/q2_forecast/ecmwf_source_pre_event_single_run.json'
    payload=json.loads(raw.read_text())
    model_z=float(payload['elevation']);source_z=limits['source_elevation_m'];lapse=limits['lapse_rate_C_per_km']
    issue=pd.Timestamp(ctx.cfg['snapshot_raw'])
    f=forecast.copy();f.index=pd.to_datetime(f.time_bjt,utc=True).dt.tz_convert('Asia/Shanghai').dt.tz_localize(None)
    f=f.loc[(f.index>=issue)&(f.index<issue+pd.Timedelta(hours=48))].copy()
    if f.empty:return f,{'status':'NO_FUTURE_HOURS'}
    f['source_temp_forecast_lapse_C']=f.temperature_2m-(source_z-model_z)*lapse/1000
    past=h.loc[h.index<issue].tail(168)
    temp_history=past.temp_2m-(source_z-limits['station_elevation_m'])*lapse/1000
    combined=pd.concat([temp_history,f.source_temp_forecast_lapse_C]).sort_index()
    future_pdd=(combined.clip(lower=0)/24).rolling('168h',min_periods=126).sum().reindex(f.index)
    rain_history=pd.concat([past.precip_mm,f.precipitation]).sort_index()
    future_rain=rain_history.rolling('24h',min_periods=18).sum().reindex(f.index)
    # Use the SAME fitted response scale as Q2's historical melt channel.
    all_source=h.temp_2m-(source_z-limits['station_elevation_m'])*lapse/1000
    pdd=(all_source.clip(lower=0)/24).rolling('168h',min_periods=126).sum()
    base=h.index.year.isin(limits['baseline_years'])&h.index.month.isin([5,6,7,8,9])
    scale=max(.05,float(pdd.loc[base].std())*.3)
    f['source_pdd7_forecast_Cday']=future_pdd
    f['hybrid_rain24_forecast_mm']=future_rain
    f['melt_gate_forecast']=physical_exceedance(future_pdd,limits['source_pdd7_warmseason_q90_Cday'],scale)
    r=future_rain;th=limits['rain_24h_nonzero_q90_mm']
    f['rain_gate_forecast']=np.where(r.notna(),np.where((r>=th)&(r>.01),physical_exceedance(r,th),0),np.nan)
    f['trigger_condition_forecast']=np.maximum(f.melt_gate_forecast,f.rain_gate_forecast)
    if limits.get('station_positive_energy_7d_q90_mmwe') is not None:
        # This is a LOWER BOUND on accumulated nonnegative energy: retain only
        # measured pre-issue contributions still inside each future window.
        # No future residual-energy observation is invented from shortwave.
        observed=h.loc[h.index<issue,'energy_balance_residual'].clip(lower=0)*3600/334000
        memory=[]
        for t in f.index:
            retained=observed.loc[(observed.index>t-pd.Timedelta(hours=168)) & (observed.index<=t)]
            memory.append(float(retained.sum()) if retained.notna().sum()>=126 else np.nan)
        f['station_energy_memory_lower_bound_mmwe']=memory
        f['energy_gate_forecast_memory_lower_bound']=physical_exceedance(
            memory,limits['station_positive_energy_7d_q90_mmwe'],limits['station_energy_gate_response_scale_mmwe'])
        f['trigger_condition_forecast']=np.maximum(f.trigger_condition_forecast,f.energy_gate_forecast_memory_lower_bound)
    from pyproj import Geod
    _,_,distance=Geod(ellps='WGS84').inv(85.52723,28.2877,payload['longitude'],payload['latitude'])
    status={'status':'ACTUAL_ARCHIVED_RUN_PROPAGATED_THROUGH_Q2',
       'n_future_hours':len(f),'model_surface_elevation_m':model_z,
       'model_grid_lonlat':[payload['longitude'],payload['latitude']],
       'requested_source_lonlat':[85.52723,28.2877],'model_grid_offset_km':float(distance/1000),
       'valid_count_by_variable':{k:int(f[k].notna().sum()) for k in FORECAST_VARS},
       'freezing_height_above_source_n':int((f.freezing_level_height.dropna()>=source_z).sum()) if f.freezing_level_height.notna().any() else None,
       'precip_future_24h_mm':float(f.precipitation.iloc[:24].sum(min_count=18)),
       'first6h_max_condition_index':float(f.trigger_condition_forecast.iloc[:7].max()),
       'temperature_adjustment':'source T = grid T - lapse*(source elevation-model surface elevation)',
       'spatial_proxy_limitation':'past rain/thermal memory from distant station, future source grid; no site-calibrated correction',
       'availability_at_issue':'run init + assumed 6h publication latency; actual historic dissemination not independently observed',
       'used_for_current_issue_sensitivity':True}
    status['energy_forecast_rule']='nonnegative seven-day energy memory lower bound from pre-issue measured contributions; future contributions unknown'
    f.to_csv(ctx.folder('q2')/'Q2_archived_forecast_trigger_conditions.csv',index_label='forecast_time_bjt')
    write_json(ctx.folder('q2')/'Q2_archived_forecast_quality.json',status)
    return f,status
