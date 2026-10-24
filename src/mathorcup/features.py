"""Shared raw-clock features, no future leakage, no event-label engineering."""
import numpy as np,pandas as pd
from .common import write_json

def make_hourly(df):
    h=df.set_index('timestamp_raw').sort_index().copy()
    def has(*cols):return all(c in h for c in cols)
    if has('temp_20m','temp_2m'):
        h['temp_gradient']=(h.temp_20m-h.temp_2m)/18
        h['temperature_inversion']=np.where(h[['temp_2m','temp_20m']].notna().all(axis=1),
                                            (h.temp_20m>h.temp_2m).astype(float),np.nan)
    if has('wind_20m','wind_2m'):
        h['wind_gradient']=(h.wind_20m-h.wind_2m)/18
        with np.errstate(divide='ignore',invalid='ignore'):
            h['wind_profile_exponent']=np.where((h.wind_20m>.2)&(h.wind_2m>.2),
                                                np.log(h.wind_20m/h.wind_2m)/np.log(10),np.nan)
        h['wind_profile_exponent']=h['wind_profile_exponent'].where(lambda x:x.between(-1.5,2))
    if has('rh_20m','rh_2m'):
        h['rh_gradient']=(h.rh_20m-h.rh_2m)/18
    if has('wd_20m','wd_2m','wind_20m','wind_2m'):
        # WD is the meteorological *from* azimuth, 0=N; vector wind components.
        for level in [2,20]:
            theta=np.deg2rad(h[f'wd_{level}m'])
            h[f'u_{level}m']=-h[f'wind_{level}m']*np.sin(theta)
            h[f'v_{level}m']=-h[f'wind_{level}m']*np.cos(theta)
        h['vector_wind_shear']=np.hypot(h.u_20m-h.u_2m,h.v_20m-h.v_2m)/18
        rawdiff=(h.wd_20m-h.wd_2m+180)%360-180
        h['wind_direction_turn']=rawdiff.abs().where((h.wind_20m>.5)&(h.wind_2m>.5))
        # Bulk Richardson number: temperature- rather than virtual-potential-temperature
        # stability proxy; only for diagnostic, not a rigorous turbulence criterion.
    if has('temp_20m','temp_2m','u_20m','u_2m','v_20m','v_2m'):
        denom=(h.u_20m-h.u_2m)**2+(h.v_20m-h.v_2m)**2
        tmean=(h.temp_20m+h.temp_2m)/2+273.15
        rb=9.81*18*(h.temp_20m-h.temp_2m)/(tmean*denom)
        h['bulk_richardson_proxy']=rb.where(denom>.25).clip(-3,3)
    if has('sw_down','sw_up','lw_down','lw_up'):
        h['net_shortwave']=h.sw_down-h.sw_up
        h['net_longwave']=h.lw_down-h.lw_up
        h['net_radiation']=h.net_shortwave+h.net_longwave
        h['albedo']=np.where((h.sw_down>80)&(h.sw_up>=0),h.sw_up/h.sw_down,np.nan)
        h['albedo']=h['albedo'].where(lambda x:x.between(0,1.1))
    for shallow,deep,width,name in [
         ('vwc_0p1','vwc_0p8',.7,'vwc_gradient'),
         ('vwc_0p1','vwc_1p6',1.5,'vwc_deep_gradient'),
         ('soilT_0p1','soilT_1p6',1.5,'soil_temperature_gradient'),
         ('soilT_0p1','soilT_0p4',.3,'soil_temperature_shallow_gradient')]:
        if has(shallow,deep):h[name]=(h[shallow]-h[deep])/width
    if has('soilT_0p1'):
        h['soil_surface_frozen']=np.where(h.soilT_0p1.notna(),(h.soilT_0p1<0).astype(float),np.nan)
    if has('temp_2m','rh_2m'):
        es=6.112*np.exp(17.67*h.temp_2m/(h.temp_2m+243.5))
        h['vpd_2m']=es*(1-h.rh_2m/100.)
    if has('sensible_heat','latent_heat'):
        h['bowen_ratio']=np.where(h.latent_heat.abs()>20,h.sensible_heat/h.latent_heat,np.nan)
        h['bowen_ratio']=h['bowen_ratio'].where(lambda x:x.between(-8,8))
        h['turbulent_total_heat']=h.sensible_heat+h.latent_heat
    if has('net_radiation','sensible_heat','latent_heat','soil_heat'):
        # This is the measured budget residual, NOT glacier melt energy, as the
        # station sits far from the source ice and all surface heat storage is unknown.
        h['surface_energy_residual']=h.net_radiation-h.sensible_heat-h.latent_heat-h.soil_heat
        h['energy_balance_residual']=h['surface_energy_residual']
    if has('vector_wind_shear'):h['vector_shear_2_20m']=h['vector_wind_shear']
    if has('rh_gradient'):h['humidity_gradient']=h['rh_gradient']
    if has('bulk_richardson_proxy'):h['bulk_richardson']=h['bulk_richardson_proxy']
    if has('soilT_0p1','soilT_0p8'):h['soil_temp_gradient']=(h.soilT_0p1-h.soilT_0p8)/.7
    soil_cols=['soilT_0p1','soilT_0p2','soilT_0p4','soilT_0p8','soilT_1p6']
    if has(*soil_cols):
        h['soil_freeze_depth_layers']=(h[soil_cols]<=0).sum(axis=1).where(h[soil_cols].notna().all(axis=1))
    if has('soil_heat','soil_heat_0p2'):h['soil_heat_gradient']=(h.soil_heat-h.soil_heat_0p2)/.1
    # Explicit aliases preserve the original 62-factor Q1 contract and units.
    if has('vapor_pressure'):h['vapour_pressure']=h['vapor_pressure']
    if has('water_density'):h['water_vapour_density']=h['water_density']
    if has('ustar','sensible_heat','temp_2m','pressure'):
        t=h.temp_2m+273.15
        rho=h.pressure*100/(287.05*t)
        L=-(rho*1005*t*h.ustar**3)/(.4*9.81*h.sensible_heat.replace(0,np.nan))
        h['z_over_L_proxy']=2/L
        h['z_over_L_proxy']=h['z_over_L_proxy'].where(
            (h.ustar>.1)&(h.sensible_heat.abs()>10)&(h.z_over_L_proxy.abs()<5))
    if has('temp_2m'):
        h['pdd_hour']=h.temp_2m.clip(lower=0)/24
        h['temp_max_24h']=h.temp_2m.rolling('24h',min_periods=18).max()
        h['pdd_7d']=h.pdd_hour.rolling('168h',min_periods=126).sum()
    if has('precip_mm'):
        p=h.precip_mm
        for hours,name in [(24,'rain_24h'),(72,'rain_72h'),(168,'rain_7d')]:
            h[name]=p.rolling(f'{hours}h',min_periods=int(np.ceil(hours*.75))).sum()
    return h

def build_features(ctx):
    out=ctx.processed/'features_hourly.csv.gz'
    if out.exists():
        metadata=ctx.processed/'features_metadata.json'
        essential={'vector_wind_shear','soil_temperature_gradient','net_shortwave',
                   'net_longwave','bowen_ratio','tau','soilT_1p6','vwc_1p6'}
        try:
            from .common import record_warning
            import json
            m=json.loads(metadata.read_text(encoding='utf-8')) if metadata.exists() else {}
            if m.get('version')=='v2.0_task_aligned' and essential.issubset(set(pd.read_csv(out,nrows=0).columns)):
                return
        except (OSError,ValueError,KeyError):
            pass
        record_warning('Stale derived-feature cache detected: rebuilding v2.0 features.')
        out.unlink(missing_ok=True)
    df=pd.read_csv(ctx.processed/'obs_hourly.csv.gz',parse_dates=['timestamp_raw'])
    h=make_hourly(df)
    h.to_csv(out,compression='gzip',index_label='timestamp_raw')
    (ctx.processed/'features_schema_v17.marker').write_text('physical candidates v1.7',encoding='utf-8')
    write_json(ctx.processed/'features_metadata.json',{'time':'Asia/Shanghai UTC+08:00 confirmed by official TPDC dataset 1',
      'forward_looking_rolling':False,
      'pdd_hour_units':'C day per hourly sample',
      'no_future_observations_in_snapshot':True,
      'version':'v2.0_task_aligned','feature_count':int(h.shape[1]),'columns':h.columns.tolist()})
