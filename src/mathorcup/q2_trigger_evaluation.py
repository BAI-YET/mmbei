"""Question 2 primary deliverable: physical-condition threshold model, temporal / spatial reliability.

The result is a reproducible *trigger-condition index*, not an empirically calibrated
probability of a catastrophic event. Morphological labels from official data 6 are
separately modelled by Q2's logistic (known historical gully type), never mixed in.
"""
from __future__ import annotations
import numpy as np, pandas as pd, rasterio, matplotlib.pyplot as plt
from pathlib import Path
from rasterio.transform import Affine
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from scipy.stats import spearmanr
from scipy.ndimage import gaussian_filter, label
from .prepare_geo import glacier_grid
from .q2_external import q2_external_status,load_q2_forecast
from .common import write_json,savefig,COLORS

EPS=1e-9

def robust_norm(z,mask=None,low=.05,high=.95):
    zz=np.asarray(z,float)
    if mask is None:mask=np.isfinite(zz)
    vals=zz[mask & np.isfinite(zz)]
    if vals.size<10:return np.zeros_like(zz,dtype=float)
    a,b=np.quantile(vals,[low,high])
    if b-a<EPS:return np.zeros_like(zz,dtype=float)
    return np.clip((zz-a)/(b-a),0,1)

def physical_exceedance(series,quantile_ref,scale=None):
    s=np.asarray(series,float)
    if scale is None:scale=max(1e-6,abs(float(quantile_ref))*.15)
    return np.clip(.5+.5*(s-quantile_ref)/(2*scale),0,1)

def dynamic_conditions(h,baseline_years=(2022,2023,2024,2025),source_z=5300,station_z=4140,lapse=6.5,quantile=.90):
    """At source elevation, apply physically transparent lapse-rate scenarios to station T.
    Nonzero rain quantile is conditional on actually raining, never q95=0 counted as rain.
    """
    d=h.copy().sort_index()
    if not {'temp_2m','precip_mm','rain_24h','pdd_7d'}.issubset(d):raise ValueError('Missing official time variables for Q2')
    t_src=d.temp_2m-(source_z-station_z)*lapse/1000.0
    src_melt_hour=np.maximum(t_src,0)/24
    src_pdd7=src_melt_hour.rolling('168h',min_periods=126).sum()
    baseline=d.index.year.isin(baseline_years)&d.index.month.isin([5,6,7,8,9])
    pdd_valid=src_pdd7.loc[baseline].dropna()
    wet24=d.rain_24h.loc[baseline].dropna()
    wet24=wet24[wet24>0.01]
    wet_hour=d.precip_mm.loc[baseline].dropna();wet_hour=wet_hour[wet_hour>0.01]
    if len(pdd_valid)<300 or len(wet24)<100 or len(wet_hour)<80:
        raise ValueError('Climate baseline too short for rainfall/melt thresholds')
    pdd_limit=float(pdd_valid.quantile(quantile))
    rain24_limit=float(wet24.quantile(quantile))
    hourly_rain_limit=float(wet_hour.quantile(quantile))
    rain_gate=np.where((d.rain_24h>=rain24_limit)&(d.rain_24h>0.01),
                       physical_exceedance(d.rain_24h,rain24_limit),0)
    rain_gate=np.where(d.rain_24h.notna(),rain_gate,np.nan)
    melt_gate=physical_exceedance(src_pdd7,pdd_limit,scale=max(.05,float(pdd_valid.std())*.3))
    melt_gate=np.where(t_src.notna()&src_pdd7.notna(),melt_gate,np.nan)
    # Gated water-supply / thermal mechanism: max avoids conflating ice-rock failure
    # with short intense rain; isolated rain on the distant station is not claimed source rain.
    d['source_temp_lapse_C']=t_src
    d['source_pdd7_lapse_Cd']=src_pdd7
    d['rain_gate']=np.asarray(rain_gate,float)
    d['melt_gate']=np.asarray(melt_gate,float)
    d['trigger_condition_index']=np.maximum(d.rain_gate,d.melt_gate)
    d['melt_condition_met']=(src_pdd7>=pdd_limit).astype('boolean').mask(t_src.isna()|src_pdd7.isna())
    d['rain_condition_met']=((d.rain_24h>=rain24_limit)&(d.rain_24h>0.01)).astype('boolean').mask(d.rain_24h.isna())
    thresholds={'baseline_years':list(baseline_years),'reference_quantile':float(quantile),
                'station_elevation_m':station_z,'source_elevation_m':float(source_z),
                'lapse_rate_C_per_km':float(lapse),
                'source_pdd7_warmseason_q90_Cday':pdd_limit,
                'rain_24h_nonzero_q90_mm':rain24_limit,
                'hourly_nonzero_precip_q90_mm_h':hourly_rain_limit,
                'zero_rain_cannot_activate_rain_trigger':True,
                'all_quantiles_are_climate_conditions_not_demonstrated_failure_strength':True}
    return d,thresholds

def velocity_strain_proxy(g,sigma_cells=1.0):
    """Surface-velocity-gradient invariant, computed in m/year/m from the
    historical vx,vy, not an instantaneous 2026 deformation measurement.

    Mask-weighted Gaussian smoothing avoids artificially treating NoData=0
    outside ice as a sharp velocity discontinuity.
    """
    mask=(g['glacier_mask'] & np.isfinite(g['vx']) & np.isfinite(g['vy']))
    denom=gaussian_filter(mask.astype(float),sigma_cells,mode='nearest')
    def smooth(field):
        num=gaussian_filter(np.where(mask,field,0.0),sigma_cells,mode='nearest')
        return np.divide(num,denom,out=np.full_like(num,np.nan),where=denom>.6)
    vx=smooth(g['vx']);vy=smooth(g['vy'])
    # NetCDF y is decreasing northing; derivatives must retain this sign.
    dvx_dy,dvx_dx=np.gradient(vx,np.asarray(g['y'],float),np.asarray(g['x'],float))
    dvy_dy,dvy_dx=np.gradient(vy,np.asarray(g['y'],float),np.asarray(g['x'],float))
    exx=dvx_dx;eyy=dvy_dy;exy=.5*(dvx_dy+dvy_dx)
    norm=np.sqrt(exx*exx+eyy*eyy+2*exy*exy)
    # Boundary cells cannot be estimated reliably from incomplete neighbours.
    norm[~mask | (denom<.85)]=np.nan
    return norm

def static_ice_susceptibility(g):
    m=(g['glacier_mask'] & np.isfinite(g['thickness']) &
       np.isfinite(g['slope']) & np.isfinite(g['v']))
    tau=900*9.81*g['thickness']*np.sin(np.deg2rad(g['slope']))/1000
    strain=velocity_strain_proxy(g)
    # Physics groups: driving shear (gravity), available ice volume, and
    # historical movement/deformation. All scores use ONLY glacier pixels.
    shear=robust_norm(tau,m)
    thick=robust_norm(g['thickness'],m)
    speed=robust_norm(g['v'],m)
    strain_normalized=robust_norm(strain,m & np.isfinite(strain))
    # Diagnostic fixed weights, not P(failure) or fitted event probabilities.
    base=.40*shear+.25*thick+.20*speed
    score=np.where(np.isfinite(strain),base+.15*strain_normalized,base/.85)
    score[~m]=np.nan
    return score,m,tau

def candidate_release_regions(g,source_score,mask,quantiles=(.8,.85,.9)):
    """Consecutive physically high-scoring ice clusters define *candidate*
    release patches and volume bounds. Ice inventory is not a measured
    collapse volume, and these are spatial scenarios rather than observations.
    """
    dx=abs(float(g['transform'][1]));dy=abs(float(g['transform'][5]))
    area=dx*dy
    records=[];masks={}
    for percentile in quantiles:
        threshold=float(np.nanquantile(source_score[mask],percentile))
        active=mask & np.isfinite(source_score)&(source_score>=threshold)
        components,n=label(active.astype('uint8'),structure=np.ones((3,3),dtype=int))
        if n<1:continue
        ids=np.arange(1,n+1)
        solid_volumes=np.bincount(components[mask],weights=g['thickness'][mask]*area,
                                 minlength=n+1)
        # Prefer a coherent high-volume patch, not isolated single pixels.
        component=int(ids[np.argmax(solid_volumes[1:])])
        region=(components==component)
        ii,jj=np.where(region)
        x=np.mean(g['x'][jj]);y=np.mean(g['y'][ii])
        vol=float((g['thickness'][region]*area).sum())
        records.append({'upper_score_quantile':float(percentile),
          'selected_score_threshold':threshold,'candidate_area_km2':float(region.sum()*area/1e6),
          'candidate_inventory_ice_volume_m3':vol,'n_pixels':int(region.sum()),
          'candidate_centroid_utm_x':float(x),'candidate_centroid_utm_y':float(y),
          'n_spatial_components':int(n)})
        masks[float(percentile)]=region
    return pd.DataFrame(records),masks


def _save_risk_tif(ctx,arr,g):
    out=ctx.folder('q2');transform=Affine.from_gdal(*g['transform'])
    f=out/'Q2_ice_source_conditional_index_30m.tif'
    with rasterio.open(f,'w',driver='GTiff',height=arr.shape[0],width=arr.shape[1],count=1,
                      dtype='float32',crs='EPSG:32645',transform=transform,nodata=-9999,
                      compress='deflate') as dst:
        dst.write(np.where(np.isfinite(arr),arr,-9999).astype('float32'),1)
    return f

def scale_reliability(arr,mask,g):
    T=Affine.from_gdal(*g['transform'])
    result=[]
    baseline=arr[mask]
    q80=float(np.quantile(baseline,.80))
    for fac in [1,2,4]:
        nh=int(np.ceil(arr.shape[0]/fac));nw=int(np.ceil(arr.shape[1]/fac))
        src=np.where(mask,arr,0).astype('float32');valid=mask.astype('float32')
        dst=np.zeros((nh,nw),np.float32);den=np.zeros_like(dst)
        TT=T*Affine.scale(fac,fac)
        reproject(src,dst,src_transform=T,src_crs='EPSG:32645',dst_transform=TT,dst_crs='EPSG:32645',resampling=Resampling.average)
        reproject(valid,den,src_transform=T,src_crs='EPSG:32645',dst_transform=TT,dst_crs='EPSG:32645',resampling=Resampling.average)
        val=dst/np.maximum(den,1e-6);val[den<.5]=np.nan
        back=np.full(arr.shape,np.nan,np.float32)
        reproject(val,back,src_transform=TT,src_crs='EPSG:32645',src_nodata=np.nan,
                  dst_transform=T,dst_crs='EPSG:32645',dst_nodata=np.nan,resampling=Resampling.nearest)
        good=mask&np.isfinite(back)
        if good.sum()<100:continue
        cor=float(spearmanr(arr[good],back[good]).statistic) if np.std(back[good])>1e-8 else np.nan
        high=back[good]>=np.quantile(back[good],.80)
        base_high=arr[good]>=np.quantile(arr[good],.80)
        overlap=float(np.logical_and(high,base_high).sum()/max(1,np.logical_or(high,base_high).sum()))
        result.append({'grid_resolution_m':fac*30,'n_valid_original_pixels':int(good.sum()),
                'spearman_vs_30m':cor,'high20percent_IoU':overlap})
    return pd.DataFrame(result)

def _plot_trigger_summary(ctx,h,lim,g,s,m,tau,spatial):
    out=ctx.folder('q2');snap=pd.Timestamp(ctx.cfg['snapshot_raw'])
    events=h.loc[(h.index>=snap-pd.Timedelta(days=28))&(h.index<snap)]
    fig,axes=plt.subplots(3,1,figsize=(10.5,8),sharex=True,layout='constrained')
    plots=[('source_pdd7_lapse_Cd','源区七日正积温 (°C·d)','source_pdd7_warmseason_q90_Cday'),
           ('rain_24h','站点24小时降水量 (mm)','rain_24h_nonzero_q90_mm'),
           ('trigger_condition_index','多通道触发条件指数',None)]
    for ax,(key,lab,th),color in zip(axes,plots,COLORS):
        ax.plot(events.index,events[key],color=color,lw=1.65)
        if th is not None:ax.axhline(lim[th],ls='--',color='#C25F4A',label='气候参考阈值')
        ax.set(ylabel=lab)
        ax.grid(axis='y',alpha=.15)
        if th is not None:ax.legend(frameon=False,fontsize=8)
    axes[0].set_title('源区热量与站点降雨的条件演变')
    axes[-1].set_xlabel('北京时间（UTC+8）')
    savefig(fig,out/'Q2_08_trigger_conditions_time_series.png')
    x=g['x'];y=g['y'];ext=[float(x.min()),float(x.max()),float(y.min()),float(y.max())]
    fig,axes=plt.subplots(1,3,figsize=(12.8,4.7),layout='constrained')
    for ax,(v,name,cmap) in zip(axes,[(g['slope'],'冰川坡度 (°)','YlOrBr'),(tau,'驱动剪应力代理 (kPa)','magma'),(s,'地形物源条件指数','viridis')]):
        img=ax.imshow(np.where(m,v,np.nan),extent=ext,origin='upper',cmap=cmap,interpolation='nearest')
        ax.set(title=name,xlabel='UTM45N 东向坐标 (m)');ax.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
        fig.colorbar(img,ax=ax,shrink=.77,pad=.03)
    axes[0].set_ylabel('UTM45N 北向坐标 (m)')
    fig.suptitle('冰川源区地形与物源条件分布',fontsize=12)
    savefig(fig,out/'Q2_09_source_spatial_condition_index.png')
    fig,ax=plt.subplots(figsize=(7.5,4.4),layout='constrained')
    ax.plot(spatial.grid_resolution_m,spatial.spearman_vs_30m,marker='o',label='空间秩相关')
    ax.plot(spatial.grid_resolution_m,spatial.high20percent_IoU,marker='s',label='高危险区IoU')
    ax.set(xlabel='空间网格分辨率 (m)',ylabel='与30m基准的一致性',ylim=(0,1.03),title='触发评估模型空间尺度可靠性')
    ax.grid(axis='y',alpha=.17);ax.legend(frameon=False)
    savefig(fig,out/'Q2_10_spatial_resolution_reliability.png')
    fig,ax=plt.subplots(figsize=(7.5,4.3),layout='constrained')
    for rate,color in [(5.,COLORS[0]),(6.5,COLORS[1]),(8.,COLORS[2])]:
        dd,_=dynamic_conditions(h,baseline_years=ctx.cfg['historic_years'],source_z=lim['source_elevation_m'],station_z=ctx.cfg.get('station_elevation_m',4140),lapse=rate)
        d=dd.loc[(dd.index>=snap-pd.Timedelta(days=21))&(dd.index<snap)]
        ax.plot(d.index,d['source_pdd7_lapse_Cd'],label=f'{rate:g}°C/km',color=color,lw=1.5)
    ax.set(title='源区热量条件垂直递减率敏感性',ylabel='源区七日正积温 (°C·d)',xlabel='北京时间（UTC+8）')
    ax.legend(frameon=False);ax.grid(axis='y',alpha=.12)
    savefig(fig,out/'Q2_11_temperature_lapse_sensitivity.png')

def execute_q2_trigger(ctx):
    out=ctx.folder('q2');extra=q2_external_status(ctx)
    h=pd.read_csv(ctx.processed/'features_hourly.csv.gz',parse_dates=['timestamp_raw']).set_index('timestamp_raw').sort_index()
    g=glacier_grid(ctx);static,mask,tau=static_ice_susceptibility(g)
    candidate_regions,candidate_masks=candidate_release_regions(g,static,mask)
    candidate_regions.to_csv(out/'Q2_source_release_candidates.csv',index=False)
    if len(candidate_regions):
        preferred=candidate_regions.iloc[(candidate_regions.upper_score_quantile-.85).abs().argmin()]
        write_json(out/'Q2_candidate_release_region.json',{
          **{k:(v.item() if hasattr(v,'item') else v) for k,v in preferred.to_dict().items()},
          'not_measured_collapse_release':True,
          'release_volume_is_total_ice_inventory_within_high_score_candidate':True,
          'historical_velocity_not_2026_speed':True})
        region=candidate_masks[float(preferred.upper_score_quantile)]
        tf=Affine.from_gdal(*g['transform'])
        with rasterio.open(out/'Q2_candidate_release_region_30m.tif','w',driver='GTiff',
               height=region.shape[0],width=region.shape[1],count=1,dtype='uint8',
               crs='EPSG:32645',transform=tf,nodata=255,compress='deflate') as ds:
            ds.write(np.where(mask,region.astype('uint8'),255).astype('uint8'),1)
    z_src=float(np.nanmedian(np.where(mask,g['dem'],np.nan)))
    d,lim=dynamic_conditions(h,baseline_years=ctx.cfg['historic_years'],
              source_z=z_src,station_z=ctx.cfg.get('station_elevation_m',4140))
    forecast,forecast_status=load_q2_forecast(ctx,ctx.cfg['snapshot_raw'])
    if forecast is not None:
        future=forecast.loc[(forecast.time_bjt>=pd.Timestamp(ctx.cfg['snapshot_raw']).tz_localize('Asia/Shanghai'))&
                    (forecast.time_bjt<pd.Timestamp(ctx.cfg['snapshot_raw']).tz_localize('Asia/Shanghai')+pd.Timedelta(days=2))]
        from .q2_external import source_forecast_conditions
        _,forecast_summary=source_forecast_conditions(ctx,h,lim,forecast)
        future.to_csv(out/'Q2_source_external_pre_event_forecast_used.csv',index=False)
        if not future.empty:
            # All points are ACTUAL archived model output with a known run initialization.
            # ECMWF grid-level 2m T is not the same as an in-situ source glacier temperature.
            fig,axes=plt.subplots(3,1,figsize=(10.4,7.8),sharex=True,layout='constrained')
            axes[0].plot(future.time_bjt,future.temperature_2m,color=COLORS[0],lw=1.8)
            axes[0].axhline(0,color='#B54F4B',ls='--',lw=1)
            axes[0].set(ylabel='模式格点2m气温 (°C)',title='a  格点气温预报')
            axes[1].plot(future.time_bjt,future.freezing_level_height,color=COLORS[1],lw=1.8)
            axes[1].axhline(z_src,color='#855DA9',ls=':',label='冰川源区代表海拔')
            axes[1].set(ylabel='零度层高程 (m)',title='b  零度层高度预报')
            if future.freezing_level_height.notna().sum()==0:
                axes[1].text(.5,.5,'该起报时次未返回有效零度层高度',transform=axes[1].transAxes,ha='center',fontsize=10)
            axes[1].legend(frameon=False,fontsize=8)
            axes[2].bar(future.time_bjt,future.precipitation,color=COLORS[2],width=.028)
            axes[2].set(ylabel='逐小时降水量 (mm)',xlabel='北京时间（UTC+8）',title='c  小时降水预报')
            for ax in axes:ax.grid(axis='y',alpha=.15)
            fig.suptitle('源区历史起报天气预报条件',fontsize=12)
            savefig(fig,out/'Q2_12_forecast_weather_conditions.png')
    else:forecast_summary={'warning':'No real downloadable archived run available in runtime; station proxy is not falsely labelled a forecast.'}
    d[['source_temp_lapse_C','source_pdd7_lapse_Cd','rain_24h','rain_gate','melt_gate',
       'trigger_condition_index','melt_condition_met','rain_condition_met']].to_csv(
          out/'Q2_station_based_dynamic_trigger_conditions.csv.gz',compression='gzip',index_label='time_bjt')
    snap=pd.Timestamp(ctx.cfg['snapshot_raw']);valid=d.loc[d.index<snap].iloc[-1]
    spatial=scale_reliability(static,mask,g)
    spatial.to_csv(out/'Q2_spatial_scale_reliability.csv',index=False)
    risk=static*float(valid['trigger_condition_index'])
    _save_risk_tif(ctx,risk,g)
    response=[]
    # Time scales defined as different issuance frequencies on the same historical series,
    # not 3 duplicated datasets or fabricated event labels.
    for step_h in [1,6,24]:
        filt=d.loc[(d.index>=snap-pd.Timedelta(days=30))&(d.index<snap)].iloc[::step_h]
        response.append({'temporal_resolution_h':step_h,'fraction_of_time_trigger_index_ge_0p5':float((filt.trigger_condition_index.dropna()>=.5).mean()),
                         'fraction_of_time_rain_gate_positive':float((filt.rain_gate.dropna()>0).mean()),
                         'n_valid_records':int(filt.trigger_condition_index.notna().sum())})
    time_table=pd.DataFrame(response);time_table.to_csv(out/'Q2_temporal_scale_reliability.csv',index=False)
    _plot_trigger_summary(ctx,d,lim,g,static,mask,tau,spatial)
    _plot_glacier_dynamics_and_candidate(ctx,g,static,mask,candidate_regions,candidate_masks)
    current={'snapshot_bjt':str(snap),'rain_24h_mm':float(valid.rain_24h),
            'rain_condition_met':bool(valid.rain_condition_met) if pd.notna(valid.rain_condition_met) else None,
            'source_pdd7_Cday':float(valid.source_pdd7_lapse_Cd),
            'melt_condition_met':bool(valid.melt_condition_met) if pd.notna(valid.melt_condition_met) else None,
            'trigger_condition_index':float(valid.trigger_condition_index),
            'spatial_index_max':float(np.nanmax(risk)),
            'spatial_index_mean':float(np.nanmean(risk))}
    report={'task':'multi-source physical threshold condition assessment',
       'official_inputs':['MET/RADM/SOIL/FLUX','historic temperature','glacier mask/thickness/slope/velocity','historical gullies','DEM'],
       'external_satellite_evidence_used':extra['satellite_summary'],
       'external_forecast_status':forecast_status,'external_forecast_summary':forecast_summary,
       'thresholds':lim,'current':current,
       'model':'conditional index = static glacier gravitational supply index x max(melt_gate, rainfall_gate)',
       'static_weights':{'gravity_shear':.40,'ice_thickness':.25,'historic_speed':.20,'velocity_strain':.15},
       'not_failure_probability':True,
       'weather_station_not_nepal_source':True,
       'release_candidate_regions':candidate_regions.to_dict('records'),'time_resolution_reliability':response,'spatial_resolution_reliability':spatial.to_dict('records')}
    write_json(out/'Q2_primary_trigger_assessment.json',report)
    print('[Q2-trigger] event rain condition',current['rain_condition_met'],'melt',current['melt_condition_met'],
          'spatial mean',round(current['spatial_index_mean'],3),'forecast',forecast_status)
    return report


def _plot_glacier_dynamics_and_candidate(ctx,g,static,mask,regions,masks):
    out=ctx.folder('q2')
    strain=velocity_strain_proxy(g)
    fig,axes=plt.subplots(1,3,figsize=(12.2,4.5),sharex=True,sharey=True,layout='constrained')
    vals=[(g['v'],'冰川历史表面速度 (m/a)','viridis'),
          (strain,'历史速度梯度 (1/年)','magma'),
          (static,'源区综合物理条件指数','YlOrRd')]
    for ax,(v,title,cmap) in zip(axes,vals):
        z=np.where(mask,v,np.nan)
        vmax=np.nanquantile(z,.98)
        im=ax.imshow(z,cmap=cmap,origin='upper',vmin=0,vmax=vmax)
        ax.set(title=title,xlabel='栅格列');fig.colorbar(im,ax=ax,shrink=.74)
    if len(regions):
        k=min(masks,key=lambda x:abs(x-.85));region=masks[k]
        axes[2].contour(region.astype(int),levels=[.5],colors=['#163C64'],linewidths=1.6)
    axes[0].set_ylabel('栅格行')
    fig.suptitle('冰川历史运动与候选失稳区空间特征')
    savefig(fig,out/'Q2_13_glacier_velocity_and_release_candidates.png')
