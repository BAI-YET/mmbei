"""Q2: source stability, critical parameters and supporting climate diagnostics.

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
    energy_limit=None;energy_scale=None
    if 'energy_balance_residual' in d:
        # Positive station energy is a thermal-support proxy, with its own
        # measured footprint. It is never relabelled as observed glacier melt.
        energy=(d.energy_balance_residual.clip(lower=0)*3600/334000).rolling(
            '168h',min_periods=126).sum()
        energy_ref=energy.loc[baseline].dropna()
        if len(energy_ref)>=300:
            energy_limit=float(energy_ref.quantile(quantile))
            energy_scale=max(.05,float(energy_ref.std())*.3)
            d['station_positive_energy_7d_mmwe']=energy
            d['energy_gate']=physical_exceedance(energy,energy_limit,energy_scale)
            d['trigger_condition_index']=np.maximum(d.trigger_condition_index,d.energy_gate)
    d['melt_condition_met']=(src_pdd7>=pdd_limit).astype('boolean').mask(t_src.isna()|src_pdd7.isna())
    d['rain_condition_met']=((d.rain_24h>=rain24_limit)&(d.rain_24h>0.01)).astype('boolean').mask(d.rain_24h.isna())
    thresholds={'baseline_years':list(baseline_years),'reference_quantile':float(quantile),
                'station_elevation_m':station_z,'source_elevation_m':float(source_z),
                'lapse_rate_C_per_km':float(lapse),
                'source_pdd7_warmseason_q90_Cday':pdd_limit,
                'rain_24h_nonzero_q90_mm':rain24_limit,
                'hourly_nonzero_precip_q90_mm_h':hourly_rain_limit,
                'station_positive_energy_7d_q90_mmwe':energy_limit,
                'station_energy_gate_response_scale_mmwe':energy_scale,
                'primary_dynamic_channels':['rain','source_positive_temperature']+(['station_energy_support'] if energy_limit is not None else []),
                'threshold_role':'historical_reference_condition',
                'station_energy_scope':'station positive residual energy equivalent; thermal support proxy',
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

def static_ice_susceptibility(g,weights=None):
    weights=np.asarray([.40,.25,.20,.15] if weights is None else weights,float)
    if weights.shape!=(4,) or not np.isfinite(weights).all() or np.any(weights<0) or not np.isclose(weights.sum(),1):
        raise ValueError('Four nonnegative susceptibility weights must sum to one')
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
    base=weights[0]*shear+weights[1]*thick+weights[2]*speed
    score=np.where(np.isfinite(strain),base+weights[3]*strain_normalized,
                   base/max(float(weights[:3].sum()),EPS))
    score[~m]=np.nan
    return score,m,tau


def source_proxy_sensitivity(ctx,h,g,static,mask,limits):
    """Recompute declared weight and elevation/lapse cases from actual source grids."""
    out=ctx.folder('q2');issue=pd.Timestamp(ctx.cfg['snapshot_raw']);rows=[]
    for q in [.1,.5,.9]:
        height=float(np.nanquantile(g['dem'][mask],q))
        for lapse in ctx.cfg['q2']['lapse_rates_C_per_km']:
            d,lim=dynamic_conditions(h,ctx.cfg['historic_years'],height,limits['station_elevation_m'],lapse)
            current=d.loc[d.index<issue].iloc[-1]
            rows.append({'glacier_elevation_quantile':q,'source_elevation_m':height,'lapse_C_per_km':lapse,
                'source_pdd7_Cday':current.source_pdd7_lapse_Cd,
                'historical_pdd7_q90_Cday':lim['source_pdd7_warmseason_q90_Cday'],
                'condition_index':current.trigger_condition_index})
    thermal=pd.DataFrame(rows);thermal.to_csv(out/'Q2_source_elevation_lapse_sensitivity.csv',index=False)
    base=static[mask];high=base>=np.quantile(base,.8);weight_rows=[]
    scenarios={'reference':[.4,.25,.2,.15],'equal':[.25]*4,
        'gravity_emphasis':[.55,.15,.15,.15],'motion_emphasis':[.25,.2,.3,.25]}
    for name,weights in scenarios.items():
        alternative,_,_=static_ice_susceptibility(g,weights);value=alternative[mask]
        top=value>=np.nanquantile(value,.8)
        weight_rows.append({'weight_experiment':name,'weights':','.join(map(str,weights)),
            'spearman_vs_reference':float(spearmanr(base,value,nan_policy='omit').statistic),
            'top20_IoU_vs_reference':float(np.logical_and(high,top).sum()/max(1,np.logical_or(high,top).sum()))})
    weight_table=pd.DataFrame(weight_rows);weight_table.to_csv(out/'Q2_static_weight_sensitivity.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.4),layout='constrained')
    for q,color in zip([.1,.5,.9],COLORS):
        d=thermal.loc[thermal.glacier_elevation_quantile.eq(q)]
        axs[0].plot(d.lapse_C_per_km,d.source_pdd7_Cday,'-o',color=color,label=f'冰川高程q{100*q:.0f}')
    axs[0].set(xlabel='温度递减率 (°C/km)',ylabel='源区7日正积温 (°C·d)',title='a  高程与递减率代理敏感性')
    axs[0].legend(frameon=False,fontsize=8)
    x=np.arange(len(weight_table))
    axs[1].bar(x-.18,weight_table.spearman_vs_reference,width=.36,label='空间秩相关',color=COLORS[0])
    axs[1].bar(x+.18,weight_table.top20_IoU_vs_reference,width=.36,label='高分区IoU',color=COLORS[2])
    axs[1].set(xticks=x,xticklabels=['参考','等权','重力优先','运动优先'],ylim=(0,1.05),title='b  固定权重的实际重算比较')
    axs[1].legend(frameon=False,fontsize=8)
    fig.suptitle('源区空间代理与静态权重敏感性')
    savefig(fig,out/'Q2_19_source_proxy_and_weight_sensitivity.png')
    return {'elevation_lapse_cases':rows,'static_weight_cases':weight_rows,
        'design_role':'declared sensitivity experiments; no disaster labels fitted'}

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
    if not ctx.cfg['q2']['initiation'].get('enabled',False):
        for pattern in ['Q2_mechanical_*','Q2_source_stability_*','Q2_20_mechanical_*']:
            for old in out.glob(pattern):
                if old.is_file():old.unlink()
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
    # Static terrain is prepared in Q2 and inherited by Q3/Q4, not inferred
    # backwards from a flow simulation. The mapped case source is a fixed
    # spatial modelling assumption, not an independently available pre-event map.
    import yaml
    from .q4_terrain_reconstruction import prepare_reconstruction_geometry
    evidence=yaml.safe_load((ctx.root/'configs/case_evidence.yaml').read_text())
    route,_,_,_,geometry=prepare_reconstruction_geometry(ctx,evidence,output_question='q2')
    with np.load(out/'Q2_reconstruction_sections.npz') as a:pixels=a['source_pixels_utm'].copy()
    study_region=np.zeros(mask.shape,bool)
    jj=np.abs(g['x'][None,:]-pixels[:,0,None]).argmin(axis=1)
    ii=np.abs(g['y'][None,:]-pixels[:,1,None]).argmin(axis=1);study_region[ii,jj]=True
    z_src=float(geometry['selected_source_mean_elevation_m'])
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
    columns=['source_temp_lapse_C','source_pdd7_lapse_Cd','rain_24h','rain_gate','melt_gate',
             'trigger_condition_index','melt_condition_met','rain_condition_met']
    columns += [c for c in ['station_positive_energy_7d_mmwe','energy_gate'] if c in d]
    d[columns].to_csv(
          out/'Q2_station_based_dynamic_trigger_conditions.csv.gz',compression='gzip',index_label='time_bjt')
    snap=pd.Timestamp(ctx.cfg['snapshot_raw']);valid=d.loc[d.index<snap].iloc[-1]
    spatial=scale_reliability(static,mask,g)
    source_sensitivity=source_proxy_sensitivity(ctx,h,g,static,mask,lim)
    spatial.to_csv(out/'Q2_spatial_scale_reliability.csv',index=False)
    sequential=sequential_trigger_assessment(ctx,h,d,geometry)
    risk=static*float(sequential['current_score'])
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
    _plot_glacier_dynamics_and_candidate(ctx,g,static,mask,candidate_regions,candidate_masks,study_region)
    current={'snapshot_bjt':str(snap),'rain_24h_mm':float(valid.rain_24h),
            'rain_condition_met':bool(valid.rain_condition_met) if pd.notna(valid.rain_condition_met) else None,
            'source_pdd7_Cday':float(valid.source_pdd7_lapse_Cd),
            'melt_condition_met':bool(valid.melt_condition_met) if pd.notna(valid.melt_condition_met) else None,
            'trigger_condition_index':float(valid.trigger_condition_index),
            'station_energy_gate':float(valid.energy_gate) if 'energy_gate' in valid else None,
            'spatial_index_max':float(np.nanmax(risk)),
            'spatial_index_mean':float(np.nanmean(risk))}
    mechanics=(execute_source_stability_experiment(ctx,h,g,candidate_masks)
               if ctx.cfg['q2']['initiation'].get('enabled',False) else None)
    report={'task':'Q1-feature-informed statistical trigger assessment with a shared case source',
       'official_inputs':['MET/RADM/SOIL/FLUX','historic temperature','glacier mask/thickness/slope/velocity','historical gullies','DEM'],
       'external_satellite_evidence_used':extra['satellite_summary'],
       'external_forecast_status':forecast_status,'external_forecast_summary':forecast_summary,
       'thresholds':lim,'current':current,
       'model':'I(t)=R(t)*(a+(1-a)*A_Q1(t)); historical statistical critical levels; fixed terrain/material source assumption',
       'sequential_model':sequential,
       'mechanics':mechanics,
       'climate_reference_model':'static source-condition index x multichannel climate-condition index; retained as supporting evidence',
       'static_weights':{'gravity_shear':.40,'ice_thickness':.25,'historic_speed':.20,'velocity_strain':.15},
       'not_failure_probability':True,
       'source_proxy_sensitivity':source_sensitivity,
       'weather_station_not_nepal_source':True,
       'release_candidate_regions':candidate_regions.to_dict('records'),'time_resolution_reliability':response,'spatial_resolution_reliability':spatial.to_dict('records')}
    write_json(out/'Q2_primary_trigger_assessment.json',report)
    print('[Q2-trigger] event rain condition',current['rain_condition_met'],'melt',current['melt_condition_met'],
          'spatial mean',round(current['spatial_index_mean'],3),'forecast',forecast_status)
    return report


def sequential_trigger_assessment(ctx,h,d,geometry):
    """Q1 -> Q2 contract. One assessment index is inherited unchanged by Q3.

    Q1 ranking supplies transparent design weights, not fitted causal effects.
    The statistical critical level is an operational modelling definition,
    not a measured universal material-failure threshold.
    """
    import hashlib
    out=ctx.folder('q2');q1=ctx.folder('q1')/'Q1_answer_mechanism_top5.csv'
    if not q1.exists():raise FileNotFoundError('Run Q1 before Q2: mechanism representatives are required')
    selected=pd.read_csv(q1);features=selected.feature.tolist()
    if not set(features).issubset(h.columns):raise ValueError('Selected Q1 factors missing from official features')
    base=h.index.year.isin(ctx.cfg['historic_years'])&h.index.month.isin([5,6,7,8,9])
    weights=selected.rarity_anomaly_score.to_numpy(float);weights/=weights.sum()
    directions=np.sign(selected.event_value.to_numpy()-selected.reference_median.to_numpy());directions[directions==0]=1
    components=pd.DataFrame(index=h.index);normalisation={}
    keys=pd.MultiIndex.from_arrays([h.index.month,h.index.hour])
    for feature,direction in zip(features,directions):
        reference=h.loc[base,feature].dropna();group=reference.groupby([reference.index.month,reference.index.hour])
        med=group.median();mad=group.apply(lambda v:1.4826*np.median(abs(v-np.median(v))))
        floor=max(.05,float(np.median(abs(reference-reference.median()))*1.4826*.2))
        centre=med.reindex(keys).to_numpy();scale=np.maximum(mad.reindex(keys).fillna(floor).to_numpy(),floor)
        components[feature]=np.clip(direction*(h[feature].to_numpy()-centre)/scale,0,4)/4
        normalisation[feature]={'direction':float(direction),'weight':float(weights[features.index(feature)]),'scale_floor':floor}
    values=components.to_numpy();known=np.isfinite(values)
    denominator=(known*weights[None,:]).sum(axis=1)
    anomaly=np.divide(np.nansum(values*weights[None,:],axis=1),denominator,
        out=np.full(len(h),np.nan),where=(known.sum(axis=1)>=int(np.ceil(.8*len(features))))&(denominator>0))
    floor=float(ctx.cfg['q2']['q1_modulation_floor'])
    if not 0<=floor<=1:raise ValueError('Q1 modulation floor must be in [0,1]')
    score=d.trigger_condition_index*(floor+(1-floor)*anomaly)
    reference=score.loc[base].dropna()
    threshold=float(reference.quantile(ctx.cfg['q2']['assessment_quantile']))
    monitor=float(reference.quantile(ctx.cfg['q2']['monitoring_quantile']))
    if threshold<=0:raise ValueError('Statistical trigger reference must be positive')
    result=pd.DataFrame({'Q1_selected_feature_anomaly':anomaly,'Q2_dynamic_water_heat':d.trigger_condition_index,
        'assessment_score':score},index=h.index)
    result.to_csv(out/'Q2_sequential_assessment_timeseries.csv.gz',index_label='time_bjt',compression='gzip')
    components.to_csv(out/'Q2_Q1_feature_components.csv.gz',index_label='time_bjt',compression='gzip')
    issue=pd.Timestamp(ctx.cfg['snapshot_raw']);current=result.loc[result.index<issue].iloc[-1]
    contract={'model':'I=R*(a+(1-a)*A_Q1)','assessment_issue_bjt':str(issue),
        'Q1_answer_sha256':hashlib.sha256(q1.read_bytes()).hexdigest(),
        'Q1_features':features,'Q1_weights':dict(zip(features,weights)),
        'Q1_definition':normalisation,'Q1_selection_role':'case-informed Q1 analysis; no claim of independent cross-event feature selection',
        'q1_modulation_floor':floor,'critical_score':threshold,'monitoring_score':monitor,
        'threshold_role':'historical statistical critical level, explicitly defined for this model',
        'current_score':float(current.assessment_score),'current_Q1_anomaly':float(current.Q1_selected_feature_anomaly),
        'current_dynamic_water_heat':float(current.Q2_dynamic_water_heat),'baseline_years':ctx.cfg['historic_years'],
        'source_geometry_file':'Q2_reconstruction_geometry.json','source_sections_file':'Q2_reconstruction_sections.npz',
        'source_route_file':'Q2_reconstruction_route.csv','source_geometry_sha256':geometry['geometry_contract_sha256'],
        'source_to_port_km':geometry['length_m']/1000,'source_area_km2':geometry['projected_source_area_m2']/1e6,
        'fixed_source_assumption':'case-mapped source, adequate material and connecting terrain; not a blind source-identification claim',
        'score_is_calibrated_disaster_probability':False,'observations_end_before_issue':True}
    write_json(out/'Q2_sequential_model_contract.json',contract)
    return contract


def mechanical_source_geometry(g,region,rho=900.0):
    """Force balance of an effective coherent ice block, not a rock-fracture model.

    Official H is vertical ice thickness. For a horizontal pixel area Ah,
    basal area Ab=Ah/cos(theta), weight=rho*g*H*Ah. Integrate W*sin(theta)
    and W*cos(theta) before dividing by total basal area. This avoids mixing
    normal thickness with vertical thickness or counting a thin edge as a block.
    The basal plane is approximated by the surface slope; its geometry is unmeasured.
    """
    valid=(np.asarray(region,bool) & np.isfinite(g['thickness']) &
           (g['thickness']>0) & np.isfinite(g['slope']) &
           (g['slope']>0) & (g['slope']<80) & np.isfinite(g['dem']))
    if valid.sum()<10:raise ValueError('Too few valid source pixels for force balance')
    theta=np.deg2rad(g['slope'][valid].astype(float))
    ah=abs(float(g['transform'][1])*float(g['transform'][5]))
    basal=ah/np.cos(theta);weight=rho*9.81*g['thickness'][valid].astype(float)*ah
    return {'horizontal_area_m2':float(valid.sum()*ah),'basal_area_m2':float(basal.sum()),
        'n_pixels':int(valid.sum()),'rho_kg_m3':float(rho),
        'source_elevation_m':float(np.mean(g['dem'][valid])),
        'mean_vertical_ice_thickness_m':float(np.mean(g['thickness'][valid])),
        'driving_stress_kPa':float(np.sum(weight*np.sin(theta))/basal.sum()/1000),
        'normal_stress_kPa':float(np.sum(weight*np.cos(theta))/basal.sum()/1000),
        'effective_slope_deg':float(np.rad2deg(np.arctan2(np.sum(weight*np.sin(theta)),np.sum(weight*np.cos(theta))))),
        'geometry_scope':'connected pre-event-data candidate, not the post-event mapped collapse outline',
        'mechanism_scope':'effective basal sliding of an ice-bearing source; rock bridges, tensile fracture and entrainment are not resolved',
        'basal_plane_assumption':'local surface slope; basal/discontinuity geometry not measured'}


def mechanical_candidate_region(g,minimum_slope_deg=30.0,minimum_pixels=10):
    """Geometric screening, not a fitted failure angle or scored ice inventory.

    Retain coherent steep ice patches; select the patch with greatest integrated
    downslope gravitational force. This avoids treating the largest flat ice
    reservoir as the detachment block. Other source screens are explicit experiments.
    """
    mask=(g['glacier_mask'] & np.isfinite(g['thickness']) & (g['thickness']>0) &
          np.isfinite(g['slope']) & (g['slope']>=minimum_slope_deg) & (g['slope']<80))
    components,n=label(mask,structure=np.ones((3,3),int))
    counts=np.bincount(components.ravel(),minlength=n+1)
    force=np.where(mask,900*9.81*g['thickness']*np.sin(np.deg2rad(g['slope'])),0)
    sums=np.bincount(components.ravel(),weights=force.ravel(),minlength=n+1)
    sums[(counts<minimum_pixels)]=0;sums[0]=0
    if not np.any(sums>0):raise ValueError('No connected steep source candidate under declared geometric screen')
    region=components==int(np.argmax(sums))
    return region


def initiation_parameter_design(config):
    """Space-filling, equal-design-weight scenarios, NOT an inferred posterior."""
    from scipy.stats import qmc
    keys=['cohesion_kPa','friction_angle_deg','storage_yield','drainage_hours',
          'degree_day_factor_mm_Cday','infiltration_fraction','lapse_C_km']
    unit=qmc.Sobol(len(keys),scramble=True,seed=int(config['seed'])).random_base2(int(config['scenario_power']))
    columns={}
    for i,key in enumerate(keys):
        lo,hi=map(float,config[key])
        if not (0<=lo<hi):raise ValueError('Invalid initiation range: '+key)
        columns[key]=(np.exp(np.log(lo)+unit[:,i]*np.log(hi/lo))
                      if key in ('storage_yield','drainage_hours') and lo>0
                      else lo+(hi-lo)*unit[:,i])
    d=pd.DataFrame(columns);d.insert(0,'scenario_id',np.arange(len(d)))
    return d


def advance_water_storage(storage_mm,input_mm_h,dt_h,drainage_h,capacity_mm=np.inf):
    """Exact linear reservoir step for constant interval forcing.

    Return storage, drained water and spill, all in mm. Capacity limits pressure
    to normal overburden; excess water spills instead of producing negative stress.
    This is a lumped drainage experiment, not a Richards infiltration solution.
    """
    s=np.asarray(storage_mm,float);q=np.asarray(input_mm_h,float);tau=np.asarray(drainage_h,float)
    if dt_h<=0 or np.any(tau<=0) or np.any(s<0) or np.any(q<0):raise ValueError('Invalid reservoir inputs')
    e=np.exp(-dt_h/tau);unbounded=s*e+q*tau*(-np.expm1(-dt_h/tau))
    result=np.minimum(unbounded,capacity_mm)
    spill=np.maximum(unbounded-result,0)
    drained=s+q*dt_h-unbounded
    return result,drained,spill


def source_safety_factor(storage_mm,parameters,geometry):
    p=parameters
    pressure_kPa=np.minimum(np.asarray(storage_mm,float)/(1000*p.storage_yield.to_numpy())*9.81,
                            geometry['normal_stress_kPa'])
    return (p.cohesion_kPa.to_numpy()+np.maximum(geometry['normal_stress_kPa']-pressure_kPa,0)*
            np.tan(np.deg2rad(p.friction_angle_deg.to_numpy())))/geometry['driving_stress_kPa']


def critical_storage(parameters,geometry):
    """Critical pore pressure and retained water. Infinite means water alone
    cannot overcome cohesion. Zero means this assumed block is unstable dry.
    """
    c=parameters.cohesion_kPa.to_numpy();mu=np.tan(np.deg2rad(parameters.friction_angle_deg.to_numpy()))
    pressure=geometry['normal_stress_kPa']-(geometry['driving_stress_kPa']-c)/mu
    pressure=np.maximum(pressure,0)
    pressure=np.where(c>geometry['driving_stress_kPa'],np.inf,pressure)
    return pressure,pressure/9.81*1000*parameters.storage_yield.to_numpy()


def source_history_states(h,issue,parameters,geometry,station_z=4140,spinup_days=30,consistency_days=7,step_h=1):
    """Use observations STRICTLY before issuance; reject inconsistent scenarios.

    A scenario predicting failure during the known-surviving last week is
    incompatible with this model and source hypothesis, not a new observed event.
    Missing forcing or clock gaps fail explicitly; no zero-rain imputation.
    """
    issue=pd.Timestamp(issue);begin=issue-pd.Timedelta(days=spinup_days)
    hours=pd.date_range(begin,issue-pd.Timedelta(hours=1),freq='h')
    observed=h.reindex(hours)[['temp_2m','precip_mm']]
    if not np.isfinite(observed.to_numpy(float)).all():raise ValueError('Incomplete initiation spin-up forcing')
    if step_h not in (1,3,6) or len(observed)%step_h:raise ValueError('Incompatible history time step')
    n=len(parameters);s=np.zeros(n);minimum=np.full(n,np.inf);curve=[]
    p=parameters;capacity=geometry['normal_stress_kPa']/9.81*1000*p.storage_yield.to_numpy()
    source_temperature=(observed.temp_2m.to_numpy()[:,None] -
        (geometry['source_elevation_m']-station_z)/1000*p.lapse_C_km.to_numpy()[None,:])
    rain=observed.precip_mm.clip(lower=0).to_numpy()[:,None]*(source_temperature>0)
    melt=np.maximum(source_temperature,0)*p.degree_day_factor_mm_Cday.to_numpy()[None,:]/24
    liquid=(rain+melt)*p.infiltration_fraction.to_numpy()[None,:]
    for k in range(0,len(hours),step_h):
        # Average rates preserve the input integral when coarsening; do not subsample.
        s,_,_=advance_water_storage(s,liquid[k:k+step_h].mean(axis=0),step_h,p.drainage_hours.to_numpy(),capacity)
        fs=source_safety_factor(s,p,geometry)
        end=hours[k+step_h-1]+pd.Timedelta(hours=1)
        if end>issue-pd.Timedelta(days=consistency_days):minimum=np.minimum(minimum,fs)
        curve.append((end,fs.copy(),s.copy()))
    return s,minimum,curve


def execute_source_stability_experiment(ctx,h,g,candidate_masks):
    """One mechanics contract, consumed unchanged by Q3."""
    cfg=ctx.cfg['q2']['initiation'];out=ctx.folder('q2');issue=pd.Timestamp(ctx.cfg['snapshot_raw'])
    region=mechanical_candidate_region(g,cfg['minimum_source_slope_deg'],cfg['minimum_source_pixels'])
    geometry=mechanical_source_geometry(g,region,cfg['ice_density_kg_m3'])
    ii,jj=np.where(region)
    geometry['centroid_utm_x']=float(np.mean(g['x'][jj]));geometry['centroid_utm_y']=float(np.mean(g['y'][ii]))
    geometry['source_selection']='connected steep patch with greatest integrated downslope force; no event-time fit'
    with rasterio.open(out/'Q2_mechanical_source_region_30m.tif','w',driver='GTiff',height=region.shape[0],width=region.shape[1],
        count=1,dtype='uint8',crs='EPSG:32645',transform=Affine.from_gdal(*g['transform']),nodata=255,compress='deflate') as ds:
        ds.write(np.where(g['glacier_mask'],region.astype('uint8'),255).astype('uint8'),1)
    parameters=initiation_parameter_design(cfg)
    s,minimum,curve=source_history_states(h,issue,parameters,geometry,ctx.cfg['station_elevation_m'],
        cfg['spinup_days'],cfg['consistency_days'])
    pressure,limit=critical_storage(parameters,geometry)
    fs=source_safety_factor(s,parameters,geometry)
    accepted=(minimum>1)&(fs>1)
    if not accepted.any():raise ValueError('No initiation scenario is consistent with source survival; inspect assumptions, do not retune to event time')
    table=parameters.copy();table['storage_at_issue_mm']=s;table['Fs_at_issue']=fs
    table['minimum_Fs_during_consistency_window']=minimum;table['survival_consistent']=accepted
    table['critical_pressure_kPa']=np.where(np.isfinite(pressure),pressure,np.nan)
    table['water_only_failure_possible']=np.isfinite(limit)
    table['critical_storage_mm']=np.where(np.isfinite(limit),limit,np.nan)
    table['additional_storage_to_failure_mm']=np.maximum(table.critical_storage_mm-s,0)
    tau=parameters.drainage_hours.to_numpy();duration=float(ctx.cfg['q3']['primary_horizon_hours'])
    required=(limit-s*np.exp(-duration/tau))/(tau*(-np.expm1(-duration/tau)))
    table['constant_effective_liquid_input_critical_mm_h']=np.where(np.isfinite(required),np.maximum(required,0),np.nan)
    table.to_csv(out/'Q2_mechanical_parameter_experiments.csv',index=False)
    curves=pd.DataFrame([{'time_bjt':t,'Fs_q025':np.quantile(f[accepted],.025),'Fs_median':np.median(f[accepted]),
        'Fs_q975':np.quantile(f[accepted],.975),'storage_median_mm':np.median(w[accepted])} for t,f,w in curve])
    curves.to_csv(out/'Q2_source_stability_history.csv',index=False)
    spatial=[]
    for screen in cfg['source_slope_experiments_deg']:
        alternative=mechanical_candidate_region(g,screen,cfg['minimum_source_pixels'])
        geo=mechanical_source_geometry(g,alternative,cfg['ice_density_kg_m3']);pc,sc=critical_storage(parameters,geo)
        other_s,other_min,_=source_history_states(h,issue,parameters,geo,ctx.cfg['station_elevation_m'],cfg['spinup_days'],cfg['consistency_days'])
        spatial.append({'minimum_screen_slope_deg':screen,**{k:geo[k] for k in ['horizontal_area_m2','driving_stress_kPa','normal_stress_kPa','source_elevation_m','effective_slope_deg']},
            'finite_critical_pressure_median_kPa':float(np.median(pc[np.isfinite(pc)])),
            'dry_unstable_scenario_fraction':float((pc<=0).mean()),'n_survival_consistent':int((other_min>1).sum())})
    pd.DataFrame(spatial).to_csv(out/'Q2_mechanical_source_geometry_experiment.csv',index=False)
    temporal=[]
    for step in (1,3,6):
        coarse,_,_=source_history_states(h,issue,parameters,geometry,ctx.cfg['station_elevation_m'],cfg['spinup_days'],cfg['consistency_days'],step)
        alternate=source_safety_factor(coarse,parameters,geometry)
        temporal.append({'forcing_aggregation_hours':step,'storage_abs_error_median_mm':float(np.median(abs(coarse-s))),
            'Fs_abs_error_q95':float(np.quantile(abs(alternate-fs),.95)),
            'at_issue_stability_disagreement_fraction':float(((alternate<=1)!=(fs<=1)).mean())})
    pd.DataFrame(temporal).to_csv(out/'Q2_mechanical_temporal_experiment.csv',index=False)
    result={'criterion':'first passage Fs<=1','geometry':geometry,
        'parameter_experiment_file':'Q2_mechanical_parameter_experiments.csv',
        'n_design_scenarios':len(table),'n_survival_consistent':int(accepted.sum()),
        'scenario_weights':'equal within declared design after pre-issue survival consistency; not calibrated probability',
        'parameter_ranges':cfg,'critical_equation':'u_crit=sigma_n-(tau_d-c)/tan(phi); S_crit=1000*Sy*u_crit/(rho_water*g)',
        'units':'stress kPa; retained water mm; pressure head m; time h',
        'forcing_equation':'q=eta*(rain_if_source_T_positive+DDF*max(source_T,0)/24)',
        'storage_equation':'dS/dt=q-S/tau; overflow at u=sigma_n',
        'energy_role':'station energy retained as independent thermal evidence; not added again to temperature-index melt',
        'source_strength_parameters_measured':False,'pore_pressure_measured':False,
        'event_time_used_to_fit_parameters':False,'post_event_source_outline_used':False,
        'parameter_evidence':'DDF is a temperature-index melt surrogate; strength, effective storage and drainage ranges are declared numerical experiments, not site measurements',
        'equation_references':['https://www.cambridge.org/core/journals/annals-of-glaciology/article/contribution-of-discontinuous-rockmass-failure-to-glacier-erosion/76A97583C001406D3F6F81015EE22BE0',
            'https://tc.copernicus.org/articles/17/211/2023/','https://tc.copernicus.org/articles/15/2295/2021/'],
        'identification':'Mechanically defined conditional initiation experiment; available observations do not identify a unique local failure strength or true event probability.',
        'spatial_experiments':spatial,'temporal_experiments':temporal}
    old_region=candidate_masks[min(candidate_masks,key=lambda x:abs(x-.85))]
    result['inventory_candidate_comparison']=mechanical_source_geometry(g,old_region,cfg['ice_density_kg_m3'])
    write_json(out/'Q2_source_stability_contract.json',result)
    fig,axs=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    valid=table.loc[table.survival_consistent & table.water_only_failure_possible]
    im=axs[0].scatter(valid.friction_angle_deg,valid.cohesion_kPa,c=valid.critical_pressure_kPa,s=12,cmap='viridis')
    axs[0].set(xlabel='有效静态摩擦角 (°)',ylabel='有效黏聚力 (kPa)',title='a  条件临界孔压；每点为公开参数试验')
    fig.colorbar(im,ax=axs[0],label='临界孔压 (kPa)')
    axs[1].plot(curves.time_bjt,curves.Fs_median,color=COLORS[0])
    axs[1].fill_between(curves.time_bjt,curves.Fs_q025,curves.Fs_q975,color=COLORS[0],alpha=.18)
    axs[1].axhline(1,color=COLORS[1],ls='--',label='极限平衡 Fs=1')
    axs[1].set(xlabel='北京时间',ylabel='源区稳定系数 Fs',title='b  灾前条件一致的参数集合；非观测置信带')
    axs[1].legend(frameon=False,fontsize=8);fig.suptitle('源区稳定性判据与关键参数临界条件')
    savefig(fig,out/'Q2_20_mechanical_critical_conditions.png')
    return result


def _plot_glacier_dynamics_and_candidate(ctx,g,static,mask,regions,masks,study_region=None):
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
        axes[2].plot([],[],color='#163C64',lw=1.6,label='原物源储量候选')
    cfg=ctx.cfg['q2']['initiation']
    steep=(mechanical_candidate_region(g,cfg['minimum_source_slope_deg'],cfg['minimum_source_pixels'])
           if study_region is None else study_region)
    axes[2].contour(steep.astype(int),levels=[.5],colors=['#DA783F'],linewidths=1.6)
    axes[2].plot([],[],color='#DA783F',lw=1.6,label='统一建模源区' if study_region is not None else '陡坡起动候选')
    axes[2].legend(frameon=False,fontsize=7,loc='lower left')
    axes[0].set_ylabel('栅格行')
    fig.suptitle('冰川历史运动与候选失稳区空间特征')
    savefig(fig,out/'Q2_13_glacier_velocity_and_release_candidates.png')
