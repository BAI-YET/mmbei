"""Use datasets (2),(3),(5),(8) without inventing nonexistent disaster labels or pixel population.

Dataset 3 is huge RAR; process only user-extracted 2015-2020 N28E085 tile when available,
never silently repurpose 2000-2020 glacier thinning as 2026 instantaneous motion.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt,rasterio,geopandas as gpd
from rasterio.warp import reproject,Resampling
from rasterio.transform import Affine
from matplotlib.colors import BoundaryNorm,ListedColormap
from .prepare_geo import glacier_grid
from .common import savefig,write_json,record_warning,find_first,COLORS


def plot_county_temperature_climate(ctx):
    fname=ctx.processed/'temp_county_daily.csv.gz'
    if not fname.exists():return
    df=pd.read_csv(fname,parse_dates=['date'])
    d=df.loc[df.date.dt.month.isin([6,7])].copy()
    county=[c for c in ['吉隆县','聂拉木县','萨嘎县','定日县','昂仁县'] if c in d]
    if not county:return
    d['year']=d.date.dt.year;yr=d.groupby('year')[county].mean()
    # Summer Jun-Jul completed for 2025 but annual 2025 is not complete
    fig,ax=plt.subplots(figsize=(10,4.4))
    for i,c in enumerate(county):
        ax.plot(yr.index,yr[c],lw=1.5,label=c,color=COLORS[i],marker='.',markersize=3)
    ax.grid(axis='y',alpha=.16);ax.set(xlabel='年份（2005–2025）',ylabel='6–7月日均气温 (°C)',
         title='吉隆县及周边县市夏季平均气温变化')
    ax.legend(frameon=False,ncol=5,fontsize=8,loc='upper left');savefig(fig,ctx.folder('q1')/'Q1_07_county_2005_2025_summer_climate.png')


def plot_landcover_and_county_population(ctx):
    out=ctx.folder('q2');demfile=find_first(ctx.raw/'unpacked'/'7','DEM.tif')
    lcfile=find_first(ctx.raw/'unpacked'/'5','*2023.tif')
    ds_g=gpd.read_file(ctx.processed/'historic_gullies.gpkg')
    pop_g=gpd.read_file(find_first(ctx.raw/'unpacked'/'8','*.shp'))
    with rasterio.open(demfile) as dem:
        factor=4;nx=int(np.ceil(dem.width/factor));ny=int(np.ceil(dem.height/factor))
        # Downsample both to exactly the same geographic grid; category labels are categorical (NEAREST)
        t=dem.transform*rasterio.Affine.scale(factor,factor)
        z=dem.read(1,out_shape=(ny,nx),resampling=Resampling.average)
        extent=[dem.bounds.left,dem.bounds.right,dem.bounds.bottom,dem.bounds.top]
        lc=np.zeros((ny,nx),dtype=np.uint8)
        with rasterio.open(lcfile) as src:
            reproject(rasterio.band(src,1),lc,src_transform=src.transform,src_crs=src.crs,
                src_nodata=0,dst_transform=t,dst_crs=dem.crs,dst_nodata=0,resampling=Resampling.nearest)
    fig,axs=plt.subplots(1,2,figsize=(11.4,5.1))
    axs[0].imshow(z,extent=extent,cmap='terrain',vmin=np.nanpercentile(z,4),vmax=np.nanpercentile(z,98))
    im=axs[1].imshow(np.ma.masked_equal(lc,0),extent=extent,cmap='tab20',vmin=1,vmax=10,
             interpolation='nearest')
    pop_g.boundary.plot(ax=axs[0],edgecolor='#293A61',lw=.9)
    pop_g.boundary.plot(ax=axs[1],edgecolor='#293A61',lw=.9)
    lon=pd.to_numeric(ds_g['x_coord'],errors='coerce');lat=pd.to_numeric(ds_g['y_coord'],errors='coerce')
    for ax in axs:
        ax.scatter(lon,lat,s=9,c='#B43036',alpha=.8,edgecolors='white',linewidths=.2)
        ax.set(xlabel='经度 (°E)',ylabel='纬度 (°N)')
    axs[0].set_title('地形高程与历史灾点空间分布')
    axs[1].set_title('土地覆盖类型及行政边界分布')
    fig.colorbar(im,ax=axs[1],shrink=.75,label='官方土地覆盖类别码（此处不臆测类别名称）')
    savefig(fig,out/'Q2_06_dem_landcover_history_overlay.png')
    total=pop_g.iloc[0].get('合计')
    write_json(out/'Q2_county_population_limit.json',{
        'population_total_census_county':int(total) if pd.notna(total) else None,
        'county_polygon_count':int(len(pop_g)),
        'no_gridded_population_available':True,
        'exposed_people_are_not_estimated':True,
        'no_pixel_level_population_reallocation':True})


def optional_glacier_thinning_2015_2020(ctx):
    """Use all four official five-year epochs when their N28E085 tiles are
    extracted. The function name is retained for pipeline backward compatibility.
    RAR extraction helper is available under scripts/; never turn historic
    elevation changes into 2026 current glacier displacement.
    """
    base=ctx.raw/'unpacked'/'3_selected'
    g=glacier_grid(ctx);sh=g['glacier_mask'].shape;tf=Affine.from_gdal(*g['transform'])
    periods=['2000-01-01_2005-01-01','2005-01-01_2010-01-01',
             '2010-01-01_2015-01-01','2015-01-01_2020-01-01']
    arrays=[];summaries=[]
    for period in periods:
        fs=list(base.rglob(f'*{period}*dhdt_UTM45_UTM45_GCS.tif'))
        if not fs:
            summaries.append({'period':period,'available':False,'valid_frac':0.0})
            continue
        dest=np.full(sh,np.nan,np.float32)
        with rasterio.open(fs[0]) as src:
            reproject(rasterio.band(src,1),dest,src_transform=src.transform,
                src_crs=src.crs,src_nodata=src.nodata,dst_transform=tf,
                dst_crs='EPSG:32645',dst_nodata=np.nan,
                resampling=Resampling.bilinear)
        valid=g['glacier_mask']&np.isfinite(dest)
        v=dest[valid];arrays.append((period,dest))
        summaries.append({'period':period,'available':True,'valid_frac':float(valid.sum()/g['glacier_mask'].sum()),
              'median_historical_elevation_change_rate_m_per_year':float(np.median(v)) if len(v) else None,
              'mean_historical_elevation_change_rate_m_per_year':float(np.mean(v)) if len(v) else None,
              'q10':float(np.quantile(v,.10)) if len(v) else None,
              'q90':float(np.quantile(v,.90)) if len(v) else None})
    out=ctx.folder('q2')
    pd.DataFrame(summaries).to_csv(out/'Q2_official_glacier_four_period_elevation_change.csv',index=False)
    if not arrays:
        write_json(out/'Q2_dhdt_optional_status.json',{'status':'NONE_EXTRACTED',
          'possible_actions':'Extract official RAR N28E085 four epochs into data/raw/unpacked/3_selected'})
        return
    fig,axes=plt.subplots(1,len(arrays),figsize=(3.3*len(arrays),4.5),
                           sharex=True,sharey=True,layout='constrained')
    if len(arrays)==1:axes=[axes]
    # Shared color scale: rate sign comparable across epochs.
    merged=np.concatenate([np.where(g['glacier_mask'],z,np.nan)[np.isfinite(z)&g['glacier_mask']]
                           for _,z in arrays]); lim=max(.1,float(np.quantile(np.abs(merged),.95)))
    for ax,(period,z) in zip(axes,arrays):
        img=ax.imshow(np.where(g['glacier_mask'],z,np.nan),cmap='RdBu',vmin=-lim,vmax=lim)
        ax.set(title=period[:4]+'–'+period[11:15],xticks=[],yticks=[])
    fig.colorbar(img,ax=axes,shrink=.7,label='高程变化率 (m/a)')
    fig.suptitle('2000—2020年冰川高程变化率')
    savefig(fig,out/'Q2_07_glacier_historic_thinning_optional.png')
    write_json(out/'Q2_dhdt_optional_status.json',{
       'status':'USED_OFFICIAL_FOUR_EPOCH_TILES' if len(arrays)==4 else 'INCOMPLETE_EPOCHS',
       'epochs_found':[p for p,_ in arrays],
       'each_epoch_is_longterm_change_not_2026_instantaneous_ice_acceleration':True})


def landcover_four_year_comparison(ctx):
    """Official dataset⑤ classes, 2020–2023, constrained by its actual valid
    footprint. In particular, it DOES NOT cover the Nepal glacier source.
    Category integer values are not assigned names without classification metadata.
    """
    out=ctx.folder('q2')
    grid=glacier_grid(ctx)
    demfile=find_first(ctx.raw/'unpacked'/'7','DEM.tif')
    with rasterio.open(demfile) as dem:
        fac=4;ny=int(np.ceil(dem.height/fac));nx=int(np.ceil(dem.width/fac))
        transform=dem.transform*rasterio.Affine.scale(fac,fac);crs=dem.crs
    # The official DEM is EPSG:4326: affine.a*affine.e is in DEGREE^2,
    # not square metres. Compute approximate cell area geodetically at study
    # latitude; never label degree-square values as m².
    if crs.is_geographic:
        from pyproj import Geod
        geod=Geod(ellps='WGS84')
        x0=transform.c+transform.a*nx/2
        y0=transform.f+transform.e*ny/2
        xs=[x0,x0+transform.a,x0+transform.a,x0]
        ys=[y0,y0,y0+transform.e,y0+transform.e]
        cell_area_m2=float(abs(geod.polygon_area_perimeter(xs,ys)[0]))
    else:
        cell_area_m2=float(abs(transform.a*transform.e-transform.b*transform.d))
    assert 100.0<cell_area_m2<1e8, 'landcover grid cell area units or projection invalid'
    layers={}
    for year in [2020,2021,2022,2023]:
        path=find_first(ctx.raw/'unpacked'/'5',f'*{year}.tif')
        target=np.zeros((ny,nx),np.uint8)
        with rasterio.open(path) as src:
            reproject(rasterio.band(src,1),target,src_transform=src.transform,src_crs=src.crs,
                src_nodata=src.nodata,dst_transform=transform,dst_crs=crs,dst_nodata=0,
                resampling=Resampling.nearest)
        layers[year]=target
    valid=np.logical_and.reduce([a>0 for a in layers.values()])
    left,right=layers[2020][valid],layers[2023][valid]
    counts=[]
    for year,array in layers.items():
        good=array>0
        count_vals=pd.Series(array[good]).value_counts()
        for k,v in count_vals.items():
            counts.append({'year':year,'category_code':int(k),'n_grid_cells':int(v),
                'grid_cell_area_m2_approx':cell_area_m2,
                'valid_percent_study_DEM':float(good.mean()*100)})
    pd.DataFrame(counts).to_csv(out/'Q2_official_landcover_four_year_class_counts.csv',index=False)
    change_mask=valid&(layers[2020]!=layers[2023])
    report={'study_DEM_covered_fraction':float((layers[2023]>0).mean()),
      'four_epoch_common_valid_fraction':float(valid.mean()),
      'changed_fraction_among_valid_cells':float(change_mask.sum()/max(1,valid.sum())),
      'source_glacier_has_landcover_coverage':False,
      'landcover_not_extrapolated_into_nepal_source':True,
      'classification_code_semantics_unverified':True,
      'cell_area_m2_at_study_centre_approx':cell_area_m2,
      'area_method':'WGS84 geodesic central-pixel approximation' if crs.is_geographic else 'affine determinant on projected CRS'}
    write_json(out/'Q2_landcover_2020_2023_change_summary.json',report)
    fig,axes=plt.subplots(1,2,figsize=(10.4,4.8),layout='constrained')
    for ax,(year,z) in zip(axes,[(2020,layers[2020]),(2023,layers[2023])]):
        im=ax.imshow(np.ma.masked_equal(z,0),cmap='tab20',vmin=1,vmax=10,interpolation='nearest')
        ax.set(title=str(year)+'年土地覆盖',xticks=[],yticks=[])
    fig.colorbar(im,ax=axes,shrink=.74,label='官方土地覆盖类别码')
    fig.suptitle('2020年与2023年土地覆盖对比')
    savefig(fig,out/'Q2_14_landcover_2020_2023_comparison.png')
    return report
