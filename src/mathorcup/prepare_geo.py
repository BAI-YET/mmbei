"""QC of the actual shapefile, DEM, NetCDF grids with provenance retained."""
import numpy as np,pandas as pd,geopandas as gpd,rasterio,h5py
from pathlib import Path
from .common import find_first,write_json

def glacier_grid(ctx):
    f=find_first(ctx.raw/'unpacked'/'4','*.nc')
    with h5py.File(f,'r') as h:
        grids={k:h[k][:] for k in ['dem','slope','thickness','v','vx','vy','glacier_mask','x','y']}
        geom=h['spatial_ref'].attrs['GeoTransform']
        if isinstance(geom,bytes):geom=geom.decode()
    # Convert mask: glacier_mask may have NaNs outside
    grids['glacier_mask']=np.isfinite(grids['glacier_mask']) & (grids['glacier_mask']>0.5)
    grids['transform']=tuple(float(q) for q in geom.split())
    return grids

def prepare_geospatial(ctx):
    out=ctx.processed/'historic_gullies.gpkg'
    if out.exists() and (ctx.processed/'glacier_summary.json').exists():return
    f=find_first(ctx.raw/'unpacked'/'6','*.shp')
    g=gpd.read_file(f,encoding='gb18030')
    if g.crs is None:raise ValueError('Historical polygons lack CRS: do not guess geographic location')
    g=g.to_crs('EPSG:4326')
    g['event_time_verified']=False
    g['debris_type']=g['类型'].astype(str)
    g['gully_id']=g.get('Code',pd.Series(np.arange(len(g)))).astype(str)
    g.to_file(out,driver='GPKG',encoding='UTF-8')
    gg=glacier_grid(ctx)
    area=30.0**2
    m=gg['glacier_mask'] & np.isfinite(gg['thickness'])
    tau=900*9.81*gg['thickness'][m]*np.sin(np.radians(gg['slope'][m]))
    summary={
        'n_gullies':int(len(g)),'types':g['debris_type'].value_counts().to_dict(),
        'historical_gully_dates_are_event_times':False,
        'ice_area_km2':float(m.sum()*area/1e6),
        'ice_total_volume_km3':float(np.sum(gg['thickness'][m])*area/1e9),
        'median_ice_thickness_m':float(np.median(gg['thickness'][m])),
        'median_slope_deg':float(np.median(gg['slope'][m])),
        'median_driving_shear_kpa':float(np.median(tau)/1000),
        'notes':'Total glacier volume is NOT avalanche release volume; glacier velocity 2014-2024 is historical.'}
    write_json(ctx.processed/'glacier_summary.json',summary)
    # DEM geometry evidence and schema
    dem=find_first(ctx.raw/'unpacked'/'7','DEM.tif')
    with rasterio.open(dem) as d:
        write_json(ctx.folder('quality')/'dem_metadata.json',{'width':d.width,'height':d.height,'crs':str(d.crs),
             'pixel_degrees':d.res,'bounds':list(d.bounds),'nodata':d.nodata,
             'source':'Official dataset 7, do not assume pre-disaster surface'})
