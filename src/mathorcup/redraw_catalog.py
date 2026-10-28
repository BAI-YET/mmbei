"""Assemble per-figure import tables, units, native spatial sources and redraw guide."""
from pathlib import Path
import csv,json,gzip,shutil,hashlib
import numpy as np
import pandas as pd
from .common import FIGURE_TITLES
from .redraw_export import _csv,_json


def _layer_name(stem,axis_order,image_order,meta):
    cmap=meta['cmap'];lab=meta.get('colorbar_label','');unit=lab
    if stem=='Q1_00_data_coverage':return '观测有效比例','比例0-1'
    if stem in ('Q1_04_multiscale_anomaly_heatmap','Q1_12_all_factor_multiscale'):return '有向MAD Z','无量纲'
    if stem=='Q1_13_all_factor_redundancy':return 'Spearman相关系数','无量纲[-1,1]'
    if stem=='Q1_05_reference_year_sensitivity':return '因素排序名次','名次'
    if stem=='Q2_05_glacier_physical_indicators':return [('冰厚','m'),('历史冰川速度','m/a'),('重力驱动剪应力代理','kPa')][axis_order-1]
    if stem=='Q2_06_dem_landcover_history_overlay':return ('地形高程','m') if axis_order==1 else ('土地覆盖类别码','整数类别；类别名称未核实')
    if stem=='Q2_07_glacier_historic_thinning_optional':return '五年期冰川高程变化率','m/a'
    if stem=='Q2_09_source_spatial_condition_index':return [('冰川坡度','degree'),('驱动剪应力代理','kPa'),('静态物源条件指数','无量纲0-1')][axis_order-1]
    if stem=='Q2_13_glacier_velocity_and_release_candidates':return [('历史冰川速度','m/a'),('历史速度梯度','1/a'),('物源条件指数','无量纲0-1')][axis_order-1]
    if stem=='Q2_14_landcover_2020_2023_comparison':return '土地覆盖类别码','整数类别；0为无效覆盖'
    if stem=='Q3_05_spatial_warning_distance':return '地形阴影底图','无量纲'
    if stem in ('Q4_00_actual_DEM_geographic_context','Q4_00B_local_true_DEM_detail'):
        return ('地形高程','m') if image_order==1 else ('地形阴影底图','无量纲')
    if stem=='Q4_01_route_geolocation_quality_gate':return '地形高程','m'
    if stem=='Q4_02_1d_runout_timespace':return '基准一维流深','m'
    if stem=='Q4_08_volume_friction_phase_diagram':return ('最大传播距离','km') if axis_order==1 else ('下游首达时间','min；空值=未到达')
    if stem=='Q4_10_impact_pressure_distribution':return ('地形阴影底图','无量纲') if image_order==1 else ('全时段峰值动压代理','kPa')
    if stem.startswith(('Q4_05','Q4_06','Q4_11')):
        if cmap in ('Greys','Greys_r','terrain'):return '地形高程底图','m'
        if cmap in ('gray','grey'):return '地形阴影底图','无量纲'
        if stem.startswith('Q4_05') and axis_order==2:return '首次达到湿润阈值的时间','min，相对源区释放'
        return '全时段最大流深' if stem.startswith('Q4_05') else '指定阶段流深','m'
    return (lab or meta['image_id'],unit or '见原图色标及变量说明')


def _add_wide_curves(folder,d):
    for ax in d['axes']:
        stem=ax['axis'];f=folder/f'{stem}_curves.csv'
        if not f.exists():continue
        c=pd.read_csv(f,keep_default_na=False,float_precision='round_trip');groups={}
        for sid,g in c.groupby('series_id',sort=False):
            key=(tuple(g.x),tuple(g.x_space),tuple(g.y_space))
            groups.setdefault(key,[]).append((sid,g))
        d['datasets']=[x for x in d['datasets'] if x['kind']!='wide_aligned_curves']
        for k,items in enumerate(groups.values(),1):
            first=items[0][1]
            header=['x','x_datetime','x_space']+[f'y_{sid}' for sid,_ in items]
            rows=[[row.x,row.x_datetime,row.x_space]+[g.iloc[i].y for _,g in items] for i,row in enumerate(first.itertuples(index=False))]
            name=f'{stem}_curves_wide_{k:02d}.csv';_csv(folder/name,header,rows)
            d['datasets'].append({'file':name,'kind':'wide_aligned_curves','rows':len(rows),'fields':header,'axis':stem,
                'note':'同一X网格的一列X、多列Y，可直接导入Origin/Excel；列名中的L编号对应figure_metadata.lines'})


def _native_geo(ctx,folder,stem,axorder,im):
    """Add authoritative georeferencing to maps originally plotted in row/column space."""
    source_grid=stem in ('Q2_07_glacier_historic_thinning_optional','Q2_13_glacier_velocity_and_release_candidates')
    landcover=stem=='Q2_14_landcover_2020_2023_comparison'
    if not source_grid and not landcover:return
    import rasterio
    from affine import Affine
    if source_grid:
        from .prepare_geo import glacier_grid
        g=glacier_grid(ctx);tr=Affine.from_gdal(*g['transform']);crs='EPSG:32645'
    else:
        from .common import find_first
        p=find_first(ctx.raw/'unpacked/7','DEM.tif')
        with rasterio.open(p) as src:tr=src.transform*Affine.scale(4,4);crs=src.crs
    data=np.load(ctx.out/'redraw_data'/im['shared_data']/'array.npz')['values']
    if data.ndim!=2:return
    name=f'{axorder}_{im["image_id"]}_native_GIS.tif'
    with rasterio.open(folder/name,'w',driver='GTiff',height=data.shape[0],width=data.shape[1],count=1,dtype='float64',
          crs=crs,transform=tr,nodata=np.nan,compress='deflate') as dest:dest.write(data,1)
    im['GeoTIFF_native']=name;im['GIS_crs_native']=str(crs)
    im['GIS_native_transform']=list(tuple(tr)[:6]);im['GIS_note']='坐标来自实际建模栅格，原图像素坐标保留在coordinates.csv'


def build_redraw_catalog(ctx):
    base=ctx.out/'redraw_data';base.mkdir(exist_ok=True)
    source_tables=base/'source_tables';source_tables.mkdir(exist_ok=True)
    source_spatial=base/'source_spatial';source_spatial.mkdir(exist_ok=True)
    for old in source_spatial.iterdir():
        if old.is_file():old.unlink()
    # Preserve all original quantitative files in the main archive; provide common imports in this annex.
    for q in ['q1','q2','q3','q4','joint','quality']:
        target=source_tables/q
        if target.exists():
            for old in target.iterdir():
                if old.is_file():old.unlink()
        target.mkdir(exist_ok=True)
        for p in (ctx.out/q).glob('*'):
            if not p.is_file() or not (p.suffix in ('.csv','.json','.npz','.tif','.gpkg') or p.name.endswith('.csv.gz')):continue
            dest=source_tables/q/p.name;dest.parent.mkdir(exist_ok=True)
            if p.name.endswith('.csv.gz'):
                dest=dest.with_suffix('');content=gzip.decompress(p.read_bytes())
                if not content.startswith(b'\xef\xbb\xbf'):content=b'\xef\xbb\xbf'+content
                dest.write_bytes(content)
            else:shutil.copy2(p,dest)
    # CSV versions of the 84 raw gully records support alternate boxplots and bins.
    import geopandas as gpd
    g=gpd.read_file(ctx.processed/'historic_gullies.gpkg');flat=pd.DataFrame(g.drop(columns=g.geometry.name));flat['geometry_WKT_EPSG4326']=g.geometry.to_wkt()
    flat.to_csv(source_tables/'Q2_84_gullies_all_attributes.csv',index=False,encoding='utf-8-sig')
    shutil.copy2(ctx.processed/'historic_gullies.gpkg',source_spatial/'historic_gullies.gpkg')
    for q in ['q2','q3']:
        for p in (ctx.out/q).glob('*'):
            if p.suffix in ('.tif','.gpkg'):shutil.copy2(p,source_spatial/p.name)
    nc=list((ctx.raw/'unpacked/4').rglob('*.nc'))
    if nc:shutil.copy2(nc[0],source_spatial/nc[0].name)
    index=[];files=[];dictionary=[];raster_index=[];warnings=[]
    originals=[p for p in ctx.out.rglob('*.png') if 'redraw_data' not in p.parts]
    active_previews={str(p.relative_to(ctx.out)) for p in originals}
    for old in (base/'previews').rglob('*.png'):
        if str(old.relative_to(base/'previews')) not in active_previews:old.unlink()
    order={'q1':1,'q2':2,'q3':3,'q4':4,'joint':5};originals.sort(key=lambda p:(order.get(p.parent.name,9),p.name))
    for number,p in enumerate(originals,1):
        folder=base/p.parent.name/p.stem;mp=folder/'figure_metadata.json'
        if not mp.exists():raise FileNotFoundError(f'Missing full plotting arrays for {p.name}')
        d=json.loads(mp.read_text());_add_wide_curves(folder,d)
        _symbolic_and_interval_notes(ctx,folder,p.stem,d)
        for k,ax in enumerate(d['axes'],1):
            if not ax['title'] and p.stem=='Q4_05_local_2d_depth_and_arrival':ax['title']='最大流深' if k==1 else '首次淹没时间'
            if not ax['title'] and p.stem=='Q4_06C_2x2_inlet_boundary_factorial':ax['title']=['侧坡+静止入流','侧坡+沿沟动量','谷底+静止入流','谷底+沿沟动量'][k-1]
            for i,im in enumerate(ax['rasters'],1):
                label,unit=_layer_name(p.stem,k,i,im);im['quantity']=label;im['value_unit']=unit
                _native_geo(ctx,folder,p.stem,ax['axis'],im)
                _display_geotiff(base,folder,im)
                raster_index.append([p.stem,ax['axis'],im['image_id'],label,unit,im['shared_data']+'/values.csv',im['coordinates'],
                    im.get('GeoTIFF_native',im.get('GeoTIFF','')),im.get('GIS_crs_native',im.get('GIS_crs','')),im['shape'][0],im['shape'][1],im['vmin'],im['vmax'],im['norm']])
            for field,meaning,unit in [('x','曲线横坐标',ax['xlabel'] or '见面板刻度/图题'),('y','曲线纵坐标',ax['ylabel'] or '见面板刻度/图题'),
                 ('x_datetime','原图横坐标时钟日期','YYYY-MM-DD HH:mm:ss；站点采用北京时间'),('y_lower/y_upper','阴影区间下限/上限',ax['ylabel'] or '同面板Y单位'),
                 ('height/width','柱形或直方图区间高度/宽度','参照XY轴；不是统一的计数单位'),('x_center/y_center','栅格像元中心坐标','原图坐标；GIS原生坐标另见GeoTIFF')]:
                dictionary.append([p.stem,ax['axis'],field,meaning,unit])
        d['source_data_location']='source_tables/<question>';d['note']='原始图值完整保存；mask、缺测、未到达与实测/情景含义均保留。'
        _json(mp,d)
        allf=sorted([f for f in folder.iterdir() if f.is_file()])
        for f in allf:
            files.append([p.stem,str(f.relative_to(base)),f.suffix,f.stat().st_size,hashlib.sha256(f.read_bytes()).hexdigest()])
        kinds=sorted(set(x['kind'] for x in d['datasets']))
        example=next((x['file'] for x in d['datasets'] if x['kind']=='wide_aligned_curves'),next((x['file'] for x in d['datasets'] if x['kind']!='tick_labels'),'figure_metadata.json'))
        index.append([number,p.parent.name.upper(),p.stem,d['title'] or FIGURE_TITLES.get(p.stem,d['axes'][0]['title']),
          str(folder.relative_to(base)),len(d['axes']),len(d['datasets']),len(allf),'；'.join(kinds),example,
          '见原样保留的PNG/SVG','完整保存'])
        warnings.extend([p.stem,x] for x in d['warnings'])
        preview=base/'previews'/p.parent.name;preview.mkdir(parents=True,exist_ok=True)
        for src in [p,p.with_suffix('.svg')]:
            if src.exists():shutil.copy2(src,preview/src.name)
    _csv(base/'00_重绘数据索引.csv',['序号','问题','图像文件名','图题','数据目录','面板数','数据集数','本图目录文件数','数据类型','建议先打开的数据表','原图','保留状态'],index)
    _csv(base/'01_字段与单位.csv',['图像','面板','字段','含义','单位或定义'],dictionary)
    _csv(base/'02_空间栅格图层.csv',['图像','面板','图层','变量','数值单位','共用矩阵','原图像元坐标表','本图GeoTIFF','坐标系','行数','列数','显示下限','显示上限','色标类型'],raster_index)
    _csv(base/'03_本图数据文件目录.csv',['图像','数据文件','格式','字节数','SHA256'],files)
    _csv(base/'04_导出说明与提示.csv',['图像','说明'],warnings)
    readme='''绘图数据保留与重绘说明

本附件覆盖63张现有图像。数值从实际绘图对象及模型输出直接导出，没有从PNG反推。

1. 先打开00_重绘数据索引.csv，按图像文件名定位数据目录。CSV采用UTF-8 BOM，Excel/Origin可直接导入。
2. 曲线优先使用A##_curves_wide_##.csv：一列X、多列Y；A##_curves.csv是完整长表。series_id/列名L编号与figure_metadata.json中的图例对应。
3. 日期曲线使用x_datetime。x是Matplotlib内部日期数值，不能直接当作Excel日期。原始站点时钟为北京时间。
4. 柱形和直方图用A##_bars.csv：x、y是矩形起点，width、height是实际绘图宽高。横向柱图数值一般是width，纵向柱图和直方图一般是height。类别标签见A##_ticks.csv。可以用source_tables里的原始数据另选分箱。
5. 置信/分位阴影用*_intervals.csv的y_lower与y_upper。多边形、等高线和盒图形状的完整顶点保存在*_paths.csv或*_patch_vertices.csv中，点序及path_code必须保留。
6. 散点用*_scatter.csv。color_value缺失表示原图直接指定颜色；RGBA列给出原图显示颜色，并不代表新的科学变量。
7. 热图和栅格数值放在shared_rasters/<编号>/values.csv。第一列是row_index，其余列为网格数值，不要把row_index当作数据列。显示掩膜display_mask为1表示原图隐藏该像元，隐藏前的数值仍保留；缺测为空，不填0。
8. 栅格坐标、单位、行序、像元中心、色标及截断上下限见02_空间栅格图层.csv与figure_metadata.json。色标截断不删除超限原始值。对数色标仍保存原始数值，并未预先取对数。
9. GeoTIFF中名称带_display者按原图掩膜显示，其他GeoTIFF保留完整原始数值。GIS图优先用GeoTIFF及source_spatial中的GeoPackage/NetCDF。UTM45N为EPSG:32645；经纬度为EPSG:4326。原图中某些轴以米存储、刻度按千米显示，不能据标签自行将原始米坐标再乘1000。
10. 同一数值矩阵在多图重复使用时只在shared_rasters保存一次，每张图均有明确引用。重绘时保留整个redraw_data目录，不要只复制某张图的子目录。
11. source_tables保留原有CSV/JSON/NPZ等；csv.gz已另解压为CSV。previews中有对应原PNG和SVG，方便比照。数组NPZ同时保留，普通制图不需要先编程解包NPZ。
12. Q4参数敏感性图中，叉号y是为展示“未到达”设置的摆放高度，并非到达时间；对应散点CSV用marker_role说明，actual_arrival_time_min保持空值。Q3_01阴影为得分的5%-95%范围，不能当作95%触发时间区间。
13. 这些数据保留原先证据含义。气象先兆覆盖率不是泥石流实测覆盖率；未到达的时间保持缺失；80百万m³二维情景与50百万m³基准分开，不将两者混用。

代码：src/mathorcup/redraw_export.py在保存图像前导出数值，redraw_catalog.py生成数据索引。可运行python scripts/export_redraw_data.py在隔离目录回放导出，并核对原结果不变。
'''
    readme=readme.replace('覆盖63张',f'覆盖{len(index)}张')
    (base/'使用说明.txt').write_text(readme,encoding='utf-8-sig')
    _json(base/'catalog_summary.json',{'figures':len(index),'panels':sum(x[5] for x in index),'plot_datasets':sum(x[6] for x in index),
       'per_figure_files':len(files),'unique_shared_rasters':len(list((base/'shared_rasters').glob('*'))),'raster_layers':len(raster_index),
       'warnings':warnings,'retained_original_results':True,'source_tables':len(list(source_tables.rglob('*.*')))})
    return index


def _display_geotiff(base,folder,im):
    import rasterio
    a=np.load(base/im['shared_data']/'array.npz')
    values=a['values'];mask=a['display_mask']
    if values.ndim!=2:return
    for field in ['GeoTIFF','GeoTIFF_native']:
        if not im.get(field):continue
        srcpath=folder/im[field];name=srcpath.stem+'_display.tif'
        with rasterio.open(srcpath) as src:profile=src.profile
        with rasterio.open(folder/name,'w',**profile) as dst:dst.write(np.where(mask,np.nan,values),1)
        im[field+'_display']=name


def _symbolic_and_interval_notes(ctx,folder,stem,d):
    if stem=='Q3_01_precursor_and_scenario_hazard':
        d['interval_note']='预报得分阴影为相似日集合5%-95%分位范围；不是触发时刻95%区间'
    if stem=='Q3_10_archived_forecast_sensitivity':
        d['interval_note']='联合气象条件得分集合2.5%-97.5%分位范围；不是灾害时间区间'
    if stem!='Q4_07_parameter_sensitivity':return
    source=pd.read_csv(ctx.out/'q4/Q4_parameter_sensitivity.csv')
    axes=d['axes'][2];counter=0
    for key,g in source.groupby('varied_parameter',sort=False):
        reached=g.reached_downstream.astype(bool)
        roles=[('actual_arrival_time',g.loc[reached])]
        if (~reached).any():roles.append(('symbolic_nonarrival',g.loc[~reached]))
        for role,subset in roles:
            counter+=1;cid=f'C{counter:03d}';p=folder/f"{axes['axis']}_{cid}_scatter.csv"
            if not p.exists():raise ValueError(f'Missing arrival scatter {p}')
            values=pd.read_csv(p,float_precision='round_trip')
            if len(values)!=len(subset):raise ValueError('Scatter count no longer matches parameter rows')
            values['varied_parameter']=key;values['parameter_cn']=g.parameter_cn.iloc[0]
            values['marker_role']=role
            values['actual_arrival_time_min']=values.y if role=='actual_arrival_time' else np.nan
            values.to_csv(p,index=False,encoding='utf-8-sig')
            axes['collections'][counter-1]['marker_role']=role
            axes['collections'][counter-1]['varied_parameter']=key
    d['arrival_marker_note']='symbolic_nonarrival的y只供重现叉号位置；真实到达时间见actual_arrival_time_min，保持空值'
