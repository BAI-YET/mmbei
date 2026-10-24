"""Export actual pre-render plotting arrays, not digitized pixels.
CSV values retain full double precision; displayed masks and axis transforms are separate.
"""
from __future__ import annotations
from pathlib import Path
import csv,json,hashlib,os,math,datetime
import numpy as np
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
from matplotlib.collections import PathCollection,PolyCollection
from matplotlib.transforms import Affine2D


def _clean(x):
    if isinstance(x,dict):return {str(k):_clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [_clean(v) for v in x]
    if isinstance(x,np.ndarray):return _clean(x.tolist())
    if isinstance(x,(np.integer,np.floating,np.bool_)):return _clean(x.item())
    if isinstance(x,float) and not math.isfinite(x):return None
    if isinstance(x,(str,int,float,bool,type(None))):return x
    return str(x)


def _json(p,data):p.write_text(json.dumps(_clean(data),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def _csv(p,header,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(header)
        for row in rows:w.writerow(['' if v is None or (isinstance(v,(float,np.floating)) and not np.isfinite(v)) else v for v in row])


def _date_axis(axis):
    conv=axis.get_converter() if hasattr(axis,'get_converter') else getattr(axis,'converter',None)
    return 'date' in type(conv).__name__.lower()


def _date(v,axis):
    if not np.isfinite(v):return ''
    tz=axis.get_units();tz=tz if isinstance(tz,datetime.tzinfo) else datetime.timezone.utc
    return mdates.num2date(float(v),tz=tz).strftime('%Y-%m-%d %H:%M:%S')


def _space(artist,ax):
    tr=artist.get_transform()
    if tr==ax.transData:return 'data','data'
    if tr==ax.transAxes:return 'axes_fraction','axes_fraction'
    if tr==ax.get_xaxis_transform():return 'data','axes_fraction'
    if tr==ax.get_yaxis_transform():return 'axes_fraction','data'
    return 'artist_transform','artist_transform'


def _main_axes(fig):return [(i,ax) for i,ax in enumerate(fig.axes,1) if getattr(ax,'_colorbar',None) is None and ax.get_label()!='<colorbar>']


def _coordinate_info(figstem,ax,extent):
    xlab=ax.get_xlabel();ylab=ax.get_ylabel();maxabs=max(abs(float(extent[0])),abs(float(extent[1])))
    if ('经度' in xlab or 'longitude' in xlab.lower()) and ('纬度' in ylab or 'latitude' in ylab.lower()):
        return 'EPSG:4326',1.0,'degrees'
    if 'UTM' in xlab or ('东向坐标' in xlab and ('北向坐标' in ylab or figstem.startswith('Q4'))):
        factor=1000. if maxabs<10000 and ('km' in xlab or '千米' in xlab) else 1.
        return 'EPSG:32645',factor,'km' if factor==1000 else 'm'
    return None,1.,'plot_axis_units'


def export_figure_data(fig,path,root=None):
    path=Path(path)
    base=Path(root or os.environ.get('MATHORCUP_REDRAW_ROOT',str(path.parent.parent/'redraw_data')))
    folder=base/path.parent.name/path.stem;folder.mkdir(parents=True,exist_ok=True)
    shared=base/'shared_rasters';shared.mkdir(exist_ok=True)
    fig.canvas.draw()
    metadata={'figure':f'{path.parent.name}/{path.name}','title':fig._suptitle.get_text() if fig._suptitle else '',
      'export_stage':'actual Matplotlib artist arrays before savefig; no pixel digitization',
      'csv_encoding':'UTF-8 with BOM','missing_values':'blank/NaN; never converted to zero',
      'time_note':'date strings follow original axis clock; naive station dates are Beijing wall time',
      'axes':[],'datasets':[],'decorative_colorbars':[],'warnings':[]}
    def add(name,kind,n,fields,axis=None,note=''):
        metadata['datasets'].append({'file':str(name),'kind':kind,'rows':int(n),'fields':fields,'axis':axis,'note':note})
    for ai,ax in _main_axes(fig):
        akey=f'A{ai:02d}';datex=_date_axis(ax.xaxis);datey=_date_axis(ax.yaxis)
        axmeta={'axis':akey,'title':ax.get_title(),'xlabel':ax.get_xlabel(),'ylabel':ax.get_ylabel(),
          'xlim':list(ax.get_xlim()),'ylim':list(ax.get_ylim()),'xscale':ax.get_xscale(),'yscale':ax.get_yscale(),
          'x_is_date':datex,'y_is_date':datey,'x_inverted':ax.xaxis_inverted(),'y_inverted':ax.yaxis_inverted(),
          'lines':[],'rasters':[],'collections':[],'patches':[]}
        metadata['axes'].append(axmeta)
        ticks=[]
        for axis,ticklabels,vals in [('x',ax.get_xticklabels(),ax.get_xticks()),('y',ax.get_yticklabels(),ax.get_yticks())]:
            for k,(tick,value) in enumerate(zip(ticklabels,vals)):ticks.append([axis,k,float(value),tick.get_text()])
        _csv(folder/f'{akey}_ticks.csv',['axis','tick_order','plot_coordinate','display_label'],ticks)
        add(f'{akey}_ticks.csv','tick_labels',len(ticks),['axis','tick_order','plot_coordinate','display_label'],akey)
        curves=[]
        for li,line in enumerate(ax.lines,1):
            x=np.ma.filled(np.ma.asarray(line.get_xdata(orig=False),dtype=float),np.nan).ravel()
            y=np.ma.filled(np.ma.asarray(line.get_ydata(orig=False),dtype=float),np.nan).ravel()
            x,y=np.broadcast_arrays(x,y)
            if not len(x):continue
            label=line.get_label();sid=f'L{li:03d}';xs,ys=_space(line,ax)
            axmeta['lines'].append({'series_id':sid,'label':label,'x_space':xs,'y_space':ys,
                'color':line.get_color(),'linestyle':line.get_linestyle(),'marker':line.get_marker(),'linewidth':line.get_linewidth(),'alpha':line.get_alpha()})
            curves.extend([sid,label,k,float(xx),float(yy),_date(xx,ax.xaxis) if datex and xs=='data' else '',_date(yy,ax.yaxis) if datey and ys=='data' else '',xs,ys] for k,(xx,yy) in enumerate(zip(x,y)))
        if curves:
            fields=['series_id','series_label','point_order','x','y','x_datetime','y_datetime','x_space','y_space']
            _csv(folder/f'{akey}_curves.csv',fields,curves);add(f'{akey}_curves.csv','lines_and_markers',len(curves),fields,akey)
        patches=[];patch_vertices=[]
        for pi,patch in enumerate(ax.patches,1):
            pid=f'P{pi:03d}';label=patch.get_label();xs,ys=_space(patch,ax)
            info={'patch_id':pid,'label':label,'type':type(patch).__name__,'facecolor':patch.get_facecolor(),'edgecolor':patch.get_edgecolor(),'alpha':patch.get_alpha(),'hatch':patch.get_hatch(),'x_space':xs,'y_space':ys}
            axmeta['patches'].append(info)
            if isinstance(patch,Rectangle):
                patches.append([pid,label,patch.get_x(),patch.get_y(),patch.get_width(),patch.get_height(),patch.get_angle(),xs,ys])
            verts=patch.get_path().vertices;codes=patch.get_path().codes
            if len(verts):
                data=ax.transData.inverted().transform(patch.get_transform().transform(verts))
                for k,(xx,yy) in enumerate(data):patch_vertices.append([pid,k,int(codes[k]) if codes is not None else '',xx,yy])
        if patches:
            fields=['patch_id','label','x','y','width','height','angle_degree','x_space','y_space']
            _csv(folder/f'{akey}_bars.csv',fields,patches);add(f'{akey}_bars.csv','bars_histogram_bins_and_spans',len(patches),fields,akey,'height is bin count/density exactly used in original plot; y+height gives histogram top')
        if patch_vertices:
            fields=['patch_id','vertex_order','path_code','x_data','y_data']
            _csv(folder/f'{akey}_patch_vertices.csv',fields,patch_vertices);add(f'{akey}_patch_vertices.csv','patch_vertices',len(patch_vertices),fields,akey)
        for ci,col in enumerate(ax.collections,1):
            cid=f'C{ci:03d}';kind=type(col).__name__;arr=col.get_array()
            cmeta={'collection_id':cid,'type':kind,'label':col.get_label(),'cmap':getattr(col.get_cmap(),'name',None),
              'clim':col.get_clim(),'alpha':col.get_alpha()};axmeta['collections'].append(cmeta)
            if isinstance(col,PathCollection):
                offsets=np.ma.filled(np.ma.asarray(col.get_offsets(),dtype=float),np.nan)
                sizes=np.asarray(col.get_sizes()).ravel();rgba=np.asarray(col.get_facecolors())
                a=np.ma.filled(np.ma.asarray(arr,dtype=float),np.nan).ravel() if arr is not None else np.array([])
                rows=[]
                for k,(xx,yy) in enumerate(offsets):
                    color=rgba[min(k,len(rgba)-1)] if len(rgba) else [np.nan]*4
                    rows.append([k,xx,yy,a[min(k,len(a)-1)] if len(a) else None,sizes[min(k,len(sizes)-1)] if len(sizes) else None,*color])
                fields=['point_order','x','y','color_value','marker_size_pt2','red','green','blue','alpha']
                name=f'{akey}_{cid}_scatter.csv';_csv(folder/name,fields,rows);add(name,'scatter_points',len(rows),fields,akey)
                continue
            paths=col.get_paths();vertices=[];bands=[]
            for po,pp in enumerate(paths):
                v=np.asarray(pp.vertices,dtype=float)
                # Most paths are in data coordinates; transform all to an explicit data space.
                try:v=ax.transData.inverted().transform(col.get_transform().transform(v))
                except (ValueError,TypeError):metadata['warnings'].append(f'{akey}/{cid}: original collection path coordinates retained')
                for k,(xx,yy) in enumerate(v):vertices.append([po,k,int(pp.codes[k]) if pp.codes is not None else '',xx,yy])
                if 'FillBetween' in kind and len(v):
                    good=v[np.isfinite(v).all(axis=1)]
                    for xx in np.unique(good[:,0]):
                        yy=good[good[:,0]==xx,1];bands.append([po,xx,float(yy.min()),float(yy.max()),_date(xx,ax.xaxis) if datex else ''])
            if vertices:
                fields=['path_id','vertex_order','path_code','x','y'];name=f'{akey}_{cid}_paths.csv'
                _csv(folder/name,fields,vertices);add(name,'contours_polygons_and_intervals',len(vertices),fields,akey)
            if bands:
                fields=['path_id','x','y_lower','y_upper','x_datetime'];name=f'{akey}_{cid}_intervals.csv'
                _csv(folder/name,fields,bands);add(name,'interval_band_lower_upper',len(bands),fields,akey)
            if arr is not None:
                a=np.ma.filled(np.ma.asarray(arr,dtype=float),np.nan).ravel();name=f'{akey}_{cid}_color_values.csv'
                _csv(folder/name,['item_order','value'],enumerate(a));add(name,'collection_color_values',len(a),['item_order','value'],akey)
            if not vertices and arr is None:metadata['warnings'].append(f'{akey}/{cid}: empty or decorative collection {kind}')
        for ii,im in enumerate(ax.images,1):
            iid=f'I{ii:03d}';masked=np.ma.asarray(im.get_array());values=np.asarray(masked.data);mask=np.ma.getmaskarray(masked)
            values=values.astype('float64');values=np.where(np.isfinite(values),values,np.nan)
            digest=hashlib.sha256(str(values.shape).encode()+values.tobytes()+np.asarray(mask,dtype='uint8').tobytes()).hexdigest()[:20]
            raster_dir=shared/digest;raster_dir.mkdir(exist_ok=True)
            if not (raster_dir/'array.npz').exists():
                np.savez_compressed(raster_dir/'array.npz',values=values,display_mask=mask)
                channels=1 if values.ndim==2 else values.shape[-1]
                for ch in range(channels):
                    matrix=values if channels==1 else values[...,ch]
                    mf=mask if mask.ndim==2 else mask[...,ch]
                    suffix='' if channels==1 else f'_ch{ch}'
                    with (raster_dir/f'values{suffix}.csv').open('w',encoding='utf-8-sig',newline='') as f:
                        w=csv.writer(f);w.writerow(['row_index']+[f'col_{j}' for j in range(matrix.shape[1])])
                        for r,row in enumerate(matrix):w.writerow([r]+['' if not np.isfinite(v) else repr(float(v)) for v in row])
                    if np.any(mf):
                        with (raster_dir/f'display_mask{suffix}.csv').open('w',encoding='utf-8-sig',newline='') as f:
                            w=csv.writer(f);w.writerow(['row_index']+[f'col_{j}' for j in range(matrix.shape[1])]);w.writerows([r]+list(row.astype(int)) for r,row in enumerate(mf))
                _json(raster_dir/'array_metadata.json',{'shape':list(values.shape),'units':'given by per-figure axis/colorbar metadata','mask_definition':'1=hidden by original plot; values still retained'})
            ext=list(map(float,im.get_extent()));nr,nc=values.shape[:2]
            xx=ext[0]+(np.arange(nc)+.5)*(ext[1]-ext[0])/nc
            yy=(ext[3]+(np.arange(nr)+.5)*(ext[2]-ext[3])/nr) if im.origin=='upper' else (ext[2]+(np.arange(nr)+.5)*(ext[3]-ext[2])/nr)
            name=f'{akey}_{iid}_coordinates.csv';n=max(nr,nc)
            _csv(folder/name,['index','x_center','y_center'],([k,float(xx[k]) if k<nc else None,float(yy[k]) if k<nr else None] for k in range(n)))
            add(name,'raster_coordinate_vectors',n,['index','x_center','y_center'],akey)
            crs,factor,unit=_coordinate_info(path.stem,ax,ext)
            rmeta={'image_id':iid,'raster_id':digest,'shared_data':f'shared_rasters/{digest}',
              'shape':list(values.shape),'extent':ext,'origin':im.origin,'coordinates':name,'cmap':im.get_cmap().name,
              'norm':type(im.norm).__name__,'vmin':im.norm.vmin,'vmax':im.norm.vmax,'alpha':im.get_alpha(),
              'interpolation':im.get_interpolation(),'raw_coordinate_unit':unit,'GIS_coordinate_factor':factor,'GIS_crs':crs}
            cbar=getattr(im,'colorbar',None)
            if cbar is not None:rmeta['colorbar_label']=cbar.ax.get_ylabel() or cbar.ax.get_xlabel()
            if crs and values.ndim==2:
                import rasterio
                from affine import Affine
                yedge=ext[3] if im.origin=='upper' else ext[2]
                step=(ext[2]-ext[3])/nr if im.origin=='upper' else (ext[3]-ext[2])/nr
                tr=Affine((ext[1]-ext[0])/nc*factor,0,ext[0]*factor,0,step*factor,yedge*factor)
                name=f'{akey}_{iid}_georeferenced.tif'
                with rasterio.open(folder/name,'w',driver='GTiff',height=nr,width=nc,count=1,dtype='float64',crs=crs,transform=tr,nodata=np.nan,compress='deflate') as ds:ds.write(values,1)
                rmeta['GeoTIFF']=name;add(name,'georeferenced_full_raster',nr*nc,['raster_values'],akey,'original numeric values; display mask is a separate file')
            axmeta['rasters'].append(rmeta)
            add(f'shared_rasters/{digest}/values.csv' if values.ndim==2 else f'shared_rasters/{digest}/values_ch0.csv','raster_numeric_matrix',nr*nc,['row_index','col_*'],akey)
        texts=[]
        for ti,t in enumerate(ax.texts,1):
            pos=t.get_position();xs,ys=_space(t,ax);texts.append([ti,pos[0],pos[1],t.get_text(),xs,ys])
        if texts:
            fields=['text_id','x','y','text','x_space','y_space'];name=f'{akey}_annotations.csv'
            _csv(folder/name,fields,texts);add(name,'numeric_annotations_and_labels',len(texts),fields,akey)
    for ax in fig.axes:
        cb=getattr(ax,'_colorbar',None)
        if cb is not None:metadata['decorative_colorbars'].append({'label':ax.get_xlabel() or ax.get_ylabel(),'ticks':cb.get_ticks(),'vmin':cb.norm.vmin,'vmax':cb.norm.vmax,'norm':type(cb.norm).__name__})
    _json(folder/'figure_metadata.json',metadata)
    return metadata
