import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mathorcup.redraw_export import export_figure_data


def test_plot_export_retains_dates_intervals_and_masked_numeric_values(tmp_path):
    fig,axs=plt.subplots(1,2)
    dates=pd.date_range('2026-08-26 09:00',periods=3,freq='h')
    axs[0].plot(dates,[1,np.nan,3],label='observed')
    axs[0].fill_between(dates,[0,1,2],[2,3,4],label='interval')
    matrix=np.ma.array([[.1,2.],[3.,4.]],mask=[[True,False],[False,False]])
    axs[1].imshow(matrix)
    result=export_figure_data(fig,Path('q1')/'example.png',tmp_path)
    folder=tmp_path/'q1/example'
    curve=pd.read_csv(folder/'A01_curves.csv')
    assert len(curve)==3 and curve.x_datetime.iloc[0]=='2026-08-26 09:00:00'
    assert pd.isna(curve.y.iloc[1])
    interval=pd.read_csv(next(folder.glob('*intervals.csv')))
    assert np.allclose(interval.y_upper-interval.y_lower,2)
    raster=result['axes'][1]['rasters'][0];a=np.load(tmp_path/raster['shared_data']/'array.npz')
    assert a['values'][0,0]==.1 and a['display_mask'][0,0]
    plt.close(fig)


def test_histogram_and_raster_coordinates_are_reconstructible(tmp_path):
    fig,axs=plt.subplots(1,2)
    axs[0].hist([0,0,1,2,3],bins=[0,1,2,4])
    axs[1].set(xlabel='UTM45N 东向坐标 (m)',ylabel='UTM45N 北向坐标 (m)')
    axs[1].imshow([[1,2],[3,4]],extent=[500000,500120,3100000,3100120],origin='upper')
    result=export_figure_data(fig,Path('q4')/'map.png',tmp_path)
    bars=pd.read_csv(tmp_path/'q4/map/A01_bars.csv')
    assert np.allclose(bars.height,[2,1,2]) and bars.height.sum()==5
    xy=pd.read_csv(tmp_path/'q4/map/A02_I001_coordinates.csv')
    assert np.allclose(xy.x_center,[500030,500090])
    assert np.allclose(xy.y_center,[3100090,3100030])
    import rasterio
    with rasterio.open(tmp_path/'q4/map/A02_I001_georeferenced.tif') as src:
        assert src.crs.to_epsg()==32645
        assert np.allclose(src.read(1),[[1,2],[3,4]])
    plt.close(fig)
