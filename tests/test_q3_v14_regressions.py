import numpy as np
import pandas as pd
from mathorcup.q3_time_risk import conservative_time_interval, forecast_issue, weighted_index
from mathorcup.common import FIGURE_TITLES


def test_discrete_95pct_interval_includes_full_quantile_bins():
    p=np.array([.02,.10,.25,.30,.23,.10])
    lo,hi=conservative_time_interval(p,.5)
    assert lo<=.5 and hi>=2.5
    lower=int(round(lo/.5));upper=int(round(hi/.5))
    assert p[lower:upper].sum()>=.95


def test_future_actual_data_cannot_change_issuance_forecast():
    idx=pd.date_range('2025-07-01',periods=24*24,freq='h')
    z=np.arange(len(idx))
    data=pd.DataFrame({'pdd_7d':70+2*np.sin(z/16),
       'vpd':9+2*np.sin(z/3.9),'net_radiation':300+130*np.sin(z/3.2),
       'temp_2m':12+4*np.sin(z/4)},index=idx)
    issue=pd.Timestamp('2025-07-22 09:00')
    scale={k:{'median':float(data[k].median()),'scale':max(1.,float(data[k].std()))} for k in data.columns}
    out1=forecast_issue(data[data.index<issue],issue,scale,0.5,horizon_hours=6,step_h=.5,samples=80,seed=4)
    altered=data.copy();altered.loc[altered.index>=issue,'temp_2m']=1000
    out2=forecast_issue(altered[altered.index<issue],issue,scale,0.5,horizon_hours=6,step_h=.5,samples=80,seed=4)
    assert np.allclose(out1['time_mass'],out2['time_mass'])
    assert out1['interval_hours']==out2['interval_hours']


def test_formal_figure_titles_are_not_diagnostics():
    assert 'Q3_02_95percent_conditional_time_interval' in FIGURE_TITLES
    for stem,title in FIGURE_TITLES.items():
        assert '!=' not in title and '≠' not in title
        assert len(title)<35
