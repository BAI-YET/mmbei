import pandas as pd,numpy as np
from mathorcup.prepare_time import _clean
from mathorcup.features import make_hourly

def test_quality_code_and_sentinel_are_excluded():
    d=pd.DataFrame({'Timestamp':['2026-08-25 01:00','2026-08-25 02:00','2026-08-25 03:00'],
      'Ta_2m(degree C)':[12,9999.9,14],'QC_Ta_2m':[0,0,8],
      'Precip(mm)':[0,9999.9,-.1]})
    x,qa=_clean(d,'MET')
    assert x.temp_2m.notna().sum()==1
    assert x.precip_mm.notna().sum()==1
    assert np.isclose(x.temp_2m.iloc[0],12)

def test_pdd_7d_rolling_only_looks_backwards():
    times=pd.date_range('2026-01-01',periods=200,freq='h')
    x=pd.DataFrame({'timestamp_raw':times,'temp_2m':[5.0]*199+[30.0],
                   'precip_mm':[0.0]*199+[15.0]})
    f=make_hourly(x)
    snap=times[-1]
    before=f.loc[f.index<snap]
    assert before.pdd_hour.max()==5/24
    assert before.precip_mm.sum()==0
    assert np.isnan(f.pdd_7d.iloc[100])
    assert f.pdd_7d.iloc[-1]>before.pdd_7d.iloc[-1]
