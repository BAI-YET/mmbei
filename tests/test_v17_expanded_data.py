"""Critical regressions for the expanded official variables and the Q1--Q4 data contracts."""
import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from mathorcup.prepare_time import _RAW_FIELDS,_clean
from mathorcup.features import make_hourly
from mathorcup.q1_event_signatures import CANDIDATES,_distinct_top5
from mathorcup.q2_trigger_evaluation import dynamic_conditions


def test_all_official_physical_fields_have_an_ingestion_mapping():
    root=Path(__file__).resolve().parents[1]/'data/raw/unpacked/1'
    required=[root/f'{group}_Gyirong_2025.csv' for group in ['MET','RADM','SOIL','FLUX']]
    if not all(path.exists() for path in required):
        pytest.skip('Official station attachments are not included in the GitHub repository')
    for group in ['MET','RADM','SOIL','FLUX']:
        header=pd.read_csv(root/f'{group}_Gyirong_2025.csv',nrows=0).columns
        raw={c for c in header if c!='Timestamp' and not c.startswith('QC_')}
        mapped=set(_RAW_FIELDS[group])
        assert raw==mapped, f'{group} unmapped={raw-mapped} invented={mapped-raw}'


def test_wind_direction_is_angular_not_linear_temperature():
    times=pd.date_range('2025-08-20',periods=2,freq='h')
    d=pd.DataFrame({'timestamp_raw':times,'wind_2m':[4.,5.],'wind_20m':[6.,7.],
      'wd_2m':[0.,360.],'wd_20m':[90.,90.], 'temp_2m':[12.,12.], 'temp_20m':[13.,13.],
      'rh_2m':[30.,30.],'rh_20m':[20.,20.]})
    h=make_hourly(d)
    assert np.isfinite(h['vector_shear_2_20m']).all()
    assert np.isclose(h['vector_shear_2_20m'].iloc[0], np.hypot(6,4)/18)


def test_radiation_components_and_soil_depth_derivations():
    t=pd.date_range('2025-08-20',periods=2,freq='h')
    d=pd.DataFrame({'timestamp_raw':t,'sw_down':[300.,300.],'sw_up':[40.,40.],
    'lw_down':[250.,250.],'lw_up':[210.,210.],
    'soilT_0p1':[1.,1.],'soilT_0p2':[2.,2.],'soilT_0p4':[3.,3.],
    'soilT_0p8':[4.,4.],'soilT_1p6':[5.,5.],
    'soil_heat':[10.,10.],'soil_heat_0p2':[8.,8.]})
    h=make_hourly(d)
    assert np.isclose(h.net_radiation.iloc[0],300)
    assert np.isclose(h.net_shortwave.iloc[0],260)
    assert np.isclose(h.net_longwave.iloc[0],40)
    assert np.isclose(h.soil_heat_gradient.iloc[0],20)


def test_expanded_q1_has_all_physical_groups_without_fetch_co2_diagnostics():
    assert len(CANDIDATES)>55
    for name in ('tau','tstar','soilT_1p6','vwc_1p6','vector_shear_2_20m','net_longwave'):
        assert name in CANDIDATES
    assert not any('fetch_' in x or x=='co2_flux' for x in CANDIDATES)


def test_q2_no_rain_never_means_rain_triggered():
    idx=pd.date_range('2022-05-01',periods=24*150,freq='h')
    rng=np.random.default_rng(12)
    rain=rng.choice([0.,0.,0.,0.,0.,1.5],size=len(idx))
    d=pd.DataFrame({'temp_2m':np.full(len(idx),12.),'precip_mm':rain,
      'rain_24h':pd.Series(rain,index=idx).rolling(24,min_periods=18).sum().values,
      'pdd_7d':np.full(len(idx),12.)},index=idx)
    h,info=dynamic_conditions(d,baseline_years=(2022,))
    zeros=h.loc[h.rain_24h.eq(0)]
    assert not zeros.rain_condition_met.any()
