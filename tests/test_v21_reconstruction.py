"""Analytical fixtures verify numerics; they are never used as case-study inputs."""
import numpy as np
import pandas as pd
import pytest
from mathorcup.q3_task_aligned import empirical_first_passage
from mathorcup.q4_cascade_reconstruction import advance_channel,section_depth,upstream_transfer,simulate_channel
from mathorcup.q4_solver import simulate_2d
from mathorcup.q4_terrain_reconstruction import priority_drainage,trace_drainage
from mathorcup.question_completion import energy_water_channels


def test_empirical_passage_uses_observed_trajectory_counts_without_a_disaster_prior():
    trajectories=np.array([[0,1,1],[0,0,1],[0,0,0],[1,1,1]],float)
    d=empirical_first_passage(trajectories,.5,.5)
    assert np.allclose(d['unconditional_bin_mass'],[.25,.25,.25])
    assert d['p_in_horizon']==.75 and d['nonoccurrence_probability']==.25
    assert np.allclose(d['conditional_bin_mass'],np.ones(3)/3)
    assert np.isclose(d['unconditional_bin_mass'].sum()+d['nonoccurrence_probability'],1)


def test_variable_width_and_bank_slopes_preserve_a_lake_at_rest():
    x=np.arange(41);bed=.2*np.sin(x/5);b=20+4*np.cos(x/7);m=.5+.2*np.sin(x/6)
    h=4-bed;A=b*h+m*h*h;Q=np.zeros_like(A);initial=A.copy()
    for _ in range(80):A,Q,out,clip=advance_channel(A,Q,bed,b,10,.02,0,1e12,bank_side_slope=m)
    assert np.allclose(A,initial,atol=1e-11)
    assert np.max(abs(Q))<1e-10 and abs(out)<1e-12 and clip==0


def test_trapezoid_area_depth_roundtrip_including_rectangular_limit():
    h=np.array([0.,.1,4.,40.]);b=np.array([20.,30.,40.,50.]);m=np.array([0.,1.,2.,3.])
    assert np.allclose(section_depth(b*h+m*h*h,b,m),h)


def test_variable_section_transport_conserves_volume_including_boundary_flux():
    n=45;x=np.arange(n);bed=-x*.05;b=30+x/10;m=np.ones(n)
    h=np.exp(-((x-10)/4)**2);A=b*h+m*h*h;Q=A*2;volume=A.sum()*10;exited=0
    for _ in range(100):
        A,Q,out,clip=advance_channel(A,Q,bed,b,10,.03,.02,2000,bank_side_slope=m);exited+=out
        assert clip<1e-10 and np.all(A>=0)
    assert np.isclose(A.sum()*10+exited,volume,rtol=1e-12)


def test_coupling_cannot_invent_inflow_volume():
    a={'history':pd.DataFrame({'time_s':[0,1,2],'outflow_cumulative_m3':[0,10,30],'outlet_velocity_m_s':[0,2,3]})}
    d=upstream_transfer(a,.6)
    assert d.cumulative_m3.tolist()==[0,6,18]
    assert d.velocity_m_s.tolist()==[0,2,3]
    with pytest.raises(ValueError):upstream_transfer(a,1.1)


def test_priority_drainage_handles_nodata_edges_and_has_no_cycles():
    z=np.array([[5,4,3,2,np.nan],[6,5,1,1,np.nan],[7,6,5,4,3],[8,7,6,5,4]],float)
    original=z.copy();filled,receivers,area=priority_drainage(z,10)
    assert np.array_equal(z,original,equal_nan=True)
    valid=np.isfinite(z)
    for rc in np.argwhere(valid):assert len(trace_drainage(receivers,tuple(rc)))<=valid.sum()
    assert np.all(filled[valid]>=z[valid])
    outlets=valid&(receivers<0)
    assert np.isclose(area[outlets].sum(),valid.sum()*100)


def test_station_energy_equivalent_has_correct_units_and_uses_only_past_data():
    idx=pd.date_range('2025-06-01',periods=200,freq='h')
    h=pd.DataFrame({'temp_2m':10.,'net_shortwave':334000/3600,
                    'energy_balance_residual':-10.,'rain_24h':0.},index=idx)
    a=energy_water_channels(h,4140,4140)
    assert a.station_net_shortwave_7d_mmwe.iloc[-1]==pytest.approx(168)
    assert a.station_positive_energy_residual_7d_mmwe.iloc[-1]==0
    h.iloc[-10:,h.columns.get_loc('net_shortwave')]=9999
    b=energy_water_channels(h,4140,4140)
    assert a.station_net_shortwave_7d_mmwe.iloc[180]==b.station_net_shortwave_7d_mmwe.iloc[180]


def test_coupled_inlet_does_not_add_an_unobserved_pressure_head():
    n=50;s=np.arange(n)*10.
    p=pd.DataFrame({'distance_m':s,'elevation_m':100-s*.01,
                    'effective_width_m':20.,'bank_side_slope':1.,'source_area_weight':0.})
    inlet=pd.DataFrame({'time_s':[0.,10.,50.],'cumulative_m3':[0.,10000.,10000.],
                         'velocity_m_s':[0.,20.,0.]})
    a=simulate_channel(p,.01,2500,total_seconds=50,inlet=inlet)
    assert a['mass']['injected_m3']==pytest.approx(10000)
    assert a['mass']['relative_mass_error']<1e-10
    assert a['mass']['maximum_wet_speed_m_s']<a['mass']['generous_energy_speed_bound_m_s']


def test_2d_injection_uses_the_recorded_cumulative_1d_volume_exactly():
    inlet=pd.DataFrame({'time_s':[0.,1.,4.],'sensor_discharge_m3_s':[0.,100.,0.],
                        'sensor_cumulative_m3':[0.,10.,10.]})
    a=simulate_2d(np.zeros((15,15)),30,(7,7),inlet,total_seconds=4,frames=3)
    assert a['mass']['input_cumulative_m3']==pytest.approx(10,abs=1e-10)
    assert a['mass']['relative_mass_error']<1e-10
