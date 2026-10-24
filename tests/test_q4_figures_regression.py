"""Regression protection: Q4 sensitivity must never produce an empty chart when arrivals are absent."""
import numpy as np
import pandas as pd
from pathlib import Path
from mathorcup.q4_solver import simulate_1d, _sensitivity
from mathorcup.common import configure_plotting


def test_nonarrival_still_has_positive_runout_and_max_depth():
    s=np.arange(0,1000,50.)
    p=pd.DataFrame({'distance_m':s,'elevation_m':3200-s*.05})
    a=simulate_1d(p,mu=.3,xi=200,width=30,volume=1000,duration=30,
                  total_seconds=70,frames=3,cfl=.25)
    assert a['mass']['first_arrival_to_downstream_end_s'] is None
    assert a['mass']['max_runout_m']>0
    assert a['mass']['max_depth_over_time_m']>0


def test_sensitivity_chart_nonempty_when_downstream_not_reached(tmp_path):
    configure_plotting()
    s=np.arange(0,2000,50.)
    p=pd.DataFrame({'distance_m':s,'elevation_m':4000-s*.025})
    class FakeContext:
        def __init__(self,root):
            self.root=root
            self.cfg={'q4':{
                'mu':.3,'xi_m_s2':200,'width_m':30.,'release_volume_m3':500.,
                'release_duration_s':20.,'total_seconds':50.,'cfl':.2,'initial_velocity_m_s':0.,
                'scenario_grid':{
                    'mus':[.2,.3,.4],'xi_m_s2':[100,200,300],
                    'volumes_m3':[250,500,750],'widths_m':[20,30,40],
                    'initial_velocities_m_s':[0,5,10]}}}
        def folder(self,name):
            d=self.root/'outputs'/name;d.mkdir(parents=True,exist_ok=True);return d
    ctx=FakeContext(tmp_path)
    df=_sensitivity(ctx,p)
    assert len(df)==18
    assert 'duration' in set(df.varied_parameter)
    assert df.runout_distance_km.notna().all()
    assert (df.runout_distance_km>=0).all()
    assert df.max_depth_m.gt(0).all()
    img=ctx.folder('q4')/'Q4_07_parameter_sensitivity.png'
    assert img.exists() and img.stat().st_size>25000
    grid=ctx.folder('q4')/'Q4_08_volume_friction_phase_diagram.png'
    assert grid.exists() and grid.stat().st_size>15000
