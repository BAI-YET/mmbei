import numpy as np
import pytest
from mathorcup.q3_time_risk import conditional_event_time,hazard_step
from mathorcup.q34_joint_arrival import joint_arrival_samples,sample_discrete_event_time

def test_first_passage_resampling_preserves_calculated_time_support():
    sampled=sample_discrete_event_time([.5,1.,1.5,2.],[.4,0.,.6,0.],200,seed=12)
    assert set(sampled)=={.5,1.5}

def test_discrete_hazard_probability_normalizes():
    h=np.array([.01,.05,.02,.005,.03])
    p,none,cond=conditional_event_time(h)
    assert np.isclose(p.sum()+none,1)
    assert np.isclose(cond.sum(),1)
    assert 0<none<1

def test_monotone_event_probability_with_intensity():
    h1=hazard_step(.0001,1,1,6)
    h2=hazard_step(.001,1,1,6)
    assert 0<h1<h2<1

def test_time_convolution_preserves_physical_delay():
    a=joint_arrival_samples([1,4,2],[12,18,6])
    assert np.allclose(a,[1.2,4.3,2.1])
    assert np.all(a>=np.array([1,4,2]))

def test_probability_mass_is_checked():
    with pytest.raises(ValueError):sample_discrete_event_time([6,12],[.2,.2],100)
