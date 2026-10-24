import numpy as np,pandas as pd
from mathorcup.common import robust_z,quantile_tail
from mathorcup.q2_spatial_mechanism import _cross_validate

def test_robust_z_is_centered():
    ref=list(range(-10,11))
    assert abs(robust_z(0,ref))<1e-12
    assert robust_z(4,ref)>robust_z(2,ref)
    assert 0<quantile_tail(4,ref)<=1

def test_spatial_block_cv_is_out_of_sample():
    rng=np.random.default_rng(44)
    n=84;y=np.r_[np.ones(43),np.zeros(41)].astype(int)
    groups=np.arange(n)%4
    X=rng.normal(size=(n,4));X[:,0]+=y*.8
    metrics,res=_cross_validate(X,y,groups,C=.5)
    probs=res['L2_Logistic_4geo']
    assert probs.shape==(84,)
    assert np.all((probs>=0)&(probs<=1))
    assert len(metrics[metrics.model=='L2_Logistic_4geo'])==4
    assert 0<=res['L2_Logistic_4geo_summary']['auc_oof']<=1
