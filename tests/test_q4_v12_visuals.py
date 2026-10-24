"""Q4 v1.2: regression guard against blank or inconsistent flow snapshots."""
from pathlib import Path
import numpy as np
import pandas as pd
from mathorcup.q4_solver import choose_inlet_on_valley
from mathorcup.common import Context

def test_valley_inlet_weights_and_dem_direction():
    # Artificial side-slope with north-west channel; source lies south-east of thalweg.
    ii,jj=np.mgrid[:17,:19]
    z=(ii+jj)*7 + 3*(ii-jj)**2
    pos=choose_inlet_on_valley(z,(8,10),flow_direction_xy=(-.7,-.7),dx=60)
    assert len(pos['cells'])>=1
    assert np.isclose(np.sum(pos['weights']),1.0,rtol=0,atol=1e-12)
    assert np.isfinite(pos['elevation_snapped_m'])
    assert pos['elevation_snapped_m']<=pos['elevation_nominal_m']+5

def test_emitted_q4_snapshot_has_ten_nonempty_event_centred_frames():
    ctx=Context(Path(__file__).resolve().parents[1])
    p=ctx.out/'q4'/'Q4_2d_snapshot_metrics.csv'
    if not p.exists():return
    d=pd.read_csv(p)
    assert len(d)==10
    assert np.all(np.diff(d.time_s.to_numpy())>0)
    assert np.all(d.wet_area_km2>0), 'Q4_06 must not waste a map panel on an empty pre-inlet frame'
    assert np.all(d.max_depth_m>0)
    assert d.time_since_2d_onset_s.iloc[0]>0

def test_q4_updated_and_legacy_both_conserve_mass():
    import json
    ctx=Context(Path(__file__).resolve().parents[1])
    p=ctx.out/'q4'/'Q4_2d_inlet_ablation.json'
    if not p.exists():return
    ans=json.loads(p.read_text())
    for k in ('legacy','valley_aligned'):
        assert ans[k]['relative_mass_error']<1e-6
        assert ans[k]['input_cumulative_m3']>0
