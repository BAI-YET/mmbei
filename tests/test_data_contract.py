from pathlib import Path
import json
import pytest
ROOT=Path(__file__).resolve().parents[1]

def test_default_settings_prohibit_unverified_claims():
    import yaml
    cfg=yaml.safe_load((ROOT/'configs'/'base.yaml').read_text(encoding='utf-8'))
    assert cfg['timezone_confirmed'] is True
    assert cfg['station_timezone'] == 'Asia/Shanghai'
    assert cfg['external_data_policy'] == 'q2_only'
    assert cfg['q3']['scenario_only'] is True
    assert cfg['q4']['route_verified'] is False

def test_optional_real_data_event_date_not_inferred_from_interpretation_time():
    path=ROOT/'data/processed/historic_gullies.gpkg'
    if not path.exists():pytest.skip('official raw datasets not mounted')
    import geopandas as gpd
    g=gpd.read_file(path)
    assert g.event_time_verified.eq(False).all()
    assert len(g)==84
