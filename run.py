#!/usr/bin/env python
"""MathorCup B 2026 end-to-end run. Default Q3/Q4 are SCENARIOS until independent evidence arrives."""
from pathlib import Path
import argparse,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context,configure_plotting,record_warning
from mathorcup.catalog import prepare_official_data

def main():
    p=argparse.ArgumentParser(description='2026 MathorCup B reproducible experiments')
    p.add_argument('--stage',choices=['all','prepare','q1','q2','q3','q4','tests'],default='all')
    p.add_argument('--fetch-q2-external',action='store_true',help='force download the real Q2 archived forecast')
    p.add_argument('--no-fetch-q2-external',action='store_true',help='disable Q2-only automatic external download when offline')
    p.add_argument('--force-rebuild',action='store_true',help='invalidate all processed official-data caches and regenerate from original files')
    p.add_argument('--official-dir',default=None,help='directory containing the official ZIP packages')
    p.add_argument('--config',default=None)
    p.add_argument('--allow-provisional-q4',action='store_true',help='legacy alias: run a physical scenario on the candidate DEM corridor')
    p.add_argument('--strict-route',action='store_true',help='require a checked route for the reference propagation calculation')
    p.add_argument('--refit-q4',action='store_true',help='repeat Q4 calibration experiments instead of reusing matching fitted parameters')
    args=p.parse_args();ctx=Context(ROOT,args.config);configure_plotting()
    ctx.cfg['refit_q4']=args.refit_q4
    # Each recomputed stage owns its output directory: never mix old and new figures.
    import shutil
    stages=['q1','q2','q3','q4','joint'] if args.stage=='all' else ([args.stage] if args.stage in ('q1','q2','q3','q4') else [])
    for stage in stages:
        target=ctx.out/stage
        if target.exists():shutil.rmtree(target)
        redraw=ctx.out/'redraw_data'/stage
        if redraw.exists():shutil.rmtree(redraw)
    # Q2 ALONE is authorised to obtain supplementary internet products. Other stages
    # reuse Q2 and official data and perform no network access.
    if not args.no_fetch_q2_external and (args.fetch_q2_external or args.stage in ('all','q2')):
        from mathorcup.q2_external import download_q2_archived_forecast,load_q2_forecast
        existing,status=load_q2_forecast(ctx)
        if args.fetch_q2_external or status!='ACTUAL_DOWNLOADED_ARCHIVED_FORECAST':
            try:
                df=download_q2_archived_forecast(ctx,timeout=12)
                print('[Q2] downloaded verified archived ECMWF run:',len(df),'forecast hours')
            except Exception as e:
                record_warning('Q2-only public weather archive could not be downloaded: '+str(e))
                record_warning('No atmospheric model hindcast is substituted or fabricated; real Sentinel-1 evidence is separately retained.')
        if args.fetch_q2_external and args.stage=='prepare':return
    if args.force_rebuild:
        for name in ['features_hourly.csv.gz','features_metadata.json','obs_hourly.csv.gz','temp_county_daily.csv.gz',
                     'historic_gullies.gpkg','glacier_summary.json']:
            (ctx.processed/name).unlink(missing_ok=True)
    if not ctx.cfg['timezone_confirmed']:
        record_warning('Original station timestamp timezone UNCONFIRMED; raw clock anchor is illustrative.')
    if args.stage not in ('tests',):
        prepare_official_data(ctx,args.official_dir)
        from mathorcup.prepare_time import prepare_timeseries
        from mathorcup.prepare_geo import prepare_geospatial
        from mathorcup.features import build_features
        prepare_timeseries(ctx)
        prepare_geospatial(ctx)
        build_features(ctx)
    if args.stage in ('all','q1'):
        from mathorcup.q1_event_signatures import execute_q1
        execute_q1(ctx)
        from mathorcup.question_completion import complete_q1
        complete_q1(ctx)
    if args.stage in ('all','q2'):
        from mathorcup.q2_spatial_mechanism import execute_q2
        execute_q2(ctx)
        from mathorcup.q2_trigger_evaluation import execute_q2_trigger
        execute_q2_trigger(ctx)
        from mathorcup.task_evidence import q2_complete_evidence
        q2_complete_evidence(ctx)
        from mathorcup.question_completion import complete_q2
        complete_q2(ctx)
    if args.stage in ('all','q3'):
        from mathorcup.q3_task_aligned import execute_q3_aligned
        execute_q3_aligned(ctx)
        from mathorcup.q3_feature_bridge import execute_q3_feature_bridge
        execute_q3_feature_bridge(ctx)
        from mathorcup.question_completion import complete_q3
        complete_q3(ctx)
    if args.stage in ('all','q4'):
        from mathorcup.q4_flowpath import execute_q4_path
        from mathorcup.q4_solver import execute_q4_solver
        route=execute_q4_path(ctx,allow_provisional=(not args.strict_route or args.allow_provisional_q4))
        execute_q4_solver(ctx,route)
        from mathorcup.q4_grid_review import execute_q4_grid_review
        execute_q4_grid_review(ctx,route)
        from mathorcup.q4_cascade_reconstruction import execute_cascade_reconstruction
        execute_cascade_reconstruction(ctx)
        from mathorcup.q3_spatial_footprint import execute_q3_spatial
        execute_q3_spatial(ctx)
        if args.stage=='q4':
            from mathorcup.publication_assets import make_publication_assets
            make_publication_assets(ctx)
    if args.stage=='all':
        from mathorcup.q34_joint_arrival import execute_joint_arrival
        execute_joint_arrival(ctx)
        from mathorcup.publication_assets import make_publication_assets
        make_publication_assets(ctx)
        from mathorcup.report import make_report
        make_report(ctx)
        from mathorcup.redraw_catalog import build_redraw_catalog
        build_redraw_catalog(ctx)
        from mathorcup.provenance import audit_project
        audit_project(ctx)
    if args.stage=='tests':
        import subprocess
        raise SystemExit(subprocess.call([sys.executable,'-m','pytest','-q','tests'],cwd=ROOT))
    print('\n[DONE] outputs:',ctx.out)

if __name__=='__main__':main()
