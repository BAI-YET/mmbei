"""Refresh publication outputs with the current Q2/Q3 contract and provenance."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context,configure_plotting
from mathorcup.publication_assets import make_publication_assets
from mathorcup.q34_joint_arrival import execute_joint_arrival
from mathorcup.q4_reconstructed_spatial import replot_warning_range_from_saved
from mathorcup.report import make_report
from mathorcup.redraw_catalog import build_redraw_catalog
from mathorcup.review_audit import build_review_verdict
from mathorcup.provenance import audit_project

if __name__=='__main__':
    ctx=Context(ROOT);configure_plotting()
    replot_warning_range_from_saved(ctx)
    make_publication_assets(ctx);execute_joint_arrival(ctx);make_report(ctx)
    build_review_verdict(ctx);build_redraw_catalog(ctx);result=audit_project(ctx)
    print('Current figures, model status and source/input hashes refreshed:',result['figure_count'])
