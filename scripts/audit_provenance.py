"""Audit the provenance of a complete local experiment."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context
from mathorcup.provenance import audit_project
if __name__=='__main__':
    result=audit_project(Context(ROOT))
    print({'figures':result['figure_count'],'input_and_source_files':len(result['files']),
       'forecast_sha_verified':result['archived_forecast_raw_sha256_verified'],
       'all_plot_inputs_recorded':result['all_figures_have_numeric_plot_inputs']})
