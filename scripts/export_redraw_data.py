"""Replay plots into an isolated directory and retain every plotting array.
Does not overwrite the user's validated numerical outputs or original images.
"""
from pathlib import Path
import sys,os,json,hashlib,shutil
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context,configure_plotting
from mathorcup.q1_event_signatures import execute_q1
from mathorcup.q2_spatial_mechanism import execute_q2
from mathorcup.q2_trigger_evaluation import execute_q2_trigger
from mathorcup.task_evidence import q2_complete_evidence
from mathorcup.q3_task_aligned import execute_q3_aligned
from mathorcup.q3_feature_bridge import execute_q3_feature_bridge
from mathorcup.q4_flowpath import execute_q4_path
from mathorcup.q4_solver import execute_q4_solver
from mathorcup.q4_grid_review import execute_q4_grid_review
from mathorcup.q3_spatial_footprint import execute_q3_spatial
from mathorcup.q34_joint_arrival import execute_joint_arrival
from mathorcup.publication_assets import make_publication_assets


def main():
    original=ROOT/'outputs';base=original/'redraw_data'
    before={str(p.relative_to(original)):hashlib.sha256(p.read_bytes()).hexdigest() for p in original.rglob('*') if p.is_file() and 'redraw_data' not in p.parts}
    base.mkdir(exist_ok=True);os.environ['MATHORCUP_REDRAW_ROOT']=str(base)
    ctx=Context(ROOT);ctx.out=ROOT/'tmp/redraw_replay_outputs';ctx.out.mkdir(parents=True,exist_ok=True);configure_plotting()
    print('Q1: export complete plotting arrays',flush=True);execute_q1(ctx)
    print('Q2: export complete plotting arrays',flush=True);execute_q2(ctx);execute_q2_trigger(ctx);q2_complete_evidence(ctx)
    print('Q3: export complete plotting arrays',flush=True);execute_q3_aligned(ctx);execute_q3_feature_bridge(ctx)
    print('Q4: export complete plotting arrays',flush=True);route=execute_q4_path(ctx,allow_provisional=True);execute_q4_solver(ctx,route);execute_q4_grid_review(ctx,route);execute_q3_spatial(ctx)
    print('Joint and publication figures',flush=True);execute_joint_arrival(ctx);make_publication_assets(ctx)
    after={name:hashlib.sha256((original/name).read_bytes()).hexdigest() for name in before}
    assert before==after,'Original published result changed during isolated replay'
    expected={str(p.relative_to(original)) for p in original.rglob('*.png') if 'redraw_data' not in p.parts}
    captured={json.loads(p.read_text())['figure'] for p in base.rglob('figure_metadata.json')}
    assert captured==expected,{'missing':sorted(expected-captured),'unexpected':sorted(captured-expected)}
    # Check numerical identity for all original CSVs recreated by the isolated run.
    import pandas as pd
    comparisons=[]
    for path in original.rglob('*.csv'):
        if 'redraw_data' in path.parts:continue
        relative=path.relative_to(original);replayed=ctx.out/relative
        if not replayed.exists():continue
        try:
            a=pd.read_csv(path);b=pd.read_csv(replayed)
            pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=1e-12,atol=1e-12)
            comparisons.append({'table':str(relative),'same':True})
        except pd.errors.EmptyDataError:
            assert path.read_bytes()==replayed.read_bytes();comparisons.append({'table':str(relative),'same':True})
    (base/'replay_validation.json').write_text(json.dumps({'original_files_unchanged':True,'original_files_checked':len(before),
       'all_figures_captured':len(captured),'expected_figures':len(expected),'numerical_table_comparisons':comparisons},ensure_ascii=False,indent=2))
    from mathorcup.redraw_catalog import build_redraw_catalog
    published=Context(ROOT);build_redraw_catalog(published)
    print(f'Complete: {len(captured)} figures. Original results unchanged. {len(comparisons)} CSV comparisons.',flush=True)

if __name__=='__main__':main()
