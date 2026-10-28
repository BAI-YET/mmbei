"""Package current verified source, actual used inputs, outputs and answer PDF."""
from pathlib import Path
import zipfile,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
DELIVERY=ROOT.parent/'output';DELIVERY.mkdir(exist_ok=True)
required=[ROOT/'docs/MathorCup_B_v2_3_1_question_answers.pdf',ROOT/'outputs/quality/reproducibility_manifest.json',ROOT/'outputs/quality/test_result.json']
if not all(p.exists() for p in required):raise FileNotFoundError('Build current report, provenance and tests before packaging')
checks=json.loads(required[-1].read_text())
if checks.get('exit_code')!=0:raise ValueError('Refusing release with failed required tests')
# Freeze current source hashes first, then copy the same quality records into
# the redraw annex. Packaging must not include an earlier provenance snapshot.
sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import Context
from mathorcup.provenance import audit_project
from mathorcup.redraw_catalog import build_redraw_catalog
ctx=Context(ROOT);audit_project(ctx);build_redraw_catalog(ctx)
files=[]
for name in ['src','scripts','tests','configs','events','assets','docs','data/processed','data/external','outputs','data/raw/unpacked']:
    for p in (ROOT/name).rglob('*'):
        if not p.is_file() or '__pycache__' in p.parts or '.pytest_cache' in p.parts:continue
        if any(part.startswith('.') for part in p.relative_to(ROOT).parts):continue
        files.append(p)
files.extend([ROOT/'README.md',ROOT/'requirements.txt',ROOT/'run.py',ROOT/'.gitignore',ROOT/'.gitattributes'])
archive=DELIVERY/'MathorCup_B_v2_3_1_complete_release.zip'
staged_archive=archive.with_suffix('.zip.partial')
with zipfile.ZipFile(staged_archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
    for p in sorted(set(files)):
        z.write(p,'MathorCup_B_v2/'+str(p.relative_to(ROOT)))
with zipfile.ZipFile(staged_archive) as z:
    if z.testzip() is not None:raise ValueError('Archive byte verification failed')
    names=z.namelist()
    assert any('data/raw/unpacked/3_selected/' in n for n in names)
    assert any(n.endswith('Q1_all_candidate_ranking.csv') for n in names)
    assert any(n.endswith('Q2_sequential_model_contract.json') for n in names)
    assert any(n.endswith('Q3_source_initiation_time_distribution.csv') for n in names)
    assert any(n.endswith('ecmwf_source_pre_event_single_run.json') for n in names)
staged_archive.replace(archive)
print(json.dumps({'archive':str(archive),'file_count':len(names),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()},ensure_ascii=False))

# A smaller download preserves all source and PNG results, with explicit scope.
light=DELIVERY/'MathorCup_B_v2_3_1_code_and_figures.zip'
staged_light=light.with_suffix('.zip.partial')
selected=[p for p in set(files) if p.relative_to(ROOT).parts[0] in ('src','scripts','tests','configs','events')]
selected += [ROOT/'README.md',ROOT/'requirements.txt',ROOT/'run.py',ROOT/'.gitignore',ROOT/'.gitattributes',required[0]]
selected += [p for q in ('q1','q2','q3','q4','joint') for p in (ROOT/'outputs'/q).glob('*.png')]
figure_count=sum(p.suffix=='.png' for p in set(selected))
with zipfile.ZipFile(staged_light,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in sorted(set(selected)):z.write(p,'MathorCup_B_v2/'+str(p.relative_to(ROOT)))
    z.writestr('MathorCup_B_v2/PACKAGE_SCOPE.txt',
        f'本包包含完整代码、配置、测试、{figure_count}张PNG与本次PDF说明。\n'
        '原始数据、处理缓存、全部数值表、SVG/GIS及逐图重绘数据在完整复现包中。\n'
        '本版按明确假设给出四问模型答案；时间区间为发生条件下的预测区间，覆盖率以题面单例为分母。\n')
with zipfile.ZipFile(staged_light) as z:
    if z.testzip() is not None:raise ValueError('Code-and-figures archive verification failed')
staged_light.replace(light)
print(json.dumps({'code_and_figures':str(light),'bytes':light.stat().st_size,'sha256':hashlib.sha256(light.read_bytes()).hexdigest()},ensure_ascii=False))
