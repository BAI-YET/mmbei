"""Package current verified source, actual used inputs, outputs and answer PDF."""
from pathlib import Path
import zipfile,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
DELIVERY=ROOT.parent/'output';DELIVERY.mkdir(exist_ok=True)
required=[ROOT/'docs/MathorCup_B_v2_四问作答与原理分析.pdf',ROOT/'outputs/quality/reproducibility_manifest.json',ROOT/'outputs/quality/test_result.json']
if not all(p.exists() for p in required):raise FileNotFoundError('Build current report, provenance and tests before packaging')
checks=json.loads(required[-1].read_text())
if checks.get('exit_code')!=0:raise ValueError('Refusing release with failed required tests')
files=[]
for name in ['src','scripts','tests','configs','events','assets','docs','data/processed','data/external','outputs','data/raw/unpacked']:
    for p in (ROOT/name).rglob('*'):
        if not p.is_file() or '__pycache__' in p.parts or '.pytest_cache' in p.parts:continue
        files.append(p)
files.extend([ROOT/'README.md',ROOT/'requirements.txt',ROOT/'run.py',ROOT/'.gitignore',ROOT/'.gitattributes'])
archive=DELIVERY/'MathorCup_B_v2_2_complete_release.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
    for p in sorted(set(files)):
        z.write(p,'MathorCup_B_v2/'+str(p.relative_to(ROOT)))
with zipfile.ZipFile(archive) as z:
    if z.testzip() is not None:raise ValueError('Archive byte verification failed')
    names=z.namelist()
    assert any('data/raw/unpacked/3_selected/' in n for n in names)
    assert any(n.endswith('Q1_all_candidate_ranking.csv') for n in names)
    assert any(n.endswith('Q3_all_models_forward_backtest.csv') for n in names)
    assert any(n.endswith('ecmwf_source_pre_event_single_run.json') for n in names)
print(json.dumps({'archive':str(archive),'file_count':len(names),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()},ensure_ascii=False))
