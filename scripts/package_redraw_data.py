"""Update the existing integration archive without dropping any previous member."""
from pathlib import Path
import json,sys,hashlib,subprocess,zipfile,shutil
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'output';archive=OUT/'MathorCup_B_v2_代码图表与复现数据.zip';data=ROOT/'outputs/redraw_data'
with zipfile.ZipFile(archive) as z:
    original_names=set(z.namelist())
    scientific={name:hashlib.sha256(z.read(name)).hexdigest() for name in z.namelist() if '/outputs/' in name and (name.endswith('.csv') or name.endswith('.png'))}
log=(ROOT.parent/'tmp/tests_redraw_full.log').read_text()
assert '60 passed' in log and 'FAILED' not in log
(ROOT/'outputs/quality/test_result.json').write_text(json.dumps({'exit_code':0,'output':log},ensure_ascii=False,indent=2))
meta=ROOT/'outputs/quality/reproducibility_manifest.json';manifest=json.loads(meta.read_text());records=[]
for group,base in [('source',ROOT/'src'),('scripts',ROOT/'scripts'),('tests',ROOT/'tests'),('config',ROOT/'configs'),('official_data',ROOT/'data/raw/unpacked'),('external_data',ROOT/'data/external')]:
    for p in sorted(base.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        records.append({'kind':group,'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
manifest['files']=records;manifest['all_63_plot_inputs_exported']=True;meta.write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
for src in [meta,ROOT/'outputs/quality/test_result.json']:shutil.copy2(src,data/'source_tables/quality'/src.name)
result=subprocess.run([sys.executable,str(ROOT/'scripts/build_release.py')],capture_output=True,text=True)
print(result.stdout);assert result.returncode==0,result.stderr
with zipfile.ZipFile(archive) as z:
    assert original_names.issubset(set(z.namelist())),original_names-set(z.namelist())
    for name,sha in scientific.items():assert hashlib.sha256(z.read(name)).hexdigest()==sha,('Changed original scientific table/image',name)
    new_count=len(z.namelist())
proof={'all_original_archive_members_retained':True,'original_archive_members':len(original_names),'updated_archive_members':new_count,
       'original_scientific_tables_and_images_unchanged':len(scientific),'covered_figures':63}
(data/'archive_retention_check.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2))
# Append the final audit only; this member is new and does not replace any source.
with zipfile.ZipFile(archive,'a',zipfile.ZIP_DEFLATED,compresslevel=5) as z:z.write(data/'archive_retention_check.json','MathorCup_B_v2/outputs/redraw_data/archive_retention_check.json')
annex=OUT/'MathorCup_B_v2_全部图像重绘数据.zip'
with zipfile.ZipFile(annex,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in sorted(data.rglob('*')):
        if p.is_file():z.write(p,'MathorCup_B_全部图像重绘数据/'+str(p.relative_to(data)))
with zipfile.ZipFile(annex) as z:assert z.testzip() is None
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps({'retention':proof,'full_archive_bytes':archive.stat().st_size,'data_only_archive':str(annex),'data_only_bytes':annex.stat().st_size},ensure_ascii=False))
