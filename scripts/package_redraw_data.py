"""Package the current verified redraw dataset, using its actual figure count."""
from pathlib import Path
import json
import zipfile
ROOT=Path(__file__).resolve().parents[1]

def main():
    data=ROOT/'outputs/redraw_data'
    summary=json.loads((data/'catalog_summary.json').read_text())
    checks=json.loads((data/'retention_checks.json').read_text())
    if checks['figure_count']!=summary['figures']:
        raise ValueError('Recheck the latest plotting inputs before packaging')
    target=ROOT.parent/'output/MathorCup_B_v2_1_全部图像重绘数据.zip'
    target.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as archive:
        for p in sorted(data.rglob('*')):
            if p.is_file():archive.write(p,'MathorCup_B_全部图像重绘数据/'+str(p.relative_to(data)))
    with zipfile.ZipFile(target) as archive:
        if archive.testzip() is not None:raise ValueError('Archive verification failed')
    print({'archive':str(target),'figures':summary['figures'],'bytes':target.stat().st_size})

if __name__=='__main__':main()
