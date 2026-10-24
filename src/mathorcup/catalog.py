"""Unpack 7 small ZIP datasets; leave the 257MB RAR optional, never duplicate it."""
import os,hashlib,zipfile
from pathlib import Path
from .common import write_json

def prepare_official_data(ctx,official_dir=None):
    candidates=[Path(official_dir)] if official_dir else [ctx.raw,ctx.root/'inputs',Path('/mnt/data')]
    unpack=ctx.raw/'unpacked';unpack.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for i in range(1,9):
        out=unpack/str(i)
        if out.exists() and any(out.rglob('*.*')):
            manifest.append({'dataset':i,'status':'existing_unpacked','path':str(out)})
            continue
        if i==3 and (unpack/'3_selected').exists():
            manifest.append({'dataset':3,'status':'four_selected_regional_periods','path':str(unpack/'3_selected'),'files':len(list((unpack/'3_selected').glob('*.tif')))})
            continue
        paths=[]
        for p in candidates:
            if p.exists():paths.extend(p.glob(f'（{i}）*.zip'))
        if not paths and i in (3,):
            manifest.append({'dataset':i,'status':'optional_rar_not_unpacked'})
            continue
        if not paths:
            raise FileNotFoundError(f'Dataset {i} ZIP not found: place official archives in {ctx.raw} or supply --official-dir')
        zpath=paths[0];out.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(zpath) as z:
            for file in z.infolist():
                target=(out/file.filename).resolve()
                if out.resolve() not in target.parents and target !=out.resolve():
                    raise ValueError('Zip path traversal blocked')
            z.extractall(out)
            manifest.append({'dataset':i,'path':str(out),'zip':zpath.name,
                             'files':len(z.infolist()),'bytes_unpacked':sum(x.file_size for x in z.infolist())})
    write_json(ctx.folder('quality')/'data_manifest.json',manifest)
    return unpack
