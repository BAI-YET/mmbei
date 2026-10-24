"""Authoritative ingestion: QC=0 when QC flag exists; never fill event-week missing hours."""
from __future__ import annotations
import pandas as pd,numpy as np
from pathlib import Path
from .common import write_json,record_warning

# Complete inventory of the 63 official non-time, non-QC variables.
# All raw physical/diagnostic fields are ingested for audit. Only independently
# interpretable, sufficiently complete measurements enter Q1 candidate ranking.
_RAW_FIELDS={
 'MET':{
  **{f'Ta_{h}m(degree C)':f'temp_{h}m' for h in [1,2,4,10,20]},
  **{f'WS_{h}m(m/s)':f'wind_{h}m' for h in [1,2,4,10,20]},
  **{f'WD_{h}m(degree)':f'wd_{h}m' for h in [1,2,4,10,20]},
  **{f'RH_{h}m(%)':f'rh_{h}m' for h in [1,2,4,10,20]},
  'Pressure(hPa)':'pressure','Precip(mm)':'precip_mm'},
 'RADM':{'Rsd(W/m2)':'sw_down','Rsu(W/m2)':'sw_up',
         'Rld(W/m2)':'lw_down','Rlu(W/m2)':'lw_up'},
 'SOIL':{
  **{f'soilT_{h}m(degree C)':f'soilT_{h.replace(".","p")}' for h in ['0.1','0.2','0.4','0.8','1.6']},
  **{f'VWC_{h}m(m3/m3)':f'vwc_{h.replace(".","p")}' for h in ['0.1','0.2','0.4','0.8','1.6']},
  'G_0.1m(W/m2)':'soil_heat','G_0.2m(W/m2)':'soil_heat_0p2'},
 'FLUX':{
  'Ux(m/s)':'flux_ux','Uy(m/s)':'flux_uy','Uz(m/s)':'flux_uz',
  'T_SONIC(deg C)':'sonic_temp','amb_e(kPa)':'ambient_vapor_pressure',
  'amb_e_sat(kPa)':'ambient_saturation_pressure','e(kPa)':'vapor_pressure',
  'e_sat(kPa)':'saturation_pressure','e_probe(kPa)':'probe_vapor_pressure',
  'e_sat_probe(kPa)':'probe_saturation_pressure','H2O_probe(g/m3)':'water_probe',
  'PA(kPa)':'flux_pressure','VPD(hPa)':'vpd','H2O_density(g/m3)':'water_density',
  'FC_mass(mg/m2*s)':'fc_mass','LE(W/m2)':'latent_heat','H(W/m2)':'sensible_heat',
  'TAU(kg/m*s2)':'tau','USTAR(m/s)':'ustar','TSTAR(deg C)':'tstar',
  'TKE(m2/s2)':'tke','FETCH_MAX(m)':'fetch_max','FETCH_90(m)':'fetch90',
  'FETCH_55(m)':'fetch55','FETCH_40(m)':'fetch40'}}

# Instrument diagnostics are retained as data-QA evidence, not as causal factors.
_DIAGNOSTIC_FIELDS={'fc_mass','fetch_max','fetch90','fetch55','fetch40',
  'ambient_saturation_pressure','probe_saturation_pressure','probe_vapor_pressure',
  'water_probe','flux_pressure'}


def _qc_column(raw):
    # Corresponding flags are case-insensitive and include field's first token before '(';
    # e.g. Ta_2m(degree C) -> QC_Ta_2m, Rsu(W/m2) -> QC_rsu.
    prefix=raw.split('(')[0]
    return 'QC_'+prefix

def _clean(df,typ):
    df.columns=[c.lstrip('\ufeff').strip() for c in df.columns]
    df['timestamp_raw']=pd.to_datetime(df['Timestamp'],errors='coerce')
    out=pd.DataFrame({'timestamp_raw':df['timestamp_raw']})
    stats=[]
    for raw,clean in _RAW_FIELDS[typ].items():
        if raw not in df.columns:
            record_warning(f'{typ} missing requested field {raw}');continue
        val=pd.to_numeric(df[raw],errors='coerce')
        val=val.mask(val.abs()>=9999)  # official 9999.9 sentinel; values this large not physical for fields selected
        qc_lookup={c.casefold():c for c in df.columns if c.casefold().startswith('qc_')}
        qc=qc_lookup.get(_qc_column(raw).casefold())
        if qc:
            flag=pd.to_numeric(df[qc],errors='coerce')
            val=val.mask(flag.ne(0))
        # Family-specific limits. Apply only unambiguous physical plausibility,
        # do not trim high events based on the reference disaster outcome.
        if clean.startswith('rh_'):val=val.where(val.between(0,100))
        if clean.startswith('wd_'):val=val.where(val.between(0,360))
        if clean.startswith('fetch_'):val=val.where(val.ge(0))
        if clean.startswith('wd_'):val=val.where(val.between(0,360))
        if clean.startswith('vwc_'):val=val.where(val.between(0,1.0))
        if clean.startswith(('temp_','soilT_')) or clean=='sonic_temp':val=val.where(val.between(-65,55))
        if clean.startswith('wind_'):val=val.where(val.between(0,90))
        if clean in ('flux_ux','flux_uy','flux_uz'):val=val.where(val.between(-90,90))
        if clean in ('pressure',):val=val.where(val.between(200,1100))
        if clean=='flux_pressure':val=val.where(val.between(20,110))
        if clean=='precip_mm':val=val.where(val.ge(0))
        if clean in ('sw_down','sw_up'):val=val.where(val.between(-15,1700))
        if clean in ('lw_down','lw_up'):val=val.where(val.between(60,800))
        if clean in ('soil_heat','soil_heat_0p2','sensible_heat','latent_heat'):
            val=val.where(val.between(-1500,2000))
        if clean in ('tau','tke','ustar','water_density','water_probe') or clean.startswith('fetch'):
            val=val.where(val.ge(0))
        if clean in ('tau',):val=val.where(val.le(60))
        if clean in ('ustar',):val=val.where(val.le(15))
        if clean in ('tke',):val=val.where(val.le(300))
        if clean in ('vpd',):val=val.where(val.between(-1,100))
        if clean in ('vapor_pressure','saturation_pressure','ambient_vapor_pressure',
                     'ambient_saturation_pressure','probe_vapor_pressure','probe_saturation_pressure'):
            val=val.where(val.between(0,15))
        if clean.startswith('fetch'):val=val.where(val.le(50000))
        out[clean]=val.values.astype(float)
        stats.append({'source':typ,'name':clean,'qc_field':qc,'valid_n':int(val.notna().sum()),
                      'missing_fraction':float(val.isna().mean()),'qc_available':bool(qc)})
    return out.dropna(subset=['timestamp_raw']),stats

def prepare_timeseries(ctx):
    out=ctx.processed/'obs_hourly.csv.gz'
    if out.exists() and (ctx.processed/'temp_county_daily.csv.gz').exists() and (ctx.processed/'obs_schema_v17.marker').exists():
        # Avoid silently reusing pre-v1.7 caches that lost 36 official raw fields.
        header=pd.read_csv(out,nrows=0).columns
        required={alias for mapping in _RAW_FIELDS.values() for alias in mapping.values()}
        inventory=ctx.folder('quality')/'complete_field_inventory_63.csv'
        if required.issubset(set(header)) and inventory.exists() and len(pd.read_csv(inventory))==63:
            return pd.read_csv(out,parse_dates=['timestamp_raw'])
        record_warning('Stale hourly cache detected: rebuilding 63-field official observation table.')
        out.unlink(missing_ok=True)
    root=ctx.raw/'unpacked'/'1'
    tables=[];audit=[]
    for typ in ['MET','RADM','SOIL','FLUX']:
        frames=[]
        for yr in range(2021,2027):
            fn=root/f'{typ}_Gyirong_{yr}.csv'
            if not fn.exists():raise FileNotFoundError(fn)
            orig=pd.read_csv(fn,low_memory=False,encoding='utf-8-sig')
            d,qa=_clean(orig,typ);qa=[dict(x,year=yr) for x in qa];audit.extend(qa)
            if d['timestamp_raw'].duplicated().any():
                record_warning(f'Raw {typ}_{yr} duplicate timestamp; retaining last')
                d=d.drop_duplicates('timestamp_raw',keep='last')
            frames.append(d)
        t=pd.concat(frames,ignore_index=True).set_index('timestamp_raw').sort_index()
        tables.append(t)
    combined=pd.concat(tables,axis=1).loc[:,lambda x:~x.columns.duplicated()]
    # Include raw rows without a valid timestamp? no. No numerical interpolation, no future fill.
    combined.to_csv(out,compression='gzip',index_label='timestamp_raw',float_format='%.8g')
    q=pd.DataFrame(audit)
    q.to_csv(ctx.folder('quality')/'qc_fields_by_year.csv',index=False)
    inventory=[]
    for group,mapping in _RAW_FIELDS.items():
        for raw,alias in mapping.items():
            sub=q[(q.source==group)&(q.name==alias)]
            inventory.append({'dataset':group,'raw_column':raw,'canonical_column':alias,
                'role':'diagnostic_only' if alias in _DIAGNOSTIC_FIELDS else 'physical_candidate',
                'valid_observations':int(sub.valid_n.sum()) if len(sub) else 0,
                'valid_years':int((sub.valid_n>0).sum()) if len(sub) else 0,
                'qc_variable_exists':bool(sub.qc_field.notna().any()) if len(sub) else False})
    pd.DataFrame(inventory).to_csv(ctx.folder('quality')/'complete_field_inventory_63.csv',index=False,encoding='utf-8-sig')
    # Annual summaries and raw clock validity / 2026 last valid observation
    q['source_variable']=q['source']+'/'+q['name']
    write_json(ctx.folder('quality')/'observation_metadata.json',{
        'clock':'Asia/Shanghai UTC+08:00 official site documentation', 'no_interpolation':True,
        'rows':int(len(combined)),'columns':int(combined.shape[1]),
        'first_valid':str(combined.dropna(how='all').index.min()),
        'last_valid':str(combined.dropna(how='all').index.max()),
        'snapshot_cut':ctx.cfg['snapshot_raw']})
    _prepare_xlsx(ctx)
    (ctx.processed/'obs_schema_v17.marker').write_text('complete physical fields v1.7',encoding='utf-8')
    return combined.reset_index()

def _prepare_xlsx(ctx):
    base=ctx.raw/'unpacked'/'2'
    dfs=[]
    for fn in sorted(base.rglob('*.xlsx')):
        try:d=pd.read_excel(fn,sheet_name=0)
        except Exception as e:
            record_warning(f'Cannot parse excel {fn.name}: {e}');continue
        if '时间' not in d.columns:continue
        d['date']=pd.to_datetime(d['时间'],errors='coerce')
        cols=[c for c in ['吉隆县','聂拉木县','萨嘎县','定日县','昂仁县'] if c in d]
        dfs.append(d[['date']+cols])
    if dfs:
        x=pd.concat(dfs,ignore_index=True).sort_values('date').drop_duplicates('date')
        x.to_csv(ctx.processed/'temp_county_daily.csv.gz',index=False,compression='gzip')
