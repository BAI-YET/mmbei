"""Readable paper figures and a complete plot-to-code index, generated from outputs."""
from pathlib import Path
import csv,json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from .common import savefig,FIGURE_TITLES,COLORS,write_json

PURPOSE={
 'Q1_11':'完整候选排序；对应Q1全部因素比较与原始前五项',
 'Q1_12':'完整多时间尺度异常；对应Q1窗口依赖',
 'Q1_13':'全部候选冗余结构；区分原始排序与机制代表项',
 'Q2_15':'全部候选模型的同一空间留出比较；非触发事件AUC',
 'Q2_16':'全部五项判别变量分布与解释',
 'Q2_17':'逐年向前阈值与不同更新间隔的监测误差',
 'Q3_08':'完整2025检验期，三组输入、区间与覆盖过程',
 'Q3_09':'预警进入、维持、续期、解除和有效期',
 'Q3_10':'实际起报天气预报进入Q2条件后的未来驱动敏感性',
 'Q4_11':'六个关键传播阶段；共同色标，来自本次二维数组',
 'Q4_12':'次级预设方量情景的真实1D通量、速度及累计量'
}
SOURCES={'q1':('q1_event_signatures.py;task_evidence.py','Q1_all_candidate_ranking.csv;Q1_all_factors_effect_and_stability.csv'),
 'q2':('q2_spatial_mechanism.py;q2_trigger_evaluation.py;q2_external.py;task_evidence.py','Q2_primary_trigger_assessment.json;Q2_all_models_comparison.csv'),
 'q3':('q3_task_aligned.py;warning_validation.py;q3_spatial_footprint.py;forecast_sensitivity.py','Q3_primary_competition_answer.json;Q3_all_models_forward_backtest.csv;Q3_warning_state_hourly.csv'),
 'q4':('q4_solver.py;q4_visuals.py;q4_grid_review.py;publication_assets.py','Q4_1d_settings_and_summary.json;Q4_2d_status.json;Q4_parameter_sensitivity.csv'),
 'joint':('q34_joint_arrival.py','Q34_status.json')}


def make_publication_assets(ctx):
    one_dimensional_grid_review(ctx)
    out=ctx.folder('q4');p=out/'Q4_2d_snapshots_for_replot.npz'
    if p.exists():
        a=np.load(p);frames=a['depth_m'];z=a['z_m'];tt=a['time_s'];maxh=a['max_depth_m']
        coeff=a['affine_coeff'];dx=abs(coeff[0]);dy=abs(coeff[4]);ex=[coeff[2]/1000,(coeff[2]+z.shape[1]*dx)/1000,(coeff[5]-z.shape[0]*dy)/1000,coeff[5]/1000]
        ids=np.round(np.linspace(0,len(frames)-1,6)).astype(int);vmax=max(.4,float(np.nanquantile(maxh[maxh>.2],.99)))
        fig,axs=plt.subplots(2,3,figsize=(12,7.4),layout='constrained');cm=plt.get_cmap('inferno').copy();cm.set_bad((1,1,1,0))
        for i,(ax,k) in enumerate(zip(axs.flat,ids)):
            ax.imshow(z,origin='upper',cmap='Greys',extent=ex,alpha=.7)
            im=ax.imshow(np.ma.masked_less_equal(frames[k],.2),origin='upper',extent=ex,cmap=cm,norm=LogNorm(.2,vmax))
            ax.set(title=f'{chr(97+i)}  {tt[k]/60:.2f} min',xlabel='东向坐标 (km)',ylabel='北向坐标 (km)');ax.tick_params(labelsize=8)
        fig.colorbar(im,ax=axs.ravel().tolist(),shrink=.75,label='等效流深 (m，共同对数色标)')
        info=json.loads((out/'Q4_2d_status.json').read_text());volume=info.get('coupled_1d_volume_m3')
        fig.suptitle(f'二维传播关键阶段：{volume/1e6:.0f}百万m³参数情景' if volume else '二维数值试验关键阶段',fontsize=13)
        savefig(fig,out/'Q4_11_six_key_frames.png')
    secondary=out/'Q4_secondary_1d_discharge.csv'
    if secondary.exists():
        h=pd.read_csv(secondary);fig,axs=plt.subplots(2,1,figsize=(10,5.8),sharex=True,layout='constrained')
        axs[0].plot(h.time_s/60,h.sensor_discharge_m3_s,color=COLORS[0]);axs[0].set(ylabel='实际1D耦合断面通量 (m³/s)')
        axs[1].plot(h.time_s/60,h.sensor_velocity_m_s,color=COLORS[1]);axs[1].set(ylabel='断面湿润单元速度 (m/s)',xlabel='源区释放后时间 (min)')
        fig.suptitle('预设上界方量情景的一维-二维实际传递')
        savefig(fig,out/'Q4_12_secondary_inflow_transfer.png')
    rows=[]
    for path in sorted(ctx.out.rglob('*.png')):
        if 'redraw_data' in path.parts:continue
        q=path.parent.name
        if q not in SOURCES:continue
        stem=path.stem;purpose=next((v for k,v in PURPOSE.items() if stem.startswith(k)),FIGURE_TITLES.get(stem,stem))
        code,data=SOURCES[q]
        if stem.startswith(('Q1_14','Q2_18','Q3_11')):
            code='question_completion.py';data={'q1':'Q1_direct_answer.json;Q1_explicit_extreme_weather.csv','q2':'Q2_expanded_mechanism_answer.json;Q2_expanded_channel_thresholds.csv','q3':'Q3_warning_lead_answer.json;Q3_proxy_backtest_2025.csv'}[q]
        elif stem.startswith('Q3_12'):
            code='q4_reconstructed_spatial.py';data='Q3_reconstructed_spatial_answer.json;Q3_reconstructed_warning_range.gpkg'
        elif q=='q4' and stem.startswith(tuple(f'Q4_{i}' for i in range(14,22))):
            code='q4_terrain_reconstruction.py;q4_cascade_reconstruction.py;q4_reconstructed_spatial.py'
            data='Q4_cascade_reconstruction.json;Q4_reconstruction_geometry.json;Q4_reconstruction_2d_summary.json'
        rows.append({'question':q.upper(),'figure':str(path.relative_to(ctx.out)),'purpose':purpose,'generating_code':'src/mathorcup/'+code,
            'numerical_evidence':data,'figure_type':'历史观测/物理条件' if q in ('q1','q2') else '气象回测/条件情景' if q=='q3' else '参数情景/数值验证'})
    with (ctx.out/'FIGURE_INDEX.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return rows


def one_dimensional_grid_review(ctx):
    """Isolate 1D numerical spacing on the SAME fixed candidate topography."""
    import hashlib
    from .q4_solver import simulate_1d
    out=ctx.folder('q4');p=pd.read_csv(out/'Q4_route_profile.csv');cfg=ctx.cfg['q4']
    fingerprint=hashlib.sha256((Path(__file__).with_name('q4_solver.py').read_text()+json.dumps(cfg,sort_keys=True)+(out/'Q4_route_profile.csv').read_text()).encode()).hexdigest()
    meta=out/'Q4_1d_grid_review_metadata.json';target=out/'Q4_1d_grid_resolution_reliability.csv'
    if meta.exists() and target.exists() and json.loads(meta.read_text()).get('fingerprint')==fingerprint:return pd.read_csv(target)
    records=[]
    for dx in [30.,60.,120.]:
        s=np.arange(0,float(p.distance_m.max()),dx)
        profile=pd.DataFrame({'distance_m':s,'elevation_m':np.interp(s,p.distance_m,p.elevation_m)})
        a=simulate_1d(profile,mu=cfg['mu'],xi=cfg['xi_m_s2'],width=cfg['width_m'],volume=cfg['release_volume_m3'],duration=cfg['release_duration_s'],
           total_seconds=cfg['total_seconds'],cfl=cfg['cfl'],frames=3,initial_velocity_m_s=cfg['initial_velocity_m_s'])
        m=a['mass'];records.append({'resolution_m':dx,'runout_km':m['max_runout_m']/1000,
           'arrival_min':m['first_arrival_to_downstream_end_s']/60 if m['first_arrival_to_downstream_end_s'] is not None else np.nan,
           'max_depth_m':m['max_depth_over_time_m'],'max_wet_velocity_m_s':m['max_speed_over_time_m_s'],
           'mass_relative_error':m['relative_absolute_mass_residual'],'topography_information':'same interpolated path; no new DEM information'})
    d=pd.DataFrame(records);d.to_csv(target,index=False)
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.3),layout='constrained')
    axs[0].plot(d.resolution_m,d.runout_km,'-o',color=COLORS[0]);axs[0].set(xlabel='一维步长 (m)',ylabel='最大传播距离 (km)')
    axs[1].plot(d.resolution_m,d.max_depth_m,'-s',color=COLORS[1]);axs[1].set(xlabel='一维步长 (m)',ylabel='最大等效流深 (m)')
    fig.suptitle('固定候选地形下的一维网格敏感性')
    savefig(fig,out/'Q4_13_1d_grid_resolution.png')
    write_json(meta,{'fingerprint':fingerprint,'interpretation':'spacing sensitivity only; not proof of asymptotic convergence','records':records})
    return d
