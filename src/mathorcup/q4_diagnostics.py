"""Event-clock, whole-route and resolution figures from accepted solver results."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .common import savefig,write_json,COLORS


def plot_review_diagnostics(ctx,upper,lower,report):
    out=ctx.folder('q4')
    residual=pd.read_csv(out/'Q4_calibration_residuals.csv')
    labels=['源区至汇入口','汇入口至口岸','源区至口岸','口岸峰值流量','峰值时断面速度']
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.7),layout='constrained')
    x=np.arange(3);d=residual.iloc[:3]
    axs[0].barh(x-.2,d.published_reconstruction_reference,height=.36,color=COLORS[0],label='文献重建参照')
    axs[0].barh(x+.2,d.computed_value,height=.36,color=COLORS[2],label='本次计算前锋')
    axs[0].set(yticks=x,yticklabels=labels[:3],xlabel='经历时间 (s)',title='a  统一源区释放起点的时间对照')
    axs[0].legend(frameon=False,fontsize=8)
    bars=axs[1].barh(labels,100*residual.relative_error,color=COLORS[1])
    axs[1].axvline(0,color='#666',lw=1);axs[1].bar_label(bars,fmt='%+.1f%%',padding=3,fontsize=8)
    axs[1].set(xlabel='相对文献重建参照偏差 (%)',title='b  全部标定量的实际偏差')
    fig.suptitle('事件时间起点与标定残差')
    savefig(fig,out/'Q4_22_time_alignment_and_residuals.png')
    plot_whole_route(ctx,upper,lower)
    plot_release_projection(ctx)
    grid=pd.read_csv(out/'Q4_reconstruction_1d_grid_experiment.csv').sort_values('dx_m')
    baseline_file=ctx.folder('quality')/'review_baseline.json'
    if not baseline_file.exists():baseline_file=ctx.root/'outputs/quality/review_baseline.json'
    baseline=None
    if baseline_file.exists():baseline=json.loads(baseline_file.read_text())
    elif (ctx.root/'configs/review_baseline.yaml').exists():
        import yaml
        baseline=yaml.safe_load((ctx.root/'configs/review_baseline.yaml').read_text())
    if baseline is not None:
        from io import StringIO
        old=pd.read_csv(StringIO(baseline['Q4_grid_csv'])).sort_values('dx_m')
        fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
        for table,name,color in [(old,'旧版一阶计算',COLORS[1]),(grid,'本次下游高阶计算',COLORS[0])]:
            axs[0].plot(table.dx_m,table.peak_discharge_m3_s/1e4,'-o',label=name,color=color)
            axs[1].plot(table.dx_m,table.arrival_s,'-o',label=name,color=color)
        axs[0].set(xlabel='一维网格 (m)',ylabel='口岸峰值流量 (万m³/s)',title='a  各自固定标定参数下的网格响应')
        axs[1].set(xlabel='一维网格 (m)',ylabel='源区释放后到达时间 (s)',title='b  各自固定标定参数下的到达时间')
        for ax in axs:ax.legend(frameon=False,fontsize=8)
        fig.suptitle('修订前后网格敏感性的完整比较')
        savefig(fig,out/'Q4_24_resolution_response_comparison.png')


def plot_whole_route(ctx,upper,lower):
    """Plan view of calculated 1D sections; this is not a whole-valley 2D inundation map."""
    out=ctx.folder('q4');route=pd.read_csv(out/'Q4_reconstruction_route.csv')
    time_end=min(float(lower['history'].time_s.max()),
                 max(480.,1.25*float(lower['history'].time_s.iloc[lower['history'].outlet_discharge_m3_s.idxmax()])))
    times=np.linspace(0,time_end,6);arrays=[]
    for answer in [upper,lower]:
        clock=np.array([c['time_s'] for c in answer['captures']])
        depth=np.array([c['h'] for c in answer['captures']])
        arrays.append(np.array([np.interp(times,clock,depth[:,k]) for k in range(depth.shape[1])]).T)
    values=np.concatenate([arrays[0],arrays[1][:,1:]],axis=1)
    if values.shape[1]!=len(route):raise ValueError('Whole-route plotting coordinates do not match the numerical meshes')
    np.savez_compressed(out/'Q4_whole_route_visualization.npz',time_s=times,depth_m=values,
        utm_x=route.utm_x.to_numpy(),utm_y=route.utm_y.to_numpy())
    write_json(out/'Q4_whole_route_visualization.json',{'source':'accepted 1D solver snapshots on the entire route',
        'time_sampling':'linear interpolation between stored computed snapshots',
        'spatial_scope':'source to port; colour is channel-section depth projected to map',
        'inundation_extent_source':'Q4_reconstruction_2d_arrays.npz is the separate local 2D calculation'})
    fig,axs=plt.subplots(2,3,figsize=(12.6,7.3),sharex=True,sharey=True,layout='constrained')
    x=route.utm_x/1000;y=route.utm_y/1000
    vmax=max(1.,float(np.quantile(values[values>.2],.99)))
    for ax,t,h in zip(axs.flat,times,values):
        ax.plot(x,y,color='#C2C6CA',lw=1,zorder=1)
        wet=h>.2
        im=ax.scatter(x[wet],y[wet],c=h[wet],s=14,cmap='turbo',vmin=0,vmax=vmax,zorder=2)
        ax.scatter(x.iloc[[0,-1]],y.iloc[[0,-1]],s=24,c=[COLORS[2],COLORS[1]],edgecolors='white',zorder=3)
        ax.set(title=f'源区释放后 {t:.0f} s',xlabel='UTM东向 (km)',ylabel='UTM北向 (km)')
        ax.set_aspect('equal',adjustable='box')
    fig.colorbar(im,ax=axs.ravel().tolist(),shrink=.8,label='沿程断面流深 (m)',extend='max')
    fig.suptitle('源区至口岸的全程传播过程投影')
    savefig(fig,out/'Q4_23_whole_route_six_stages.png')


def plot_release_projection(ctx):
    out=ctx.folder('q4');p=pd.read_csv(out/'Q4_release_pixel_projection.csv')
    route=pd.read_csv(out/'Q4_reconstruction_route.csv')
    geometry=json.loads((out/'Q4_reconstruction_geometry.json').read_text())
    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    image=axs[0].scatter(p.source_utm_x/1000,p.source_utm_y/1000,c=p.source_elevation_m,
        s=12,cmap='terrain',label='官方源区网格')
    local=route.loc[route.distance_m<=p.route_distance_m.max()+100]
    axs[0].plot(local.utm_x/1000,local.utm_y/1000,color=COLORS[0],lw=1.3,label='实际地形路径')
    axs[0].set(xlabel='UTM东向 (km)',ylabel='UTM北向 (km)',title='a  源区真实像元与路径的空间位置')
    axs[0].legend(frameon=False,fontsize=8);fig.colorbar(image,ax=axs[0],label='源区像元高程 (m)',shrink=.85)
    axs[1].scatter(p.source_elevation_m,p.route_elevation_m,s=11,alpha=.55,color=COLORS[0])
    lo=min(p.source_elevation_m.min(),p.route_elevation_m.min());hi=max(p.source_elevation_m.max(),p.route_elevation_m.max())
    axs[1].plot([lo,hi],[lo,hi],color='#777',ls='--',label='相同高程')
    axs[1].set(xlabel='官方源区像元高程 (m)',ylabel='投影路径高程 (m)',title='b  初始释放高程与势能一致性')
    axs[1].text(.04,.96,f"平均绝对高程差 {geometry['release_projection_mean_abs_elevation_error_m']:.1f} m\n释放物质平均床高 {geometry['assigned_release_mean_bed_elevation_m']:.1f} m",
        transform=axs[1].transAxes,va='top',fontsize=9)
    axs[1].legend(frameon=False,fontsize=8,loc='lower right')
    fig.suptitle('源区三维投影与初始释放势能核查')
    savefig(fig,out/'Q4_25_release_geometry_and_elevation.png')
