"""Two-reach reconstruction: measured-DEM geometry and published event constraints.

No synthetic inlet hydrograph is generated. The downstream input is the recorded
finite-volume discharge of the upstream avalanche, with the published mobile fraction.
"""
from __future__ import annotations
import json
import hashlib
import itertools
import inspect
import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from .q4_solver import G, EPS, rusanov_1d, flux_1d, simulate_2d
from .q4_flowpath import local_dem
from .q4_terrain_reconstruction import prepare_reconstruction_geometry
from .common import COLORS, savefig, write_json


def section_depth(area,width,side):
    """Invert A=b*h+m*h², using a stable expression also valid for m=0."""
    return 2*np.asarray(area)/(width+np.sqrt(width*width+4*side*np.asarray(area)))


def section_pressure(depth,width,side):
    return G*(.5*width*depth*depth+side*depth**3/3)


def advance_channel(A,Q,bed,width,dx,dt,mu,xi,qin=0.,vin=0.,bank_side_slope=None,return_interface_flux=False):
    """Conservative variable-width hydrostatic reconstruction, including width pressure.

    Hydrostatic face corrections use each cell's own width. Therefore a lake at
    rest remains at rest even where both the bed and channel width change.
    """
    side=np.zeros_like(width) if bank_side_slope is None else bank_side_slope
    h=section_depth(A,width,side);u=Q/np.maximum(A,EPS);n=len(h)
    crest=np.maximum(bed[:-1],bed[1:])
    hl=np.maximum(h[:-1]+bed[:-1]-crest,0);hr=np.maximum(h[1:]+bed[1:]-crest,0)
    facewidth=np.minimum(width[:-1],width[1:])
    faceside=np.minimum(side[:-1],side[1:])
    al=facewidth*hl+faceside*hl**2;ar=facewidth*hr+faceside*hr**2
    ql=al*u[:-1];qr=ar*u[1:]
    pl=section_pressure(hl,facewidth,faceside);pr=section_pressure(hr,facewidth,faceside)
    fl=np.stack([ql,ql*u[:-1]+pl]);fr=np.stack([qr,qr*u[1:]+pr])
    wave=np.maximum(abs(u[:-1])+np.sqrt(G*al/np.maximum(facewidth+2*faceside*hl,EPS)),
                    abs(u[1:])+np.sqrt(G*ar/np.maximum(facewidth+2*faceside*hr,EPS)))
    face=.5*(fl+fr)-.5*wave[None,:]*np.stack([ar-al,qr-ql])
    outgoing=np.zeros((2,n));incoming=np.zeros((2,n))
    outgoing[:,:-1]=face;incoming[:,1:]=face
    pressure=section_pressure(h,width,side)
    outgoing[1,:-1]+=pressure[:-1]-pl
    incoming[1,1:]+=pressure[1:]-pr
    outgoing[:,-1]=[max(float(Q[-1]),0),max(float(Q[-1]),0)*u[-1]+pressure[-1]]
    incoming[1,0]=pressure[0]  # closed upstream face at source
    if qin>0:
        # Transfer the computed volume and advective momentum only. Prescribing
        # an additional stage from Q/u can inject unobserved pressure energy at
        # the interface, especially when the upstream discharge drops rapidly.
        incoming[:,0]=[qin,qin*vin+pressure[0]]
    updated=np.stack([A,Q])-dt/dx*(outgoing-incoming)
    if np.min(updated[0])<-1e-6:raise FloatingPointError('Negative channel area')
    clip=float(-np.minimum(updated[0],0).sum()*dx)
    area=np.maximum(updated[0],0);velocity=updated[1]/np.maximum(area,EPS)
    hn=section_depth(area,width,side);cos=np.cos(np.arctan(np.gradient(bed,dx)))
    # Backward-Euler quadratic drag is dissipative without reversing velocity.
    speed=np.maximum(np.abs(velocity)-dt*mu*G*cos,0)
    hydraulic_depth=area/np.maximum(width+2*side*hn,EPS)
    bdt=dt*G/(xi*np.maximum(hydraulic_depth,.01))
    speed=2*speed/(1+np.sqrt(1+4*bdt*speed))
    discharge=area*np.sign(velocity)*speed
    discharge[hn<1e-6]=0
    result=(area,discharge,float(outgoing[0,-1]*dt),clip)
    return result+(face[0].copy(),) if return_interface_flux else result


def simulate_channel(profile,mu,xi,total_seconds=1000,cfl=.3,initial_volume_m3=0.,
                     inlet=None,frames=12,sensor_upstream_m=2500):
    p=profile.reset_index(drop=True);s=p.distance_m.to_numpy(float);s=s-s[0]
    z=p.elevation_m.to_numpy(float);width=p.effective_width_m.to_numpy(float)
    side=p.bank_side_slope.to_numpy(float) if 'bank_side_slope' in p else np.zeros(len(p))
    dx=float(np.median(np.diff(s)))
    if np.max(abs(np.diff(s)-dx))>1e-7:raise ValueError('Channel grid must be uniform')
    w=p.source_area_weight.to_numpy(float) if 'source_area_weight' in p else np.zeros(len(p))
    if initial_volume_m3 and w.sum()<=0:raise ValueError('Initial volume needs a measured release footprint')
    A=initial_volume_m3*w/max(w.sum(),EPS)/dx;Q=np.zeros(len(p))
    sensor=min(len(p)-2,max(1,int(np.argmin(abs(s-(s[-1]-sensor_upstream_m))))))
    if inlet is not None:
        times=inlet.time_s.to_numpy(float);cumulative=inlet.cumulative_m3.to_numpy(float)
        speeds=inlet.velocity_m_s.to_numpy(float)
        if np.any(np.diff(times)<=0) or np.any(np.diff(cumulative)<-1e-6):raise ValueError('Invalid upstream cumulative hydrograph')
        rates=np.diff(cumulative)/np.diff(times)
    else:times=cumulative=speeds=rates=None
    t=0.;outvol=0.;inputvol=float(initial_volume_m3);clip=0.;first=None;first_sensor=None
    sensor_cumulative=0.;sensor_flux=0.;sensor_reverse=0.
    transfer_velocity=0.
    captures=[];frame_times=np.linspace(0,total_seconds,frames);next_frame=0;rows=[]
    max_depth=section_depth(A,width,side);max_speed=np.zeros(len(A));peak_port=0.;peak_v=0.;steps=0
    arrival_field=np.where(max_depth>.2,0.,np.nan)
    def record():
        h=section_depth(A,width,side);u=Q/np.maximum(A,EPS)
        rows.append({'time_s':t,'outflow_cumulative_m3':outvol,'input_cumulative_m3':inputvol,
          'stored_m3':float(A.sum()*dx),'outlet_discharge_m3_s':max(float(Q[-1]),0.),
          'outlet_velocity_m_s':max(float(u[-1]),0.),'outlet_depth_m':float(h[-1]),
          'outlet_transfer_velocity_m_s':transfer_velocity,
          'sensor_discharge_m3_s':max(sensor_flux,0.),'sensor_velocity_m_s':max(float(u[sensor]),0.),
          'sensor_cumulative_m3':sensor_cumulative,'sensor_reverse_cumulative_m3':sensor_reverse,
          'sensor_depth_m':float(h[sensor]),'max_depth_m':float(h.max())})
    record()
    while t<total_seconds-1e-8:
        h=section_depth(A,width,side);u=Q/np.maximum(A,EPS)
        wave=float(np.max(abs(u)+np.sqrt(G*A/np.maximum(width+2*side*h,EPS))))
        dt=min(1.,cfl*dx/max(wave,.5),total_seconds-t)
        qin=vin=0.
        if inlet is not None:
            k=min(len(rates)-1,max(0,int(np.searchsorted(times,t+1e-8,side='right')-1)))
            if t<times[-1]:
                if times[k+1]>t+1e-8:dt=min(dt,float(times[k+1]-t))
                qin=max(float(rates[k]),0.);vin=float(speeds[k+1])
                if qin>0:
                    ab=qin/max(vin,.01);hb=section_depth(ab,width[0],side[0])
                    dt=min(dt,cfl*dx/(max(vin,.01)+np.sqrt(G*ab/(width[0]+2*side[0]*hb))))
            before=float(np.interp(t,times,cumulative));after=float(np.interp(t+dt,times,cumulative))
            qin=max((after-before)/dt,0.)
        transfer_velocity=max(float(Q[-1]/max(A[-1],EPS)),0.)
        A,Q,dv,loss,flux=advance_channel(A,Q,z,width,dx,dt,mu,xi,qin,vin,bank_side_slope=side,return_interface_flux=True)
        sensor_flux=float(flux[sensor]);sensor_cumulative+=max(sensor_flux,0.)*dt
        sensor_reverse+=max(-sensor_flux,0.)*dt
        t+=dt;steps+=1;outvol+=dv;inputvol+=qin*dt;clip+=loss
        h=section_depth(A,width,side);u=Q/np.maximum(A,EPS);max_depth=np.maximum(max_depth,h)
        newly=(h>.2)&np.isnan(arrival_field);arrival_field[newly]=t
        max_speed=np.maximum(max_speed,np.where(h>.2,np.abs(u),0))
        if first is None and h[-1]>.2:first=t
        if first_sensor is None and h[sensor]>.2:first_sensor=t
        if Q[-1]>peak_port:peak_port=float(Q[-1]);peak_v=float(u[-1])
        while next_frame<len(frame_times) and t>=frame_times[next_frame]-1e-8:
            captures.append({'time_s':t,'h':h.copy(),'q':(Q/width).copy()});next_frame+=1
        record()
        if steps>150000:raise RuntimeError('Channel CFL step limit exceeded')
    stored=float(A.sum()*dx)
    mass={'injected_m3':inputvol,'outflow_m3':outvol,'stored_m3':stored,
      'relative_mass_error':abs(stored+outvol-inputvol)/max(inputvol,1),
      'negative_clip_m3':clip,'arrival_s':first,'sensor_arrival_s':first_sensor,
      'peak_discharge_m3_s':peak_port,'velocity_at_peak_discharge_m_s':peak_v,
      'maximum_depth_m':float(max_depth.max()),'maximum_wet_speed_m_s':float(max_speed.max()),
      'steps':steps,'dx_m':dx,'sensor_distance_m':float(s[sensor])}
    inlet_speed_max=float(np.max(speeds)) if speeds is not None else 0.
    bound=np.sqrt(inlet_speed_max**2+2*G*(float(np.max(z)-np.min(z))+2*mass['maximum_depth_m']))
    mass['generous_energy_speed_bound_m_s']=float(bound)
    if mass['maximum_wet_speed_m_s']>bound*1.1:
        raise FloatingPointError('Computed wet-flow speed exceeds the terrain/input energy bound')
    return {'profile':p,'s':s,'z':z,'width':width,'history':pd.DataFrame(rows),
            'captures':captures,'mass':mass,'max_depth':max_depth,'max_speed':max_speed,
            'sensor_idx':sensor,'A':A,'Q':Q,'arrival_field_s':arrival_field}


def upstream_transfer(answer,mobile_fraction):
    if not 0<=mobile_fraction<=1:raise ValueError('Mobile fraction cannot create upstream material')
    h=answer['history']
    v=h.outlet_transfer_velocity_m_s if 'outlet_transfer_velocity_m_s' in h else h.outlet_velocity_m_s
    return pd.DataFrame({'time_s':h.time_s,'cumulative_m3':h.outflow_cumulative_m3*mobile_fraction,
                         'velocity_m_s':v})


def fit_cascade(ctx,profile,junction,evidence):
    model_path=ctx.root/'configs/calibrated_q4.yaml'
    if model_path.exists() and not ctx.cfg.get('refit_q4',False):
        model=yaml.safe_load(model_path.read_text())
        if model.get('calibration_basis_sha256')==calibration_fingerprint(profile,evidence):
            return replay_calibrated_cascade(ctx,profile,junction,evidence,model)
    values={k:v['value'] for k,v in evidence['published_values'].items()};cfg=evidence['reconstruction']
    first_profile=profile.iloc[:junction+1].copy();second_profile=profile.iloc[junction:].copy()
    candidates=[];stage1=[];best1=None;loss1=np.inf
    for mu,xi in itertools.product(cfg['stage1_mu_candidates'],cfg['stage1_xi_candidates']):
        a=simulate_channel(first_profile,mu,xi,total_seconds=cfg['total_seconds'],cfl=cfg['cfl'],
                           initial_volume_m3=values['release_volume_m3'],frames=12)
        arrival=a['mass']['arrival_s']
        passed=a['mass']['outflow_m3'];mobile=values['mobile_debris_volume_m3']
        # The literature's 30+100 s describes the transition interval, rather
        # than an independently timed first 0.2 m wetting front. Use it as an
        # upper bound, and fit the observed retained/mobile volume partition.
        if arrival and passed>=mobile:
            error=float(np.log(passed/mobile)**2+max(0,np.log(arrival/values['confluence_elapsed_s']))**2)
        else:error=np.inf
        stage1.append({'stage':1,'mu':mu,'xi':xi,'arrival_s':arrival,'objective':error,
                       'passed_to_confluence_m3':a['mass']['outflow_m3']})
        if error<loss1:best1=(mu,xi,a);loss1=error
    if best1 is None:raise RuntimeError('No upstream parameter set reaches the confluence')
    mu1,xi1,upper=best1
    # Upstream stopping already accounts for part of the 67 Mm³ deposit. Only
    # the additional retained fraction is removed at the transfer; do not
    # apply 36/103 again and double-count the simulated upstream deposition.
    fraction=values['mobile_debris_volume_m3']/upper['mass']['outflow_m3']
    inlet=upstream_transfer(upper,fraction)
    target_time=values['confluence_to_port_elapsed_s'];target_q=values['port_peak_discharge_m3_s']
    target_v=values['port_cross_section_mean_velocity_m_s'];best2=None;loss2=np.inf
    for mu,xi in itertools.product(cfg['stage2_mu_candidates'],cfg['stage2_xi_candidates']):
        a=simulate_channel(second_profile,mu,xi,total_seconds=cfg['total_seconds'],cfl=cfg['cfl'],inlet=inlet,frames=12)
        m=a['mass'];transit=m['arrival_s']-upper['mass']['arrival_s'] if m['arrival_s'] else None
        if transit and transit>0 and m['peak_discharge_m3_s']>0:
            terms=[np.log(transit/target_time),np.log(m['peak_discharge_m3_s']/target_q),
                   np.log(m['velocity_at_peak_discharge_m_s']/target_v)]
            error=float(np.mean(np.square(terms)))
        else:error=float('inf')
        candidates.append({'stage':2,'mu':mu,'xi':xi,'arrival_s':m['arrival_s'],'transit_s':transit,
             'peak_discharge_m3_s':m['peak_discharge_m3_s'],'velocity_m_s':m['velocity_at_peak_discharge_m_s'],
             'objective':error,'mass_error':m['relative_mass_error']})
        if error<loss2:best2=(mu,xi,a);loss2=error
    if best2 is None:raise RuntimeError('No downstream calibration member reaches the port')
    # The effective capacity of a trapezoidal section is an estimated model
    # parameter. Q_peak/u_peak supplies an independent reconstructed section
    # area; fitting a capacity factor avoids calling DEM idealisation exact.
    initial=best2;best_scale=1.;evaluations=0
    def objective(x):
        nonlocal best2,loss2,best_scale,evaluations
        mu,xi,scale=np.exp(x)
        trial_profile=second_profile.copy()
        trial_profile['effective_width_m']*=scale
        trial_profile['bank_side_slope']*=scale
        a=simulate_channel(trial_profile,float(mu),float(xi),total_seconds=cfg['total_seconds'],
                           cfl=cfg['cfl'],inlet=inlet,frames=12)
        m=a['mass'];transit=m['arrival_s']-upper['mass']['arrival_s'] if m['arrival_s'] else None
        if transit and transit>0 and m['peak_discharge_m3_s']>0:
            error=float(np.mean(np.square([np.log(transit/target_time),
                np.log(m['peak_discharge_m3_s']/target_q),np.log(m['velocity_at_peak_discharge_m_s']/target_v)])))
        else:error=100.
        candidates.append({'stage':2,'mu':mu,'xi':xi,'section_capacity_scale':scale,
             'arrival_s':m['arrival_s'],'transit_s':transit,'peak_discharge_m3_s':m['peak_discharge_m3_s'],
             'velocity_m_s':m['velocity_at_peak_discharge_m_s'],'objective':error,'mass_error':m['relative_mass_error']})
        evaluations+=1
        if error<loss2:best2=(float(mu),float(xi),a);best_scale=float(scale);loss2=error
        if evaluations%10==0:print('[Q4 calibration]',evaluations,'best relative-log objective',loss2,flush=True)
        return error
    x0=np.log([initial[0],initial[1],1.])
    simplex=np.array([x0,x0+[.15,0,0],x0+[0,.25,0],x0+[0,0,.35]])
    optimize=minimize(objective,x0,method='Nelder-Mead',
       bounds=np.log(np.array([[.003,.12],[500,15000],[.5,3.]])),
       options={'maxfev':48,'xatol':.015,'fatol':.0005,'initial_simplex':simplex})
    mu2,xi2,lower=best2
    return save_cascade_results(ctx,profile,junction,evidence,upper,lower,
        (mu1,xi1),(mu2,xi2,best_scale),(loss1,loss2),stage1+candidates)


def calibration_fingerprint(profile,evidence):
    p=profile.copy()
    for name in p:p[name]=p[name].round(12 if name=='source_area_weight' else 6)
    code=''.join(inspect.getsource(f) for f in [section_depth,section_pressure,advance_channel,simulate_channel,upstream_transfer])
    return hashlib.sha256((p.to_csv(index=False)+yaml.safe_dump(evidence,sort_keys=True)+code).encode()).hexdigest()


def save_cascade_results(ctx,profile,junction,evidence,upper,lower,first_params,second_params,losses,trials):
    mu1,xi1=first_params;mu2,xi2,best_scale=second_params;loss1,loss2=losses
    values={k:v['value'] for k,v in evidence['published_values'].items()}
    fraction=values['mobile_debris_volume_m3']/upper['mass']['outflow_m3']
    inlet=upstream_transfer(upper,fraction);out=ctx.folder('q4')
    pd.DataFrame(trials).to_csv(out/'Q4_reconstruction_calibration_trials.csv',index=False)
    if len(trials)>10:pd.DataFrame(trials).to_csv(out/'Q4_reconstruction_training_trials.csv',index=False)
    else:pd.DataFrame(trials).to_csv(out/'Q4_reconstruction_sensitivity.csv',index=False)
    upper['history'].to_csv(out/'Q4_reconstruction_upper_history.csv',index=False)
    lower['history'].to_csv(out/'Q4_reconstruction_lower_history.csv',index=False)
    inlet.to_csv(out/'Q4_reconstruction_transfer.csv',index=False)
    for name,a in [('upper',upper),('lower',lower)]:
        np.savez_compressed(out/f'Q4_reconstruction_{name}_arrays.npz',distance_m=a['s'],
          time_s=[c['time_s'] for c in a['captures']],depth_m=[c['h'] for c in a['captures']],
          elevation_m=a['z'],width_m=a['width'],max_depth_m=a['max_depth'])
    retained=upper['mass']['outflow_m3']*(1-fraction)
    total=upper['mass']['stored_m3']+retained+lower['mass']['stored_m3']+lower['mass']['outflow_m3']
    report={'role':'post-event reconstruction calibrated against published event evidence',
      'calibration_reference':evidence['references']['case_study'],
      'published_constraints':values,'geometry_file':'Q4_reconstruction_geometry.json',
      'release_initial_velocity_m_s':0,'release_spatial_distribution':'official glacier pixels near mapped source, projected to route',
      'stage1_parameters':{'mu':mu1,'xi':xi1},
      'stage2_parameters':{'mu':mu2,'xi':xi2,'section_capacity_scale':best_scale},
      'capacity_scale_basis':'fit effective cross-section to reconstructed peak area Q/u and propagation constraints; not a surveyed width',
      'upstream':upper['mass'],'downstream':lower['mass'],
      'mobile_fraction_of_computed_upper_outflow':fraction,
      'published_mobile_fraction_of_release':values['mobile_debris_volume_m3']/values['release_volume_m3'],
      'mobile_input_derived_m3':float(inlet.cumulative_m3.iloc[-1]),
      'deposit_partition_m3':retained,'cascade_volume_residual':abs(total-values['release_volume_m3'])/values['release_volume_m3'],
      'calibration_objective':{'stage1_deposition_and_transition_constraint':loss1,'stage2_mean_squared_relative_log_error':loss2},
      'downstream_transit_s':lower['mass']['arrival_s']-upper['mass']['arrival_s'],
      'tongue_elevation_crossing_s':float(upper['arrival_field_s'][int(np.argmin(abs(upper['z']-values['glacier_impact_elevation_m'])))]),
      'input_hydrograph':'recorded upstream finite-volume flux times published mobile fraction',
      'calibration_is_independent_validation':False}
    write_json(out/'Q4_cascade_reconstruction.json',report)
    model={'calibration_basis_sha256':calibration_fingerprint(profile,evidence),
       'calibrated_from':evidence['references']['case_study']['url'],
       'parameter_role':'fitted effective parameters, with official DEM and published reconstruction constraints',
       'stage1':report['stage1_parameters'],'stage2':report['stage2_parameters'],
       'calibration_objective':report['calibration_objective']}
    if ctx.cfg.get('write_calibration_model',True):
        (ctx.root/'configs/calibrated_q4.yaml').write_text(yaml.safe_dump(model,sort_keys=False),encoding='utf-8')
    return upper,lower,report,pd.DataFrame(trials).loc[lambda d:d.stage.eq(2)]


def replay_calibrated_cascade(ctx,profile,junction,evidence,model):
    cfg=evidence['reconstruction'];values={k:v['value'] for k,v in evidence['published_values'].items()}
    first=model['stage1'];second=model['stage2'];scale=second['section_capacity_scale']
    upper=simulate_channel(profile.iloc[:junction+1],first['mu'],first['xi'],
        total_seconds=cfg['total_seconds'],cfl=cfg['cfl'],initial_volume_m3=values['release_volume_m3'])
    inlet=upstream_transfer(upper,values['mobile_debris_volume_m3']/upper['mass']['outflow_m3'])
    downstream=profile.iloc[junction:].copy()
    downstream['effective_width_m']*=scale;downstream['bank_side_slope']*=scale
    lower=simulate_channel(downstream,second['mu'],second['xi'],total_seconds=cfg['total_seconds'],cfl=cfg['cfl'],inlet=inlet)
    trials=[]
    for fm,fx in [(1,1),(.8,.8),(1,.8),(1.2,.8),(.8,1.2),(1,1.2),(1.2,1.2)]:
        a=lower if (fm,fx)==(1,1) else simulate_channel(downstream,second['mu']*fm,second['xi']*fx,
                   total_seconds=cfg['total_seconds'],cfl=cfg['cfl'],inlet=inlet)
        m=a['mass'];transit=m['arrival_s']-upper['mass']['arrival_s'] if m['arrival_s'] else None
        terms=([np.log(transit/values['confluence_to_port_elapsed_s']),
                np.log(m['peak_discharge_m3_s']/values['port_peak_discharge_m3_s']),
                np.log(m['velocity_at_peak_discharge_m_s']/values['port_cross_section_mean_velocity_m_s'])]
                if transit and m['peak_discharge_m3_s'] else [np.inf])
        trials.append({'stage':2,'mu':second['mu']*fm,'xi':second['xi']*fx,'section_capacity_scale':scale,
           'transit_s':transit,'arrival_s':m['arrival_s'],'peak_discharge_m3_s':m['peak_discharge_m3_s'],
           'velocity_m_s':m['velocity_at_peak_discharge_m_s'],'objective':float(np.mean(np.square(terms))),
           'mass_error':m['relative_mass_error']})
    loss1=float(np.log(upper['mass']['outflow_m3']/values['mobile_debris_volume_m3'])**2)
    return save_cascade_results(ctx,profile,junction,evidence,upper,lower,
       (first['mu'],first['xi']),(second['mu'],second['xi'],scale),(loss1,trials[0]['objective']),trials)


def plot_reconstruction(ctx,profile,junction,z,tr,geometry,upper,lower,report,trials):
    out=ctx.folder('q4');extent=[tr.c,tr.c+z.shape[1]*tr.a,tr.f+z.shape[0]*tr.e,tr.f]
    fig,axs=plt.subplots(1,3,figsize=(15.0,4.7),layout='constrained')
    axs[0].imshow(z,extent=extent,cmap='terrain')
    axs[0].plot(profile.utm_x,profile.utm_y,color='#B94135',lw=1.7)
    p=profile.iloc[[0,junction,-1]]
    axs[0].scatter(p.utm_x,p.utm_y,c=COLORS[:3],s=45,edgecolors='white')
    for r,label in zip(p.itertuples(),['释放区','河流汇入口','口岸河段']):axs[0].annotate(label,(r.utm_x,r.utm_y),xytext=(4,5),textcoords='offset points',fontsize=8)
    axs[0].set(xlabel='UTM东向 (m)',ylabel='UTM北向 (m)',title='a  DEM汇流约束路径')
    axs[1].plot(profile.distance_m/1000,profile.elevation_m,color=COLORS[0])
    axs[1].axvline(profile.distance_m.iloc[junction]/1000,ls='--',color=COLORS[1]);axs[1].set(xlabel='沿程距离 (km)',ylabel='高程 (m)',title='b  两阶段地形纵剖面')
    axs[2].plot(profile.distance_m/1000,profile.effective_width_m,color=COLORS[2])
    axs[2].set(xlabel='沿程距离 (km)',ylabel='DEM等效断面宽度 (m)',title='c  地形变化与释放区宽度')
    fig.suptitle('源区—汇入口—吉隆口岸的地形与断面')
    savefig(fig,out/'Q4_14_drainage_route_and_cross_sections.png')
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for a,name,col in [(upper,'上游冰岩崩',COLORS[0]),(lower,'下游泥石流',COLORS[1])]:
        h=a['history'];axs[0,0].plot(h.time_s,h.outlet_discharge_m3_s,label=name,color=col)
    axs[0,0].set(xlabel='源区释放后 (s)',ylabel='断面流量 (m³/s)',title='a  动力学求解的通量传递');axs[0,0].legend(frameon=False)
    h=lower['history'];axs[0,1].plot(h.time_s,h.outlet_velocity_m_s,color=COLORS[1])
    axs[0,1].axhline(report['published_constraints']['port_cross_section_mean_velocity_m_s'],ls='--',color='#666',label='文献断面平均速度')
    axs[0,1].set(xlabel='源区释放后 (s)',ylabel='口岸断面流速 (m/s)',title='b  口岸速度过程');axs[0,1].legend(frameon=False,fontsize=8)
    for a,ax,title in [(upper,axs[1,0],'c  冰岩崩流深演进'),(lower,axs[1,1],'d  河道泥石流流深演进')]:
        depth=np.array([c['h'] for c in a['captures']]);tt=[c['time_s'] for c in a['captures']]
        im=ax.imshow(depth,aspect='auto',origin='lower',extent=[0,a['s'][-1]/1000,tt[0],tt[-1]],cmap='turbo',vmin=0,vmax=max(1,float(np.quantile(depth,.99))))
        ax.set(xlabel='本河段距离 (km)',ylabel='源区释放后 (s)',title=title);fig.colorbar(im,ax=ax,label='流深 (m)',shrink=.75)
    fig.suptitle('分阶段冰岩崩—泥石流的通量、流速与流深')
    savefig(fig,out/'Q4_15_calibrated_cascade_process.png')
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    obs=report['published_constraints'];sim=[report['downstream_transit_s'],lower['mass']['peak_discharge_m3_s'],lower['mass']['velocity_at_peak_discharge_m_s']]
    truth=[obs['confluence_to_port_elapsed_s'],obs['port_peak_discharge_m3_s'],obs['port_cross_section_mean_velocity_m_s']]
    for ax,a,b,title,unit in zip(axs,sim,truth,['河道传播时间','口岸峰值流量','峰值流量对应速度'],['s','m³/s','m/s']):
        bars=ax.bar(['文献约束','本次计算'],[b,a],color=[COLORS[0],COLORS[1]])
        ax.bar_label(bars,labels=[f'{b:,.1f}',f'{a:,.1f}'],padding=3,fontsize=9)
        ax.set(ylabel=unit,title=title,ylim=(0,max(a,b)*1.2));ax.grid(axis='y',alpha=.14)
    fig.suptitle('事后标定的定量对照')
    savefig(fig,out/'Q4_16_event_constraint_comparison.png')
    fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    grid=trials.loc[trials.section_capacity_scale.isna()] if 'section_capacity_scale' in trials and trials.section_capacity_scale.isna().any() else trials
    for xi,g in grid.groupby('xi'):
        axs[0].plot(g.mu,g.transit_s,'-o',label=f'ξ={xi:g}')
        axs[1].plot(g.mu,g.peak_discharge_m3_s,'-o',label=f'ξ={xi:g}')
    axs[0].axhline(truth[0],color='#777',ls='--');axs[1].axhline(truth[1],color='#777',ls='--')
    axs[0].set(xlabel='库仑摩擦 μ',ylabel='河道传播时间 (s)',title='a  传播时间响应')
    axs[1].set(xlabel='库仑摩擦 μ',ylabel='口岸峰值流量 (m³/s)',title='b  峰值流量响应')
    for ax in axs:ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
    fig.suptitle('同一真实输入下的标定参数试验')
    savefig(fig,out/'Q4_17_calibration_parameter_experiments.png')


def execute_cascade_reconstruction(ctx):
    evidence=yaml.safe_load((ctx.root/'configs/case_evidence.yaml').read_text())
    profile,junction,z,tr,geometry=prepare_reconstruction_geometry(ctx,evidence)
    upper,lower,report,trials=fit_cascade(ctx,profile,junction,evidence)
    plot_reconstruction(ctx,profile,junction,z,tr,geometry,upper,lower,report,trials)
    from .q4_reconstructed_spatial import execute_reconstructed_spatial
    execute_reconstructed_spatial(ctx,lower,report,evidence)
    reconstruction_grid_experiment(ctx,upper,lower,report,evidence)
    print('[Q4 reconstruction]',json.dumps(report,ensure_ascii=False),flush=True)
    return upper,lower,report


def reconstruction_grid_experiment(ctx,upper,lower,report,evidence):
    cfg=evidence['reconstruction'];p=lower['profile'];out=ctx.folder('q4')
    fraction=report['mobile_input_derived_m3']/upper['mass']['outflow_m3'];inlet=upstream_transfer(upper,fraction)
    rows=[]
    for dx in [cfg['solver_dx_m']/2,cfg['solver_dx_m'],cfg['solver_dx_m']*2]:
        if dx==cfg['solver_dx_m']:a=lower
        else:
            s=np.arange(float(p.distance_m.iloc[0]),float(p.distance_m.iloc[-1]),dx)
            q=pd.DataFrame({name:np.interp(s,p.distance_m,p[name]) for name in p if name!='distance_m'})
            q.insert(0,'distance_m',s)
            a=simulate_channel(q,report['stage2_parameters']['mu'],report['stage2_parameters']['xi'],
                total_seconds=cfg['total_seconds'],cfl=cfg['cfl'],inlet=inlet)
        rows.append({'dx_m':dx,'arrival_s':a['mass']['arrival_s'],
           'peak_discharge_m3_s':a['mass']['peak_discharge_m3_s'],
           'peak_depth_m':a['mass']['maximum_depth_m'],'mass_error':a['mass']['relative_mass_error']})
    pd.DataFrame(rows).to_csv(out/'Q4_reconstruction_1d_grid_experiment.csv',index=False)
    sections=np.load(out/'Q4_reconstruction_sections.npz');offsets=sections['offsets_m'];z=sections['elevation_m']
    geometry=json.loads((out/'Q4_reconstruction_geometry.json').read_text())
    fig,axs=plt.subplots(1,3,figsize=(14,4.2),layout='constrained')
    for k,label in [(0,'源区'),(geometry['confluence_index'],'汇入口'),(len(z)-1,'口岸')]:
        axs[0].plot(offsets,z[k]-np.nanmin(z[k]),label=label)
    axs[0].set(xlabel='横向距中心 (m)',ylabel='相对最低高程 (m)',title='a  官方DEM实际横断面');axs[0].legend(frameon=False,fontsize=8)
    axs[1].plot([r['dx_m'] for r in rows],[r['arrival_s'] for r in rows],'-o',color=COLORS[0]);axs[1].set(xlabel='数值步长 (m)',ylabel='源区释放后到达时间 (s)',title='b  一维到达时间的网格响应')
    axs[2].plot([r['dx_m'] for r in rows],[r['peak_discharge_m3_s'] for r in rows],'-o',color=COLORS[1]);axs[2].set(xlabel='数值步长 (m)',ylabel='口岸峰值流量 (m³/s)',title='c  一维流量的网格响应')
    fig.suptitle('原始横断面与固定参数的数值分辨率试验')
    savefig(fig,out/'Q4_21_cross_sections_and_1d_grid_experiment.png')
    reconstruction_deposition_geometry(ctx)


def reconstruction_deposition_geometry(ctx):
    """DEM fill area per standard 100 m of channel, not an invented lake outline."""
    out=ctx.folder('q4');geometry=json.loads((out/'Q4_reconstruction_geometry.json').read_text())
    report=json.loads((out/'Q4_cascade_reconstruction.json').read_text())
    sections=np.load(out/'Q4_reconstruction_sections.npz');x=sections['offsets_m']
    section=sections['elevation_m'][geometry['confluence_index']];bed=float(np.nanmin(section))
    ds=float(np.median(np.diff(x)));rows=[]
    for height in [5.,10.,20.]:
        area=float(np.nansum(np.maximum(bed+height-section,0))*ds)
        rows.append({'relative_fill_height_m':height,'fill_cross_section_area_m2':area,
          'volume_per_100m_channel_m3':area*100,
          'fraction_of_model_confluence_deposit':area*100/report['deposit_partition_m3']})
    write_json(out/'Q4_confluence_deposition_geometry.json',{'source':'official DEM transverse section at model transition',
       'model_confluence_deposit_m3':report['deposit_partition_m3'],'unit_channel_length_m':100,
       'interpretation':'geometric occupancy per unit reach; the reporting length is not a measured dam length',
       'results':rows})
