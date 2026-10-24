"""Q4: conservative 1D shallow mixture (Rusanov + explicit gravity/friction) and optional local 2D
finite-volume (Rusanov) spread. Explicit scenario-only provenance and mass-balance audit.

The PDE is an effective continuum; it does not model multi-phase ice->debris->flood conversion.
Routes from DEM cost-path are provisional and must not be cited as actual avalanche route.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np,pandas as pd,rasterio,matplotlib.pyplot as plt
from matplotlib import animation
from rasterio.transform import rowcol
from scipy.ndimage import gaussian_filter, distance_transform_edt
from .q4_visuals import plot_local_2d, plot_2d_inlet_diagnostic, plot_geographic_reference, plot_inlet_factorial
from .q4_flowpath import local_dem,UTM_CRS
from .common import COLORS,savefig,write_json,record_warning

G=9.81
EPS=1e-8

def analytic_voellmy_velocity(t,theta_deg,mu,xi,h):
    """Uniform-slope, uniform-depth analytic solution from rest (when a>0).
    u(t) = sqrt(a/b) tanh(sqrt(a*b)t), a=g(sin theta-mu cos theta), b=g/(xi*h).
    """
    theta=np.deg2rad(theta_deg);a=G*(np.sin(theta)-mu*np.cos(theta));b=G/(xi*h)
    if a<=0 or b<=0:raise ValueError('Closed-form case requires a>0,b>0')
    return np.sqrt(a/b)*np.tanh(np.sqrt(a*b)*np.asarray(t))

def _friction_velocity(v,h,dt,mu,xi,cos_theta):
    speed=np.abs(v)
    friction=mu*G*cos_theta + G*speed*speed/(xi*np.maximum(h,0.01))
    magnitude=np.maximum(speed-dt*friction,0)
    return np.sign(v)*magnitude

def flux_1d(h,q):
    vel=q/np.maximum(h,EPS)
    return np.stack((q,q*vel+0.5*G*h*h),axis=0)

def rusanov_1d(L,R):
    hL,qL=L;hR,qR=R
    uL=qL/np.maximum(hL,EPS);uR=qR/np.maximum(hR,EPS)
    a=np.maximum(np.abs(uL)+np.sqrt(G*np.maximum(hL,0)),np.abs(uR)+np.sqrt(G*np.maximum(hR,0)))
    return 0.5*(flux_1d(hL,qL)+flux_1d(hR,qR))-0.5*a*(R-L)

def hydrostatic_flux_1d(L,R,zL,zR):
    """Audusse hydrostatic reconstruction: preserve eta=h+z at rest."""
    bed=np.maximum(zL,zR)
    hL=np.maximum(0,L[0]+zL-bed);hR=np.maximum(0,R[0]+zR-bed)
    LL=np.stack([hL,L[1]*hL/np.maximum(L[0],EPS)])
    RR=np.stack([hR,R[1]*hR/np.maximum(R[0],EPS)])
    F=rusanov_1d(LL,RR);FL=F.copy();FR=F.copy()
    FL[1]+=.5*G*(L[0]**2-hL**2)
    FR[1]+=.5*G*(R[0]**2-hR**2)
    return FL,FR


def advance_1d(h,q,sin_theta,cos_theta,dx,dt,mu,xi,release_m3_s,width,initial_velocity_m_s=0.0,bed=None):
    """Update volume per unit width and momentum; boundary is wall at source and free outflow at sink.
    Released material is injected into first cell with specified mean inlet velocity; then accelerated by gravity.
    Return new h,q,outflow_m3 and clipped_negative_volume_m3.
    """
    n=len(h);U=np.stack([h,q]);fluxes=np.zeros((2,n+1))
    fluxes[:,1:n]=rusanov_1d(U[:,:-1],U[:,1:])
    # reflective wall at upstream boundary -> zero mass flux
    fluxes[0,0]=0;fluxes[1,0]=.5*G*h[0]**2
    # right boundary transmissive for positive outgoing velocity, no inward injection
    fluxes[:,n]=flux_1d(h[-1],max(float(q[-1]),0))
    if bed is None:
        Unew=U-(dt/dx)*(fluxes[:,1:]-fluxes[:,:-1])
    else:
        left,right=hydrostatic_flux_1d(U[:,:-1],U[:,1:],bed[:-1],bed[1:])
        outgoing=fluxes[:,1:].copy();incoming=fluxes[:,:-1].copy()
        outgoing[:,:-1]=left;incoming[:,1:]=right
        Unew=U-(dt/dx)*(outgoing-incoming)
    hn=Unew[0];qn=Unew[1]
    negative=np.minimum(hn,0)
    lost_m3=float(-np.sum(negative)*dx*width)
    if np.nanmin(hn)<-1e-4:raise FloatingPointError(f'negative depth beyond tolerance: {np.nanmin(hn)}')
    hn=np.maximum(hn,0)
    # slope gravity drive and dissipative friction; the friction step cannot reverse velocity
    raw_v=(qn+(dt*G*hn*sin_theta if bed is None else 0))/np.maximum(hn,EPS)
    vf=_friction_velocity(raw_v,hn,dt,mu,xi,cos_theta)
    qn=hn*vf
    # Material source contributes both volume and momentum (the case u0=0 recovers original model).
    dh0=release_m3_s*dt/(width*dx)
    hn[0]+=dh0
    qn[0]+=dh0*initial_velocity_m_s
    if not np.isfinite(hn).all() or not np.isfinite(qn).all():raise FloatingPointError('nonfinite finite-volume solver state')
    qout=max(0,float(fluxes[0,n]))*width
    return hn,qn,qout*dt,lost_m3

def simulate_1d(profile,mu=.12,xi=400,width=30,volume=6e5,duration=300,
                total_seconds=2400,cfl=.35,frames=12,sensor_upstream_m=2500,initial_velocity_m_s=0.0):
    prof=profile.reset_index(drop=True)
    s=prof.distance_m.to_numpy(float);z=prof.elevation_m.to_numpy(float)
    dx=float(np.median(np.diff(s)))
    if np.max(np.abs(np.diff(s)-dx))>1e-5:
        # Last partial cell common: resample to strict uniform spacing
        sn=np.arange(0,s[-1],dx);zn=np.interp(sn,s,z);s,z=sn,zn
    n=len(s)
    if n<15:raise ValueError('Real slope profile too short for finite volume 1D solver')
    grade=-np.gradient(z,dx)
    theta=np.arctan(grade)
    sin_t=np.sin(theta);cos_t=np.cos(theta)
    h=np.zeros(n);q=np.zeros(n)
    t=0.0;input_vol=0.0;output_vol=0.0;numeric_loss=0.0
    steps=0;frame_t=np.linspace(0,total_seconds,frames)
    runout_max_m=0.;max_h_seen=0.;max_v_seen=0.;max_q_sensor=0.
    captures=[];samples=[{'time_s':0.,'mass_domain_m3':0.,'released_cum_m3':0.,'outflow_cum_m3':0.,'numerical_negative_loss_m3':0.,'sensor_discharge_m3_s':0.,'outlet_discharge_m3_s':0.,'sensor_velocity_m_s':0.,'sensor_depth_m':0.,'max_depth_m':0.,'max_velocity_m_s':0.}];next_frame=0
    sens_idx=int(np.argmin(np.abs(s-(s[-1]-sensor_upstream_m))))
    arrival_time_sensor=None;arrival_time_outlet=None
    while t<total_seconds-1e-9:
        u=q/np.maximum(h,EPS)
        c=np.sqrt(G*np.maximum(h,0))
        maxwave=float(np.max(np.abs(u)+c))
        dt=min(1.0,float(cfl*dx/max(maxwave,.5)),total_seconds-t)
        if dt<1e-5:raise RuntimeError('CFL time step too small - check slope, release volume and bathymetry')
        qin=volume/duration if t<duration else 0.0
        if t+dt>duration and t<duration:dt=min(dt,duration-t)
        h,q,dvout,lost=advance_1d(h,q,sin_t,cos_t,dx,dt,mu,xi,qin,width,initial_velocity_m_s=initial_velocity_m_s,bed=z)
        t+=dt;steps+=1;input_vol+=qin*dt;output_vol+=dvout;numeric_loss+=lost
        wet=np.flatnonzero(h>.20)
        if len(wet):runout_max_m=max(runout_max_m,float(s[wet[-1]]))
        max_h_seen=max(max_h_seen,float(np.max(h)))
        max_v_seen=max(max_v_seen,float(np.max(np.where(h>.2,np.abs(q/np.maximum(h,EPS)),0.))))
        sensor_flux,_=hydrostatic_flux_1d(np.array([h[sens_idx],q[sens_idx]]),
            np.array([h[sens_idx+1],q[sens_idx+1]]),z[sens_idx],z[sens_idx+1])
        sensor_q=float(max(0,sensor_flux[0])*width)
        port_q=float(max(0,q[-1])*width)
        max_q_sensor=max(max_q_sensor,sensor_q)
        if arrival_time_sensor is None and h[sens_idx]>.2:arrival_time_sensor=t
        if arrival_time_outlet is None and h[-1]>.2:arrival_time_outlet=t
        if True:  # retain every accepted solver step for conservative 1D -> 2D forcing
            mass=float(h.sum()*dx*width)
            samples.append({'time_s':t,'mass_domain_m3':mass,'released_cum_m3':input_vol,
                'outflow_cum_m3':output_vol,'numerical_negative_loss_m3':numeric_loss,
                'sensor_discharge_m3_s':sensor_q,'outlet_discharge_m3_s':port_q,
                'sensor_velocity_m_s':float(q[sens_idx]/max(h[sens_idx],EPS)) if h[sens_idx]>.2 else 0.,'sensor_depth_m':float(h[sens_idx]),
                'max_depth_m':float(h.max()),'max_velocity_m_s':float(np.max(np.where(h>.2,np.abs(q/np.maximum(h,EPS)),0.)))})
        while next_frame<len(frame_t) and t>=frame_t[next_frame]-1e-6:
            captures.append({'time_s':float(t),'h':h.copy(),'q':q.copy()})
            next_frame+=1
        if steps>200000:raise RuntimeError('More than 200,000 time steps. Stop rather than return an unverified result')
    while next_frame<len(frame_t):
        captures.append({'time_s':float(t),'h':h.copy(),'q':q.copy()});next_frame+=1
    final_vol=float(h.sum()*dx*width)
    residual=final_vol+output_vol-input_vol
    residual_ratio=abs(residual)/max(input_vol,1)
    mass={'injected_m3':float(input_vol),'outflow_m3':float(output_vol),
      'material_in_domain_m3':final_vol,'relative_absolute_mass_residual':residual_ratio,
      'numerical_negative_clip_m3':numeric_loss,'n_steps':steps,'dx_m':dx,
      'first_arrival_to_coupling_sensor_s':arrival_time_sensor,
      'first_arrival_to_downstream_end_s':arrival_time_outlet,
      'coupling_sensor_distance_m':float(s[sens_idx]),
      'max_runout_m':float(runout_max_m),
      'max_runout_fraction_of_route':float(runout_max_m/max(float(s[-1]),1)),
      'max_depth_over_time_m':float(max_h_seen),
      'max_depth_to_assumed_width_ratio':float(max_h_seen/width),
      'max_adjacent_bed_step_m':float(np.max(np.abs(np.diff(z)))),
      'coarse_steep_dry_front_caution':True,
      'velocity_statistics_wet_threshold_m':.2,
      'depth_averaging_geometry_caution':bool(max_h_seen/width>.3),
      'max_speed_over_time_m_s':float(max_v_seen),
      'max_coupling_discharge_m3_s':float(max_q_sensor),
      'initial_velocity_m_s':float(initial_velocity_m_s),
      'model_scope':'competition simulation conditional on stated assumptions; not a uniquely identified physical reconstruction'}
    return {'s':s,'z':z,'captures':captures,'history':pd.DataFrame(samples),
            'mass':mass,'sensor_idx':sens_idx,'last_h':h,'last_q':q,'profile':prof}

def _save_1d(ctx,ans,report_tag='SCENARIO'):
    out=ctx.folder('q4');s,z=ans['s'],ans['z'];capture=ans['captures']
    ans['history'].to_csv(out/'Q4_1d_discharge_and_mass_timeseries.csv',index=False)
    pd.DataFrame([ans['mass']]).to_csv(out/'Q4_1d_mass_balance.csv',index=False)
    depth=np.array([a['h'] for a in capture]);frames=np.array([a['time_s'] for a in capture]);
    np.savez_compressed(out/'Q4_1d_frames_arrays.npz',distance_m=s,time_s=frames,depth_m=depth,
                        elevation_m=z)
    fig,axes=plt.subplots(2,1,figsize=(10,7.2),gridspec_kw={'height_ratios':[1,1.2]})
    axes[0].plot(s/1000,z,color='#63767F',lw=1.5)
    axes[0].set(ylabel='高程 (m)',title='传播路径地形纵剖面');axes[0].grid(True,alpha=.2)
    im=axes[1].imshow(depth,aspect='auto',extent=[0,s[-1]/1000,frames[-1]/60,frames[0]/60],
                     cmap='turbo',vmin=0,vmax=max(np.percentile(depth,99),.1))
    axes[1].invert_yaxis();axes[1].set(xlabel='沿程距离 (km)',ylabel='模拟启动后分钟数',title='一维流深时空分布')
    fig.colorbar(im,ax=axes[1],label='流深 (m)',shrink=.85)
    fig.suptitle('一维泥石流传播过程及地形纵剖面',fontsize=12)
    savefig(fig,out/'Q4_02_1d_runout_timespace.png')
    hist=ans['history']
    fig,axs=plt.subplots(1,2,figsize=(11,4.3))
    axs[0].plot(hist.time_s/60,hist.sensor_discharge_m3_s,color=COLORS[0],lw=1.6,label='二维耦合断面处')
    axs[0].plot(hist.time_s/60,hist.outlet_discharge_m3_s,color=COLORS[1],lw=1.6,label='下游出口')
    axs[0].set(xlabel='模拟时间 (min)',ylabel='体积流量 (m³/s)',title='耦合断面流量过程');axs[0].legend(frameon=False)
    axs[1].plot(hist.time_s/60,hist.released_cum_m3,label='累计释放量',color=COLORS[0])
    axs[1].plot(hist.time_s/60,hist.mass_domain_m3+hist.outflow_cum_m3,
                label='域内量 + 累计出流',color=COLORS[2],ls='--')
    axs[1].set(xlabel='模拟时间 (min)',ylabel='累计体积 (m³)',title='流体体积收支变化')
    axs[1].legend(frameon=False)
    for ax in axs:ax.grid(True,alpha=.16)
    savefig(fig,out/'Q4_03_1d_hydrograph_and_mass_conservation.png')
    # 12 sequential snapshots: an explicit visualizable process, not merely static map
    fig,axs=plt.subplots(3,4,figsize=(13,8.2),sharex=True,sharey=True)
    for ax,record in zip(axs.flat,capture):
        ax.plot(s/1000,record['h'],color=COLORS[0],lw=1.6)
        ax.fill_between(s/1000,record['h'],color=COLORS[0],alpha=.2)
        ax.set_title(f'{record["time_s"]/60:.1f} min',fontsize=9);ax.grid(True,alpha=.15)
    for ax in axs[-1]:ax.set_xlabel('距离 (km)')
    for ax in axs[:,0]:ax.set_ylabel('流深 (m)')
    fig.suptitle('一维泥石流流深演进过程',fontsize=12)
    savefig(fig,out/'Q4_04_1d_12frames.png')
    # GIF is for supplementary display only. Keep PNG/SVG as paper-ready main figures.
    try:
        fig,ax=plt.subplots(figsize=(8,3.7));line,=ax.plot([],[],color=COLORS[0],lw=2)
        ax.set(xlim=(0,s[-1]/1000),ylim=(0,max(1,np.max(depth)*1.1)),xlabel='沿程公里',ylabel='流深 (m)',title='一维泥石流流深动态变化')
        def upd(k):line.set_data(s/1000,depth[k]);ax.set_title(f'一维流深演进 · t={frames[k]/60:.1f} min');return (line,)
        ani=animation.FuncAnimation(fig,upd,frames=len(frames),blit=False,interval=300)
        ani.save(out/'Q4_1d_runout_animation.gif',writer='pillow',fps=3);plt.close(fig)
    except Exception as e:record_warning(f'GIF skipped: {e}')

def flux_x(U):
    h=U[0];qx=U[1];qy=U[2];ux=qx/np.maximum(h,EPS)
    return np.stack([qx,qx*ux+0.5*G*h*h,qy*ux])

def flux_y(U):
    h=U[0];qx=U[1];qy=U[2];uy=qy/np.maximum(h,EPS)
    return np.stack([qy,qx*uy,qy*uy+0.5*G*h*h])

def _rusanov_2d(L,R,axis):
    fun=flux_x if axis=='x' else flux_y
    hL,hR=L[0],R[0]
    uL=L[1 if axis=='x' else 2]/np.maximum(hL,EPS)
    uR=R[1 if axis=='x' else 2]/np.maximum(hR,EPS)
    a=np.maximum(np.abs(uL)+np.sqrt(G*np.maximum(hL,0)),np.abs(uR)+np.sqrt(G*np.maximum(hR,0)))
    return .5*(fun(L)+fun(R))-.5*a*(R-L)

def hydrostatic_flux_2d(L,R,zL,zR,axis):
    bed=np.maximum(zL,zR)
    hl=np.maximum(0,L[0]+zL-bed);hr=np.maximum(0,R[0]+zR-bed)
    LL=L*np.divide(hl,np.maximum(L[0],EPS))[None,...]
    RR=R*np.divide(hr,np.maximum(R[0],EPS))[None,...]
    LL[0]=hl;RR[0]=hr
    F=_rusanov_2d(LL,RR,axis);FL=F.copy();FR=F.copy()
    k=1 if axis=='x' else 2
    FL[k]+=.5*G*(L[0]**2-hl**2)
    FR[k]+=.5*G*(R[0]**2-hr**2)
    return FL,FR


def _fill_nan_elevation(z):
    mask=~np.isfinite(z)
    if not mask.any():return z
    _,ids=distance_transform_edt(mask,return_indices=True)
    return z[tuple(ids)]

def choose_inlet_on_valley(z, source_rc, flow_direction_xy=None, dx=60.,
                           mode='valley_cross_section', halfwidth_cells=5):
    """Select inlet cells from *measured DEM*, not an unconstrained uphill 3x3 patch.

    With flow direction, search across the local *normal* cross-section (not along the
    river) for the lowest defensible valley cell. Source snapping is a model choice,
    recorded for scrutiny, never silently treated as a measured stream centreline.
    We inject only in the centre pixel and adjacent cross-section pixels that lie
    within 40 m of its elevation; weights are elevation-based and sum to exactly 1.
    """
    nr,nc=z.shape
    r,c=map(int,source_rc)
    if flow_direction_xy is None or np.linalg.norm(flow_direction_xy)<1e-9:
        direction=np.array([-1.,-1.])/np.sqrt(2)
    else:
        direction=np.array(flow_direction_xy,dtype=float)
        direction/=np.linalg.norm(direction)
    normal=np.array([-direction[1],direction[0]]) # in east,north coordinates
    cand=[]
    ks=np.arange(-halfwidth_cells,halfwidth_cells+1) if mode=='valley_cross_section' else np.array([0])
    for k in ks:
        rr=int(round(r-k*normal[1]));cc=int(round(c+k*normal[0]))
        if 1<=rr<nr-1 and 1<=cc<nc-1 and np.isfinite(z[rr,cc]):
            # weak lateral penalty avoids snapping far away for negligible elevation difference
            objective=float(z[rr,cc])+.05*abs(k)*dx
            cand.append((objective,rr,cc,int(k)))
    if not cand:raise ValueError('Cannot place valid inlet on DEM')
    _,ar,ac,k=min(cand)
    # A one-cell inlet would generate an unrealistically deep point source.  Instead
    # distribute over the 3 lowest neighbouring cells near the channel thalweg.
    # This does NOT send material into uphill bank pixels without any justification.
    candidates=[]
    for rr in range(ar-1,ar+2):
        for cc in range(ac-1,ac+2):
            if 0<=rr<nr and 0<=cc<nc and np.isfinite(z[rr,cc]):
                e=float(z[rr,cc]);dist=np.hypot(rr-ar,cc-ac)
                candidates.append((e+12*dist,rr,cc,e))
    candidates.sort()
    eligible=[v for v in candidates if v[3]<=float(z[ar,ac])+45]
    if not eligible:eligible=[(float(z[ar,ac]),ar,ac,float(z[ar,ac]))]
    selected=eligible[:3]
    cells=np.array([[item[1],item[2]] for item in selected],dtype=int)
    # Allocate slightly more to lower cells without inventing water on the steep bank.
    heights=np.array([item[3] for item in selected],dtype=float)
    weight=np.exp(-np.maximum(heights-heights.min(),0)/35.)
    weight=weight/weight.sum()
    return {'nominal_rc':(r,c),'snapped_rc':(ar,ac),'cross_section_offset_cells':int(k),
            'cells':cells,'weights':weight,'elevation_nominal_m':float(z[r,c]),
            'elevation_snapped_m':float(z[ar,ac]),'flow_direction_xy':direction,
            'mode':mode,'search_halfwidth_m':float(halfwidth_cells*dx)}

def advance_2d(U,z,dx,dt,mu,xi,source_i,source_j,inject_m3,velocity_xy=(0,0),
               inlet_cells=None,inlet_weights=None):
    ny,nx=z.shape
    # Separate left/right bed-pressure corrections balance still water on nonflat DEM.
    fx=np.zeros((3,ny,nx+1));fy=np.zeros((3,ny+1,nx))
    fl,fr=hydrostatic_flux_2d(U[:,:,:-1],U[:,:,1:],z[:,:-1],z[:,1:],'x')
    ft,fb=hydrostatic_flux_2d(U[:,:-1,:],U[:,1:,:],z[:-1,:],z[1:,:],'y')
    # Outflow-only boundaries; closed normal faces retain hydrostatic pressure.
    for col,allowed in [(0,U[1,:,0]<0),(-1,U[1,:,-1]>0)]:
        f=flux_x(U[:,:,col]);wall=np.zeros_like(f);wall[1]=.5*G*U[0,:,col]**2
        fx[:,:,col]=np.where(allowed[None,:],f,wall)
    for row,allowed in [(0,U[2,0,:]<0),(-1,U[2,-1,:]>0)]:
        f=flux_y(U[:,row,:]);wall=np.zeros_like(f);wall[2]=.5*G*U[0,row,:]**2
        fy[:,row,:]=np.where(allowed[None,:],f,wall)
    xout=fx[:,:,1:].copy();xin=fx[:,:,:-1].copy()
    yout=fy[:,1:,:].copy();yin=fy[:,:-1,:].copy()
    xout[:,:,:-1]=fl;xin[:,:,1:]=fr
    yout[:,:-1,:]=ft;yin[:,1:,:]=fb
    fx[:,:,1:nx]=fl;fy[:,1:ny,:]=ft  # mass flux is identical on both sides
    Unew=U-(dt/dx)*(xout-xin+yout-yin)
    h=np.maximum(Unew[0],0)
    negative_cut=float(-np.minimum(Unew[0],0).sum()*dx*dx)
    if np.nanmin(Unew[0])<-1e-3:raise FloatingPointError('2D negative depth too large')
    u=Unew[1]/np.maximum(h,EPS)
    v=Unew[2]/np.maximum(h,EPS)
    speed=np.hypot(u,v)
    friction=mu*G+G*speed**2/(xi*np.maximum(h,.05))
    reduced=np.maximum(0,speed-dt*friction)
    factor=np.divide(reduced,speed,out=np.zeros_like(reduced),where=speed>1e-7)
    u*=factor;v*=factor
    qx=h*u;qy=h*v
    dry=h<1e-5
    qx[dry]=0;qy[dry]=0
    if inlet_cells is None:
        # Compatibility/reference baseline: original 3x3 injection used in v1.1.
        ii0=max(0,source_i-1);ii1=min(ny,source_i+2)
        jj0=max(0,source_j-1);jj1=min(nx,source_j+2)
        coords=np.array([(i,j) for i in range(ii0,ii1) for j in range(jj0,jj1)])
        weights=np.ones(len(coords),dtype=float)/len(coords)
    else:
        coords=np.asarray(inlet_cells,dtype=int)
        weights=np.asarray(inlet_weights,dtype=float)
        if not np.isclose(weights.sum(),1,rtol=0,atol=1e-12):
            raise ValueError('inlet weights must sum to one')
    for (i,j),w in zip(coords,weights):
        dh=inject_m3*w/(dx*dx)
        h[i,j]+=dh
        qx[i,j]+=dh*velocity_xy[0]
        # UTM northing increases upwards; raster flux_y increases with row (south).
        qy[i,j]-=dh*velocity_xy[1]
    net_flux_out=dt*dx*(fx[0,:,-1].sum()-fx[0,:,0].sum()+fy[0,-1,:].sum()-fy[0,0,:].sum())
    return np.stack([h,qx,qy]),net_flux_out,negative_cut

def simulate_2d(z,dx,source_rc,inlet_curve,total_seconds=2400,mu=.12,xi=400,
                cfl=.24,frames=10,inlet_mode='valley_cross_section',flow_direction_xy=None,
                density_kg_m3=1800.0):
    """2D simulation on real DEM; time zero = 1D release, NOT 2D inlet activation.

    If source is delayed, select snapshots relative to FIRST 2D inlet entry rather
    than wasting 2 of 10 panels on entirely dry terrain. Never rescale each frame.
    """
    z=_fill_nan_elevation(z.astype(np.float64))
    if min(z.shape)<15:raise ValueError('ROI not suitable for 2D dynamics')
    ny,nx=z.shape;U=np.zeros((3,ny,nx),dtype=np.float64)
    t=0.;steps=0;injected=0.;net_out=0.;cut=0.
    times=np.array(inlet_curve.time_s,float)
    flow=np.maximum(np.array(inlet_curve.sensor_discharge_m3_s,float),0)
    if 'sensor_velocity_m_s' in inlet_curve:
        source_velocity=np.maximum(np.array(inlet_curve.sensor_velocity_m_s,float),0)
    else: source_velocity=np.zeros_like(flow)
    qthreshold=max(.005*float(np.max(flow)),1.)
    if not np.any(flow>qthreshold):
        first_inlet=float('nan');frame_times=np.linspace(0,total_seconds,frames)
    else:
        first_inlet=float(times[np.argmax(flow>qthreshold)])
        # Key phases of flood advance; absolute time remains the source-release clock.
        # Use event-driven frames rather than 0,3.3,6.7,10,... which hide the onset.
        offsets=np.array([12,28,48,78,118,178,278,418,618,898],dtype=float)
        frame_times=first_inlet+np.interp(np.linspace(0,len(offsets)-1,frames),np.arange(len(offsets)),offsets)
        if frame_times[-1]>total_seconds:
            frame_times=np.linspace(min(first_inlet+5,total_seconds),total_seconds,frames)
    capt=[];next_f=0;wet_onset_s=None
    max_depth=np.zeros((ny,nx));arrivals=np.full((ny,nx),np.nan)
    # Proxy mechanical dynamic pressure q = 0.5*rho*|u|^2. This is not the
    # structure-specific impact pressure, which needs an impact coefficient.
    if density_kg_m3<=0 or not np.isfinite(density_kg_m3):
        raise ValueError('density_kg_m3 must be positive and finite')
    max_dynamic_pressure=np.zeros((ny,nx),dtype=float)
    max_velocity=np.zeros((ny,nx),dtype=float)
    if inlet_mode not in ('legacy_3x3','legacy_directional','valley_zero_momentum','valley_cross_section'):
        raise ValueError('Unknown inlet mode: '+str(inlet_mode))
    legacy_position=inlet_mode in ('legacy_3x3','legacy_directional')
    uses_momentum=inlet_mode in ('legacy_directional','valley_cross_section')
    placement=choose_inlet_on_valley(z,source_rc,flow_direction_xy,dx,
         mode='legacy_3x3' if legacy_position else 'valley_cross_section',halfwidth_cells=5)
    inlet_cells=None if legacy_position else placement['cells']
    inlet_weights=None if legacy_position else placement['weights']
    uvec=placement['flow_direction_xy']
    hist=[];last_report=-1
    ii,jj=np.mgrid[:ny,:nx]
    # projection along downstream path (column/east, -row/north)
    dr=-(ii-placement['snapped_rc'][0])*dx
    dc=(jj-placement['snapped_rc'][1])*dx
    projected=dc*uvec[0]+dr*uvec[1]
    while t<total_seconds-1e-8:
        h=U[0];u=U[1]/np.maximum(h,EPS);v=U[2]/np.maximum(h,EPS)
        maxspeed=float(np.nanmax(np.hypot(u,v)+np.sqrt(G*np.maximum(h,0))))
        dt=min(1.1,cfl*dx/max(maxspeed,.2),total_seconds-t)
        if dt<1e-5:raise FloatingPointError('2D dt collapsed')
        knot=np.searchsorted(times,t+1e-9,side='right')
        if knot<len(times):dt=min(dt,float(times[knot]-t))
        q0=float(np.interp(t,times,flow,left=0,right=0));q1=float(np.interp(t+dt,times,flow,left=0,right=0))
        qin=.5*(q0+q1)
        speedin=float(np.interp(t+.5*dt,times,source_velocity,left=0,right=0))
        velocity_xy=tuple(uvec*speedin) if uses_momentum else (0,0)
        U,diver,negative=advance_2d(U,z,dx,dt,mu,xi,source_rc[0],source_rc[1],qin*dt,
                                     velocity_xy=velocity_xy,inlet_cells=inlet_cells,
                                     inlet_weights=inlet_weights)
        t+=dt;steps+=1;injected+=qin*dt;net_out+=diver;cut+=negative
        max_depth=np.maximum(max_depth,U[0]);new=(U[0]>.2)&np.isnan(arrivals);arrivals[new]=t
        wet=U[0]>.2
        if wet_onset_s is None and wet.any():
            wet_onset_s=t
            frame_times=frame_times+max(0.,t+1.-frame_times[0])
        velocity=np.hypot(U[1],U[2])/np.maximum(U[0],EPS)
        velocity=np.where(wet,velocity,0.0)
        max_velocity=np.maximum(max_velocity,velocity)
        max_dynamic_pressure=np.maximum(max_dynamic_pressure,
                density_kg_m3*velocity**2/(2*1000.0))  # kPa
        while wet_onset_s is not None and next_f<len(frame_times) and t>=frame_times[next_f]-1e-6:
            capt.append({'time_s':t,'h':U[0].copy(),'cum_injected_m3':injected,
                         'wet_area_km2':np.count_nonzero(U[0]>.2)*dx*dx/1e6})
            next_f+=1
        if int(t/10)>last_report:
            last_report=int(t/10)
            upper=float(U[0][projected < -2*dx].sum()*dx**2)
            stored=float(U[0].sum()*dx**2)
            hist.append({'time_s':t,'inlet_m3_s':qin,'input_cumulative_m3':injected,
                         'stored_m3':stored,'wet_area_km2':np.count_nonzero(U[0]>.2)*dx**2/1e6,
                         'max_depth_m':float(U[0].max()),
                         'upstream_lobe_volume_fraction':upper/max(stored,1.)})
        if steps>150000:raise RuntimeError('2D solver did not complete within safe step count')
    while next_f<len(frame_times):
        capt.append({'time_s':t,'h':U[0].copy(),'cum_injected_m3':injected,
                     'wet_area_km2':np.count_nonzero(U[0]>.2)*dx*dx/1e6})
        next_f+=1
    final_v=float(np.sum(U[0])*dx*dx)
    err=(final_v+net_out-injected)
    mass={'input_cumulative_m3':float(injected),'signed_boundary_outflow_m3':float(net_out),
          'stored_final_m3':final_v,'relative_mass_error':abs(err)/max(1,injected),
          'clipped_negative_m3':float(cut),'n_steps':int(steps),
          'peak_depth_m':float(np.max(max_depth)),
          'wet_area_h_gt_0p2_km2':float(np.sum(max_depth>.2)*dx**2/1e6)}
    mass['max_dynamic_pressure_kPa']=float(np.max(max_dynamic_pressure))
    mass['max_wet_cell_velocity_m_s']=float(np.max(max_velocity))
    mass['dynamic_pressure_coefficient']=0.5
    mass['assumed_density_kg_m3_for_pressure']=float(density_kg_m3)
    serial_place={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in placement.items()}
    return {'z':z,'captures':capt,'max_depth':max_depth,'arrival':arrivals,
            'max_dynamic_pressure_kPa':max_dynamic_pressure,'max_wet_velocity_m_s':max_velocity,
            'mass':mass,'U':U,'timeline':pd.DataFrame(hist),'inlet':serial_place,
            'first_inlet_s':first_inlet,'first_wet_s':wet_onset_s,'frame_times_s':frame_times.tolist(),
            'dx_m':dx,'inlet_mode':inlet_mode}

def _plot_dynamic_pressure(ctx,ans,transform):
    """Raster of maximum dynamic-pressure proxy, mapped on official local DEM.
    Only cells that became wet (h>0.2 m) appear; saturation is percentile based.
    """
    z=ans['z'];depth=ans['max_depth'];pressure=ans['max_dynamic_pressure_kPa']
    from scipy.ndimage import gaussian_filter
    from matplotlib.colors import LightSource
    from matplotlib import cm
    from matplotlib.ticker import MaxNLocator
    surface=gaussian_filter(z,1)
    hill=LightSource(azdeg=315,altdeg=45).hillshade(surface,vert_exag=1.1,
                        dx=float(transform.a),dy=float(abs(transform.e)))
    ext=[transform.c,transform.c+z.shape[1]*transform.a,
         transform.f+z.shape[0]*transform.e,transform.f]
    fig,ax=plt.subplots(figsize=(9.6,6.2),layout='constrained')
    ax.imshow(hill,extent=ext,origin='upper',cmap='Greys',alpha=.8)
    pressureplot=np.ma.masked_where(depth<=.2,pressure)
    exposed=pressure[depth>.2]
    high=float(np.quantile(exposed,.98)) if len(exposed) else 1.0
    im=ax.imshow(pressureplot,extent=ext,origin='upper',cmap='magma',vmin=0,vmax=max(.01,high),alpha=.84)
    cb=fig.colorbar(im,ax=ax,shrink=.80,fraction=.04,pad=.02)
    cb.set_label('峰值动压代理 (kPa)',fontsize=10)
    ax.set(xlabel='UTM 东向坐标 (m)',ylabel='UTM 北向坐标 (m)')
    ax.set_aspect('equal')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    fig.suptitle('二维流动峰值动压分布',fontsize=13)
    savefig(fig,ctx.folder('q4')/'Q4_10_impact_pressure_distribution.png')
    write_json(ctx.folder('q4')/'Q4_dynamic_pressure_diagnostic.json',{
       'equation':'q = 0.5 * rho * velocity^2',
       'pressure_units':'kPa','wet_threshold_m':.2,
       'density_kg_m3':ans['mass']['assumed_density_kg_m3_for_pressure'],
       'peak_dynamic_pressure_kPa':ans['mass']['max_dynamic_pressure_kPa'],
       'interpretation':'depth-averaged dynamic-pressure proxy, not actual structural impact load',
       'plot_vmax_percentile_98_kPa':high})

def _plot_2d(ctx,ans,transform,scenario_label='ONE_WAY_1D_TO_2D_SCENARIO'):
    plot_local_2d(ctx,ans,transform,scenario_label)


def _sensitivity(ctx,profile):
    """Competition-grade OAT sweep: runout and depth are defined EVEN if no downstream arrival.

    Unlike v1.0, the figure never uses downstream arrival as its sole y-data.
    A second 2D heatmap maps reachability and predicted arrival under volume/friction scenarios.
    """
    cfg=ctx.cfg['q4'];out=ctx.folder('q4')
    baseline=dict(mu=cfg['mu'],xi=cfg['xi_m_s2'],width=cfg['width_m'],
                  volume=cfg['release_volume_m3'],duration=cfg['release_duration_s'],
                  total_seconds=cfg['total_seconds'],cfl=cfg['cfl'],frames=3,
                  initial_velocity_m_s=cfg.get('initial_velocity_m_s',0.0))
    sweep=cfg.get('scenario_grid',{})
    variations=[('mu',sweep.get('mus',[baseline['mu']*.7,baseline['mu'],baseline['mu']*1.5])),
                ('xi',sweep.get('xi_m_s2',[baseline['xi']*.6,baseline['xi'],baseline['xi']*1.4])),
                ('volume',sweep.get('volumes_m3',[baseline['volume']*.6,baseline['volume'],baseline['volume']*1.5])),
                ('width',sweep.get('widths_m',[baseline['width']*.7,baseline['width'],baseline['width']*1.5])),
                ('duration',sweep.get('release_durations_s',[baseline['duration']*.75,baseline['duration'],baseline['duration']*2.])),
                ('initial_velocity_m_s',sweep.get('initial_velocities_m_s',[
                    max(0,baseline['initial_velocity_m_s']*.7),baseline['initial_velocity_m_s'],baseline['initial_velocity_m_s']*1.3]))]
    names={'mu':'库仑摩擦 μ','xi':'湍流阻力 ξ','volume':'初始释放方量 V',
           'width':'沟道等效宽度 W','initial_velocity_m_s':'初始速度 u0','duration':'释放历时 T'}
    rows=[];sim_cache={}
    for key,vals in variations:
        for v in vals:
            v=float(v);kw=baseline.copy();kw[key]=v
            k=tuple(sorted((kk,round(float(vv),7)) for kk,vv in kw.items()))
            if k not in sim_cache:sim_cache[k]=simulate_1d(profile,**kw)
            a=sim_cache[k];m=a['mass']
            import hashlib,json
            scenario_id=hashlib.sha256(json.dumps({kk:round(float(vv),7) for kk,vv in kw.items()},sort_keys=True).encode()).hexdigest()[:12]
            rows.append({'scenario_id':scenario_id,'varied_parameter':key,'parameter_cn':names[key],
                'value':v,'baseline_value':baseline[key],
                'relative_change_percent':100*(v/baseline[key]-1) if baseline[key] else np.nan,
                'runout_distance_km':m['max_runout_m']/1000.,
                'runout_fraction':m['max_runout_fraction_of_route'],
                'arrival_time_min':m['first_arrival_to_downstream_end_s']/60. if m['first_arrival_to_downstream_end_s'] is not None else np.nan,
                'coupling_arrival_time_min':m['first_arrival_to_coupling_sensor_s']/60. if m['first_arrival_to_coupling_sensor_s'] is not None else np.nan,
                'reached_downstream':bool(m['first_arrival_to_downstream_end_s'] is not None),
                'max_depth_m':m['max_depth_over_time_m'],
                'max_velocity_m_s':m['max_speed_over_time_m_s'],
                'relative_mass_residual':m['relative_absolute_mass_residual']})
    df=pd.DataFrame(rows);df.to_csv(out/'Q4_parameter_sensitivity.csv',index=False)
    fig,axs=plt.subplots(1,3,figsize=(14.1,4.8),constrained_layout=True)
    for i,(key,gr) in enumerate(df.groupby('varied_parameter',sort=False)):
        x=gr.relative_change_percent.to_numpy(float)
        axs[0].plot(x,gr.runout_distance_km,marker='o',lw=1.7,color=COLORS[i%len(COLORS)],label=names[key])
        axs[1].plot(x,gr.max_depth_m,marker='o',lw=1.7,color=COLORS[i%len(COLORS)],label=names[key])
        # Preserve missing arrival visually: no-arrival is marked above as a cross, not omitted.
        ok=gr.reached_downstream.astype(bool).to_numpy()
        axs[2].scatter(x[ok],gr.arrival_time_min.to_numpy(float)[ok],color=COLORS[i%len(COLORS)],s=48,marker='o')
        if (~ok).any():
            ymax=max(10.,float(np.nanmax(df.arrival_time_min)) if df.arrival_time_min.notna().any() else 10.)
            axs[2].scatter(x[~ok],np.full(np.sum(~ok),ymax*1.15),color=COLORS[i%len(COLORS)],s=64,marker='x')
    L=float(profile.distance_m.max())/1000
    axs[0].axhline(L,color='#333333',lw=1.2,ls='--',label='下游断面距离')
    axs[0].set(ylabel='最大传播距离 (km)',title='a  最大传播距离',xlabel='参数相对基线变化 (%)')
    axs[1].set(ylabel='全时段最大等效流深 (m)',title='b  最大流深',xlabel='参数相对基线变化 (%)')
    axs[2].set(ylabel='下游首次到达时间 (min)',title='c  下游首次到达时间',xlabel='参数相对基线变化 (%)')
    for ax in axs:ax.axvline(0,color='#8C919A',lw=1,ls=':');ax.grid(True,alpha=.16)
    axs[0].legend(frameon=False,fontsize=7,loc='best')
    fig.suptitle('泥石流传播模型参数敏感性分析',fontsize=13)
    savefig(fig,out/'Q4_07_parameter_sensitivity.png')

    # Two-parameter envelope, useful for identifying reachability threshold rather than guessing a single parameter.
    volumes=np.array(sorted(set(float(x) for x in sweep.get('volumes_m3',[baseline['volume']*.6,baseline['volume'],baseline['volume']*1.4]))))
    mus=np.array(sorted(set(float(x) for x in sweep.get('mus',[baseline['mu']*.6,baseline['mu'],baseline['mu']*1.4]))))
    travel=np.full((len(mus),len(volumes)),np.nan);runout=np.zeros_like(travel)
    extra=[]
    for r,mu in enumerate(mus):
        for col,vol in enumerate(volumes):
            kw=baseline.copy();kw.update(mu=float(mu),volume=float(vol))
            key=tuple(sorted((kk,round(float(vv),7)) for kk,vv in kw.items()))
            if key not in sim_cache:sim_cache[key]=simulate_1d(profile,**kw)
            m=sim_cache[key]['mass'];runout[r,col]=m['max_runout_m']/1000
            if m['first_arrival_to_downstream_end_s'] is not None:travel[r,col]=m['first_arrival_to_downstream_end_s']/60
            extra.append({'mu':mu,'volume_m3':vol,'runout_km':runout[r,col],
                'arrival_time_min':travel[r,col], 'reached':bool(np.isfinite(travel[r,col]))})
    pd.DataFrame(extra).to_csv(out/'Q4_volume_mu_phase_diagram.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(11.4,4.7),constrained_layout=True)
    im=axs[0].imshow(runout,cmap='viridis',origin='lower',aspect='auto',vmin=0,vmax=L)
    im2=axs[1].imshow(np.ma.masked_invalid(travel),cmap='plasma',origin='lower',aspect='auto')
    for r in range(len(mus)):
        for c in range(len(volumes)):
            axs[0].text(c,r,f'{runout[r,c]:.1f}',ha='center',va='center',fontsize=10,color='white' if runout[r,c]<L*.6 else 'black')
            axs[1].text(c,r,f'{travel[r,c]:.1f}' if np.isfinite(travel[r,c]) else '未到达',ha='center',va='center',fontsize=9,
                       color='black' if not np.isfinite(travel[r,c]) else 'white')
    for ax in axs:
        ax.set_xticks(range(len(volumes)),[f'{v/1e6:.0f}' for v in volumes]);ax.set_yticks(range(len(mus)),[f'{x:.3f}' for x in mus])
        ax.set(xlabel='释放量 (百万 m³)',ylabel='库仑摩擦系数 μ')
    axs[0].set_title('最大传播距离 (km)');axs[1].set_title('下游首次到达 (min)')
    fig.colorbar(im,ax=axs[0],label='km',fraction=.044,pad=.04)
    if np.isfinite(travel).any():fig.colorbar(im2,ax=axs[1],label='min',fraction=.044,pad=.04)
    fig.suptitle('物源体积与摩擦系数对传播的影响',fontsize=13)
    savefig(fig,out/'Q4_08_volume_friction_phase_diagram.png')
    return df

def _synthetic_benchmark(ctx):
    # Verifiable numerical benchmark even without route evidence
    s=np.arange(0,3000,50.)
    z=5000-np.tan(np.deg2rad(22))*s
    p=pd.DataFrame({'distance_m':s,'elevation_m':z})
    run=simulate_1d(p,mu=.12,xi=400,width=30,volume=10000,duration=50,total_seconds=550,frames=12)
    _save_1d(ctx,run,'SYNTHETIC_SLOPE_BENCHMARK_NO_REAL_EVENT')
    write_json(ctx.folder('q4')/'Q4_solver_execution_status.json',{'status':'SYNTHETIC_ONLY','reason':'Route not verified; realistic reconstruction prohibited'})
    return run

def execute_q4_solver(ctx,profile):
    if profile is None:
        record_warning('Q4 runs only synthetic slope benchmark. Use --allow-provisional-q4 to show unverified DEM scenario, or provide verified_route.csv.')
        return _synthetic_benchmark(ctx)
    cfg=ctx.cfg['q4'];out=ctx.folder('q4')
    kws=dict(mu=cfg['mu'],xi=cfg['xi_m_s2'],width=cfg['width_m'],volume=cfg['release_volume_m3'],
             duration=cfg['release_duration_s'],total_seconds=cfg['total_seconds'],
             cfl=cfg['cfl'],frames=cfg['output_frames'],
             initial_velocity_m_s=cfg.get('initial_velocity_m_s',0.0))
    # Crosswalk to Q2-derived *ice inventory* scenario. Do not equate this
    # inventory volume to debris-flow volume without a solid/liquid mixture model.
    q2path=ctx.folder('q2')/'Q2_source_release_candidates.csv'
    if q2path.exists():
        q2release=pd.read_csv(q2path)
        bounds=sorted(q2release.candidate_inventory_ice_volume_m3.astype(float).tolist())
        write_json(out/'Q4_source_volume_data_crosswalk.json',{
           'Q2_physical_candidate_ice_inventory_m3_sorted':bounds,
           'Q4_effective_mixture_volume_scenario_m3':float(kws['volume']),
           'scenario_ratio_to_mid_inventory':float(kws['volume']/np.median(bounds)),
           'quantity_type_warning':'ice inventory != actual release != downstream water/solid mixture',
           'Q4_release_parameters_are_assessed_with_sensitivity':True})
    answer=simulate_1d(profile,**kws)
    plot_geographic_reference(ctx,profile)
    _save_1d(ctx,answer,'物理参数参考情景·DEM候选通道')
    write_json(out/'Q4_1d_settings_and_summary.json',{'parameters':kws,'mass':answer['mass'],
       'path_verified':bool(cfg['route_verified']),
       'scenario_only':not bool(cfg['route_verified']),
       'source_ice_volume_is_not_release_volume':True,
       'spatial_discretization':'Audusse hydrostatic reconstruction + Rusanov finite-volume; lake-at-rest tested; steep coarse dry-front limitation recorded',
       'literature_scale_total_avalanche_m3':103000000,
       'literature_scale_mobile_debris_m3':36000000})
    # A 2D local downstream computation is only useful if 1D flux actually enters the coupling area.
    hydro=answer['history']
    has_actual_1d_inflow=(hydro.sensor_discharge_m3_s.max()>=0.01)
    primary_1d_coupled=bool(has_actual_1d_inflow)
    secondary=None
    if not has_actual_1d_inflow:
        # Predetermined upper-volume sensitivity, not a fitted claim of the real event.
        secondary_settings=kws.copy();secondary_settings['volume']=max(cfg['scenario_grid']['volumes_m3'])
        secondary=simulate_1d(profile,**secondary_settings)
        if secondary['history'].sensor_discharge_m3_s.max()>=.01:
            hydro=secondary['history'];has_actual_1d_inflow=True
            write_json(out/'Q4_secondary_coupling_case.json',{'role':'prespecified upper-volume sensitivity, separate from primary',
                 'parameters':secondary_settings,'mass':secondary['mass'],'primary_case_reaches':False})
            hydro.to_csv(out/'Q4_secondary_1d_discharge.csv',index=False)
    if has_actual_1d_inflow:
        sens=profile.iloc[min(answer['sensor_idx'],len(profile)-1)]
        dest=profile.iloc[-1]
        hydro_input=hydro
        scenario_label='ONE_WAY_1D_TO_2D_SCENARIO' if primary_1d_coupled else 'UPPER_VOLUME_ONE_WAY_SENSITIVITY_SCENARIO'
        state='PROVISIONAL_ONE_WAY_COUPLED_SCENARIO' if primary_1d_coupled else 'SENSITIVITY_ONE_WAY_COUPLED_SCENARIO'
    else:
        # Non-arrival is valid, hence 2D must use an independent scenario forcing.
        # Offer an explicit *standalone* 2D numerical test
        # on the same official DEM with an externally specified Gaussian hydrograph.
        record_warning('1D flood did NOT reach sensor. 2D figures are a separate SYNTHETIC-INFLOW TEST, NOT coupled reconstruction.')
        sens=profile.iloc[max(0,len(profile)-35)]
        dest=profile.iloc[-1]
        tgrid=np.arange(0,cfg['total_seconds']+10,10)
        qsynthetic=600*np.exp(-0.5*((tgrid-360)/180)**2) # illustrative 2D solver QA forcing, not inferred from data
        hydro_input=pd.DataFrame({'time_s':tgrid,'sensor_discharge_m3_s':qsynthetic})
        hydro_input.to_csv(out/'Q4_STANDALONE_SYNTHETIC_2d_inflow.csv',index=False)
        scenario_label='STANDALONE_SYNTHETIC_INFLOW_TEST_NOT_CONNECTED'
        state='STANDALONE_SYNTHETIC_HYDROGRAPH_2D_TEST'
    z,tr,crs=local_dem(ctx,[(sens.utm_x,sens.utm_y),(dest.utm_x,dest.utm_y)],buffer_m=500,resolution_m=cfg['grid_step_m'])
    ri,ci=rowcol(tr,float(sens.utm_x),float(sens.utm_y))
    inlet_profile_idx=answer['sensor_idx'] if has_actual_1d_inflow else max(0,len(profile)-35)
    ahead=profile.iloc[min(inlet_profile_idx+3,len(profile)-1)]
    pdir=np.array([ahead.utm_x-sens.utm_x,ahead.utm_y-sens.utm_y])
    pdir=pdir/max(np.linalg.norm(pdir),1e-10)
    hydro_input.to_csv(out/'Q4_actual_2d_inflow_hydrograph.csv',index=False)
    write_json(out/'Q4_2d_inlet_geometry.json',{'sensor_utm_xy':[float(sens.utm_x),float(sens.utm_y)],
       'target_utm_xy':[float(dest.utm_x),float(dest.utm_y)],'flow_direction_xy':pdir.tolist(),
       'coupled_to_primary_1d':primary_1d_coupled})
    try:
        two=simulate_2d(z,cfg['grid_step_m'],(int(ri),int(ci)),hydro_input,
             total_seconds=cfg['total_seconds'],mu=cfg['mu'],xi=cfg['xi_m_s2'],cfl=.22,frames=10,
             inlet_mode='valley_cross_section',flow_direction_xy=pdir,
             density_kg_m3=cfg.get('density_kg_m3',1800.0))
        _plot_2d(ctx,two,tr,scenario_label=scenario_label)
        _plot_dynamic_pressure(ctx,two,tr)
        if cfg.get('run_inlet_ablation',True):
            old=simulate_2d(z,cfg['grid_step_m'],(int(ri),int(ci)),hydro_input,
                 total_seconds=cfg['total_seconds'],mu=cfg['mu'],xi=cfg['xi_m_s2'],
                 cfl=.22,frames=10,inlet_mode='legacy_3x3',flow_direction_xy=pdir)
            plot_2d_inlet_diagnostic(ctx,two,old,tr)
            old_dir=simulate_2d(z,cfg['grid_step_m'],(int(ri),int(ci)),hydro_input,
                 total_seconds=cfg['total_seconds'],mu=cfg['mu'],xi=cfg['xi_m_s2'],
                 cfl=.22,frames=10,inlet_mode='legacy_directional',flow_direction_xy=pdir)
            valley_novel=simulate_2d(z,cfg['grid_step_m'],(int(ri),int(ci)),hydro_input,
                 total_seconds=cfg['total_seconds'],mu=cfg['mu'],xi=cfg['xi_m_s2'],
                 cfl=.22,frames=10,inlet_mode='valley_zero_momentum',flow_direction_xy=pdir)
            four={'侧坡+静止入流':old,'侧坡+沿沟动量':old_dir,
                  '谷底+静止入流':valley_novel,'谷底+沿沟动量':two}
            plot_inlet_factorial(ctx,four,tr)
            pd.DataFrame([{'condition':key,'peak_depth_m':value['mass']['peak_depth_m'],
                'affected_area_km2':value['mass']['wet_area_h_gt_0p2_km2'],
                'mass_error':value['mass']['relative_mass_error'],
                'inlet_offset_cells':value['inlet']['cross_section_offset_cells']}
                for key,value in four.items()]).to_csv(out/'Q4_2d_boundary_factorial.csv',index=False)
            write_json(out/'Q4_2d_inlet_ablation.json',{'legacy':old['mass'],
                 'valley_aligned':two['mass'],'legacy_inlet':old['inlet'],
                 'valley_aligned_inlet':two['inlet'],
                 'purpose':'numerical inlet-location ablation to diagnose artificial bifurcation'})
        write_json(out/'Q4_2d_status.json',{'status':state,'mass':two['mass'],
             'inlet':two['inlet'],'first_inlet_s':two['first_inlet_s'],
             'chosen_frames_s':two['frame_times_s'],
             'inlet_boundary_condition':'valley_cross_section (not source-cell 3x3)',
             '1d_inflow_reached_coupling':bool(has_actual_1d_inflow),'primary_1d_coupled':primary_1d_coupled,
             'coupled_1d_volume_m3':float(kws['volume']) if primary_1d_coupled else (float(secondary_settings['volume']) if has_actual_1d_inflow else None),'verified_route':bool(cfg['route_verified']),
             'assumptions':'shallow effective continuum, stated DEM inlet; input momentum from 1D only if coupled',
             'inflow_source':('primary 1D hydrograph' if primary_1d_coupled else 'prespecified upper-volume 1D sensitivity hydrograph') if has_actual_1d_inflow else 'illustrative synthetic Gaussian inflow: no calibration, not a real observation'})
    except Exception as e:
        record_warning(f'2D solver did not converge: {e}. No image may be claimed as a result.')
        write_json(out/'Q4_2d_status.json',{'status':'NUMERICAL_FAILURE','error':str(e)})
    df=_sensitivity(ctx,profile)
    print('[Q4] 1D mass residual',answer['mass']['relative_absolute_mass_residual'],
          'downstream arrival_s',answer['mass']['first_arrival_to_downstream_end_s'])
    return answer
