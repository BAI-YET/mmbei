"""Well-balanced MUSCL-HLL / SSP-RK2 update on measured variable sections.

The old first-order Rusanov method remains an explicit comparison. Reconstruction
uses water-surface elevation and velocity; both preserve a stationary lake.
"""
from __future__ import annotations
import numpy as np
from .q4_solver import G, EPS


def limited_slope(values, theta=1.3):
    x=np.asarray(values,float);s=np.zeros_like(x)
    left=x[1:-1]-x[:-2];right=x[2:]-x[1:-1]
    centre=.5*(left+right)
    same=(left*right)>0
    s[1:-1]=np.where(same,np.sign(centre)*np.minimum.reduce(
        [theta*np.abs(left),np.abs(centre),theta*np.abs(right)]),0.)
    return s


def muscl_euler(A,Q,bed,width,dx,dt,mu,xi,qin=0.,vin=0.,bank_side_slope=None):
    from .q4_cascade_reconstruction import section_depth, section_pressure
    side=np.zeros_like(width) if bank_side_slope is None else bank_side_slope
    h=section_depth(A,width,side);u=Q/np.maximum(A,EPS);n=len(A)
    se=limited_slope(bed+h);su=limited_slope(u)
    wet=h>1e-6;interior=wet.copy()
    interior[1:]&=wet[:-1];interior[:-1]&=wet[1:]
    se=np.where(interior,np.clip(se,-2*h,2*h),0.)
    su=np.where(interior,su,0.)
    hp=np.maximum(h+.5*se,0.);hm=np.maximum(h-.5*se,0.)
    crest=np.maximum(bed[:-1],bed[1:])
    hl=np.maximum(bed[:-1]+hp[:-1]-crest,0.)
    hr=np.maximum(bed[1:]+hm[1:]-crest,0.)
    b=np.minimum(width[:-1],width[1:]);m=np.minimum(side[:-1],side[1:])
    al=b*hl+m*hl*hl;ar=b*hr+m*hr*hr
    ul=np.where(hl>1e-8,u[:-1]+.5*su[:-1],0.)
    ur=np.where(hr>1e-8,u[1:]-.5*su[1:],0.)
    ql=al*ul;qr=ar*ur
    pl=section_pressure(hl,b,m);pr=section_pressure(hr,b,m)
    fl=np.stack([ql,ql*ul+pl]);fr=np.stack([qr,qr*ur+pr])
    cl=np.sqrt(G*al/np.maximum(b+2*m*hl,EPS))
    cr=np.sqrt(G*ar/np.maximum(b+2*m*hr,EPS))
    sl=np.minimum(np.minimum(ul-cl,ur-cr),0.)
    sr=np.maximum(np.maximum(ul+cl,ur+cr),0.)
    face=(sr[None,:]*fl-sl[None,:]*fr+(sl*sr)[None,:]*
          np.stack([ar-al,qr-ql]))/np.maximum(sr-sl,EPS)[None,:]
    outgoing=np.zeros((2,n));incoming=np.zeros((2,n))
    outgoing[:,:-1]=face;incoming[:,1:]=face
    # Use reconstructed OWN-face pressures. Using the centre pressure on both
    # sides would cancel real pressure gradients in second-order reconstruction.
    outgoing[1,:-1]+=section_pressure(hp[:-1],width[:-1],side[:-1])-pl
    incoming[1,1:]+=section_pressure(hm[1:],width[1:],side[1:])-pr
    pressure=section_pressure(h,width,side)
    qout=max(float(Q[-1]),0.)
    outgoing[:,-1]=[qout,qout*u[-1]+pressure[-1]]
    incoming[1,0]=pressure[0]
    if qin>0:incoming[:,0]=[qin,qin*vin+pressure[0]]
    updated=np.stack([A,Q])-dt/dx*(outgoing-incoming)
    if np.min(updated[0])<-1e-6:raise FloatingPointError('Negative MUSCL channel area; reduce CFL')
    clip=float(-np.minimum(updated[0],0.).sum()*dx)
    area=np.maximum(updated[0],0.);velocity=updated[1]/np.maximum(area,EPS)
    hn=section_depth(area,width,side);cos=np.cos(np.arctan(np.gradient(bed,dx)))
    speed=np.maximum(abs(velocity)-dt*mu*G*cos,0.)
    depth=area/np.maximum(width+2*side*hn,EPS)
    drag=dt*G/(xi*np.maximum(depth,.01))
    speed=2*speed/(1+np.sqrt(1+4*drag*speed))
    discharge=area*np.sign(velocity)*speed;discharge[hn<1e-6]=0.
    return area,discharge,float(qout*dt),clip,face[0].copy()


def advance_channel_muscl(A,Q,bed,width,dx,dt,mu,xi,qin=0.,vin=0.,bank_side_slope=None,
                         return_interface_flux=False):
    args=(bed,width,dx,dt,mu,xi,qin,vin,bank_side_slope)
    a1,q1,v1,c1,f1=muscl_euler(A,Q,*args)
    a2,q2,v2,c2,f2=muscl_euler(a1,q1,*args)
    result=(.5*(A+a2),.5*(Q+q2),.5*(v1+v2),.5*(c1+c2))
    return result+(.5*(f1+f2),) if return_interface_flux else result
