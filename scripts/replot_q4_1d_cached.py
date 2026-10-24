"""Regenerate Q4 one-dimensional runout plots from cached finite-volume output.
No physical parameters or numerical solution are changed here.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import configure_plotting,savefig,COLORS
configure_plotting()
q=ROOT/'outputs'/'q4'
a=np.load(q/'Q4_1d_frames_arrays.npz')
s=a['distance_m'];z=a['elevation_m'];frames=a['time_s'];depth=a['depth_m']
fig,axes=plt.subplots(2,1,figsize=(10.4,7.0),gridspec_kw={'height_ratios':[1,1.15]})
axes[0].plot(s/1000,z,color='#63767F',lw=1.6)
axes[0].set(ylabel='高程 (m)',title='a  传播路径地形纵剖面')
axes[0].grid(axis='y',alpha=.16)
im=axes[1].imshow(depth,aspect='auto',extent=[0,s[-1]/1000,frames[-1]/60,frames[0]/60],
                  cmap='turbo',vmin=0,vmax=max(np.percentile(depth,99),.1))
axes[1].invert_yaxis()
axes[1].set(xlabel='沿程距离 (km)',ylabel='源区启动后时间 (min)',title='b  一维流深时空分布')
cb=fig.colorbar(im,ax=axes[1],label='流深 (m)',shrink=.84)
fig.suptitle('一维泥石流传播及地形纵剖面',fontsize=13,y=.965)
fig.subplots_adjust(left=.10,right=.94,top=.87,bottom=.10,hspace=.32)
savefig(fig,q/'Q4_02_1d_runout_timespace.png')
print('Regenerated cached Q4_02 figure (PNG and SVG).')
