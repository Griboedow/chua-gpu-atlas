"""Plot the saved source-control mismatch; no integrations or fitted parameters."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from microfluidic import ROOT

data=np.load(ROOT/'independent-source-control.npz')
t=.25*data['times'];states=data['dop853_states'];other=data['independent_radau_states']
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for ax,col,label in [(axs[0,0],2,'Transverse flow Q3 (µL/s per mm depth)'),
                      (axs[0,1],3,'Vortex size (µm)')]:
    scale=1 if col==2 else 500/12
    ax.plot(t,states[:,col]*scale,label='Reduced ODE / DOP853',lw=1.7)
    ax.plot(t[::130],other[::130,col]*scale,'o',ms=3,label='Five-flow reconstruction / Radau')
    ax.set(xlabel='Physical time (s)',ylabel=label,xlim=(0,.5))
    ax.grid(alpha=.2)
axs[0,0].legend(fontsize=8)
flow=np.linspace(0,4500,1000);xi=(5.5**3*(.12+8e-6*flow)-25)/5.5
axs[1,0].plot(flow,xi,color='#8e44ad')
axs[1,0].axhline(0,color='black',lw=.8)
axs[1,0].axvspan(0,900,color='#f1c40f',alpha=.2,label='Q1 ≤ plotted total flow (<900)')
axs[1,0].scatter([states[-1,0]],[(5.5**3*(.12+8e-6*states[-1,0])-25)/5.5],color='#c0392b',label='Our 250 Pa equilibrium')
axs[1,0].set(xlabel='Inlet flow Q1 (µL/s per mm depth)',ylabel='ξ from printed S1, S2 and S13')
axs[1,0].legend(fontsize=8);axs[1,0].grid(alpha=.2)
eig=data['equilibrium_eigenvalues']
axs[1,1].scatter(eig.real,eig.imag,s=55,color='#c0392b')
axs[1,1].axvline(0,color='black',lw=.8)
axs[1,1].set(xlabel='Real part (per source time)',ylabel='Imaginary part (per source time)',
             title='Four eigenvalues of the stable equilibrium',xlim=(-60,5))
axs[1,1].grid(alpha=.2)
fig.suptitle('Literal 2020 microfluidic closure: published periodic control not recovered\n'
             '250 Pa at each inlet; all fitted coefficients unchanged',fontsize=13)
path=Path(__file__).resolve().parents[1]/'docs/assets/microfluidic-controls.png'
fig.savefig(path,dpi=180);print(path)
