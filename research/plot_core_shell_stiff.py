"""Draw the saved implicit periodic-orbit control, without new integration."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core_shell_gel import ROOT

data=np.load(ROOT/'stiff-followup-2882.npz')
t=data['cycle_times'];s=data['cycle_states']
fig=plt.figure(figsize=(11,4.5),layout='constrained')
ax=fig.add_subplot(121)
ax.plot(t,s[:,2],label='Proton ratio x',color='#1565c0')
ax.plot(t,s[:,3],label='Polymer fraction φ',color='#d35400')
ax.set(xlabel='Scaled time within one period',ylabel='State value',title='Independent BDF orbit closure')
ax.legend();ax.grid(alpha=.2)
bx=fig.add_subplot(122,projection='3d')
bx.plot(s[:,0],s[:,2],s[:,3],color='#1565c0',lw=1.2)
bx.scatter(*s[0,[0,2,3]],color='#d35400',s=30,label='Section φ = 0.3')
bx.set(xlabel='Scaled glucose G',ylabel='Proton ratio x',zlabel='Polymer fraction φ',
       title='Projection of the four-state stable cycle')
bx.view_init(elev=25,azim=-55)
fig.suptitle('Core–shell gel: previously unresolved endpoint 2882\n'
             'Period 0.8998653920; dominant transverse multiplier 1.71 × 10⁻⁷',fontsize=12)
path=Path(__file__).resolve().parents[1]/'docs/assets/core-shell-stiff-cycle.png'
fig.savefig(path,dpi=180);print(path)
