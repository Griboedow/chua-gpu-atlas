"""Plot saved saddle-departure evidence; no numerical integration."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'results/light-operators/saddle-departure'
ASSETS = ROOT/'docs/assets'


def main():
    report=json.loads((DATA/'probe.json').read_text())
    rows=report['rows']; colors={-1:'#c14634',1:'#176c9c'}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(15,4.7),layout='constrained')
    ax=fig.add_subplot(131); ar=fig.add_subplot(132); ap=fig.add_subplot(133,projection='3d')
    for row in rows:
        d=np.load(DATA/row['array_file']); color=colors[row['sign']]
        label=f"{'Irregular side' if row['sign']<0 else 'Cycle side'}, epsilon={row['amplitude']:g}"
        ax.semilogy(d['section_times'],np.maximum(d['section_distance_to_stable'],1e-14),
                    color=color,alpha=.7,lw=1,ls='-' if row['amplitude']==1e-4 else '--',label=label)
        mask=d['section_times']>100; ct=d['section_times'][mask]; periods=np.diff(ct)
        ar.scatter(periods[:-1],periods[1:],s=8,color=color,alpha=.4)
        if row['amplitude']==1e-5 and row['sign']<0:
            y=d['states'][d['times']>300]
            ap.scatter(y[:,0],y[:,3],y[:,2],color=color,s=3,alpha=.45,label='Irregular samples')
    stable=np.load(ROOT/'results/light-operators/cycle-arclength/control-initial.npz')['trajectory']
    ap.plot(stable[:,0],stable[:,3],stable[:,2],color=colors[1],lw=1.7,label='Periodic cycle')
    ax.set(xlabel='Time from saddle perturbation',ylabel='Distance to stable cycle on section',
           title='Two signs, two perturbation sizes')
    ax.legend(fontsize=8,loc='lower left'); ax.grid(alpha=.2)
    ar.set(xlabel='Return interval n',ylabel='Return interval n+1',title='Returns after time 100')
    ar.grid(alpha=.2)
    ap.set(xlabel='Displacement minus',ylabel='Displacement plus',zlabel='Temperature minus',
           title='Physical projection of six-state flow')
    ap.view_init(elev=23,azim=35); ap.legend(fontsize=8)
    fig.suptitle('Two observed departure outcomes near a saddle cycle at optical width 0.04\n'
                 'Finite-time evidence only; assumed smooth illumination law',fontsize=13)
    fig.savefig(ASSETS/'optical_saddle_departure.png',dpi=170)
    plt.close(fig)
    fig=plt.figure(figsize=(7.2,6),layout='constrained')
    ax=fig.add_subplot(111,projection='3d')
    for row in rows:
        if row['amplitude']!=1e-5 or row['sign']>0: continue
        d=np.load(DATA/row['array_file']);y=d['states'][d['times']>300]
        ax.scatter(y[:,0],y[:,3],y[:,2],color=colors[row['sign']],s=4,alpha=.5,
                label='Irregular finite-time trajectory samples')
    ax.plot(stable[:,0],stable[:,3],stable[:,2],color=colors[1],lw=1.7,label='Stable periodic cycle')
    ax.set(xlabel='Displacement minus',ylabel='Displacement plus',zlabel='Temperature minus')
    ax.legend(loc='upper left',fontsize=8)
    fig.suptitle('Two departures from the same optical saddle cycle\n'
                 'Three physical observables of a six-state model',fontsize=12)
    def frame(i):
        angle=2*np.pi*i/72
        ax.view_init(elev=22+14*np.sin(angle),azim=360*i/72,roll=6*np.sin(2*angle))
        return []
    movie=FuncAnimation(fig,frame,frames=72,interval=85)
    movie.save(ASSETS/'optical_saddle_departure_rotation.gif',writer=PillowWriter(fps=12),dpi=95)
    plt.close(fig)


if __name__=='__main__':main()
