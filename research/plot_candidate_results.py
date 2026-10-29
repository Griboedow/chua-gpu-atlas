"""Original figures from this project's numerical results, not source-paper art."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if __name__=='__main__':
    root=Path('results/candidate-research');out=Path('docs/assets');out.mkdir(parents=True,exist_ok=True)
    a=np.load(root/'scout/scan.npz');b=np.load(root/'refined/scan.npz')
    ids=np.load(root/'refined/parent_indices.npy');control=np.load(root/'periodic_reference.npz')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(12,8),layout='constrained')
    ax=fig.add_subplot(2,2,1);mask=control['t']<8
    ax.plot(control['t'][mask],control['state'][3,mask],color='#087e8b',lw=1.4)
    ax.set(xlabel='Time after the reference run',ylabel='Voltage V',title='A. Recovered periodic control')
    ax=fig.add_subplot(2,2,2,projection='3d')
    xyz=control['state'][:3,:2000];ax.plot(*xyz,color='#087e8b',lw=1.1)
    ax.set(xlabel='X',ylabel='Y',zlabel='Z',title='B. Closed orbit in Bloch coordinates')
    ax=fig.add_subplot(2,2,3);xx=np.arange(1,len(ids)+1)
    ax.vlines(xx,b['lle'],a['lle'][ids],color='#bfc7cc',lw=2)
    ax.scatter(xx,a['lle'][ids],label='Observation 400; transient 200',color='#cc6133',zorder=3)
    ax.scatter(xx,b['lle'],label='Observation 8,000; transient 5,000',color='#087e8b',marker='s',zorder=3)
    ax.axhline(0,color='black',lw=.7);ax.set(xlabel='First-survey candidate',ylabel='Mean finite-time largest exponent',title='C. Ten apparent candidates fail the long check',xticks=xx)
    ax.set_ylim(-.0005,.014)
    ax.legend(fontsize=8,loc='upper left')
    ax=fig.add_subplot(2,2,4)
    distributions={name:np.load(root/name/'scan.npz')['lle'].reshape(-1,3).max(axis=1) for name in ['scout','geometry-scout']}
    bins=np.linspace(min(v.min() for v in distributions.values()),.03,100)
    for name,label,color in [('scout','Fixed x0: 4,096 tuples','#087e8b'),('geometry-scout','Varied x0: 8,192 tuples','#cc6133')]:
        data=np.load(root/name/'scan.npz');lle=data['lle'].reshape(-1,3).max(axis=1)
        ax.hist(lle,bins=bins,histtype='step',lw=1.7,label=label,color=color)
    ax.set(xlabel='Largest mean estimate among three starts',ylabel='Tuples per bin (log scale)',yscale='log',title='D. Finite-time estimates, not chaos labels')
    ax.axvline(0,color='black',lw=.7);ax.legend(fontsize=8)
    fig.suptitle('Beyond Chua: quantum-neuron search and numerical controls',fontsize=16)
    fig.savefig(out/'candidate_research.png',dpi=170);plt.close(fig)
    print(out/'candidate_research.png')
