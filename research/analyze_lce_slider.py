"""Local parameter sensitivity and original figures for the slider candidate."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from lce_slider import SliderGPU, affine_flow, next_crossing


ROOT=Path('results/lce-slider')


def robustness():
    data=np.load(ROOT/'wide/scan.npz');p=data['parameters'][2222]
    q=p[2]/(1+p[2]);rng=np.random.default_rng(2026102922)
    parameters=[];starts=[];labels=[]
    for level in [0.,.001,.01,.03,.1]:
        for sample in range(128):
            candidate=p*rng.uniform(1-level,1+level,4)
            qq=candidate[2]/(1+candidate[2])
            for seed in range(4):
                parameters.append(candidate)
                starts.append([np.sqrt(candidate[0])*qq*10**rng.uniform(-1,1),qq*rng.uniform(.01,.99)])
                labels.append([level,sample,seed])
    result=SliderGPU().scan(np.array(parameters),np.array(starts),transient=2000,cycles=8000,subdivisions=256,
                           progress=lambda f,t:print(f'{f:.1%}, {t:.1f}s',flush=True))
    labels=np.array(labels);result['labels']=labels
    path=ROOT/'local';path.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path/'scan.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    good=(result['status']==0)&(result['lle']>1e-4)&(result['lle_windows'].min(axis=1)>1e-5)
    groups=[]
    for level in np.unique(labels[:,0]):
        mask=labels[:,0]==level
        groups.append(dict(relative_jitter=float(level),trajectories=int(mask.sum()),retained=int(good[mask].sum()),
                           valid=int((result['status'][mask]==0).sum()),median_lle=float(np.nanmedian(result['lle'][mask]))))
    meta={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    meta.update(seed=2026102922,base_parameters=p.tolist(),groups=groups,
                perturbation='Independent uniform relative jitter of k, b, I, d; four starts per sampled tuple',
                caveat='Finite sampled ensembles, not confidence bounds or a proof of an open chaotic region')
    (path/'summary.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(groups),flush=True)


def figures():
    data=np.load(ROOT/'analytic-2222/trace.npz');p=data['parameters']
    section=data['section'];s=section[-1,1:3].copy();dense=[];clock=0.
    for _ in range(100):
        for on in [True,False]:
            end,dt,j=next_crossing(s,p,on,256)
            t=np.linspace(0,dt,65,endpoint=False)
            states=affine_flow(np.r_[0.,s],p,on,t).T
            dense.extend(np.column_stack((clock+t,states,np.full(len(t),on))))
            clock+=dt;s=end
    dense=np.array(dense)
    # C and g only specify an illustrative back-transformation, not a calibration.
    C=.2;g=1.2;h=g/p[0]-C*p[3]
    x=h+C*dense[:,1];elastic=x+C*dense[:,3]
    geometry=dict(C=C,g=g,h=float(h),minimum_x=float(x.min()),maximum_x=float(x.max()),
                  minimum_elastic_strain=float(elastic.min()),maximum_elastic_strain=float(elastic.max()),
                  min_phi=float(dense[:,3].min()),max_phi=float(dense[:,3].max()),
                  caveat='Illustrative nondimensional reconstruction on 100 late cycles; sampled extrema, no physical calibration')
    path=ROOT/'geometry-check';path.mkdir(exist_ok=True)
    (path/'summary.json').write_text(json.dumps(geometry,indent=2))
    np.savez_compressed(path/'trace.npz',trajectory=dense,parameters=p)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(13,9),layout='constrained')
    fig.suptitle('Light-switching slider: chaos candidate beyond the original numerical parameter range',fontsize=16)
    ax=fig.add_subplot(221,projection='3d')
    ax.plot(dense[:4000,1]*1e3,dense[:4000,2],dense[:4000,3]*1e3,lw=.55,color='#157f9b',alpha=.8)
    ax.set(xlabel='Displacement u × 1000',ylabel='Velocity v',zlabel='Cis fraction φ × 1000',title='Published three-state equations')
    ax=fig.add_subplot(222)
    ax.scatter(section[:-1,1],section[1:,1],s=2,alpha=.35,color='#e26936')
    ax.set(xlabel='Upward-crossing speed vₙ',ylabel='Next crossing speed vₙ₊₁',title='Event-defined return map (8,000 cycles)')
    ax=fig.add_subplot(223)
    for name,label,color in [('cpu-2222','Independent DOP853','#a55b9d'),('analytic-2222','Analytic segments','#157f9b')]:
        trace=np.load(ROOT/name/'trace.npz');hist=trace['history']
        ax.plot(hist[:,0],hist[:,1],label=label,color=color,lw=1.4)
    ax.axhline(0,color='black',lw=.7);ax.legend()
    ax.set(xlabel='Elapsed time after transient',ylabel='Largest transverse Lyapunov exponent',title='Positive stretching with both solvers')
    ax=fig.add_subplot(224);local=json.loads((ROOT/'local/summary.json').read_text())['groups']
    ax.bar([str(round(100*r['relative_jitter'],2))+'%' for r in local],[r['retained']/r['trajectories'] for r in local],color='#269b72')
    ax.set(xlabel='Independent relative parameter jitter',ylabel='Retained fraction of sampled starts',ylim=(0,1.08),title='Finite ensemble sensitivity, 512 starts per level')
    for i,r in enumerate(local):ax.text(i,r['retained']/r['trajectories']+.02,str(r['retained'])+'/512',ha='center',fontsize=9)
    fig.savefig('docs/assets/lce_slider_evidence.png',dpi=180);plt.close(fig)
    print(json.dumps(geometry,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['all','robustness','figures'],default='all')
    args=parser.parse_args()
    if args.mode in ['all','robustness']:robustness()
    if args.mode in ['all','figures']:figures()
