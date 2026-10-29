"""GPU return-map survey of the published light-switching slider model."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from lce_slider import SliderGPU


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--power',type=int,default=12)
    parser.add_argument('--seeds',type=int,default=2)
    parser.add_argument('--wide',action='store_true')
    parser.add_argument('--refine',type=Path)
    parser.add_argument('--top',type=int,default=16)
    parser.add_argument('--transient',type=int,default=200)
    parser.add_argument('--cycles',type=int,default=800)
    parser.add_argument('--subdivisions',type=int,default=64)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    parent=None
    if args.refine:
        data=np.load(args.refine/'scan.npz')
        meta=json.loads((args.refine/'summary.json').read_text())
        valid=np.flatnonzero(data['status']==0)
        if not len(valid):raise SystemExit('No valid trajectories to refine')
        top=valid[np.argsort(data['lle'][valid])[-args.top:]] if args.top else np.array([],int)
        parent=np.unique(np.r_[meta['candidate_indices'],top]).astype(int)
        p,starts=data['parameters'][parent],data['starts'][parent]
        design=dict(parent=str(args.refine),parent_indices=parent.tolist())
    else:
        sobol=qmc.Sobol(4,scramble=True,seed=2026102920).random_base2(args.power)
        if args.wide:
            k=10**(-2.3+5.3*sobol[:,0]);b=2*np.sqrt(k)*10**(-3+2.875*sobol[:,1])
        else:
            k=10**(-1+2*sobol[:,0]);b=10**(np.log10(.002)+np.log10(.5/.002)*sobol[:,1])
        intensity=10**(np.log10(.005)+2*sobol[:,2]);q=intensity/(1+intensity)
        d=q*(.005+.99*sobol[:,3])
        p=np.repeat(np.column_stack((k,b,intensity,d)),args.seeds,axis=0)
        rng=np.random.default_rng(2026102921);qq=p[:,2]/(1+p[:,2])
        starts=np.column_stack((np.sqrt(p[:,0])*qq*10**rng.uniform(-1,1,len(p)),qq*rng.uniform(.01,.99,len(p))))
        p[0]=[5.8,.1,.1,(1.2/5.8-.2)/.2];starts[0]=[.3,.03]
        design=dict(seed=2026102920,samples=2**args.power,seeds=args.seeds,wide=args.wide,
                    k_range=[10**-2.3,1000] if args.wide else [.1,10],
                    damping_range='damping ratio .001 to .75' if args.wide else [.002,.5],
                    intensity_range=[.005,.5],threshold_fraction_range=[.005,.995],
                    scope='Expanded mathematical parameters' if args.wide else 'k, damping, intensity within the paper table; threshold chosen between the two branch equilibria')
    print(f'Starting {len(p)} slider trajectories',flush=True)
    result=SliderGPU().scan(p,starts,args.transient,args.cycles,args.subdivisions,
                            progress=lambda f,t:print(f'{f:.1%}, {t:.1f}s',flush=True))
    path=args.output or Path('results/lce-slider')/('refined' if args.refine else ('wide' if args.wide else 'scout'))
    path.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path/'scan.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    good=result['status']==0
    candidates=good&(result['lle']>1e-4)&(np.min(result['lle_windows'],axis=1)>1e-5)
    summary={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    codes,counts=np.unique(result['status'],return_counts=True)
    summary.update(design=design,trajectories=len(p),status_counts=dict(zip(map(str,codes),map(int,counts))),
                   candidate_indices=np.flatnonzero(candidates).tolist(),
                   max_lle=float(np.nanmax(result['lle'])) if good.any() else None,
                   criterion='Transverse return-map exponent per elapsed time > 1e-4, all four windows > 1e-5, status 0',
                   interpretation='Finite-time event-map screening; a periodic orbit has a negative transverse exponent, with the neutral flow exponent omitted.',
                   status_key={'0':'valid','1':'no return within bracketing horizon','2':'near grazing / invalid crossing direction','3':'nonfinite segment','4':'degenerate tangent'})
    (path/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(dict(output=str(path),status_counts=summary['status_counts'],retained=int(candidates.sum()),
                         max_lle=summary['max_lle'],wall_seconds=result['wall_seconds'])),flush=True)


if __name__=='__main__':main()
