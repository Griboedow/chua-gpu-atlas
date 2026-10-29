"""Search small noiseless synchronematic clusters at published coefficients."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from synchronematic import SynchronematicGPU,geometry,flow_and_jac,COEFFICIENTS,REPOSITORY_NUMPY_COEFFICIENTS

def save(path,result,parents=None):
    path.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path/'scan.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    valid=result['status']==0;windows=result['lle_windows']
    mask=valid&(result['lle']>.002)&(windows.min(axis=1)>.001)&(result['stats'][:,-1]>1e-5)
    summary={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    summary.update(trajectories=len(valid),valid=int(valid.sum()),candidate_indices=np.flatnonzero(mask).tolist(),max_lle=float(np.nanmax(result['lle'])),
                   criterion='mean > .002 and all four windows > .001; late flow speed > 1e-5; status 0',
                   interpretation='Finite-time candidates in the specified deterministic Fourier model; actual coefficients are retained in the arrays. Not proof of novelty or physical chaos.')
    if parents is not None:summary['parent_indices']=np.asarray(parents).tolist()
    (path/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({'output':str(path),'trajectories':len(valid),'retained':int(mask.sum()),'max_lle':summary['max_lle'],'timing':result['timing']}),flush=True)

def validate():
    p=np.array([[0.,0.],[1.,0.],[.4,.9]]);p/=min(np.linalg.norm(p[i]-p[j]) for i in range(3) for j in range(i))
    pairs=geometry(p);s=np.random.default_rng(93).uniform(-2,2,6)
    flow,jac=flow_and_jac(s,pairs);epsilon=1e-6;eye=np.eye(6)
    numerical=np.column_stack([(flow_and_jac(s+epsilon*e,pairs)[0]-flow_and_jac(s-epsilon*e,pairs)[0])/(2*epsilon) for e in eye])
    reference=solve_ivp(lambda t,x:flow_and_jac(x,pairs)[0],(0,1),s,method='DOP853',rtol=1e-12,atol=1e-14)
    gpu=SynchronematicGPU(3);errors=[]
    for dt in [.01,.005,.0025]:
        a=gpu.scan(p[None],s[None],dt=dt,transient=0,duration=1)
        difference=np.angle(np.exp(1j*(a['final_states'][0]-reference.y[:,-1])))
        errors.append(float(np.max(abs(difference))))
    r={'jacobian_error':float(np.max(abs(jac-numerical))),'state_errors_vs_DOP853':errors,'dt':[.01,.005,.0025],
       'global_phase_null_error':float(np.max(abs(jac@np.r_[np.ones(3),np.zeros(3)]))),
       'backend':gpu.info,'passed':bool(np.max(abs(jac-numerical))<1e-7 and errors[-1]<1e-6 and errors[-1]<errors[0]/100)}
    path=Path('results/synchronematic');path.mkdir(exist_ok=True,parents=True);(path/'validation.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2),flush=True)
    if not r['passed']:raise RuntimeError('Validation failed')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,default=3);parser.add_argument('--geometries',type=int,default=128)
    parser.add_argument('--seeds',type=int,default=8);parser.add_argument('--dt',type=float,default=.01)
    parser.add_argument('--transient',type=float,default=200);parser.add_argument('--duration',type=float,default=400)
    parser.add_argument('--refine',type=Path);parser.add_argument('--top',type=int,default=0);parser.add_argument('--output',type=Path)
    parser.add_argument('--variant',choices=['article','repository'],default='article')
    parser.add_argument('--validate',action='store_true');args=parser.parse_args()
    if args.validate:validate();raise SystemExit()
    parents=None
    if args.refine:
        a=np.load(args.refine/'scan.npz');meta=json.loads((args.refine/'summary.json').read_text())
        parents=np.array(meta['candidate_indices'],int)
        if args.top:parents=np.unique(np.r_[parents,np.argsort(a['lle'])[-args.top:]])
        positions=a['positions'][parents];starts=a['starts'][parents];args.n=positions.shape[1]
        coefficients=a['coefficients'][parents] if 'coefficients' in a else COEFFICIENTS
        if not len(starts):raise SystemExit('No selected trajectories')
    else:
        coefficients=COEFFICIENTS if args.variant=='article' else REPOSITORY_NUMPY_COEFFICIENTS
        rng=np.random.default_rng(20261029+args.n);positions=[]
        for k in range(args.geometries):
            if k==0:
                angle=np.arange(args.n)*2*np.pi/args.n;p=np.column_stack([np.cos(angle),np.sin(angle)])
            else:
                while True:
                    p=rng.uniform(-1,1,(args.n,2));minimum=min(np.linalg.norm(p[i]-p[j]) for i in range(args.n) for j in range(i))
                    if minimum>.25:break
            minimum=min(np.linalg.norm(p[i]-p[j]) for i in range(args.n) for j in range(i));positions.extend([p/minimum]*args.seeds)
        positions=np.array(positions);starts=rng.uniform(-np.pi,np.pi,(len(positions),2*args.n))
    gpu=SynchronematicGPU(args.n);print('Starting',len(starts),'clusters, N =',args.n,flush=True)
    a=gpu.scan(positions,starts,dt=args.dt,transient=args.transient,duration=args.duration,coefficients=coefficients,
               progress=lambda f,s:print(f'{f:.1%}, {s:.1f}s',flush=True))
    a['coefficient_source']='inherited from saved ensemble' if args.refine else args.variant
    suffix='refined' if args.refine else ('scout' if args.variant=='article' else 'repository-scout')
    save(args.output or Path(f'results/synchronematic/n{args.n}-{suffix}'),a,parents)
