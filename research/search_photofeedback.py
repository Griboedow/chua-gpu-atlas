"""Mechanism exploration, explicitly not hardware reproduction or novelty proof."""
from pathlib import Path
import argparse,json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import qmc
from photofeedback import PhotoFeedbackGPU,parameter_at_equilibrium,stability,rhs

def validate(gpu):
    cases=[parameter_at_equilibrium(.8,.18,.8,2.36,.2),parameter_at_equilibrium(1.,.18,.02,0.,.2)]
    params=np.array([v[0] for v in cases]);starts=np.array([v[1]+[.01,.01,0] for v in cases])
    p,s=cases[0];tan=np.array([.3,.2,-.1]);eps=1e-6
    err=float(np.max(abs(rhs(0,np.r_[s,tan],p)[3:]-(rhs(0,s+eps*tan,p)-rhs(0,s-eps*tan,p))/(2*eps))))
    result=gpu.scan(params,starts,transient=0,duration=4,rtol=1e-10,atol=1e-12)
    diffs=[]
    for p,s,end in zip(params,starts,result['final_states']):
        sol=solve_ivp(lambda t,x:rhs(t,x,p),(0,4),s,method='DOP853',rtol=1e-12,atol=1e-14)
        if not sol.success:raise RuntimeError(sol.message)
        diffs.append(float(np.max(abs(sol.y[:,-1]-end))))
    r={'jacobian_error':err,'gpu_vs_cpu_state_errors':diffs,'status':result['status'].tolist(),'passed':bool(err<1e-7 and max(diffs)<1e-7 and (result['status']==0).all()),'model':'Phenomenological photoresponse extension, not calibrated to the experiment'}
    path=Path('results/photofeedback');path.mkdir(exist_ok=True,parents=True);(path/'validation.json').write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
    if not r['passed']:raise RuntimeError('Validation failed')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--validate',action='store_true');parser.add_argument('--power',type=int,default=12)
    parser.add_argument('--refine',type=Path);parser.add_argument('--output',type=Path,default=Path('results/photofeedback/scout'))
    parser.add_argument('--transient',type=float,default=200);parser.add_argument('--duration',type=float,default=400);parser.add_argument('--rtol',type=float,default=1e-7)
    args=parser.parse_args();gpu=PhotoFeedbackGPU()
    if args.validate:validate(gpu);raise SystemExit()
    if args.refine:
        old=np.load(args.refine/'scan.npz');meta=json.loads((args.refine/'summary.json').read_text());ids=meta['candidate_indices'];p=old['parameters'][ids];s=old['starts'][ids]
        if not len(p):raise SystemExit('No candidates to refine')
    else:
        u=qmc.Sobol(5,scramble=True,seed=20261029).random_base2(args.power);p=[];s=[];eigenvalues=[]
        for row in u:
            x=.55+.7*row[0];b=.025*16**row[1];tau=.01*1000**row[2];k=.2*50**row[3];width=.06*6**row[4]
            pp,ss=parameter_at_equilibrium(x,b,tau,k,width);ev,_=stability(pp)
            for start in [ss+np.array([.005,.001,0]),np.array([.65,.3,.1]),np.array([1.4,-.3,.8])]:p.append(pp);s.append(start);eigenvalues.append(ev)
        p=np.array(p);s=np.array(s)
    print('Scanning',len(p),'phenomenological model trajectories',flush=True)
    result=gpu.scan(p,s,transient=args.transient,duration=args.duration,rtol=args.rtol,atol=args.rtol/1000,
                    progress=lambda f,t:print(f'{f:.1%}, {t:.1f}s',flush=True))
    args.output.mkdir(exist_ok=True,parents=True);arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray)}
    if not args.refine:arrays['equilibrium_eigenvalues']=np.array(eigenvalues)
    np.savez_compressed(args.output/'scan.npz',**arrays)
    good=(result['status']==0)&(result['lle']>.005)&(result['lle_windows'].min(axis=1)>.002)&((result['stats'][:,7]-result['stats'][:,6])>1e-5)
    meta={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    meta.update(trajectories=len(p),status_counts={str(i):int((result['status']==i).sum()) for i in np.unique(result['status'])},candidate_indices=np.flatnonzero(good).tolist(),max_lle=float(np.nanmax(result['lle'])),
                parameter_columns=['a','b','tau','k','width','equilibrium_height','unused','unused'],
                criterion='status 0, mean LLE > .005, all four windows > .002, late height span > 1e-5',
                interpretation='A deliberately phenomenological extension. A candidate here is not evidence of chaos in the published physical device.')
    if args.refine:meta.update(parent_source=str(args.refine),parent_indices=ids)
    (args.output/'summary.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
