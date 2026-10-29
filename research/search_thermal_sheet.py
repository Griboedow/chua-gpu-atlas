"""Survey the literal sheet ODE family, with explicit state-domain flags."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import qmc

from thermomechanical_sheet import equilibrium,field_jacobian
from thermal_sheet_gpu import ThermalSheetGPU


def save(result,path,extra=None):
    path.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(path/'scan.npz',**arrays)
    good=(result['status']==0)&(result['lle']>.005)&(result['lle_windows'].min(axis=1)>.002)&((result['stats'][:,7]-result['stats'][:,6])>1e-6)
    admissible=(result['domain_flags'][:,1]==0)
    meta={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    meta.update(trajectories=len(good),candidate_indices=np.flatnonzero(good).tolist(),
                admissible_candidate_indices=np.flatnonzero(good&admissible).tolist(),
                admissible_completed=int(((result['status']==0)&admissible).sum()),
                status_counts={str(i):int((result['status']==i).sum()) for i in np.unique(result['status'])},
                domain_counts={str(i):int((result['domain_flags'][:,1]==i).sum()) for i in np.unique(result['domain_flags'][:,1])},
                max_lle=float(np.nanmax(result['lle'])) if np.isfinite(result['lle']).any() else None,
                criterion='completed, mean LLE > .005, all four windows > .002, late X span > 1e-6; admissibility additionally requires no observation domain flags',
                interpretation='Exploration of the printed ODE family, including source-reported cubic damping; no calibrated device or novelty claim.')
    if extra:meta.update(extra)
    (path/'summary.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps({k:meta[k] for k in ['trajectories','status_counts','admissible_completed','candidate_indices','admissible_candidate_indices','max_lle','timing']}),flush=True)


def validate(gpu):
    p=np.array([[10,5,.2,.5,.1],[560,74,.237,.8,10],[2,.1,.1,3,1]])
    starts=np.array([equilibrium(pp)+[.01,0,-.01] for pp in p])
    result=gpu.scan(p,starts,transient=0,duration=2,rtol=1e-10,atol=1e-12)
    differences=[];lle_errors=[]
    for pp,start,end,lle in zip(p,starts,result['final_states'],result['lle']):
        def rhs(t,state):
            f,j=field_jacobian(state[:3],pp)
            return np.r_[f,j@state[3:]]
        initial=np.r_[start,np.array([1,2,3])/np.sqrt(14)]
        sol=solve_ivp(rhs,(0,2),initial,method='DOP853',rtol=1e-12,atol=1e-14)
        if not sol.success:raise RuntimeError(sol.message)
        differences.append(float(np.max(abs(sol.y[:3,-1]-end))))
        lle_errors.append(float(abs(np.log(np.linalg.norm(sol.y[3:,-1]))/2-lle)))
    passed=bool(max(differences)<1e-7 and max(lle_errors)<1e-7 and (result['status']==0).all())
    output=dict(gpu_vs_dop853_state_errors=differences,finite_time_lle_errors=lle_errors,passed=passed,backend=gpu.info)
    Path('results/thermomechanical-sheet/validation.json').write_text(json.dumps(output,indent=2));print(json.dumps(output),flush=True)
    if not passed:raise RuntimeError('GPU state/tangent validation failed')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--validate',action='store_true')
    parser.add_argument('--power',type=int,default=12);parser.add_argument('--seed',type=int,default=2026102943)
    parser.add_argument('--refine',type=Path);parser.add_argument('--top',type=int,default=24)
    parser.add_argument('--indices',type=int,nargs='+',help='Explicit source indices for a follow-up, preserving selection provenance')
    parser.add_argument('--transient',type=float,default=300);parser.add_argument('--duration',type=float,default=600)
    parser.add_argument('--rtol',type=float,default=1e-7);parser.add_argument('--step-factor',type=float,default=1.)
    parser.add_argument('--max-attempts',type=int,default=3000000);parser.add_argument('--max-wall-seconds',type=float,default=600)
    parser.add_argument('--output',type=Path,default=Path('results/thermomechanical-sheet/scout'))
    args=parser.parse_args();gpu=ThermalSheetGPU()
    if args.validate:validate(gpu);return
    if args.indices is not None and args.refine is None:parser.error('--indices requires --refine')
    if args.refine:
        old=np.load(args.refine/'scan.npz');valid=np.flatnonzero((old['status']==0)&np.isfinite(old['lle']))
        ids=np.asarray(args.indices,int) if args.indices is not None else valid[np.argsort(-old['lle'][valid])[:args.top]]
        if len(ids)==0 or len(np.unique(ids))!=len(ids) or not np.isin(ids,valid).all():
            raise ValueError('Follow-up indices must be unique completed source trajectories')
        p=old['parameters'][ids];starts=old['starts'][ids]
        extra=dict(parent_source=str(args.refine),parent_indices=ids.tolist())
    else:
        rows=qmc.Sobol(5,scramble=True,seed=args.seed).random_base2(args.power)
        p=[];starts=[];labels=[];diagnostics=[]
        for sample,row in enumerate(rows):
            K=10**(-1+5*row[0]);G=K*10**(-3+4*row[1]);S=10**(-2+4*row[2])
            eq=equilibrium([K,G,1,S,0]);ke=K*(1-S*eq[2]);F=G+K*S*eq[0]
            critical=2*F/(ke+1+np.sqrt((ke+1)**2+4*F))
            B=critical*10**(-1.5+2*row[3]);D=10**(-3+5*row[4])/(K*eq[0]**2)
            pp=np.array([K,G,B,S,D])
            amplitude=min(eq[0],eq[2])
            for seed,offset in enumerate([.001,-.1]):
                p.append(pp);starts.append(eq+offset*amplitude*np.array([1,0,-1]))
                labels.append([sample,seed]);diagnostics.append([critical,ke,F])
        p=np.array(p);starts=np.array(starts)
        extra=dict(seed=args.seed,sampling='Sobol: log10 K[-1,4], G/K[-3,1], S[-2,2], B/B_H[-1.5,.5], D*K*Xeq^2[-3,2]; two equilibrium perturbations',
                   source_relation='Freely explores source equation parameters; the source does not provide an unambiguous baseline cubic damping calibration')
    print(f'Scanning {len(p)} thermal sheet trajectories',flush=True)
    result=gpu.scan(p,starts,transient=args.transient,duration=args.duration,rtol=args.rtol,
                    atol=args.rtol/1000,step_factor=args.step_factor,
                    max_attempts=args.max_attempts,max_wall_seconds=args.max_wall_seconds,
                    progress=lambda f,t:print(f'{f:.1%}, {t:.1f}s',flush=True))
    if not args.refine:result.update(labels=np.array(labels),linear_diagnostics=np.array(diagnostics))
    save(result,args.output,extra)


if __name__=='__main__':main()
