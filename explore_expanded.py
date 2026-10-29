"""Expand the previous survey fourfold in alpha and nonlinear-slope depth.

Samples beta relative to the analytic outer-equilibrium stability scale, so
the upper beta search limit grows quadratically with alpha. Old data remain
unchanged. Bulk integration uses the adaptive FP64 OpenCL backend.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.stats import qmc
from adaptive_gpu import AdaptiveGPU
from explore_extent import starting_states
from run import save_json, progress_printer

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs'/'expanded'

def write_run(path,parameters,result,**extra):
    arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(path,parameters=parameters,**arrays,**extra)
    save_json(path.with_suffix('.json'),{k:v for k,v in result.items() if not isinstance(v,np.ndarray)})

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--power',type=int,default=14)
    parser.add_argument('--output',type=Path,default=OUT)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if not 5<=args.power<=18:parser.error('power must be 5–18')
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    bounds=np.array([[2,.01,.003],[1024,1.3,28.]])
    sample=qmc.Sobol(3,scramble=True,seed=20261030).random_base2(args.power)
    a,q,depth=np.exp(np.log(bounds[0])+sample*np.log(bounds[1]/bounds[0])).T
    beta=q*(5*a/7+10*a*a/49)
    parameters=np.c_[a,beta,-1-depth,np.full(len(a),-5/7)]
    order=np.argsort(np.maximum(a*(2/7),np.sqrt(beta)))
    parameters=parameters[order];q=q[order]
    gpu=AdaptiveGPU();started=time.perf_counter()
    if args.resume and (out/'scout.npz').exists():
        scout=dict(np.load(out/'scout.npz'))
        if not np.array_equal(scout['parameters'],parameters):raise ValueError('Saved scout uses different sampling')
    else:
        result=gpu.scan(parameters,initial_states=starting_states(parameters),transient=100,duration=200,
                        rtol=1e-6,atol=1e-9,progress=progress_printer('Expanded scout'))
        write_run(out/'scout.npz',parameters,result,beta_over_critical=q,sobol_index=order)
        scout=dict(np.load(out/'scout.npz'))
    indices=np.flatnonzero((scout['classes']==4).any(axis=1))
    print(f'{len(indices)} scout candidates; refining tolerance and observation',flush=True)
    if args.resume and (out/'verified.npz').exists():
        fine=dict(np.load(out/'verified.npz'))
        if not np.array_equal(fine['indices'],indices):raise ValueError('Saved verification has different candidates')
    else:
        if not len(indices):raise RuntimeError('No scout candidates; inspect scout data before continuing')
        result=gpu.scan(parameters[indices],initial_states=starting_states(parameters[indices]),transient=200,duration=600,
                        rtol=1e-8,atol=1e-11,progress=progress_printer('Expanded verification'))
        write_run(out/'verified.npz',parameters[indices],result,indices=indices)
        fine=dict(np.load(out/'verified.npz'))
    keep=(fine['classes']==4).any(axis=1);accepted=fine['parameters'][keep]
    scout_meta=json.loads((out/'scout.json').read_text());fine_meta=json.loads((out/'verified.json').read_text())
    original=np.load(ROOT/'outputs/extent/verified.npz')
    original_p=original['parameters'][(original['classes']==4).any(axis=1)]
    summary={'sampling':'Scrambled Sobol, uniform in log(alpha), log(beta/beta_critical), log(-1-m0)',
        'seed':20261030,'points':len(parameters),'starts_per_point':3,'fixed_m1':-5/7,
        'ranges':{'alpha':[2,1024],'beta_over_critical':[.01,1.3],'depth_minus_one_minus_m0':[.003,28]},
        'beta_critical':'5*alpha/7 + 10*alpha**2/49','theoretical_beta_max':1.3*(5*1024/7+10*1024**2/49),
        'escape_cutoff':1e6,'scout':scout_meta,'verification':fine_meta,
        'scout_candidates':len(indices),'verified_candidates':int(keep.sum()),
        'numerical_failures_scout':int((scout['classes']==5).sum()),'numerical_failures_verification':int((fine['classes']==5).sum()),
        'candidate_parameter_min':accepted.min(axis=0).tolist(),'candidate_parameter_max':accepted.max(axis=0).tolist(),
        'previous_candidate_min':original_p.min(axis=0).tolist(),'previous_candidate_max':original_p.max(axis=0).tolist(),
        'alpha_above_256':int((accepted[:,0]>256).sum()),'m0_below_minus_8':int((accepted[:,2]<-8).sum()),
        'alpha_above_768':int((accepted[:,0]>768).sum()),'m0_below_minus_22':int((accepted[:,2]<-22).sum()),
        'outside_previous_search':int(((accepted[:,0]>256)|(accepted[:,2]<-8)|(accepted[:,1]/accepted[:,0]>64)).sum()),
        'device':gpu.info,'current_invocation_seconds':time.perf_counter()-started,
        'scope':'Sparse finite-time screen. Adaptive tolerances and escape cutoff differ from the old survey. '
                'No global boundedness, connectedness or rigorous attractor-existence conclusion. '
                'Candidate fractions depend on the changed sampling measure; do not compare them as volumes.'}
    save_json(out/'summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
