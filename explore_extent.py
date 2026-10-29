"""Coarse multiscale exploration of the Chua parameter region, not an infinity proof."""
from __future__ import annotations
import argparse
import time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from chua_gpu import GPU
from run import save_json,progress_printer

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs'/'extent'

def starting_states(p):
    xstar=(p[:,3]-p[:,2])/(1+p[:,3])
    starts=np.zeros((len(p),3,3))
    starts[:,0,0]=.7
    starts[:,1,0]=.9*xstar;starts[:,1,1]=.02;starts[:,1,2]=-.9*xstar
    starts[:,2,0]=1.1*xstar;starts[:,2,1]=-.02;starts[:,2,2]=-1.1*xstar
    return starts

def scan_groups(gpu,p,transient,duration,refine=0,label='scout'):
    # Bound the time step against inner/outer rates and the LC frequency.
    rate=np.maximum.reduce([p[:,0]*np.abs(1+p[:,2]),p[:,0]*np.abs(1+p[:,3]),np.sqrt(p[:,1]),np.ones(len(p))])
    groups=np.maximum(0,np.ceil(np.log2(.01*rate/.15))).astype(int)+refine
    output={"classes":np.empty((len(p),3),np.uint8),"lle":np.empty((len(p),3)),
            "lle_windows":np.empty((len(p),3,4)),"stats":np.empty((len(p),3,8)),
            "final_states":np.empty((len(p),3,3)),"dt":.01/2**groups}
    records=[]
    for group in np.unique(groups):
        idx=np.flatnonzero(groups==group);dt=.01/2**int(group)
        print(f'{label}: {len(idx)} points, dt={dt:.9g}',flush=True)
        result=gpu.scan(p[idx],initial_states=starting_states(p[idx]),dt=dt,transient=transient,
                        duration=duration,escape=1000,progress=progress_printer(f'{label}/dt{group}'))
        for key in ['classes','lle','lle_windows','stats','final_states']:
            output[key][idx]=result[key]
        records.append({'points':len(idx),'dt':dt,'timing':result['timing']})
    return output,records

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--power',type=int,default=15,help='Number of Sobol points = 2**power')
    parser.add_argument('--alpha-max',type=float,default=256)
    parser.add_argument('--depth-max',type=float,default=7)
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    if args.alpha_max<=2 or args.depth_max<=.003 or not 5<=args.power<=20:
        parser.error('invalid exploration limits')
    sample=qmc.Sobol(3,scramble=True,seed=20261029).random_base2(args.power)
    low=np.log([2,.25,.003]);high=np.log([args.alpha_max,64,args.depth_max])
    transformed=np.exp(low+sample*(high-low))
    alpha,ratio,depth=transformed.T
    p=np.column_stack([alpha,alpha*ratio,-1-depth,np.full(len(alpha),-5/7)])
    gpu=GPU();start=time.perf_counter()
    coarse,coarse_records=scan_groups(gpu,p,100,200,label='wide scout')
    np.savez_compressed(out/'scout.npz',parameters=p,ratio=ratio,depth=depth,**coarse)
    candidates=(coarse['classes']==4).any(axis=1)
    idx=np.flatnonzero(candidates)
    print(f'{len(idx)} scout candidates; independently refining them',flush=True)
    fine,fine_records=scan_groups(gpu,p[idx],300,800,refine=1,label='candidate verification') if len(idx) else ({},[])
    np.savez_compressed(out/'verified.npz',indices=idx,parameters=p[idx],**fine)
    keep=(fine['classes']==4).any(axis=1) if len(idx) else np.zeros(0,bool)
    accepted=p[idx][keep]
    # Recheck extremal examples for much longer than the scouting pass.
    extremes=[]
    if len(accepted):
        for column in range(3):
            extremes.extend([int(np.argmin(accepted[:,column])),int(np.argmax(accepted[:,column]))])
        extremes.extend(np.linspace(0,len(accepted)-1,min(12,len(accepted)),dtype=int).tolist())
        extremes=np.unique(extremes)
        long,long_records=scan_groups(gpu,accepted[extremes],600,1600,refine=2,label='long extremal checks')
        np.savez_compressed(out/'extremes.npz',parameters=accepted[extremes],**long)
    else:long_records=[]
    summary={'sampling':'Scrambled Sobol, uniform in log(alpha), log(beta/alpha), log(-1-m0)',
        'seed':20261029,'points':len(p),'starts_per_point':3,'fixed_m1':-5/7,
        'ranges':{'alpha':[2,args.alpha_max],'beta_over_alpha':[.25,64],'depth_minus_one_minus_m0':[.003,args.depth_max]},
        'initial_states':'(0.7,0,0), (0.9*x*,0.02,-0.9*x*), (1.1*x*,-0.02,-1.1*x*), x*=(m1-m0)/(1+m1)',
        'escape_cutoff':1000,'scout_transient':100,'scout_observation':200,
        'verification_transient':300,'verification_observation':800,
        'scout_candidates':int(candidates.sum()),'verified_candidates':int(keep.sum()),
        'verified_parameter_min':accepted.min(axis=0).tolist() if len(accepted) else None,
        'verified_parameter_max':accepted.max(axis=0).tolist() if len(accepted) else None,
        'coarse_groups':coarse_records,'verification_groups':fine_records,'long_check_groups':long_records,
        'device':gpu.info,'wall_seconds':time.perf_counter()-start,
        'scope':'Sparse finite-domain evidence. No proof of global closure, unboundedness, connectedness, or absence in empty regions. '
                'Only m0<-1 is explored in this continuation of the double-scroll family. '
                'Raw positive-node fractions depend on the log sampling density and are not parameter-space volumes.'}
    save_json(out/'summary.json',summary)
    print('Verified extent',summary['verified_parameter_min'],summary['verified_parameter_max'],flush=True)
    print('Total seconds',summary['wall_seconds'],flush=True)

if __name__=='__main__':main()
