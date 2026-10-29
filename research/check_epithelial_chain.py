"""Independent adaptive DOP853 checks of saved cell-train trajectories."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from epithelial_chain import field_jacobian


def check(source,index,out,transient,duration,rtol):
    data=np.load(source/'scan.npz');p=data['parameters'][index];y=data['final_states'][index].copy();d=len(y);n=(d+1)//3
    v=np.sin(np.arange(d)+1);v/=np.linalg.norm(v);z=np.r_[y,v]
    logs=np.zeros(4);mins=np.full(3,np.inf);maxstate=0.;records=[];wall=time.perf_counter();calls=0
    def rhs(t,z):
        f,j=field_jacobian(z[:d],p);return np.r_[f,j@z[d:]]
    for window,length in [(-1,transient)]+[(i,duration/4) for i in range(4)]:
        chunks=int(np.ceil(length));interval=length/chunks if chunks else 0
        for i in range(chunks):
            sol=solve_ivp(rhs,(0,interval),z,method='DOP853',rtol=rtol,atol=rtol/100,max_step=.1)
            if not sol.success:raise RuntimeError(sol.message)
            z=sol.y[:,-1];norm=np.linalg.norm(z[d:]);z[d:]/=norm;calls+=sol.nfev
            if window>=0:
                logs[window]+=np.log(norm)
                u=np.vstack([np.zeros(sol.y.shape[1]),sol.y[:n-1],np.zeros(sol.y.shape[1])])
                mins=np.minimum(mins,[np.min(1+np.diff(u,axis=0)),np.min(1-sol.y[n-1:2*n-1]),np.min(1+sol.y[2*n-1:d])])
                maxstate=max(maxstate,float(np.max(abs(sol.y[:d]))))
                if window==3:records.append(z[:d].copy())
        print('CPU window',window,'seconds',round(time.perf_counter()-wall,1),flush=True)
    out.mkdir(parents=True,exist_ok=True)
    result=dict(source=str(source),source_index=index,parameters=p.tolist(),cells=n,method='DOP853 analytic tangent, renormalized each <=1 unit',
                initialization='Saved GPU final state, followed by additional transient',transient=transient,duration=duration,rtol=rtol,atol=rtol/100,max_step=.1,
                lle=float(logs.sum()/duration),lle_windows=(logs/(duration/4)).tolist(),domain_minima=mins.tolist(),max_abs_state=maxstate,
                domain_columns=['cell length','effective rest length','relative concentration'],wall_seconds=time.perf_counter()-wall,rhs_calls=calls)
    (out/'summary.json').write_text(json.dumps(result,indent=2));np.savez_compressed(out/'trajectory.npz',final_state=z[:d],late_samples=np.array(records),parameters=p)
    print(result,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,default=Path('results/epithelial-chain/n16-scout'));ap.add_argument('--index',type=int,default=340)
    ap.add_argument('--output',type=Path,default=Path('results/epithelial-chain/cpu-340'));ap.add_argument('--transient',type=float,default=500);ap.add_argument('--duration',type=float,default=1000);ap.add_argument('--rtol',type=float,default=1e-11)
    args=ap.parse_args();check(args.source,args.index,args.output,args.transient,args.duration,args.rtol)
