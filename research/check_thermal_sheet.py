"""Independent SciPy full-spectrum check of selected thermal sheet trajectories."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from thermomechanical_sheet import field_jacobian


def check(parameters, start, transient, duration, rtol, max_step, qr_interval, method):
    p=np.asarray(parameters,float);state=np.r_[start,np.eye(3).ravel(),0.]
    logs=np.zeros((4,3));divergences=np.zeros(4);bounds=np.array([np.inf,-np.inf,np.inf,-np.inf,np.inf])
    flags=np.zeros(2,dtype=int);nfev=0;wall=time.perf_counter()
    def rhs(t,z):
        f,j=field_jacobian(z[:3],p)
        return np.r_[f,(j@z[3:12].reshape(3,3)).ravel(),np.trace(j)]
    def escaped(t,z):return 1e4-np.max(abs(z[:3]))
    escaped.terminal=True;escaped.direction=-1
    for window,length in [(-1,transient)]+[(i,duration/4) for i in range(4)]:
        t=0.
        while t<length:
            dt=min(qr_interval,length-t)
            sol=solve_ivp(rhs,(0,dt),state,method=method,rtol=rtol,atol=rtol/100,
                          max_step=max_step,events=escaped)
            nfev+=sol.nfev
            if not sol.success or sol.t[-1]<dt:
                raise RuntimeError(f'Integration incomplete: {sol.message}; end {sol.t[-1]} / {dt}')
            state=sol.y[:,-1].copy();q,r=np.linalg.qr(state[3:12].reshape(3,3))
            diagonal=np.diag(r)
            if np.any(abs(diagonal)<np.finfo(float).tiny):raise RuntimeError('Singular tangent frame')
            q*=np.where(diagonal<0,-1,1)[None,:]
            if window>=0:
                logs[window]+=np.log(abs(diagonal));divergences[window]+=state[-1]
                x,y=sol.y[0],sol.y[2]
                bounds=np.array([min(bounds[0],x.min()),max(bounds[1],x.max()),
                                 min(bounds[2],y.min()),max(bounds[3],y.max()),
                                 min(bounds[4],(1-p[3]*y).min())])
            x,y=sol.y[0],sol.y[2]
            flag=int((x < -1e-8).any())|int((y < -1e-8).any())*2|int((y > 1+1e-8).any())*4|int((1-p[3]*y < -1e-8).any())*8
            flags[int(window>=0)]|=flag
            state[3:12]=q.ravel();state[-1]=0.;t+=dt
            if length-t < 1e-12*max(1,length):t=length
        print(f'Completed window {window}; time {time.perf_counter()-wall:.1f}s',flush=True)
    spectrum=logs.sum(axis=0)/duration
    return dict(parameters=p.tolist(),start=np.asarray(start).tolist(),spectrum=spectrum.tolist(),
                spectrum_windows=(logs/(duration/4)).tolist(),mean_divergence=float(divergences.sum()/duration),
                divergence_identity_error=float(abs(spectrum.sum()-divergences.sum()/duration)),
                final_state=state[:3].tolist(),observation_bounds=bounds.tolist(),
                bound_columns=['min_X','max_X','min_Y','max_Y','min_stiffness_ratio'],
                domain_flags=flags.tolist(),domain_flag_columns=['transient','observation'],
                settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol/100,
                              max_step=max_step,qr_interval=qr_interval,method=method),
                nfev=nfev,wall_seconds=time.perf_counter()-wall,
                caveat='Finite-time numerical evidence; domain flags sample accepted solver nodes, not an interval bound.')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path);parser.add_argument('--index',type=int,required=True)
    parser.add_argument('--transient',type=float,default=600);parser.add_argument('--duration',type=float,default=1600)
    parser.add_argument('--rtol',type=float,default=1e-10);parser.add_argument('--max-step',type=float,default=.02)
    parser.add_argument('--qr-interval',type=float,default=.1);parser.add_argument('--method',choices=['DOP853','Radau'],default='DOP853')
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    old=np.load(args.source/'scan.npz');p=old['parameters'][args.index];start=old['starts'][args.index]
    result=check(p,start,args.transient,args.duration,args.rtol,args.max_step,args.qr_interval,args.method)
    result.update(source=str(args.source),source_index=args.index)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
