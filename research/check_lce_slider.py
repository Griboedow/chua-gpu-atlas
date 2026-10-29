"""Independent DOP853 event returns and two-exponent transverse spectrum."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from lce_slider import flow, next_crossing


def dop853_crossing(section, p, on, rtol, normalized=False):
    k,b,intensity,d=p
    A=np.array([[0.,1.,0.],[-k,-b,-k],[0.,0.,-1-intensity*on]])
    basis=np.array([[0.,0.],[1.,0.],[0.,1.]])
    initial=np.r_[0.,section,basis.ravel()]
    def rhs(t,y):return np.r_[flow(y[:3],p,on,normalized),(A@y[3:].reshape(3,2)).ravel()]
    def event(t,y):return y[0]
    event.direction=-1 if on else 1
    event.terminal=True
    period=2*np.pi/np.sqrt(k-b*b/4)
    result=solve_ivp(rhs,(0,128*period),initial,method='DOP853',rtol=rtol,atol=rtol*.01,
                     events=event,max_step=period/64)
    if not result.success or not len(result.t_events[0]):raise RuntimeError('DOP853 return failed')
    t=result.t_events[0][0];end=result.y_events[0][0];f=flow(end[:3],p,on,normalized)
    if t<=0 or abs(f[0])<1e-10:raise RuntimeError('Nontransverse DOP853 return')
    derivative=end[3:].reshape(3,2)
    derivative-=np.outer(f,derivative[0]/f[0])
    return end[1:3],t,derivative[1:]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=Path('results/lce-slider/wide'))
    parser.add_argument('--index',type=int,default=2222)
    parser.add_argument('--method',choices=['analytic','DOP853'],default='DOP853')
    parser.add_argument('--transient',type=int,default=300)
    parser.add_argument('--cycles',type=int,default=1200)
    parser.add_argument('--rtol',type=float,default=1e-11)
    parser.add_argument('--subdivisions',type=int,default=256)
    parser.add_argument('--normalized',action='store_true')
    parser.add_argument('--output',type=Path,default=Path('results/lce-slider/cpu'))
    args=parser.parse_args()
    if args.cycles%4:parser.error('Cycles must be divisible by four')
    data=np.load(args.source/'scan.npz');p=data['parameters'][args.index];s=data['starts'][args.index].copy()
    Q=np.eye(2);sums=np.zeros((4,2));durations=np.zeros(4);divergence=np.zeros(4)
    traces=[];history=[];wall=time.perf_counter();last=wall
    for cycle in range(args.transient+args.cycles):
        jac=np.eye(2);elapsed=0.;div=0.;vstart=s[0]
        for on in [True,False]:
            s,dt,j=(dop853_crossing(s,p,on,args.rtol,args.normalized) if args.method=='DOP853'
                    else next_crossing(s,p,on,args.subdivisions,args.normalized))
            jac=j@jac;elapsed+=dt;div-=dt*(p[1]+1+p[2]*on)
        Q,R=np.linalg.qr(jac@Q);logs=np.log(abs(np.diag(R)))
        # Map Jacobian includes the velocity ratio between the two sections.
        exact_logdet=div+np.log(abs(vstart/s[0]))
        if cycle>=args.transient:
            window=min(3,4*(cycle-args.transient)//args.cycles)
            sums[window]+=logs;durations[window]+=elapsed;divergence[window]+=exact_logdet
            traces.append([durations.sum(),*s,elapsed])
            if (cycle-args.transient)%20==19:history.append([durations.sum(),*(sums.sum(axis=0)/durations.sum())])
        if time.perf_counter()-last>20:
            print(f'{args.method} {(cycle+1)/(args.transient+args.cycles):.1%}, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
    spectrum=sums.sum(axis=0)/durations.sum()
    result=dict(source=str(args.source),index=args.index,parameters=p.tolist(),start=data['starts'][args.index].tolist(),
                method=args.method,rtol=args.rtol,subdivisions=args.subdivisions,normalized_photo=args.normalized,transient_cycles=args.transient,
                observation_cycles=args.cycles,elapsed_observation_time=float(durations.sum()),
                transverse_spectrum=spectrum.tolist(),window_spectra=(sums/durations[:,None]).tolist(),
                sum_rule_error=float(abs(spectrum.sum()-divergence.sum()/durations.sum())),
                wall_seconds=time.perf_counter()-wall,
                interpretation='Two exponents on the event-defined return section per elapsed time; autonomous neutral flow direction omitted.')
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'spectrum.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(args.output/'trace.npz',section=np.array(traces),history=np.array(history),parameters=p)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
