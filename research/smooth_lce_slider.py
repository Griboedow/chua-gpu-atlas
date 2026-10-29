"""Sensitivity to a finite switching width, explicitly changing the ideal model."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit


def smooth_field(state, p, width):
    k,b,intensity,d=p;u,v,phi=state
    light=expit(u/width);slope=light*(1-light)/width
    f=np.array([v,k*(d-u-phi)-b*v,intensity*light*(1-phi)-phi])
    jac=np.array([[0.,1.,0.],[-k,-b,-k],[intensity*slope*(1-phi),0.,-1-intensity*light]])
    return f,jac


def spectrum(p,start,width,transient,duration,rtol,max_step=.01):
    interval=.1;steps=round((transient+duration)/interval);burn=round(transient/interval)
    state=np.r_[start,np.eye(3).ravel(),0.];sums=np.zeros((4,3));div=np.zeros(4);history=[]
    wall=time.perf_counter();last=wall
    def rhs(t,y):
        f,j=smooth_field(y[:3],p,width)
        return np.r_[f,(j@y[3:12].reshape(3,3)).ravel(),np.trace(j)]
    for step in range(steps):
        state[-1]=0.
        sol=solve_ivp(rhs,(0,interval),state,method='DOP853',rtol=rtol,atol=rtol*.01,max_step=max_step)
        if not sol.success:raise RuntimeError(sol.message)
        state=sol.y[:,-1];Q,R=np.linalg.qr(state[3:12].reshape(3,3));state[3:12]=Q.ravel()
        if step>=burn:
            window=min(3,4*(step-burn)//(steps-burn));sums[window]+=np.log(abs(np.diag(R)));div[window]+=state[-1]
            if (step-burn)%100==99:history.append([(step-burn+1)*interval,*list(sums.sum(axis=0)/((step-burn+1)*interval))])
        if time.perf_counter()-last>20:
            print(f'width={width:.3g}: {(step+1)/steps:.1%}, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
    exponents=sums.sum(axis=0)/duration
    result=dict(width=float(width),transient=transient,duration=duration,rtol=rtol,max_step=max_step,spectrum=exponents.tolist(),
                window_spectra=(sums/(duration/4)).tolist(),divergence_identity_error=float(abs(exponents.sum()-div.sum()/duration)),
                wall_seconds=time.perf_counter()-wall,final_state=state[:3].tolist())
    return result,np.array(history)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fractions',nargs='+',type=float,default=[1e-6,1e-5,1e-4,1e-3,.01])
    parser.add_argument('--transient',type=int,default=200);parser.add_argument('--duration',type=int,default=800)
    parser.add_argument('--rtol',type=float,default=1e-10);parser.add_argument('--output',type=Path,default=Path('results/lce-slider/smooth'))
    parser.add_argument('--max-step',type=float,default=.01)
    args=parser.parse_args();data=np.load('results/lce-slider/analytic-2222/trace.npz')
    p=data['parameters'];start=np.r_[0.,data['section'][-1,1:3]];q=p[2]/(1+p[2])
    args.output.mkdir(parents=True,exist_ok=True);results=[]
    # Verify the new smooth Jacobian at a transition-scale state before integration.
    width=q*1e-3;sample=np.array([width*.3,.001,.0003]);f,j=smooth_field(sample,p,width);eps=width*1e-4
    numerical=np.column_stack([(smooth_field(sample+eps*e,p,width)[0]-smooth_field(sample-eps*e,p,width)[0])/(2*eps) for e in np.eye(3)])
    error=float(np.max(abs(j-numerical)))
    if error>1e-4:raise RuntimeError(f'Smooth Jacobian check failed: {error}')
    for fraction in args.fractions:
        result,history=spectrum(p,start,q*fraction,args.transient,args.duration,args.rtol,args.max_step)
        result['width_over_q']=fraction;results.append(result)
        tag=f'{fraction:.6g}' if fraction not in [1e-6,1e-5,1e-4,1e-3,.01] else f'{fraction:.0e}'
        (args.output/f'width-{tag}.json').write_text(json.dumps(result,indent=2))
        np.savez_compressed(args.output/f'width-{tag}.npz',history=history,parameters=p,start=start)
        summary=dict(parameters=p.tolist(),start=start.tolist(),changes_model=True,
                     smoothing='H(u) replaced by logistic(u/width)',jacobian_error=error,results=results,
                     interpretation='Finite-width model sensitivity, not a calibrated contact or optical response')
        (args.output/'summary.json').write_text(json.dumps(summary,indent=2))
        print(json.dumps(result),flush=True)


if __name__=='__main__':main()
