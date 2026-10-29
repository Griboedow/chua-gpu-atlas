"""Independent event-located DOP853/Radau and full-spectrum QR check.

The ODE is integrated numerically between switches. Only the model vector
field and the physical saltation rule are shared with the analytic solver.
The stored three-dimensional views of this six-state model are projections.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from coupled_operators import flow, initial_modes


def integrate(p, start, transient=100., observation=400., method='DOP853',rtol=1e-11,qr_interval=.05,save_prefix=None,width=0.):
    k,b,dm,dp=p
    if transient<0 or observation<=0 or rtol<=0 or qr_interval<=0 or width<0:
        raise ValueError('Invalid integration settings')
    a=np.array([[0,1,0],[-k,-b,k],[0,0,-1.]])
    jac=np.zeros((6,6));jac[:3,:3]=a;jac[3:,3:]=a
    y=start.copy();q=np.eye(6);modes=initial_modes(y,p)
    qr_interval=min(qr_interval,2/b)
    time_now=0.; event_count=0; min_speed=np.inf; min_separation=np.inf; last_event=-np.inf
    sums=np.zeros((4,6));durations=np.zeros(4);states=[];times=[];sections=[];event_records=[]
    minimum=y.copy();maximum=y.copy();wall=time.perf_counter();last_print=wall;max_sum_error=0.
    for window,duration in [(-1,transient)]+[(i,observation/4) for i in range(4)]:
        target=time_now+duration
        while time_now < target-1e-11:
            step=min(qr_interval,target-time_now)
            def rhs(t,z):
                if width<=0:return np.r_[flow(z[:6],p,modes),(jac@z[6:].reshape(6,6)).ravel()]
                sm=expit((z[3]-dp)/width);sp=expit((dm-z[0])/width)
                f=flow(z[:6],p,[sm,sp]);j=jac.copy()
                j[2,3]=sm*(1-sm)/width;j[5,0]=-sp*(1-sp)/width
                return np.r_[f,(j@z[6:].reshape(6,6)).ravel()]
            def event_m(t,z):return z[0]-dm
            def event_p(t,z):return z[3]-dp
            event_m.terminal=True;event_p.terminal=True
            event_m.direction=-(1-2*modes[1]);event_p.direction=-(2*modes[0]-1)
            sol=solve_ivp(rhs,(0,step),np.r_[y,q.ravel()],method=method,rtol=rtol,atol=rtol*.01,
                          max_step=min(step,2*np.pi/np.sqrt(k-b*b/4)/32),events=[event_m,event_p] if width<=0 else None)
            if not sol.success:raise RuntimeError(sol.message)
            elapsed=sol.t[-1]
            if elapsed<1e-12:raise RuntimeError('Unresolved zero-time event')
            y=sol.y[:6,-1];matrix=sol.y[6:,-1].reshape(6,6)
            event_indices=[i for i,ev in enumerate(sol.t_events) if len(ev)] if width<=0 else []
            if len(event_indices)>1:raise RuntimeError('Simultaneous events require separate analysis')
            if event_indices:
                which=event_indices[0];j=3*which;speed=y[j+1]
                if abs(speed)<1e-9:raise RuntimeError('Grazing event')
                y[j]=p[2+which];changed=modes.copy();changed[1-which]=1-changed[1-which]
                jump=flow(y,p,changed)-flow(y,p,modes)
                matrix+=np.outer(jump,matrix[j]/speed)
                modes=changed;event_count+=1
                min_speed=min(min_speed,abs(speed));min_separation=min(min_separation,time_now+elapsed-last_event);last_event=time_now+elapsed
                if window>=0:
                    event_records.append([time_now+elapsed,which,speed])
                    if which==0 and speed>0:sections.append(y.copy())
            q,r=np.linalg.qr(matrix);logs=np.log(abs(np.diag(r)))
            max_sum_error=max(max_sum_error,abs(logs.sum()+2*(b+1)*elapsed))
            if window>=0:
                sums[window]+=logs;durations[window]+=elapsed
                states.append(y.copy());times.append(time_now+elapsed)
            minimum=np.minimum(minimum,y);maximum=np.maximum(maximum,y);time_now+=elapsed
            if time.perf_counter()-last_print>30:
                print(f't={time_now:.3f}/{transient+observation}, events={event_count}, {time.perf_counter()-wall:.1f}s',flush=True);last_print=time.perf_counter()
    spectrum=sums.sum(axis=0)/durations.sum();windows=sums/durations[:,None]
    out=dict(parameters=p.tolist(),initial_state=start.tolist(),final_state=y.tolist(),spectrum=spectrum.tolist(),spectrum_windows=windows.tolist(),
             divergence=-2*(b+1),spectrum_sum_error=float(abs(spectrum.sum()+2*(b+1))),maximum_local_volume_error=max_sum_error,
             event_count=event_count,minimum_crossing_speed=min_speed if event_count else None,minimum_event_interval=min_separation if event_count else None,
             minimum_sampled_state=minimum.tolist(),maximum_sampled_state=maximum.tolist(),
             settings=dict(transient=transient,observation=observation,method=method,rtol=rtol,qr_interval=qr_interval,logistic_width=width),
             model='Continuous heat with ideal relays' if width<=0 else 'Explicit modified model: logistic optical transitions of common normalized width',
             poincare_section='x_minus=delta_minus, v_minus>0; 5D section (full 6 coordinates saved)' if width<=0 else 'Not sampled in this smooth-control run',
             wall_seconds=time.perf_counter()-wall)
    if save_prefix:
        path=Path(save_prefix);path.parent.mkdir(parents=True,exist_ok=True)
        path.with_suffix('.json').write_text(json.dumps(out,indent=2))
        np.savez_compressed(path.with_suffix('.npz'),times=times,states=states,sections=sections,event_records=event_records,spectrum_windows=windows)
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',default='results/light-operators/coupled-scout.npz');parser.add_argument('--index',type=int,default=3054);parser.add_argument('--transient',type=float,default=100);parser.add_argument('--observation',type=float,default=400);parser.add_argument('--method',default='DOP853');parser.add_argument('--rtol',type=float,default=1e-11);parser.add_argument('--qr-interval',type=float,default=.05);parser.add_argument('--width',type=float,default=0.);parser.add_argument('--output',default='results/light-operators/coupled-cpu')
    parser.add_argument('--periodic-control',action='store_true')
    args=parser.parse_args()
    if args.periodic_control:
        source='results/light-operators/coupled-structure.json';d=json.loads(Path(source).read_text())['source_parameter_periodic_control']
        p=np.array(d['parameters']);start=np.array(d['section_state']);args.transient=0.;args.observation=d['period'];args.source=source;args.index=None
    else:
        data=np.load(args.source);p=data['parameters'][args.index];start=data['final_states'][args.index]
        if p.shape==(5,):
            if args.width and args.width!=p[4]:raise ValueError('Width disagrees with the selected smooth-source row')
            args.width=float(p[4]);p=p[:4]
    result=integrate(p,start,args.transient,args.observation,args.method,args.rtol,args.qr_interval,args.output,args.width)
    result.update(source_file=args.source,source_index=args.index)
    if args.periodic_control:result['periodic_closure_error']=float(max(abs(np.array(result['final_state'])-start)))
    Path(args.output).with_suffix('.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
