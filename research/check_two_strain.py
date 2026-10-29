"""Independent CPU full physical-coordinate spectra and Poincare traces."""
import argparse
import time
import numpy as np
from scipy.integrate import solve_ivp
from two_strain import ROOT, parameters, initial, encode, decode, log_rhs, balanced_jac, log_jac, rhs_jac, save


def check(p, q0, transient, duration, rtol, max_step, qr_interval, method, output, coordinates='log'):
    wall=time.perf_counter()
    q=q0.copy();basis=np.linalg.qr(np.random.default_rng(2026102971).normal(size=(5,5)))[0]
    logs=np.zeros((4,5));traces=np.zeros(4);physical_traces=np.zeros(4);post=None
    states=[];log_states=[];times=[];sections=[];section_times=[];convergence=[]
    minimums=np.full(5,np.inf);max_population=0.;t=0.;blocks=0
    section_s=(p[3]+p[0])/p[1]
    def field(t,z):
        j=log_jac(z[:5],p) if coordinates=='log' else balanced_jac(z[:5],p)
        return np.r_[log_rhs(z[:5],p),(j@z[5:30].reshape(5,5)).ravel(),np.trace(j),np.trace(rhs_jac(decode(z[:5],p),p)[1])]
    def event(t,z):return z[0]-section_s
    event.direction=1
    while t<transient+duration:
        window=-1 if t<transient else min(3,int((t-transient)/(duration/4)))
        boundary=transient if window<0 else transient+(window+1)*duration/4
        target=min(transient+duration,t+qr_interval,boundary)
        z=np.r_[q,basis.ravel(),0.,0.]
        sol=solve_ivp(field,(t,target),z,method=method,rtol=rtol,atol=rtol*.01,max_step=max_step,events=event)
        if not sol.success:raise RuntimeError(sol.message)
        q=sol.y[:5,-1];basis,upper=np.linalg.qr(sol.y[5:30,-1].reshape(5,5))
        diag=np.diag(upper)
        if np.any(diag==0):raise RuntimeError('Singular numerical QR')
        basis*=np.sign(diag)[None,:]
        if window>=0:
            logs[window]+=np.log(abs(diag));traces[window]+=sol.y[30,-1];physical_traces[window]+=sol.y[31,-1]
        y=decode(q,p);minimums=np.minimum(minimums,y);max_population=max(max_population,float(y.sum()))
        if min(y)<-1e-9 or y.sum()>1+1e-9:raise RuntimeError('Population simplex failure')
        t=target;blocks+=1
        if t==transient:post=q.copy();print(dict(transient_complete=t,wall_seconds=time.perf_counter()-wall),flush=True)
        if window>=0:
            if blocks%2==0:times.append(t-transient);states.append(y);log_states.append(q.copy())
            if len(sol.t_events[0]):
                section_times.extend(sol.t_events[0]-transient);sections.extend([decode(z[:5],p) for z in sol.y_events[0]])
            if t==boundary:
                convergence.append((logs.sum(axis=0)/(t-transient)).tolist())
                print(dict(elapsed=t-transient,spectrum=convergence[-1],wall_seconds=time.perf_counter()-wall),flush=True)
    spectrum=logs.sum(axis=0)/duration;windows=logs/(duration/4)
    if post is None:post=q0.copy()
    coordinate_boundary=(q[1]+2*q[2]-post[1]-2*post[2])/duration
    out=dict(parameters=p,initial_log_state=q0,post_transient_log_state=post,final_log_state=q,
        spectrum=spectrum,windows=windows,trace_integrals=traces,physical_trace_integrals=physical_traces,convergence=np.array(convergence),
        times=np.array(times),states=np.array(states),log_states=np.array(log_states),sections=np.array(sections),section_times=np.array(section_times),
        settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol*.01,max_step=max_step,qr_interval=qr_interval,method=method,coordinates=coordinates),
        leading_exponent=float(spectrum[0]),sum_rule_error=float(abs(spectrum.sum()-traces.sum()/duration)),
        min_sampled_compartments=minimums.tolist(),max_sampled_population_sum=max_population,
        wall_seconds=time.perf_counter()-wall,section='S=(mu+gamma1)/beta1, upward crossing; four-dimensional section shown through projections',
        tangent_coordinates=coordinates,
        trace_coordinate_boundary=float(coordinate_boundary),
        trace_coordinate_identity_error=float(abs((physical_traces.sum()-traces.sum())/duration-coordinate_boundary)) if coordinates=='log' else None,
        coordinate_caveat='Log full spectrum is explicitly a different numerical metric. Equality of asymptotic exponents requires a bounded or subexponentially conditioned transformation on the limiting orbit; finite-time spectra can differ. Original-coordinate leading GPU exponent is a separate check.',
        scope='Finite-time numerical diagnosis; no claim of novelty or biologically calibrated chaos')
    save(output,out)
    print(dict(output=str(output),spectrum=spectrum.tolist(),sum_rule_error=out['sum_rule_error'],wall_seconds=out['wall_seconds']),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--index',type=int);parser.add_argument('--r2',type=float,default=1.58)
    parser.add_argument('--start',choices=['invasion','coexistence'],default='coexistence');parser.add_argument('--source',default=str(ROOT/'slice-scout.npz'))
    parser.add_argument('--transient',type=float,default=20000);parser.add_argument('--duration',type=float,default=80000)
    parser.add_argument('--rtol',type=float,default=1e-11);parser.add_argument('--max-step',type=float,default=.5);parser.add_argument('--qr',type=float,default=5.)
    parser.add_argument('--method',choices=['DOP853','Radau'],default='DOP853');parser.add_argument('--coordinates',choices=['log','physical'],default='log');parser.add_argument('--output',required=True);args=parser.parse_args()
    if args.index is not None:
        data=np.load(args.source);p=data['parameters'][args.index];q=data['final_log_states'][args.index]
    else:p=parameters(args.r2);q=encode(initial(p,args.start),p)
    check(p,q,args.transient,args.duration,args.rtol,args.max_step,args.qr,args.method,ROOT/args.output,args.coordinates)
