"""Independent CPU full-spectrum and closed-orbit checks for gel candidates."""
import argparse
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from core_shell_gel import rhs_jac,ROOT,save


def cycle(file,index,output=None):
    old=np.load(file);p=old['parameters'][index];y=old['final_states'][index].copy()
    def event(t,y):return y[3]-.3
    event.direction=1
    pre=solve_ivp(lambda t,y:rhs_jac(y,p)[0],(0,100),y,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01,events=event)
    if not pre.success or len(pre.t_events[0])<3:raise RuntimeError('No regular crossing seed')
    seed=pre.y_events[0][-1,:3]
    event.terminal=True
    def return_map(v):
        initial=np.r_[v,.3]
        def aug(t,z):
            f,j=rhs_jac(z[:4],p);return np.r_[f,(j@z[4:].reshape(4,4)).ravel()]
        z0=np.r_[initial,np.eye(4).ravel()]
        skip=solve_ivp(aug,(0,1e-6),z0,method='DOP853',rtol=1e-12,atol=1e-14)
        sol=solve_ivp(aug,(1e-6,100),skip.y[:,-1],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01,events=event)
        if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('Return map failed')
        z=sol.y_events[0][0];f=rhs_jac(z[:4],p)[0];phi=z[4:].reshape(4,4)
        dp=(phi-np.outer(f,phi[3])/f[3])[:3,:3]
        return z[:3],dp,float(sol.t_events[0][0])
    fit=root(lambda v:return_map(v)[0]-v,seed,jac=lambda v:return_map(v)[1]-np.eye(3),tol=1e-10)
    v,dp,period=return_map(fit.x);initial=np.r_[fit.x,.3]
    closure=float(max(abs(v-fit.x)));eig=np.linalg.eigvals(dp)
    direct=solve_ivp(lambda t,y:rhs_jac(y,p)[0],(0,period),initial,method='Radau',jac=lambda t,y:rhs_jac(y,p)[1],rtol=1e-11,atol=1e-13,max_step=.01,dense_output=True)
    times=np.linspace(0,period,2001);states=direct.sol(times).T
    out=dict(parameters=p,section_state=initial,multipliers=eig,times=times,states=states,period=period,section_closure=closure,radau_closure=float(max(abs(direct.y[:,-1]-initial))),spectral_radius=float(max(abs(eig))),source=Path(file).name,index=index,method='Newton Poincare return map with DOP853 variational derivative; independent Radau closed-orbit check',rtol=1e-12,atol=1e-14,max_step=.01,section='phi=0.3, upward',scope='Small multipliers below floating-point resolution are not assigned accurate relative values')
    if closure>1e-8 or out['radau_closure']>1e-8:raise RuntimeError('Closed-orbit validation failed')
    save(ROOT/(output or f'{Path(file).stem}-cycle-{index}'),out);print({k:out[k] for k in ['period','section_closure','radau_closure','spectral_radius']},flush=True)


def check(file,index,transient=2000.,duration=4000.,max_step=.05,rtol=1e-11,qr_interval=.1,output=None):
    old=np.load(file);p=old['parameters'][index];y=old['final_states'][index].copy()
    rng=np.random.default_rng(2026103004);q,_=np.linalg.qr(rng.normal(size=(4,4)))
    def aug(t,z):
        f,j=rhs_jac(z[:4],p);return np.r_[f,(j@z[4:20].reshape(4,4)).ravel(),np.trace(j)]
    sums=np.zeros((4,4));div=np.zeros(4);ranges=np.array([[np.inf,-np.inf]]*4);trace=[];times=[]
    wall=time.perf_counter();t=0.;last=wall
    while t<transient+duration:
        window=-1 if t<transient else min(3,int((t-transient)/(duration/4)))
        edge=transient if window<0 else transient+(window+1)*duration/4
        dt=min(qr_interval,edge-t)
        if dt<1e-10:t=edge;continue
        sol=solve_ivp(aug,(0,dt),np.r_[y,q.ravel(),0.],method='DOP853',rtol=rtol,atol=rtol*.001,max_step=max_step)
        if not sol.success:raise RuntimeError(sol.message)
        y=sol.y[:4,-1];q,r=np.linalg.qr(sol.y[4:20,-1].reshape(4,4));t+=dt
        if y[2]<=0 or not 0<y[3]<1 or np.any(y[:2]<0):raise RuntimeError('Physical domain exit')
        if window>=0:
            sums[window]+=np.log(abs(np.diag(r)));div[window]+=sol.y[20,-1]
            ranges[:,0]=np.minimum(ranges[:,0],sol.y[:4].min(axis=1));ranges[:,1]=np.maximum(ranges[:,1],sol.y[:4].max(axis=1))
            trace.append(y.copy());times.append(t)
        if time.perf_counter()-last>20:print(f'CPU {t/(transient+duration):.1%}',flush=True);last=time.perf_counter()
    windows=sums/(duration/4);spectrum=windows.mean(axis=0)
    out=dict(parameters=p,initial_state=old['final_states'][index],final_state=y,spectrum=spectrum,spectrum_windows=windows,divergence_windows=div/(duration/4),ranges=ranges,times=np.array(times),states=np.array(trace),
        source=Path(file).name,index=index,transient=transient,duration=duration,max_step=max_step,rtol=rtol,atol=rtol*.001,qr_interval=qr_interval,method='DOP853 and QR, full tangent burn-in, constant scaled physical concentrations',
        sum_divergence_error=float(abs(spectrum.sum()-div.sum()/duration)),wall_seconds=time.perf_counter()-wall)
    save(ROOT/(output or f'{Path(file).stem}-cpu-{index}'),out);print({'spectrum':spectrum.tolist(),'divergence_error':out['sum_divergence_error']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',required=True);parser.add_argument('--index',type=int,required=True);parser.add_argument('--transient',type=float,default=2000);parser.add_argument('--duration',type=float,default=4000);parser.add_argument('--max-step',type=float,default=.05);parser.add_argument('--rtol',type=float,default=1e-11);parser.add_argument('--qr-interval',type=float,default=.1);parser.add_argument('--output');parser.add_argument('--cycle',action='store_true');a=parser.parse_args()
    if a.cycle:cycle(a.source,a.index,a.output)
    else:check(a.source,a.index,a.transient,a.duration,a.max_step,a.rtol,a.qr_interval,a.output)
