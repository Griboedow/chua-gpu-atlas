"""Independent source-cycle closure on the exact equal-outflow mass plane."""
import argparse
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from urease import ROOT,SCALE,parameters,rate,save


def field_jac(y,p):
    s,d,h=y;v,vs,vh=rate(s,h,p)
    a=p[0]*SCALE;r=p[1];ks=p[5];kh=p[6];k=p[7];se=p[9]/SCALE;he=p[10]/SCALE
    if k!=ks:raise ValueError('Requires k=kS')
    f=np.array([-v+ks*(se-s),2*v-k*d+(kh-k)*h-kh*he,
                -a*(d+h)*h+r*(2*se-2*s-d-h)+kh*(he-h)])
    j=np.array([[-vs-ks,0,-vh],[2*vs,-k,2*vh+kh-k],[-2*r,-a*h-r,-a*(d+2*h)-r-kh]])
    return f,j


def main(source=None,index=0):
    if source:
        record=np.load(source);p=record['parameters'][index];ref=record['final_states'][index]
    else:
        p=parameters();ref=np.load(ROOT/'reference.npz')['states'][-1]
    y=np.array([ref[0],ref[2]-ref[1],ref[1]])
    level=.001
    def event(t,y):return y[2]-level
    event.direction=1
    horizon=max(3000.,6/p[5])
    pre=solve_ivp(lambda t,y:field_jac(y,p)[0],(0,horizon),y,method='Radau',jac=lambda t,y:field_jac(y,p)[1],rtol=1e-11,atol=1e-14,max_step=2.,events=event)
    if not pre.success or len(pre.t_events[0])<2:raise RuntimeError('No section seed')
    seed=pre.y_events[0][-1,:2]
    event.terminal=True
    def mapping(v):
        initial=np.r_[v,level]
        def augmented(t,z):
            f,j=field_jac(z[:3],p);return np.r_[f,(j@z[3:].reshape(3,3)).ravel()]
        z=np.r_[initial,np.eye(3).ravel()]
        skip=solve_ivp(augmented,(0,1e-5),z,method='Radau',rtol=1e-11,atol=1e-14)
        sol=solve_ivp(augmented,(1e-5,horizon),skip.y[:,-1],method='Radau',rtol=1e-11,atol=1e-14,max_step=2.,events=event)
        if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('No return')
        end=sol.y_events[0][0];f=field_jac(end[:3],p)[0];phi=end[3:].reshape(3,3)
        deriv=(phi-np.outer(f,phi[2])/f[2])[:2,:2]
        return end[:2],deriv,float(sol.t_events[0][0])
    fit=root(lambda v:mapping(v)[0]-v,seed,jac=lambda v:mapping(v)[1]-np.eye(2),tol=1e-9)
    end,deriv,period=mapping(fit.x);initial=np.r_[fit.x,level]
    check=solve_ivp(lambda t,y:field_jac(y,p)[0],(0,period),initial,method='BDF',jac=lambda t,y:field_jac(y,p)[1],rtol=1e-12,atol=1e-15,max_step=.2,dense_output=True)
    closure=float(np.max(abs(end-fit.x)));independent=float(np.max(abs(check.y[:,-1]-initial)))
    if not check.success or closure>1e-7 or independent>1e-6:raise RuntimeError('Cycle closure failed')
    times=np.linspace(0,period,4001);z=check.sol(times).T
    physical=np.column_stack([z[:,0],z[:,2],z[:,1]+z[:,2],2*p[9]/SCALE-2*z[:,0]-z[:,1]-z[:,2]])
    multipliers=np.linalg.eigvals(deriv);normal=float(np.exp(-p[5]*period))
    out=dict(parameters=p,section_state_transport=initial,times=times,states=physical,section_multipliers=multipliers,
        period_seconds=period,period_minutes=period/60,section_closure=closure,independent_bdf_closure=independent,
        exact_normal_multiplier=normal,maximum_transverse_multiplier=float(max(normal,max(abs(multipliers)))),
        method='Newton section map with implicit Radau variational derivative on exact N plane; independent BDF orbit',
        section='H=1e-7 M upward',state_scale_molar=SCALE,time_unit='s',source_period_minutes=[17.78,17.80] if source is None else None,
        source=Path(source).name if source else 'reference.npz',source_index=index if source else None,
        caveat='Unresolved extremely small multipliers have no claimed relative accuracy; numerical closure is not a rigorous orbit-existence proof')
    name=f'{Path(source).stem}-cycle-{index}' if source else 'source-cycle'
    save(ROOT/name,out);print({k:v for k,v in out.items() if not isinstance(v,np.ndarray)},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source');parser.add_argument('--index',type=int,default=0);a=parser.parse_args();main(a.source,a.index)
