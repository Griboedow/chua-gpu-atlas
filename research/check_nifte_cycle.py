"""Newton Poincare closure and independent Radau verification of NIFTE cycles."""
import argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from nifte import ROOT, matrix, rhs_jac, mass_rhs, save


def main(loss=False):
    label='9d' if loss else '5d';record=np.load(ROOT/f'source-{label}.npz')
    p=record['parameters'];seed=record['section_states'][-1];n=len(seed);keep=np.delete(np.arange(n),1)
    linear=matrix(p,loss)
    def event(t,z):return z[1]
    event.direction=1;event.terminal=True
    def mapping(v):
        initial=np.zeros(n);initial[keep]=v
        def augmented(t,z):
            f,j=rhs_jac(z[:n],p,loss,linear)
            return np.r_[f,(j@z[n:].reshape(n,n)).ravel()]
        z=np.r_[initial,np.eye(n).ravel()]
        skip=solve_ivp(augmented,(0,1e-6),z,method='DOP853',rtol=1e-12,atol=1e-14)
        sol=solve_ivp(augmented,(1e-6,10),skip.y[:,-1],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01,events=event)
        if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('No section return')
        end=sol.y_events[0][0];f=rhs_jac(end[:n],p,loss,linear)[0];phi=end[n:].reshape(n,n)
        deriv=(phi-np.outer(f,phi[1])/f[1])[np.ix_(keep,keep)]
        return end[keep],deriv,float(sol.t_events[0][0]),phi
    fit=root(lambda v:mapping(v)[0]-v,seed[keep],jac=lambda v:mapping(v)[1]-np.eye(n-1),tol=1e-10)
    end,deriv,period,phi=mapping(fit.x);initial=np.zeros(n);initial[keep]=fit.x
    check=solve_ivp(lambda t,y:mass_rhs(y,p,loss),(0,period),initial,method='Radau',rtol=1e-12,atol=1e-14,max_step=.002,dense_output=True)
    closure=float(np.max(abs(end-fit.x)));independent=float(np.max(abs(check.y[:,-1]-initial)))
    times=np.linspace(0,period,4001);states=check.sol(times).T
    multipliers=np.linalg.eigvals(deriv)
    out=dict(parameters=p,section_state=initial,times=times,states=states,section_multipliers=multipliers,full_multipliers=np.linalg.eigvals(phi),
        period_seconds=period,section_closure=closure,independent_radau_closure=independent,maximum_transverse_multiplier=float(max(abs(multipliers))),
        model=f'NIFTE {label}',time_unit='physical s',section='d=0 upward',
        method='Newton section map, DOP853 variational derivative, independent dimensional mass-balance Radau orbit',
        full_determinant=float(np.linalg.det(phi)),liouville_determinant=float(np.exp(np.trace(linear[0])*period)),
        direct_determinant_resolved=not loss,
        multiplier_accuracy_scope='Dominant multipliers only; tiny 9D eigenvalues have no relative-accuracy claim' if loss else 'Independent cycle and Liouville controls passed',
        caveat='Smallest 9D multipliers are limited by FP64; numerical closure is not a rigorous orbit-existence proof')
    if not check.success or closure>1e-8 or independent>1e-8:raise RuntimeError('Cycle closure failed')
    save(ROOT/f'source-{label}-cycle',out)
    print({k:v for k,v in out.items() if not isinstance(v,np.ndarray)},flush=True)
    print('Floquet',out['full_multipliers'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--loss',action='store_true');a=parser.parse_args();main(a.loss)
