"""Recover the published unstable five-state NIFTE cycle by CPU shooting.

Uses the literal nifte.rhs_jac model and checks the final orbit independently
with nifte.mass_rhs/Radau. The eigenvector seeds are our choices, not source
initial conditions. No physical parameter is fitted.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from nifte import COLUMNS, mass_rhs, matrix, parameters, rhs_jac

OUT=Path(__file__).resolve().parents[1]/'results/nifte/source-5d-unstable-cycle'
SCALE=np.array([.4,.005,.02,.2,.7])


class Shooting:
    def __init__(self,p,deadline,rtol=2e-10,atol=2e-12,max_step=.02):
        self.p=p;self.linear=matrix(p);self.deadline=deadline
        self.rtol=rtol;self.atol=atol;self.max_step=max_step
        self.cached=None;self.calls=0

    def flow(self,x0,period,dense=False):
        if time.perf_counter()>self.deadline:
            raise TimeoutError('CPU shooting wall-time budget reached')
        if not .3<period<2.0:
            raise ValueError('Newton period left the chosen search window')
        def fun(t,z):
            f,j=rhs_jac(z[:5],self.p,linear=self.linear)
            return np.r_[f,(j@z[5:].reshape(5,5)).ravel()]
        sol=solve_ivp(fun,(0,period),np.r_[x0,np.eye(5).ravel()],
                      method='DOP853',rtol=self.rtol,atol=self.atol,
                      max_step=self.max_step,dense_output=dense)
        if not sol.success:raise RuntimeError(sol.message)
        return sol

    def evaluate(self,z):
        if self.cached is not None and np.array_equal(z,self.cached[0]):
            return self.cached[1:]
        x0=z[:5]*SCALE;period=.88*np.exp(z[5])
        sol=self.flow(x0,period)
        x1=sol.y[:5,-1];mon=sol.y[5:,-1].reshape(5,5)
        f1=rhs_jac(x1,self.p,linear=self.linear)[0]
        residual=np.r_[(x1-x0)/SCALE,x0[1]/SCALE[1]]
        jac=np.zeros((6,6))
        jac[:5,:5]=(mon-np.eye(5))*SCALE[None,:]/SCALE[:,None]
        jac[:5,5]=period*f1/SCALE;jac[5,1]=1
        self.cached=(z.copy(),residual,jac);self.calls+=1
        return residual,jac

    def fun(self,z):return self.evaluate(z)[0]
    def jac(self,z):return self.evaluate(z)[1]


def complex_pairs(values):
    return [[float(z.real),float(z.imag)] for z in values]


def main(max_wall_seconds=480.):
    start=time.perf_counter();deadline=start+max_wall_seconds
    p=parameters(True);_,j0=rhs_jac(np.zeros(5),p)
    eig,vec=np.linalg.eig(j0)
    index=int(np.argmax(eig.imag));lam=eig[index];mode=vec[:,index]
    phase=np.pi/2-np.angle(mode[1]);base=(mode*np.exp(1j*phase)).real
    if rhs_jac(base,p)[0][1]<0:base=-base
    period0=2*np.pi/lam.imag
    attempts=[];found=None
    # A bounded deterministic list, centered on the high-frequency linear mode.
    for amp in [1.,2.,.5,3.,1.5,.25,4.,.75]:
        solver=Shooting(p,deadline)
        x0=amp*base;guess=np.r_[x0/SCALE,np.log(period0/.88)]
        attempt={'amplitude_factor':amp,'seed_state':x0.tolist(),'seed_period_s':period0}
        try:
            result=root(solver.fun,guess,jac=solver.jac,method='hybr',
                        options={'xtol':1e-10,'maxfev':100})
            xs=result.x[:5]*SCALE;period=.88*np.exp(result.x[5])
            residual=float(np.max(abs(solver.fun(result.x))))
            sol=solver.flow(xs,period,True)
            states=sol.sol(np.linspace(0,period,1001))[:5].T
            amplitude=float(np.max(np.ptp(states,axis=0)/SCALE))
            attempt.update(success=bool(result.success),message=str(result.message),
                           scaled_residual=residual,period_s=float(period),
                           scaled_peak_to_peak=amplitude,integrations=solver.calls)
            # Reject the trivial equilibrium and the lower-frequency source cycle.
            if residual<1e-8 and amplitude>.05 and .7<period<1.05:
                found=(result.x.copy(),attempt.copy());attempts.append(attempt);break
        except (TimeoutError,ValueError,RuntimeError) as exc:
            attempt['error']=str(exc)
        attempts.append(attempt)
        print(json.dumps(attempt),flush=True)
        if time.perf_counter()>deadline:break

    record={'model':'NIFTE NTP 5D','status':'not_recovered','attempts':attempts,
            'parameters':dict(zip(COLUMNS,p.tolist())),
            'source':'Wang et al. 2019, DOI 10.1186/s42480-019-0019-3, Table3',
            'source_target_period_s':.88,
            'source_target_multipliers':[[1.,0.],[.0348,0.],[.7958,0.],[-.3848,1.2078],[-.3848,-1.2078]],
            'chosen_seed_scope':'Our deterministic scaled high-frequency origin eigenvector; source initial state unavailable.',
            'origin_eigenvalue':complex_pairs([lam])[0],
            'time_unit':'physical seconds','no_parameters_fitted':True,
            'nontriviality_test':'scaled peak-to-peak > 0.05 and 0.7 < T < 1.05 s',
            'scientific_scope':'Recovery of a published unstable periodic orbit; no chaos or novelty claim.'}
    arrays={'parameters':p,'origin_eigenvalues':eig}
    if found is not None:
        guess,accepted=found
        fine=Shooting(p,deadline,rtol=3e-13,atol=3e-15,max_step=.005)
        result=root(fine.fun,guess,jac=fine.jac,method='hybr',options={'xtol':2e-11,'maxfev':50})
        xs=result.x[:5]*SCALE;period=float(.88*np.exp(result.x[5]))
        sol=fine.flow(xs,period,True)
        times=np.linspace(0,period,4001);states=sol.sol(times)[:5].T
        mon=sol.y[5:,-1].reshape(5,5);floquet=np.linalg.eigvals(mon)
        order=np.lexsort((floquet.imag,floquet.real));floquet=floquet[order]
        independent=solve_ivp(lambda t,y:mass_rhs(y,p),(0,period),xs,
                              method='Radau',rtol=2e-12,atol=2e-14,
                              max_step=.003,dense_output=True)
        ref_states=independent.sol(times).T
        def independent_variational(t,z):
            x=z[:5]
            eps=1e-25
            jac=np.column_stack([mass_rhs(x.astype(complex)+1j*eps*np.eye(5)[k],p).imag/eps for k in range(5)])
            return np.r_[mass_rhs(x,p),(jac@z[5:].reshape(5,5)).ravel()]
        independent_mon=solve_ivp(independent_variational,(0,period),np.r_[xs,np.eye(5).ravel()],
                                  method='Radau',rtol=2e-12,atol=2e-14,max_step=.003)
        ref_mon=independent_mon.y[5:,-1].reshape(5,5)
        ref_floquet=np.linalg.eigvals(ref_mon)
        ref_floquet=ref_floquet[np.lexsort((ref_floquet.imag,ref_floquet.real))]
        closure=float(np.max(abs(sol.y[:5,-1]-xs)))
        independent_closure=float(np.max(abs(independent.y[:,-1]-xs)))
        cross_error=float(np.max(abs(ref_states-states)))
        trivial=int(np.argmin(abs(floquet-1)))
        nontrivial=np.delete(floquet,trivial)
        # Discard endpoint neighborhoods: roundoff in the shooting section can
        # otherwise count the initial and final crossing as two distinct ones.
        signs=states[:,1]
        interior=(times[:-1]>period*1e-5)&(times[1:]<period*(1-1e-5))
        positive_crossings=int(np.count_nonzero(interior&(signs[:-1]<0)&(signs[1:]>=0)))
        negative_crossings=int(np.count_nonzero(interior&(signs[:-1]>0)&(signs[1:]<=0)))
        half_state=sol.sol(period/2)[:5]
        half_closure=float(np.max(abs(half_state-xs)))
        half_antisymmetry=float(np.max(abs(half_state+xs)))
        multiplier_cross_error=float(np.max(abs(floquet-ref_floquet)))
        determinant_expected=float(np.exp(np.trace(matrix(p)[0])*period))
        valid=bool(result.success and independent.success and independent_mon.success and closure<1e-9
                   and independent_closure<1e-9 and cross_error<1e-8
                   and abs(floquet[trivial]-1)<1e-7 and max(abs(nontrivial))>1.01
                   and multiplier_cross_error<1e-8 and positive_crossings==0
                   and negative_crossings==1 and half_closure>.05)
        record.update(status='recovered_and_independently_checked' if valid else 'candidate_validation_incomplete',
                      period_s=period,section_state=xs.tolist(),
                      floquet_multipliers=complex_pairs(floquet),
                      nontrivial_multiplier_moduli=abs(nontrivial).tolist(),
                      closure_max_abs=closure,independent_Radau_closure_max_abs=independent_closure,
                      independent_trajectory_max_abs_difference=cross_error,
                      independent_floquet_multipliers=complex_pairs(ref_floquet),
                      independent_multiplier_max_abs_difference=multiplier_cross_error,
                      independent_variational_method='Radau with complex-step derivative of independent dimensional mass_rhs',
                      determinant_from_monodromy=float(np.linalg.det(mon)),
                      determinant_from_constant_divergence=determinant_expected,
                      positive_d_crossings_in_open_period=positive_crossings,
                      negative_d_crossings_in_open_period=negative_crossings,
                      half_period_closure_max_abs=half_closure,
                      half_period_antisymmetry_max_abs=half_antisymmetry,
                      fundamental_period_check='One descending d=0 crossing and no ascending crossing in the open period; half-period state differs and equals minus the initial state.',
                      peak_to_peak=np.ptp(states,axis=0).tolist(),
                      shooting_settings={'method':'DOP853','rtol':fine.rtol,'atol':fine.atol,'max_step':fine.max_step},
                      independent_settings={'method':'Radau','rtol':2e-12,'atol':2e-14,'max_step':.003},
                      final_root_message=str(result.message),passed=valid)
        source=np.array([-.3848-1.2078j,-.3848+1.2078j,.0348,.7958,1.])
        record['published_multiplier_max_abs_difference']=float(np.max(abs(floquet-source)))
        record['source_comparison_scope']='Period and instability agree with the published control, but multipliers differ slightly from Table3. Cause unresolved; rounded source components were used without fitting.'
        arrays.update(section_state=xs,times=times,states=states,independent_states=ref_states,
                      monodromy=mon,floquet_multipliers=floquet,
                      independent_monodromy=ref_mon,independent_floquet_multipliers=ref_floquet)
    record['wall_seconds']=time.perf_counter()-start
    OUT.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(str(OUT)+'.npz',**arrays)
    record['array_shapes']={k:list(v.shape) for k,v in arrays.items()}
    Path(str(OUT)+'.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k not in ['parameters','attempts','array_shapes']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--max-wall-seconds',type=float,default=480.)
    args=parser.parse_args();main(args.max_wall_seconds)
