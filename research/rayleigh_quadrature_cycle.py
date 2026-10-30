"""Numerical shooting control for the weak-coupling quadrature saddle cycle.

Exact 4D quadratic drag; no attraction or novelty claim. Canonical sigma.
"""
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from rayleigh_pair import ROOT, field_jac, modal_energy, save


def main():
    started = time.perf_counter()
    deadline = started+140
    p = np.array([.0331323264448661, .002])
    mu, kappa = p
    amplitude = 3*np.pi/8
    omega = np.sqrt(1+kappa)
    ratio = kappa/mu
    # Original sum coupling: positive x2 at x1=0,v1>0 requires A1<A2.
    seeds = [np.array([amplitude*(1-ratio), amplitude*(1+ratio)/omega, 0., 2*np.pi/omega]),
             np.array([amplitude, amplitude/omega, 0., 2*np.pi/omega]),
             np.array([amplitude*(1+ratio), -amplitude*(1-ratio)/omega, 0., 2*np.pi/omega])]
    attempts = []
    accepted = None

    def augmented(t, z):
        f, j = field_jac(z[:4], p)
        out = np.empty(21)
        out[:4] = f;out[4:20] = (j@z[4:20].reshape(4,4)).ravel();out[20] = np.trace(j)
        return out

    for seed in seeds:
        cache = {}

        def shooting(v):
            if time.perf_counter() > deadline:
                raise TimeoutError('Bounded shooting budget reached')
            key = v.tobytes()
            if key in cache:
                return cache[key]
            if not 3 < v[3] < 12 or np.linalg.norm(v[:3]) > 10:
                raise RuntimeError('Newton iterate left the prescribed fundamental-period neighborhood')
            initial = np.r_[0., v[:3]]
            sol = solve_ivp(augmented, (0, v[3]), np.r_[initial,np.eye(4).ravel(),0.],
                            method='DOP853',rtol=2e-13,atol=2e-15,max_step=.015)
            if not sol.success:
                raise RuntimeError(sol.message)
            end = sol.y[:,-1];phi = end[4:20].reshape(4,4)
            derivative = np.column_stack(((phi-np.eye(4))[:,1:],field_jac(end[:4],p)[0]))
            cache[key] = end[:4]-initial,derivative,end
            return cache[key]

        try:
            fit = root(lambda v:shooting(v)[0],seed,jac=lambda v:shooting(v)[1],tol=1e-11,options={'maxfev':45})
            residue,_,end = shooting(fit.x)
            initial = np.r_[0.,fit.x[:3]]
            sym_distance = np.linalg.norm(initial[:2]-initial[2:])
            anti_distance = np.linalg.norm(initial[:2]+initial[2:])
            attempts.append({'seed':seed.tolist(),'solution':fit.x.tolist(),'solver_message':str(fit.message),
                             'closure':float(max(abs(residue))),'symmetric_distance':float(sym_distance),
                             'antisymmetric_distance':float(anti_distance)})
            if max(abs(residue)) < 1e-10 and min(sym_distance,anti_distance) > .1:
                accepted = initial,float(fit.x[3]),end
                break
        except (TimeoutError,RuntimeError) as exc:
            attempts.append({'seed':seed.tolist(),'failure':str(exc)})
    target = ROOT/'quadrature-cycle'
    if accepted is None:
        save(target,dict(status='No off-modal cycle verified',parameters=p,attempts=attempts,wall_seconds=time.perf_counter()-started))
        print('No off-modal cycle verified',attempts,flush=True)
        return
    initial,period,end = accepted
    phi = end[4:20].reshape(4,4)
    flow = field_jac(initial,p)[0]
    section_derivative = (phi-np.outer(flow,phi[0])/flow[0])[1:,1:]
    multipliers = np.linalg.eigvals(section_derivative)

    def modal_rhs(t,y):
        x,v,z,w=y;s,a=(x+z)/2,(x-z)/2
        gv,gw=(1-abs(v))*v,(1-abs(w))*w
        sa=-(1+2*kappa)*s+mu*(gv+gw)/2;aa=-a+mu*(gv-gw)/2
        return np.array([v,sa+aa,w,sa-aa])

    check = solve_ivp(modal_rhs,(0,period),initial,method='Radau',rtol=2e-12,atol=2e-14,max_step=.002,dense_output=True)
    times = np.linspace(0,period,6001);states=check.sol(times).T
    energies=modal_energy(states,p);fractions=np.min(energies,axis=1)/np.sum(energies,axis=1)
    angles=np.arctan2(-states[:,3]/omega,states[:,2])-np.arctan2(-states[:,1]/omega,states[:,0])
    full=np.linalg.eigvals(phi);flow_error=float(np.linalg.norm(phi@flow-flow)/(1+np.linalg.norm(flow)))
    data=dict(status='Numerically closed off-modal saddle cycle' if max(abs(multipliers))>1 else 'Numerically closed off-modal cycle',
              parameters=p,time_unit='canonical sigma',state_columns=['x1','v1','x2','v2'],attempts=attempts,
              initial_state=initial,period_sigma=period,newton_closure=float(max(abs(end[:4]-initial))),
              independent_modal_radau_closure=float(max(abs(check.y[:,-1]-initial))),
              section_multipliers=multipliers,full_multipliers=full,
              floquet_exponents=np.sort(np.r_[0.,np.log(abs(multipliers))/period])[::-1],
              floquet_exponents_parent_scaled_time=np.sort(np.r_[0.,np.log(abs(multipliers))/period])[::-1]/mu,
              maximum_transverse_multiplier=float(max(abs(multipliers))),
              monodromy=phi,section_derivative=section_derivative,flow_identity_scaled_error=flow_error,
              integrated_divergence=float(end[-1]),determinant=float(np.linalg.det(phi)),
              liouville_determinant=float(np.exp(end[-1])),
              times=times,states=states,modal_energies=energies,
              mean_minimum_modal_energy_fraction=float(np.mean(fractions)),
              kinematic_phase_difference_circular_mean=float(np.angle(np.mean(np.exp(1j*angles)))),
              phase_scope='Kinematic atan2 angles using the approximate common frequency; not exact asymptotic phases.',
              script='research/rayleigh_quadrature_cycle.py',wall_seconds=time.perf_counter()-started,
              scope='Unstable periodic-orbit control; not an attractor or novelty claim. Attribution to a scout requires a separate phase/proximity check.')
    save(target,data)
    print({k:v for k,v in data.items() if not isinstance(v,np.ndarray)},flush=True)
    print('Multipliers',multipliers,'Floquet exponents',data['floquet_exponents'],flush=True)


if __name__=='__main__':main()
