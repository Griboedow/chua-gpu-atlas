"""Reciprocity audit and conditional DC energy exclusion for a moving coil.

The reciprocal example is deliberately illustrative, not a fit or replacement
for the Dinavahi-Harris apparatus. Source equations are audited separately.
"""
import json
import argparse
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit


def magnetic(x):
    t=np.tanh(x);s=1-t*t
    return 1+.3*t,.3*s,.2*t,.2*s


def energy(state,i0=.8):
    x,v,current=state[:3];L,Lx,psi,psix=magnetic(x)
    return .5*v*v+.5*x*x+.5*L*current*current-i0*(L*current+psi)


def field(state,R,i0=.8,b=.3):
    x,v,current=state[:3];L,Lx,psi,psix=magnetic(x)
    force=.5*Lx*current*current+psix*current
    di=(R*(i0-current)-(Lx*current+psix)*v)/L
    return np.array([v,force-x-b*v,di])


def main():
    rng=np.random.default_rng(2026102945);identity_errors=[];gradient_errors=[]
    for _ in range(128):
        z=rng.uniform(-2,2,3);R=10**rng.uniform(-1,1);i0=.8;b=.3
        x,v,current=z;L,Lx,psi,psix=magnetic(x);f=field(z,R)
        gradient=np.array([x+.5*Lx*current**2-i0*(Lx*current+psix),v,L*(current-i0)])
        dissipation=-b*v*v-R*(current-i0)**2
        identity_errors.append(float(abs(gradient@f-dissipation)))
        h=1e-6
        finite=np.array([(energy(z+h*e)-energy(z-h*e))/(2*h) for e in np.eye(3)])
        gradient_errors.append(float(np.max(abs(finite-gradient))))
    source=[]
    for x in [.022,.03,.0475]:
        s=expit(250*(x-.03));mu=8-7*s;mu_x=-1750*s*(1-s)
        current=.02;c=.15/8
        source.append(dict(height_m=x,current_A=current,mu=mu,
                           force_current_derivative=7.6e-6*mu/x**4,
                           flux_position_derivative=c*mu_x*current,
                           scale_note='c=0.15/8 is an illustrative positive inductance scale; the opposite-sign mismatch does not depend on its magnitude.'))
    histories=[];summaries=[]
    for R in [.1,1.,10.]:
        initial=np.array([-1.5,1.,-.5]);W0=energy(initial)
        def rhs(t,z):return np.r_[field(z,R),.3*z[1]**2+R*(z[2]-.8)**2]
        times=np.linspace(0,200,4001)
        sol=solve_ivp(rhs,(0,200),np.r_[initial,0.],method='DOP853',rtol=1e-11,atol=1e-13,
                      max_step=.05,t_eval=times)
        if not sol.success:raise RuntimeError(sol.message)
        W=np.array([energy(z) for z in sol.y.T]);balance=W-W0+sol.y[3]
        histories.append(dict(R=R,time=times,state=sol.y[:3].T,W=W,dissipation=sol.y[3]))
        summaries.append(dict(R=R,initial=initial.tolist(),final_state=sol.y[:3,-1].tolist(),
                              maximum_energy_balance_error=float(np.max(abs(balance))),
                              maximum_sampled_energy_increase=float(np.max(np.diff(W))),
                              final_field_norm=float(np.linalg.norm(field(sol.y[:3,-1],R)))))
    result=dict(source_doi='10.1038/s41598-025-12529-y',source_equations=[2,3,5,9,13],
                reciprocity_condition='partial F/partial i = partial Phi/partial x for a reversible magnetic coenergy model',
                source_audit=source,source_conclusion='The printed Phi=L(x)i and F=a(x)i have opposite cross-derivative signs for positive current and decreasing permeability. They cannot jointly define a complete reciprocal magnetic transducer under these assumptions.',
                conditional_identity='W=kinetic+U+H(x,Phi)-i0*Phi; dW/dt=-b*v^2-R*(i-i0)^2, i0=Vs/R',
                illustrative_model='L=1+0.3*tanh(x), psi=0.2*tanh(x), U=x^2/2, m=1,b=0.3,i0=0.8; Phi=L*i+psi; F=Lprime*i^2/2+psiprime*i',
                is_device_calibration=False,algebra_samples=128,
                max_energy_derivative_identity_error=max(identity_errors),max_energy_gradient_fd_error=max(gradient_errors),
                trajectories=summaries,solver=dict(method='DOP853',duration=200,rtol=1e-11,atol=1e-13,max_step=.05),
                scope='Excludes non-equilibrium recurrent compact dynamics for the specified reciprocal fixed-DC fixed-resistance system; does not exclude optical feedback, other energy ports, hysteresis or a singular model.')
    if max(identity_errors)>1e-10 or max(gradient_errors)>1e-7 or max(a['maximum_energy_balance_error'] for a in summaries)>1e-8:
        raise RuntimeError('Energy audit failed')
    out=Path('results/meminductor-energy');out.mkdir(parents=True,exist_ok=True)
    (out/'audit.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(out/'illustrative-traces.npz',**{f'case{n}_{k}':v for n,h in enumerate(histories) for k,v in h.items()})
    print(json.dumps(result,indent=2))


def audit_impacts():
    """Check the energy jump only; this is not a calibrated hybrid device model."""
    rng=np.random.default_rng(2026102968);errors=[];increases=[]
    for _ in range(256):
        z=rng.uniform(-2,2,3);restitution=rng.uniform(0,1)
        after=z.copy();after[1]=-restitution*z[1]
        # At fixed x, flux continuity implies current continuity in this example.
        jump=energy(after)-energy(z)
        expected=-.5*(1-restitution**2)*z[1]**2
        errors.append(abs(jump-expected));increases.append(jump)
    result=dict(samples=256,seed=2026102968,
        assumptions='Fixed wall, continuous position and flux, v_after=-e*v_before, 0<=e<=1; m=1 in numeric example',
        identity='Delta W=-m*(1-e^2)*v_before^2/2',
        maximum_identity_error=float(max(errors)),maximum_energy_jump=float(max(increases)),
        is_device_calibration=False,
        scope='Passive impacts cannot restore the energy lost by the reciprocal fixed-DC flow; no claim about sticking, Zeno solutions or measured electrical impulses')
    result['passed']=bool(max(errors)<1e-12 and max(increases)<=1e-14)
    Path('results/meminductor-energy/passive-impact-audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    if not result['passed']:raise RuntimeError('Impact identity failed')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--impacts',action='store_true');args=parser.parse_args()
    if args.impacts:audit_impacts()
    else:main()
