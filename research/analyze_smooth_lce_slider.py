"""Equilibria, numerical periodic orbits and switching-width sensitivity.

This analyzes the explicitly modified logistic-switch model, not a calibrated
material response. Orbit closure and Floquet estimates are numerical evidence.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root
from scipy.special import expit

from smooth_lce_slider import smooth_field


BASE = Path('results/lce-slider')


def equilibrium(p, width):
    k, b, intensity, d = p
    q = intensity / (1 + intensity)
    def balance(u):
        light = expit(u / width)
        return u + intensity * light / (1 + intensity * light) - d
    u = brentq(balance, d-q, d, xtol=1e-18)
    state = np.array([u, 0., d-u])
    light = expit(u / width)
    a = 1 + intensity * light
    c = intensity * (1-state[2]) * light * (1-light) / width
    # det(lambda I-J)=(lambda+a)(lambda^2+b lambda+k)+k c.
    margin = b*(k+a*b+a*a)-k*c
    eig = np.linalg.eigvals(smooth_field(state, p, width)[1])
    return state, margin, eig, np.sqrt(k+a*b)


def periodic_orbit(p, width, start):
    k = p[0]; q = p[2]/(1+p[2]); scale = np.array([q, q*np.sqrt(k), q])
    eq = equilibrium(p, width)[0] / scale
    def field(t, y):
        return smooth_field(y*scale, p, width)[0] / scale
    def crossing(t, y):
        return y[0] - eq[0]
    crossing.direction = 1
    tail = solve_ivp(field, (0, 40), start/scale, method='DOP853',
                     rtol=1e-12, atol=1e-14, max_step=.005, events=crossing)
    if not tail.success:
        raise RuntimeError(tail.message)
    times = tail.t_events[0]; states = tail.y_events[0]
    if len(times) < 12:
        raise RuntimeError('Too few transverse crossings for the periodic control')
    candidates = [(lag, np.linalg.norm(states[-1]-states[-1-lag])) for lag in range(1,9)]
    valid = [(lag, err) for lag, err in candidates if err < 1e-7]
    if not valid:
        raise RuntimeError(f'No short periodic recurrence: {candidates}')
    lag, error = valid[0]
    period = times[-1]-times[-1-lag]
    initial = np.r_[states[-1], period]
    cached = {}
    def flow_and_tangent(t, y):
        f, j = smooth_field(y[:3]*scale, p, width)
        j = j*scale[None,:]/scale[:,None]
        return np.r_[f/scale, (j@y[3:].reshape(3,3)).ravel()]
    def evaluate(z):
        if 'z' in cached and np.array_equal(z, cached['z']):
            return cached['residual'], cached['jac'], cached['M']
        if z[3] <= 0:
            raise RuntimeError('Orbit solve attempted a nonpositive period')
        sol = solve_ivp(flow_and_tangent, (0,z[3]), np.r_[z[:3],np.eye(3).ravel()],
                        method='DOP853', rtol=1e-12, atol=1e-14, max_step=.005)
        if not sol.success:
            raise RuntimeError(sol.message)
        end = sol.y[:3,-1]; M = sol.y[3:,-1].reshape(3,3)
        residual = np.r_[end-z[:3], z[0]-eq[0]]
        jac = np.zeros((4,4));jac[:3,:3]=M-np.eye(3);jac[:3,3]=field(0,end);jac[3,0]=1
        cached.update(z=z.copy(),residual=residual,jac=jac,M=M)
        return residual,jac,M
    solved = root(lambda z:evaluate(z)[0], initial, jac=lambda z:evaluate(z)[1], tol=1e-10)
    residual, _, M = evaluate(solved.x)
    closure = float(np.linalg.norm(residual))
    if closure > 1e-8 or solved.x[3] < .5*period:
        raise RuntimeError(f'Unresolved nontrivial orbit: {solved.message}, {closure}')
    multipliers = np.linalg.eigvals(M)
    multipliers = multipliers[np.argsort(-abs(multipliers))]
    trace_t = np.linspace(0,solved.x[3],1001)
    trace = solve_ivp(field,(0,solved.x[3]),solved.x[:3],method='DOP853',rtol=1e-12,
                      atol=1e-14,max_step=.005,t_eval=trace_t)
    if not trace.success:
        raise RuntimeError(trace.message)
    result = dict(width=float(width),width_over_q=float(width/q),period=float(solved.x[3]),
                  section_returns_per_period=lag,initial_recurrence_error=float(error),
                  scaled_closure_error=closure,state=(solved.x[:3]*scale).tolist(),
                  floquet_multipliers=[dict(real=float(z.real),imag=float(z.imag),absolute=float(abs(z))) for z in multipliers],
                  method='DOP853 shooting with analytic variational matrix and equilibrium-displacement phase section',
                  rtol=1e-12,max_step=.005,interpretation='Numerical periodic orbit and stability evidence, not a validated proof')
    return result, np.column_stack((trace_t,trace.y.T*scale))


def main():
    summary = json.loads((BASE/'smooth/summary.json').read_text())
    p = np.array(summary['parameters']); q = p[2]/(1+p[2])
    path = BASE/'smooth-controls'; path.mkdir(exist_ok=True)
    ratios = np.geomspace(1e-6,.1,501)
    margins = np.array([equilibrium(p,q*r)[1] for r in ratios])
    roots = []
    for i in np.where(margins[:-1]*margins[1:]<0)[0]:
        ratio = brentq(lambda r:equilibrium(p,q*r)[1],ratios[i],ratios[i+1],xtol=1e-14)
        state,margin,eig,omega = equilibrium(p,q*ratio)
        roots.append(dict(width_over_q=ratio,width=float(q*ratio),state=state.tolist(),
                          angular_frequency=float(omega),linear_period=float(2*np.pi/omega),
                          eigenvalues=[dict(real=float(z.real),imag=float(z.imag)) for z in eig],
                          interpretation='Linear Hopf candidate; nonlinear criticality and a nondegeneracy proof are not supplied'))
    controls = []
    for entry in summary['results']:
        ratio = entry['width_over_q']
        if ratio in [1e-4,1e-3]:
            result,trace = periodic_orbit(p,q*ratio,np.array(entry['final_state']))
            controls.append(result)
            tag = f'{ratio:.0e}'
            np.savez_compressed(path/f'orbit-{tag}.npz',trace=trace,parameters=p)
            print(json.dumps(result),flush=True)
    state, margin, eig, _ = equilibrium(p,q*.01)
    steady = dict(width_over_q=.01,state=state.tolist(),eigenvalue_real_parts=eig.real.tolist(),
                  max_eigenvalue_imaginary=float(abs(eig.imag).max()),routh_hurwitz_margin=float(margin),
                  endpoint_error=float(np.linalg.norm(np.array(summary['results'][-1]['final_state'])-state)))
    output = dict(parameters=p.tolist(),changes_model=True,linear_hopf_candidates=roots,
                  periodic_orbits=controls,steady_control=steady,
                  stability_polynomial='(lambda+a)(lambda^2+b lambda+k)+k c; a=1+I S, c=I(1-phi) Sprime',
                  stable_iff='b*(k+a*b+a*a)>k*c (all coefficients positive)',
                  caveat='A local equilibrium criterion does not rule out coexisting remote attractors')
    (path/'summary.json').write_text(json.dumps(output,indent=2))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(2,2,figsize=(12,8),layout='constrained')
    fig.suptitle('Finite switching width changes the slider dynamics',fontsize=15)
    x = [r['width_over_q'] for r in summary['results']]
    vals = np.array([r['spectrum'] for r in summary['results']])
    for j,color in enumerate(['#cc5e29','#1593a5','#67539c']):
        axes[0,0].plot(x,vals[:,j],'o-',color=color,label=f'Exponent {j+1}')
    axes[0,0].set(xscale='log',xlabel='Switching scale w/q',ylabel='Exponent / dimensionless time',title='Five finite-time spectra, duration 800')
    axes[0,0].axhline(0,color='gray',lw=.8);axes[0,0].legend()
    for ratio,color in [(1e-4,'#1593a5'),(1e-3,'#cc5e29')]:
        trace = np.load(path/f'orbit-{ratio:.0e}.npz')['trace']
        axes[0,1].plot(trace[:,1]/q,trace[:,2]/(q*np.sqrt(p[0])),color=color,label=f'w/q={ratio:g}')
    axes[0,1].set(xlabel='u/q',ylabel='v/(q sqrt(k))',title='Closed attracting periodic orbits');axes[0,1].legend()
    eqmax = [np.linalg.eigvals(smooth_field(equilibrium(p,q*r)[0],p,q*r)[1]).real.max() for r in ratios]
    axes[1,0].plot(ratios,eqmax,color='#67539c');axes[1,0].axhline(0,color='gray',lw=.8)
    for hopf in roots:axes[1,0].axvline(hopf['width_over_q'],ls='--',color='black',lw=.8)
    axes[1,0].set(xscale='log',xlabel='Switching scale w/q',ylabel='Largest real part of eigenvalues',title='Unique equilibrium: local stability')
    for ratio,color in [(1e-6,'#1593a5'),(1e-5,'#cc5e29')]:
        hist = np.load(BASE/f'smooth/width-{ratio:.0e}.npz')['history']
        axes[1,1].plot(hist[:,0],hist[:,1],color=color,label=f'w/q={ratio:g}')
    refined = BASE/'smooth-refined/width-1e-05.npz'
    if refined.exists():
        hist=np.load(refined)['history'];axes[1,1].plot(hist[:,0],hist[:,1],color='black',label='1e-5, tighter / longer')
    axes[1,1].set(xlabel='Observation time',ylabel='Leading finite-time exponent',title='Chaos survives a narrow smooth transition');axes[1,1].legend()
    fig.savefig('docs/assets/lce_slider_switch_width.png',dpi=170);plt.close(fig)
    print(json.dumps(output,indent=2))


if __name__ == '__main__':
    main()
