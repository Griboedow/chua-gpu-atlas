"""Bounded implicit follow-up of the one unresolved gel scout, source row 2882.

Exact four-state source equations and saved parameters are retained. Only new
stiff-followup-2882 files are written; original GPU results remain unchanged.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from core_shell_gel import ROOT, COLUMNS, HC, K1, K2, KA, PMIN, PMAX, rhs_jac, save


def response_second(x, n):
    a, b, k = K2/HC, HC/K1, KA/HC
    d = x+a+b*x*x
    u = x/d
    up = (a-b*x*x)/d**2
    upp = -2*b*x/d**2-2*(a-b*x*x)*(1+2*b*x)/d**3
    xn, kn = x**n, k**n
    g = PMIN+(PMAX-PMIN)*xn/(xn+kn)
    gp = (PMAX-PMIN)*n*kn*x**(n-1)/(xn+kn)**2
    gpp = (PMAX-PMIN)*n*kn*((n-1)*x**(n-2)/(xn+kn)**2-2*n*x**(2*n-2)/(xn+kn)**3)
    return upp*g+2*up*gp+u*gpp


def jacobian_derivatives(s, p):
    from core_shell_gel import response
    G, U, x, phi = s
    A, C, gamma, rho, B, eps, lam, n = p
    f, fp = response(x, n)
    fpp = response_second(x, n)
    h = np.zeros((4, 4, 4))  # derivative state, field row, Jacobian column
    h[3, 0, 3] = lam*B*B*A*np.exp(-B*phi)
    h[2, 1, 1] = -rho*fp; h[1, 1, 2] = -rho*fp; h[2, 1, 2] = -rho*fpp*U
    h[2, 2, 1] = -fp; h[1, 2, 2] = -fp; h[2, 2, 2] = -fpp*U
    h[3, 3, 3] = -150*(phi-.3)/eps
    return h


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--duration', type=float, default=200.)
    ap.add_argument('--wall-limit', type=float, default=540.)
    args = ap.parse_args()
    began = time.perf_counter()
    deadline = began+args.wall_limit
    source = ROOT/'scout-refined.npz'
    old = np.load(source)
    row = int(np.flatnonzero(old['source_indices'] == 2882)[0])
    p, initial = old['parameters'][row].copy(), old['final_states'][row].copy()
    elapsed, y = 0., initial.copy()
    ts, ys, st, sy = [], [], [], []
    nfev, njev, nlu = 0, 0, 0
    status = 'completed requested state continuation'

    def wall_check():
        if time.perf_counter() > deadline:
            raise TimeoutError('wall budget exhausted; completed segments retained, classification unresolved')

    def fun(t, z):
        wall_check()
        return rhs_jac(z, p)[0]

    def jac(t, z):
        return rhs_jac(z, p)[1]

    def section(t, z):
        return z[3]-.3
    section.direction = 1
    try:
        while elapsed < args.duration:
            length = min(20., args.duration-elapsed)
            sol = solve_ivp(fun, (0, length), y, method='Radau', jac=jac,
                            rtol=1e-10, atol=1e-13, max_step=.1, events=section, dense_output=True)
            if not sol.success:
                raise RuntimeError(sol.message)
            tt = np.linspace(0, length, int(np.ceil(length/.01))+1)[1:]
            states = sol.sol(tt).T
            # Include adaptive solver nodes as well as plotted samples.
            if np.min(states[:, :3]) <= 0 or np.min(sol.y[:3]) <= 0 or \
                    min(states[:, 3].min(), sol.y[3].min()) <= 0 or max(states[:, 3].max(), sol.y[3].max()) >= 1:
                raise RuntimeError('physical domain exit; no clipping applied')
            ts.append(elapsed+tt); ys.append(states)
            st.extend(elapsed+sol.t_events[0]); sy.extend(sol.y_events[0])
            y = sol.y[:, -1]
            elapsed += length
            nfev += sol.nfev; njev += sol.njev; nlu += sol.nlu
            print(dict(elapsed=elapsed, section_returns=len(st), wall_seconds=time.perf_counter()-began), flush=True)
            if elapsed >= 100 and len(sy) >= 20:
                ss = np.asarray(sy)
                diff = np.max(abs(ss[-12:]-ss[-13:-1])/(1+abs(ss[-12:])))
                if diff < 1e-8:
                    status = 'stopped after >=100 units with 12 repeat section returns; proceeding to periodic closure'
                    break
    except (TimeoutError, RuntimeError, ValueError) as exc:
        status = str(exc)
    times = np.concatenate(ts) if ts else np.empty(0)
    states = np.concatenate(ys) if ys else np.empty((0, 4))
    section_times, section_states = np.asarray(st), np.asarray(sy).reshape(-1, 4)
    result = dict(source_row=2882, refined_row=row, parameters=p, parameter_columns=COLUMNS,
                  initial_saved_endpoint=initial, final_state=y,
                  source_file='results/core-shell-gel/scout-refined.npz',
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  previous_GPU_status=int(old['status'][row]), previous_GPU_elapsed=float(old['elapsed'][row]),
                  original_GPU_interval_still_incomplete=True,
                  state_columns=['G', 'U', 'x', 'phi'], time_unit='time scaled by proton relaxation; not calibrated seconds',
                  continuation_status=status, requested_continuation=args.duration, actual_continuation=elapsed,
                  continuation_solver=dict(method='Radau', analytic_jacobian=True, rtol=1e-10, atol=1e-13, max_step=.1),
                  times=times, states=states, section_times=section_times, section_states=section_states,
                  section='phi=.3 upward', nfev=nfev, njev=njev, nlu=nlu,
                  observed_minima=states.min(axis=0) if len(states) else np.empty(0),
                  observed_maxima=states.max(axis=0) if len(states) else np.empty(0),
                  interpretation='Unresolved until periodic-orbit controls complete; this one-case follow-up does not classify other step-budget cases.',
                  wall_seconds=time.perf_counter()-began)
    path = ROOT/'stiff-followup-2882'
    save(path, result)
    try:
        if elapsed < 100 or len(st) < 20:
            raise RuntimeError('insufficient continuation or section returns for periodic closure')
        scales = np.maximum(states.max(axis=0), 1.)
        order, differences = None, {}
        for k in range(1, 9):
            err = float(np.max(abs((section_states[-12:]-section_states[-12-k:-k])/scales)))
            differences[str(k)] = err
            if order is None and err < 1e-7:
                order = k
        result['section_return_scaled_differences'] = differences
        if order is None:
            raise RuntimeError('no periodic section return through order eight; retain unresolved classification')
        seed = section_states[-1]
        Tguess = float(np.median(section_times[-12:]-section_times[-12-order:-order]))
        h = jacobian_derivatives(seed, p)
        fd = np.stack([(rhs_jac(seed+np.eye(4)[k]*1e-6, p)[1]-rhs_jac(seed-np.eye(4)[k]*1e-6, p)[1])/2e-6 for k in range(4)])
        herr = float(np.max(abs(h-fd)/(1+abs(h))))
        result['augmented_Hessian_scaled_finite_difference_error'] = herr
        if herr > 1e-6:
            raise RuntimeError('augmented Jacobian derivative validation failed')
        cache = {}

        def mapping(v):
            wall_check()
            key = v.tobytes()
            if key in cache:
                return cache[key]
            w0 = np.r_[v[:3], .3/scales[3]]
            T = float(v[3])
            if not 0 < T < 200:
                raise RuntimeError('shooting proposed inadmissible period')

            def aug(t, a):
                wall_check()
                f, j = rhs_jac(a[:4]*scales, p)
                jw = j*scales[None, :]/scales[:, None]
                return np.r_[f/scales, (jw@a[4:20].reshape(4, 4)).ravel(), np.trace(j)]

            def augjac(t, a):
                s = a[:4]*scales
                j = rhs_jac(s, p)[1]
                jw = j*scales[None, :]/scales[:, None]
                hess = jacobian_derivatives(s, p)
                out = np.zeros((21, 21))
                out[:4, :4] = jw
                out[4:20, 4:20] = np.kron(jw, np.eye(4))
                phi = a[4:20].reshape(4, 4)
                for k in range(4):
                    hw = hess[k]*scales[k]*scales[None, :]/scales[:, None]
                    out[4:20, k] = (hw@phi).ravel()
                    out[20, k] = np.trace(hw)
                return out

            sol = solve_ivp(aug, (0, T), np.r_[w0, np.eye(4).ravel(), 0.], method='Radau',
                            jac=augjac, rtol=2e-11, atol=2e-13, max_step=.03)
            if not sol.success:
                raise RuntimeError(sol.message)
            end = sol.y[:4, -1]
            phi = sol.y[4:20, -1].reshape(4, 4)
            derivative = np.column_stack(((phi-np.eye(4))[:, :3], rhs_jac(end*scales, p)[0]/scales))
            cache[key] = end-w0, derivative, phi, float(sol.y[-1, -1])
            return cache[key]

        fitted = root(lambda v: mapping(v)[0], np.r_[seed[:3]/scales[:3], Tguess],
                      jac=lambda v: mapping(v)[1], tol=1e-10, options={'maxfev': 20})
        residue, _, monodromy, trace = mapping(fitted.x)
        start = np.r_[fitted.x[:3]*scales[:3], .3]
        period = float(fitted.x[3])
        other = solve_ivp(fun, (0, period), start, method='BDF', jac=jac,
                          rtol=1e-11, atol=1e-13, max_step=.01, dense_output=True)
        if not other.success:
            raise RuntimeError(other.message)
        raw_newton, raw_bdf = residue*scales, other.y[:, -1]-start
        local_scale = np.maximum(1., abs(start))
        closure = float(np.max(abs(raw_newton)/local_scale))
        independent = float(np.max(abs(raw_bdf)/local_scale))
        flow = rhs_jac(start, p)[0]/scales
        section_map = (monodromy-np.outer(flow, monodromy[3])/flow[3])[:3, :3]
        multipliers = np.linalg.eigvals(section_map)
        floor = max(1e-8, 500*np.finfo(float).eps*np.linalg.norm(section_map, 2))
        resolved = abs(multipliers) > floor
        ct = np.linspace(0, period, 8001)
        cs = other.sol(ct).T
        result.update(cycle_status='numerically closed cycle' if max(closure, independent) < 1e-7 else 'periodic closure unresolved',
                      cycle_initial_state=start, cycle_period=period, cycle_times=ct, cycle_states=cs,
                      cycle_section_return_order=order, cycle_coordinate_scales=scales,
                      cycle_newton_raw_closure=raw_newton, cycle_BDF_raw_closure=raw_bdf,
                      cycle_newton_section_local_closure=closure, cycle_BDF_section_local_closure=independent,
                      cycle_section_local_scale=local_scale, cycle_scaled_monodromy=monodromy,
                      cycle_section_map=section_map, cycle_multipliers=multipliers,
                      cycle_multiplier_resolution_floor=floor, cycle_multiplier_resolved_mask=resolved,
                      cycle_resolved_transverse_exponents=np.log(abs(multipliers[resolved]))/period,
                      cycle_transverse_spectral_radius=float(max(abs(multipliers))),
                      cycle_section_map_two_norm=float(np.linalg.norm(section_map, 2)),
                      cycle_integrated_divergence=trace,
                      cycle_flow_identity_error=float(np.linalg.norm(monodromy@flow-flow)/(1+np.linalg.norm(flow))),
                      cycle_solver=dict(method='Radau exact full augmented Jacobian and Newton shooting', rtol=2e-11, atol=2e-13, max_step=.03),
                      independent_cycle_solver=dict(method='BDF analytic Jacobian', rtol=1e-11, atol=1e-13, max_step=.01),
                      tiny_multiplier_caveat='Below the stated absolute resolution floor, individual multipliers are unresolved. No direct determinant or accurate small-mode exponent claim.',
                      wall_seconds=time.perf_counter()-began)
        if result['cycle_status'] == 'numerically closed cycle' and result['cycle_transverse_spectral_radius'] < 1:
            result['interpretation'] = 'The saved unresolved source2882 endpoint approaches a numerically closed stable periodic orbit. No chaotic attractor supported for this trajectory. Original long GPU interval remains incomplete; other budget cases remain unclassified.'
    except (TimeoutError, RuntimeError, ValueError) as exc:
        result['cycle_status'] = str(exc)
        result['wall_seconds'] = time.perf_counter()-began
    save(path, result)
    print({k: v for k, v in result.items() if k in ['cycle_status', 'cycle_period', 'cycle_newton_section_local_closure',
          'cycle_BDF_section_local_closure', 'cycle_transverse_spectral_radius', 'cycle_section_map_two_norm', 'wall_seconds', 'interpretation']}, flush=True)


if __name__ == '__main__':
    main()
