"""Bounded implicit CPU follow-up of the unfinished extreme HHS source row 3009.

Radau integrates the unchanged fixed-reference single-cell equations. Large D is
not an experimentally calibrated regime. A budget failure remains unresolved.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from hhs import ROOT, field_jac, original_field, save


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-row', type=int, default=3009)
    ap.add_argument('--input', type=Path, default=ROOT/'scout-refined.npz')
    ap.add_argument('--duration', type=float, default=300.)
    ap.add_argument('--wall-limit', type=float, default=225.)
    args = ap.parse_args()
    started = time.perf_counter()
    deadline = started+args.wall_limit
    data = np.load(args.input)
    row = int(np.flatnonzero(data['source_indices'] == args.source_row)[0])
    p = data['parameters'][row].copy()
    initial = data['final_states'][row].copy()
    endpoint = initial.copy()
    elapsed, status = 0., 'completed'
    times, states, st, sy = [], [], [], []
    nfev, njev, nlu = 0, 0, 0
    level = .4

    def check_wall():
        if time.perf_counter() > deadline:
            raise TimeoutError('CPU wall budget exhausted; completed implicit segments retained')

    def rhs(t, z):
        check_wall()
        return field_jac(z, p)[0]

    def jac(t, z):
        return field_jac(z, p)[1]

    def section(t, z):
        return z[0]-level
    section.direction = 1
    try:
        while elapsed < args.duration:
            length = min(25., args.duration-elapsed)
            sol = solve_ivp(rhs, (0, length), endpoint, method='Radau', jac=jac,
                            rtol=1e-10, atol=1e-13, max_step=.5, dense_output=True, events=section)
            if not sol.success:
                raise RuntimeError(sol.message)
            grid = np.linspace(0, length, int(np.ceil(length/.02))+1)[1:]
            times.append(elapsed+grid); states.append(sol.sol(grid).T)
            st.extend(elapsed+sol.t_events[0]); sy.extend(sol.y_events[0])
            endpoint = sol.y[:, -1]
            elapsed += length
            nfev += sol.nfev; njev += sol.njev; nlu += sol.nlu
            print(dict(elapsed=elapsed, section_returns=len(st),
                       wall_seconds=time.perf_counter()-started), flush=True)
    except (TimeoutError, RuntimeError) as exc:
        status = str(exc)
    tt = np.concatenate(times) if times else np.empty(0)
    zz = np.concatenate(states) if states else np.empty((0, 5))
    st, sy = np.asarray(st), np.asarray(sy).reshape(-1, 5)
    scales = np.maximum(np.max(abs(zz), axis=0), [1., 1., 1., 1., 1.]) if len(zz) else np.ones(5)
    record = dict(status=status, source_row=args.source_row, refined_row=row, parameters=p,
                  parameter_columns=['mu', 'tauD', 'tauL'], state_columns=['M', 'E', 'D', 'y', 'x'],
                  initial_saved_endpoint=initial, endpoint=endpoint,
                  source_elapsed_before_CPU=float(data['elapsed'][row]),
                  input_file=args.input.resolve().relative_to(ROOT.parent.parent).as_posix(),
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  requested_duration=args.duration, actual_duration=elapsed,
                  method='Radau with analytic state Jacobian', rtol=1e-10, atol=1e-13, max_step=.5,
                  times=tt, states=zz, section_times=st, section_states=sy,
                  section_definition='M=.4 upward; time begins at saved incomplete GPU endpoint',
                  state_minima=zz.min(axis=0) if len(zz) else np.empty(0),
                  state_maxima=zz.max(axis=0) if len(zz) else np.empty(0),
                  minimum_length_at_alpha_point1=float((1-.1*zz[:, 3:]).min()) if len(zz) else None,
                  nfev=nfev, njev=njev, nlu=nlu, wall_seconds=time.perf_counter()-started,
                  time_unit='source model time (tauR=1), not calibrated seconds',
                  physical_scope='Extreme mathematical parameters and huge D; not an experimentally calibrated regime.',
                  interpretation='State-only implicit continuation; no chaos exclusion without further diagnostics.')
    save(ROOT/f'cpu-row{args.source_row}-stiff-continuation', record)
    cycle = dict(status='periodic closure not attempted', parameters=p, source_row=args.source_row)
    try:
        if status != 'completed' or len(st) < 16:
            raise RuntimeError('Insufficient completed continuation or section returns')
        order = None
        differences = {}
        for k in range(1, 9):
            if len(st) < 8+k:
                continue
            error = float(np.max(abs((sy[-8:]-sy[-8-k:-k])/scales)))
            differences[str(k)] = error
            if error < 1e-6 and order is None:
                order = k
        cycle['section_return_scaled_differences'] = differences
        if order is None:
            raise RuntimeError('No repeating section through return order eight; classification unresolved')
        seed = sy[-1]
        period_guess = float(np.median(st[-8:]-st[-8-order:-order]))
        cache = {}

        def mapping(v):
            check_wall()
            key = v.tobytes()
            if key in cache:
                return cache[key]
            w0 = np.r_[level/scales[0], v[:4]]
            T = float(v[4])
            if not 0 < T < 500:
                raise RuntimeError('Newton proposed inadmissible period')

            def aug_rhs(t, q):
                check_wall()
                f, j = field_jac(q[:5]*scales, p)
                jw = j*scales[None, :]/scales[:, None]
                return np.r_[f/scales, (jw@q[5:30].reshape(5, 5)).ravel(), np.trace(j)]

            def aug_jac(t, q):
                # Exact augmented Jacobian avoids finite-difference overflow
                # in strongly contracted, effectively zero tangent columns.
                z = q[:5]*scales
                m, e, d, y, x = z
                mu, td, tl = p
                j = field_jac(z, p)[1]
                jw = j*scales[None, :]/scales[:, None]
                h = np.zeros((5, 5, 5))  # derivative index, field row, field column
                h[0, 0, 0] = 2*e; h[1, 0, 0] = 2*m; h[0, 0, 1] = 2*m
                h[0, 1, 0] = -2*e; h[1, 1, 0] = -2*m; h[0, 1, 1] = -2*m
                h[2, 1, 1] = -1.; h[1, 1, 2] = -1.
                h[4, 2, 2] = mu/td; h[2, 2, 4] = mu/td
                out = np.zeros((31, 31))
                out[:5, :5] = jw
                out[5:30, 5:30] = np.kron(jw, np.eye(5))
                phi = q[5:30].reshape(5, 5)
                for k in range(5):
                    hw = h[k]*scales[k]*scales[None, :]/scales[:, None]
                    out[5:30, k] = (hw@phi).ravel()
                    out[30, k] = np.trace(hw)
                return out

            sol = solve_ivp(aug_rhs, (0, T), np.r_[w0, np.eye(5).ravel(), 0.],
                            method='Radau', jac=aug_jac, rtol=2e-11, atol=2e-13, max_step=.2)
            if not sol.success:
                raise RuntimeError(sol.message)
            endpoint_w = sol.y[:5, -1]
            phi = sol.y[5:30, -1].reshape(5, 5)
            velocity_w = field_jac(endpoint_w*scales, p)[0]/scales
            derivative = np.column_stack(((phi-np.eye(5))[:, 1:], velocity_w))
            cache[key] = endpoint_w-w0, derivative, phi, float(sol.y[-1, -1])
            return cache[key]

        guess = np.r_[seed[1:]/scales[1:], period_guess]
        fit = root(lambda v: mapping(v)[0], guess, jac=lambda v: mapping(v)[1],
                   tol=2e-10, options={'maxfev': 20})
        residue, _, phi, divergence = mapping(fit.x)
        start = np.r_[level, fit.x[:4]*scales[1:]]
        period = float(fit.x[4])
        alpha = .1
        start_original = np.r_[start[:3], 1-alpha*start[3:]]

        def original_rhs(t, z):
            check_wall()
            return original_field(z, p, alpha)

        def original_jac(t, z):
            normalized = np.r_[z[:3], (1-z[3:])/alpha]
            transform = np.array([1., 1., 1., -alpha, -alpha])
            return field_jac(normalized, p)[1]*transform[:, None]/transform[None, :]

        # Independent coordinates and a second implicit integrator.
        other = solve_ivp(original_rhs, (0, period), start_original, method='BDF',
                          jac=original_jac, rtol=2e-11, atol=2e-13, max_step=.05, dense_output=True)
        if not other.success:
            raise RuntimeError(other.message)
        end_normalized = np.r_[other.y[:3, -1], (1-other.y[3:, -1])/alpha]
        independent = float(np.max(abs((end_normalized-start)/scales)))
        independent_componentwise = float(np.max(abs(end_normalized-start)/(1+abs(start))))
        newton_componentwise = float(np.max(abs(residue*scales)/(1+abs(start))))
        # A section-local check prevents the enormous D peak from hiding error
        # when D is small at the chosen return section.
        state_control = solve_ivp(rhs, (0, period), start, method='Radau', jac=jac,
                                  rtol=2e-11, atol=2e-13, max_step=.2)
        if not state_control.success:
            raise RuntimeError(state_control.message)
        state_residual = state_control.y[:, -1]-start
        local_scale = np.maximum(1., abs(start))
        state_local_error = float(np.max(abs(state_residual)/local_scale))
        bdf_local_error = float(np.max(abs(end_normalized-start)/local_scale))
        flow = field_jac(start, p)[0]/scales
        section_map = (phi-np.outer(flow, phi[0])/flow[0])[1:, 1:]
        multipliers = np.linalg.eigvals(section_map)
        resolution = max(1e-8, 500*np.finfo(float).eps*np.linalg.norm(section_map, 2))
        resolved = abs(multipliers) > resolution
        ct = np.linspace(0, period, 8001)
        physical = other.sol(ct).T
        normalized = np.column_stack((physical[:, :3], (1-physical[:, 3:])/alpha))
        closure = float(max(abs(residue)))
        cycle.update(status='numerically closed cycle' if max(closure, independent, independent_componentwise, newton_componentwise, state_local_error, bdf_local_error) < 1e-7 else 'periodic closure unresolved',
                     period=period, section_return_order=order, initial_state=start,
                     coordinate_scales=scales, scaled_newton_closure_inf=closure,
                     independent_original_BDF_scaled_closure_inf=independent,
                     newton_componentwise_relative_closure_inf=newton_componentwise,
                     independent_original_BDF_componentwise_relative_closure_inf=independent_componentwise,
                     componentwise_relative_scale_definition='1+abs(initial_state); additional max(1,abs(initial_state)) checks below',
                     state_only_Radau_closure_vector=state_residual,
                     original_BDF_closure_vector=end_normalized-start,
                     state_only_Radau_section_local_closure_inf=state_local_error,
                     original_BDF_section_local_closure_inf=bdf_local_error,
                     section_local_scale=local_scale,
                     implicit_solvers_successful=['Radau state continuation with analytic Jacobian',
                         'Radau scaled shooting with exact augmented Hessian Jacobian',
                         'BDF original coordinates with analytic transformed Jacobian',
                         'Radau state-only section-local closure check'],
                     warning_context='An initial successful finite-difference augmented-Jacobian run produced SciPy numerical-Jacobian overflow warnings in strongly contracted tangent columns. The recorded repeat uses the exact augmented Jacobian and analytic original-coordinate BDF Jacobian and completed without those warnings.',
                     alpha_for_original_coordinate_check=alpha,
                     section_multipliers=multipliers, full_multipliers=np.linalg.eigvals(phi),
                     scaled_monodromy=phi, multiplier_resolved_mask=resolved,
                     multiplier_absolute_resolution_floor=resolution,
                     resolved_transverse_exponents=np.log(abs(multipliers[resolved]))/period,
                     maximum_transverse_multiplier=float(max(abs(multipliers))),
                     flow_identity_scaled_error=float(np.linalg.norm(phi@flow-flow)/(1+np.linalg.norm(flow))),
                     integrated_divergence=divergence, times=ct, states=normalized,
                     original_coordinate_states=physical, wall_seconds=time.perf_counter()-started,
                     caveat='Scaled coordinates condition shooting only; exact field unchanged. Tiny multipliers unresolved. No direct determinant/all-exponent certification. Not an experimental validity claim.')
        if cycle['status'] == 'numerically closed cycle' and cycle['maximum_transverse_multiplier'] < 1:
            record['interpretation'] = 'The saved endpoint approaches a numerically closed stable cycle under implicit continuation. The original 10000-unit GPU interval remains incomplete; this local control does not complete that interval or exclude other attractors.'
            record['cycle_record'] = f'cpu-row{args.source_row}-cycle-check.json'
            save(ROOT/f'cpu-row{args.source_row}-stiff-continuation', record)
    except (TimeoutError, RuntimeError) as exc:
        cycle['status'] = str(exc)
        cycle['wall_seconds'] = time.perf_counter()-started
    save(ROOT/f'cpu-row{args.source_row}-cycle-check', cycle)
    print({k: v for k, v in cycle.items() if not isinstance(v, np.ndarray)}, flush=True)


if __name__ == '__main__':
    main()
