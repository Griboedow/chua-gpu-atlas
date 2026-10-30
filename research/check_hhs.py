"""Independent, wall-bounded full-spectrum CPU check of an exact HHS GPU scout.

The field is the fixed-reference, unsaturated, five-dimensional single-cell S2.
No continuum saturation or alternative reference concentration is introduced.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from hhs import ROOT, MREF, field_jac, original_field, save


def cycle_check(p, st, sy, deadline):
    if len(st) < 24:
        return {"status": "insufficient section returns for periodic closure"}
    order = None
    comparisons = {}
    for k in range(1, 9):
        err = float(np.max(abs(sy[-12:]-sy[-12-k:-k])))
        comparisons[str(k)] = err
        if err < 1e-7 and order is None:
            order = k
    if order is None:
        return {"status": "no repeated section orbit through return order eight",
                "section_return_max_absolute_differences": comparisons}
    period = float(np.median(st[-12:]-st[-12-order:-order]))
    cache = {}

    def mapping(v):
        if time.perf_counter() > deadline:
            raise TimeoutError("wall limit during periodic closure")
        key = v.tobytes()
        if key in cache:
            return cache[key]
        start, T = np.r_[MREF, v[:4]], float(v[4])
        if not 0 < T < 500:
            raise RuntimeError("periodic solver proposed inadmissible period")

        def augmented(t, z):
            if time.perf_counter() > deadline:
                raise TimeoutError("wall limit during periodic closure")
            f, j = field_jac(z[:5], p)
            return np.r_[f, (j@z[5:30].reshape(5, 5)).ravel(), np.trace(j)]

        sol = solve_ivp(augmented, (0, T), np.r_[start, np.eye(5).ravel(), 0.],
                        method="DOP853", rtol=1e-12, atol=1e-14, max_step=.015)
        if not sol.success:
            raise RuntimeError(sol.message)
        endpoint = sol.y[:5, -1]
        phi = sol.y[5:30, -1].reshape(5, 5)
        derivative = np.column_stack(((phi-np.eye(5))[:, 1:], field_jac(endpoint, p)[0]))
        cache[key] = endpoint-start, derivative, phi, float(sol.y[30, -1])
        return cache[key]

    fit = root(lambda v: mapping(v)[0], np.r_[sy[-1, 1:], period],
               jac=lambda v: mapping(v)[1], tol=1e-10, options={"maxfev": 25})
    residue, _, phi, integrated_trace = mapping(fit.x)
    initial, period = np.r_[MREF, fit.x[:4]], float(fit.x[4])
    alpha = .1
    original = np.r_[initial[:3], 1-alpha*initial[3:]]

    def original_rhs(t, z):
        if time.perf_counter() > deadline:
            raise TimeoutError("wall limit during independent original-coordinate Radau closure")
        return original_field(z, p, alpha)

    other = solve_ivp(original_rhs, (0, period), original, method="Radau",
                      rtol=1e-12, atol=1e-14, max_step=.01, dense_output=True)
    if not other.success:
        raise RuntimeError(other.message)
    closure = float(max(abs(residue)))
    check = float(max(abs(other.y[:, -1]-original)))
    flow = field_jac(initial, p)[0]
    section_map = (phi-np.outer(flow, phi[0])/flow[0])[1:, 1:]
    multipliers = np.linalg.eigvals(section_map)
    # Direct monodromy cannot resolve the extremely contracting filter mode.
    resolution = max(1e-12, 200*np.finfo(float).eps*np.linalg.norm(section_map, 2))
    resolved = abs(multipliers) > resolution
    times = np.linspace(0, period, 4001)
    physical = other.sol(times).T
    normalized = np.column_stack((physical[:, :3], (1-physical[:, 3:])/alpha))
    return {"status": "numerically closed cycle" if max(closure, check) < 1e-8 else "cycle closure unresolved",
            "period": period, "section_return_order": order, "initial_state": initial,
            "newton_closure_inf": closure, "independent_original_radau_closure_inf": check,
            "alpha_for_original_coordinate_check": alpha, "section_multipliers": multipliers,
            "full_multipliers": np.linalg.eigvals(phi), "monodromy": phi,
            "multiplier_absolute_resolution_floor": resolution, "multiplier_resolved_mask": resolved,
            "resolved_transverse_exponents": np.log(abs(multipliers[resolved]))/period,
            "maximum_transverse_multiplier": float(max(abs(multipliers))),
            "flow_identity_scaled_error": float(np.linalg.norm(phi@flow-flow)/(1+np.linalg.norm(flow))),
            "integrated_divergence": integrated_trace,
            "times": times, "states": normalized, "original_coordinate_states": physical,
            "minimum_length_at_alpha_point1": float(physical[:, 3:].min()),
            "caveat": "Numerical periodic-orbit evidence. Tiny direct monodromy multipliers are unresolved; no determinant or all-exponent accuracy claim from them."}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--row', type=int, default=1445)
    ap.add_argument('--input', type=Path, default=ROOT/'scout.npz')
    ap.add_argument('--warm', type=float, default=1000.)
    ap.add_argument('--observe', type=float, default=4000.)
    ap.add_argument('--qr', type=float, default=.5)
    ap.add_argument('--wall-limit', type=float, default=285.)
    args = ap.parse_args()
    began = time.perf_counter()
    deadline = began+args.wall_limit
    data = np.load(args.input)
    p, initial = data['parameters'][args.row].copy(), data['final_states'][args.row].copy()
    nw, no = int(round(args.warm/args.qr)), int(round(args.observe/args.qr))
    if nw*args.qr != args.warm or no*args.qr != args.observe or no % 4:
        raise ValueError('Durations must be QR multiples and observation must contain four equal windows')
    y, q, elapsed = initial.copy(), np.eye(5), 0.
    logs, flows, traces, qt, qy, tt, ty, st, sy = ([] for _ in range(9))
    propagation_error = 0.
    status, equilibrium = 'completed', None

    def augmented(t, z):
        if time.perf_counter() > deadline:
            raise TimeoutError('CPU wall limit reached; partial completed QR segments saved')
        f, j = field_jac(z[:5], p)
        out = np.empty(36)
        out[:5], out[5:30] = f, (j@z[5:30].reshape(5, 5)).ravel()
        out[30:35], out[35] = j@z[30:35], np.trace(j)
        return out

    def section(t, z):
        return z[0]-MREF
    section.direction = 1
    try:
        for i in range(nw+no):
            if time.perf_counter() > deadline:
                raise TimeoutError('CPU wall limit reached; partial completed QR segments saved')
            flow, _ = field_jac(y, p)
            speed = np.linalg.norm(flow)
            if speed < 1e-12:
                fit = root(lambda z: field_jac(z, p)[0], y,
                           jac=lambda z: field_jac(z, p)[1], tol=1e-11)
                ev = np.linalg.eigvals(field_jac(fit.x, p)[1])
                equilibrium = {"state": fit.x.tolist(), "residual_inf": float(max(abs(field_jac(fit.x, p)[0]))),
                               "eigenvalues_real_imag": [[float(v.real), float(v.imag)] for v in ev],
                               "physical_chemical_domain": bool(np.all(fit.x[:3] > 0))}
                status = 'flow below 1e-12; refined equilibrium instead of normalizing vanishing flow'
                break
            sol = solve_ivp(augmented, (0, args.qr), np.r_[y, q.ravel(), flow/speed, 0.],
                            method='DOP853', rtol=1e-11, atol=1e-13, max_step=.025, events=section)
            if not sol.success:
                raise RuntimeError(sol.message)
            z = sol.y[:, -1]
            y = z[:5]
            q, r = np.linalg.qr(z[5:30].reshape(5, 5))
            q *= np.where(np.diag(r) < 0, -1., 1.)
            if i >= nw:
                logs.append(np.log(abs(np.diag(r))))
                flows.append(float(np.log(np.linalg.norm(z[30:35]))))
                traces.append(float(z[35]))
                expected = field_jac(y, p)[0]/speed
                propagation_error = max(propagation_error, float(np.linalg.norm(z[30:35]-expected)/(1+np.linalg.norm(expected))))
                qt.append(elapsed+args.qr-args.warm); qy.append(y.copy())
                tt.append(elapsed-args.warm+sol.t[1:]); ty.append(sol.y[:5, 1:].T.copy())
                if len(sol.t_events[0]):
                    st.extend(elapsed-args.warm+sol.t_events[0]); sy.extend(sol.y_events[0][:, :5])
            elapsed += args.qr
            if (i+1) % int(round(1000/args.qr)) == 0:
                print(dict(elapsed=elapsed, wall_seconds=time.perf_counter()-began,
                           spectrum=(np.sum(logs, axis=0)/(elapsed-args.warm)).tolist() if logs else None), flush=True)
    except (TimeoutError, RuntimeError) as exc:
        status = str(exc)
    duration = len(logs)*args.qr
    array = np.asarray(logs).reshape(-1, 5)
    spectrum = array.sum(axis=0)/duration if duration else np.full(5, np.nan)
    windows, fw, tw = [], [], []
    width = no//4
    for k in range(4):
        lo, hi = k*width, (k+1)*width
        if len(logs) < hi:
            break
        windows.append(array[lo:hi].sum(axis=0)/(width*args.qr))
        fw.append(sum(flows[lo:hi])/(width*args.qr))
        tw.append(sum(traces[lo:hi])/(width*args.qr))
    sections, section_times = np.asarray(sy).reshape(-1, 5), np.asarray(st)
    trajectory = np.concatenate(ty) if ty else np.empty((0, 5))
    result = dict(status=status, parameters=p, parameter_columns=['mu', 'tauD', 'tauL'],
                  state_columns=['M', 'E', 'D', 'y', 'x'], initial_saved_endpoint=initial, endpoint=y,
                  row=args.row, input_file=args.input.resolve().relative_to(ROOT.parent.parent).as_posix(),
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  time_unit='source model time, not calibrated seconds',
                  requested_warm=args.warm, actual_warm=min(elapsed, args.warm),
                  requested_observation=args.observe, actual_observation=duration,
                  qr_interval=args.qr, method='DOP853 full analytic variational system', rtol=1e-11,
                  atol=1e-13, max_step=.025, full_spectrum=spectrum, spectrum_windows=np.asarray(windows),
                  neutral_flow_windows=np.asarray(fw), divergence_windows=np.asarray(tw),
                  mean_divergence=float(sum(traces)/duration) if duration else None,
                  spectrum_sum_minus_mean_divergence=float(spectrum.sum()-sum(traces)/duration) if duration else None,
                  maximum_leading_minus_flow_window_rate=float(np.max(abs(np.asarray(windows)[:, 0]-fw))) if windows else None,
                  flow_propagation_max_scaled_error=propagation_error,
                  qr_times=np.asarray(qt), qr_states=np.asarray(qy), qr_log_diagonals=array,
                  qr_flow_logs=np.asarray(flows), qr_trace_integrals=np.asarray(traces),
                  trajectory_times=np.concatenate(tt) if tt else np.empty(0), trajectory_states=trajectory,
                  section_times=section_times, section_states=sections,
                  section_definition='M=MREF upward; observation-relative time after warmup',
                  chemical_state_minima=trajectory[:, :3].min(axis=0) if len(trajectory) else np.empty(0),
                  minimum_length_at_alpha_point1=float((1-.1*trajectory[:, 3:]).min()) if len(trajectory) else None,
                  equilibrium_diagnostic=equilibrium, wall_seconds=time.perf_counter()-began,
                  interpretation='Finite-time diagnostics only; a positive leading exponent alone does not prove an attractor or novelty.')
    output = ROOT/f'cpu-row{args.row}-spectrum'
    save(output, result)
    print(dict(status=status, spectrum=spectrum.tolist(), windows=np.asarray(windows).tolist(),
               flow_windows=fw, trace_error=result['spectrum_sum_minus_mean_divergence'],
               section_count=len(st), wall_seconds=result['wall_seconds']), flush=True)
    if status == 'completed':
        try:
            cycle = cycle_check(p, section_times, sections, deadline)
        except (TimeoutError, RuntimeError) as exc:
            cycle = {'status': str(exc)}
        cycle.update(parameters=p, row=args.row, time_unit='source model time, not calibrated seconds')
        save(ROOT/f'cpu-row{args.row}-cycle-check', cycle)
        if cycle.get('status') == 'numerically closed cycle' and cycle['maximum_transverse_multiplier'] < 1:
            result['interpretation'] = 'Saved endpoint approaches a numerically closed stable cycle; no persistent chaos supported by this audit.'
            result['cycle_record'] = f'cpu-row{args.row}-cycle-check.json'
            save(output, result)
        print({k: v for k, v in cycle.items() if not isinstance(v, np.ndarray)}, flush=True)


if __name__ == '__main__':
    main()
