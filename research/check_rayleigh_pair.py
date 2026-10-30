"""Bounded independent CPU audit of an exact four-state Rayleigh-pair scout.

Canonical sigma, original canonical coordinates (x1,v1,x2,v2). Full warmed QR
spectrum, integrated divergence, normalized-flow control, trajectory and section.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from rayleigh_pair import ROOT, field_jac, modal_energy, save


def try_cycle(p, section_times, section_states, deadline):
    if len(section_times) < 24:
        return {"status": "insufficient section returns for cycle test"}
    order = None
    for m in range(1, 9):
        spread = float(np.max(abs(section_states[-12:]-section_states[-12-m:-m])))
        if spread < 1e-7:
            order = m
            break
    if order is None:
        return {"status": "no repeated section orbit through period eight", "tested_return_orders": list(range(1, 9))}
    seed = section_states[-1].copy()
    guess_period = float(np.median(section_times[-12:]-section_times[-12-order:-order]))
    cache = {}

    def mapping(v):
        if time.perf_counter() > deadline:
            raise TimeoutError("Wall cap before cycle closure")
        key = v.tobytes()
        if key in cache:
            return cache[key]
        initial = np.r_[0., v[:3]]
        period = v[3]
        if period <= 0:
            raise RuntimeError("Newton proposed nonpositive period")

        def aug(t, z):
            f, j = field_jac(z[:4], p)
            return np.r_[f, (j@z[4:].reshape(4, 4)).ravel()]

        sol = solve_ivp(aug, (0, period), np.r_[initial, np.eye(4).ravel()], method="DOP853",
                        rtol=1e-12, atol=1e-14, max_step=.015)
        end, phi = sol.y[:4, -1], sol.y[4:, -1].reshape(4, 4)
        derivative = np.column_stack(((phi-np.eye(4))[:, 1:], field_jac(end, p)[0]))
        cache[key] = end-initial, derivative, phi
        return cache[key]

    fit = root(lambda v: mapping(v)[0], np.r_[seed[1:], guess_period],
               jac=lambda v: mapping(v)[1], tol=1e-10, options={"maxfev": 25})
    residue, _, phi = mapping(fit.x)
    initial, period = np.r_[0., fit.x[:3]], float(fit.x[3])

    def modal_rhs(t, y):
        x, v, z, w = y
        mu, kappa = p
        s, a = (x+z)/2, (x-z)/2
        gv, gw = (1-abs(v))*v, (1-abs(w))*w
        sa, aa = -(1+2*kappa)*s+mu*(gv+gw)/2, -a+mu*(gv-gw)/2
        return [v, sa+aa, w, sa-aa]

    check = solve_ivp(modal_rhs, (0, period), initial, method="Radau", rtol=1e-12, atol=1e-14, max_step=.005, dense_output=True)
    closure = float(max(abs(residue)))
    independent = float(max(abs(check.y[:, -1]-initial)))
    flow = field_jac(initial, p)[0]
    section_derivative = (phi-np.outer(flow, phi[0])/flow[0])[1:, 1:]
    multipliers = np.linalg.eigvals(section_derivative)
    times = np.linspace(0, period, 4001)
    return {"status": "numerically closed cycle" if closure < 1e-8 and independent < 1e-8 else "cycle closure unresolved",
            "section_return_order": order, "period_sigma": period, "initial_state": initial,
            "newton_closure": closure, "independent_modal_radau_closure": independent,
            "section_multipliers": multipliers, "full_multipliers": np.linalg.eigvals(phi),
            "maximum_transverse_multiplier": float(max(abs(multipliers))),
            "floquet_exponents": np.sort(np.r_[0., np.log(abs(multipliers))/period])[::-1],
            "times": times, "states": check.sol(times).T,
            "caveat": "Numerical closure and Floquet check, not a rigorous orbit proof."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--row", type=int, default=1865)
    parser.add_argument("--input", type=Path, default=ROOT/"scout.npz")
    parser.add_argument("--warm", type=float, default=1000)
    parser.add_argument("--observe", type=float, default=4000)
    parser.add_argument("--qr", type=float, default=.5)
    parser.add_argument("--wall-limit", type=float, default=260)
    args = parser.parse_args()
    started = time.perf_counter()
    deadline = started+args.wall_limit
    record = np.load(args.input)
    p = record['parameters'][args.row].copy()
    initial = record['final_states'][args.row].copy()
    y, q = initial.copy(), np.eye(4)
    warm_steps = int(round(args.warm/args.qr))
    observe_steps = int(round(args.observe/args.qr))
    if warm_steps*args.qr != args.warm or observe_steps*args.qr != args.observe or observe_steps%4:
        raise ValueError("Durations must be QR multiples, with four equal observation windows")
    logs, flow_logs, traces, qr_states, qr_times, trajectory_t, trajectory_y, section_t, section_y = ([] for _ in range(9))
    flow_error = 0.
    t = 0.

    def augmented(tau, z):
        f, j = field_jac(z[:4], p)
        out = np.empty(25)
        out[:4] = f
        out[4:20] = (j@z[4:20].reshape(4, 4)).ravel()
        out[20:24] = j@z[20:24]
        out[24] = np.trace(j)
        return out

    def section(tau, z):
        return z[0]

    section.direction = 1
    termination = "completed"
    try:
        for i in range(warm_steps+observe_steps):
            if time.perf_counter() > deadline:
                raise TimeoutError("CPU wall budget exhausted")
            flow = field_jac(y, p)[0]
            speed = np.linalg.norm(flow)
            if speed < 1e-14:
                raise RuntimeError("Flow control undefined near equilibrium")
            z0 = np.r_[y, q.ravel(), flow/speed, 0.]
            sol = solve_ivp(augmented, (0, args.qr), z0, method="DOP853", rtol=1e-11,
                            atol=1e-13, max_step=.025, events=section)
            if not sol.success:
                raise RuntimeError(sol.message)
            end = sol.y[:, -1]
            y = end[:4]
            q, r = np.linalg.qr(end[4:20].reshape(4, 4))
            q *= np.where(np.diag(r) < 0, -1., 1.)
            if i >= warm_steps:
                logs.append(np.log(abs(np.diag(r))))
                flow_logs.append(float(np.log(np.linalg.norm(end[20:24]))))
                traces.append(float(end[24]))
                expected = field_jac(y, p)[0]/speed
                flow_error = max(flow_error, float(np.linalg.norm(end[20:24]-expected)/(1+np.linalg.norm(expected))))
                qr_states.append(y.copy());qr_times.append(t+args.qr-args.warm)
                trajectory_t.append(t-args.warm+sol.t[1:]);trajectory_y.append(sol.y[:4, 1:].T.copy())
                if len(sol.t_events[0]):
                    section_t.extend(t-args.warm+sol.t_events[0]);section_y.extend(sol.y_events[0][:, :4])
            t += args.qr
            if (i+1)%int(round(1000/args.qr)) == 0:
                print({"sigma": t, "wall_seconds": time.perf_counter()-started,
                       "spectrum": (np.sum(logs, axis=0)/(t-args.warm)).tolist() if t > args.warm else None}, flush=True)
    except (TimeoutError, RuntimeError) as exc:
        termination = str(exc)
    observed = len(logs)*args.qr
    log_array = np.asarray(logs).reshape(-1, 4)
    spectrum = log_array.sum(axis=0)/observed if observed else np.full(4, np.nan)
    windows, flow_windows, trace_windows = [], [], []
    window_steps = observe_steps//4
    for i in range(4):
        lo, hi = i*window_steps, (i+1)*window_steps
        if len(logs) < hi:
            break
        duration = window_steps*args.qr
        windows.append(log_array[lo:hi].sum(axis=0)/duration)
        flow_windows.append(sum(flow_logs[lo:hi])/duration)
        trace_windows.append(sum(traces[lo:hi])/duration)
    sections = np.asarray(section_y).reshape(-1, 4)
    section_times = np.asarray(section_t)
    input_label = args.input.resolve().relative_to(ROOT.parent.parent).as_posix()
    output = dict(parameters=p, parameter_columns=['mu', 'kappa'], state_columns=['x1', 'v1', 'x2', 'v2'],
                  initial_saved_endpoint=initial, endpoint=y, row=args.row, input_file=input_label,
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(), time_unit='canonical sigma',
                  warm_sigma=args.warm, actual_warm_sigma=min(t,args.warm), requested_observation_sigma=args.observe,
                  actual_observation_sigma=observed, qr_interval_sigma=args.qr, method='DOP853 full analytic variational system',
                  rtol=1e-11, atol=1e-13, max_step=.025, termination=termination,
                  spectrum=spectrum, spectrum_windows=np.asarray(windows), flow_neutral_windows=np.asarray(flow_windows),
                  divergence_windows=np.asarray(trace_windows), integrated_divergence=float(sum(traces)),
                  mean_divergence=float(sum(traces)/observed) if observed else None,
                  spectrum_sum_minus_mean_divergence=float(spectrum.sum()-sum(traces)/observed) if observed else None,
                  flow_propagation_max_scaled_error=flow_error, qr_times=np.asarray(qr_times), qr_states=np.asarray(qr_states),
                  qr_log_diagonals=log_array, qr_flow_logs=np.asarray(flow_logs), qr_trace_integrals=np.asarray(traces),
                  trajectory_times=np.concatenate(trajectory_t) if trajectory_t else np.empty(0),
                  trajectory_states=np.concatenate(trajectory_y) if trajectory_y else np.empty((0, 4)),
                  section_times=section_times, section_states=sections,
                  section_definition='x1=0 upward; times relative to observation start after warmup',
                  modal_energies=modal_energy(np.asarray(qr_states),p) if qr_states else np.empty((0,2)),
                  physical_scope='Canonical translated/scaled mechanical coordinates; no physical LCE operating-domain certificate.',
                  wall_seconds=time.perf_counter()-started,
                  interpretation='Finite-time spectrum and section evidence only; positive leading exponent does not by itself prove an attractor or novelty.')
    save(ROOT/f'cpu-row{args.row}-spectrum',output)
    print({"termination":termination,"spectrum":spectrum.tolist(),"trace_error":output['spectrum_sum_minus_mean_divergence'],
           "flow_windows":flow_windows,"section_count":len(section_t),"wall_seconds":output['wall_seconds']},flush=True)
    if termination == 'completed':
        try:
            cycle=try_cycle(p,section_times,sections,deadline)
        except (TimeoutError,RuntimeError) as exc:
            cycle={"status":str(exc)}
        cycle.update(parameters=p,row=args.row,time_unit='canonical sigma',input_file=input_label)
        if cycle.get('status') == 'numerically closed cycle':
            energies = modal_energy(cycle['states'], p)
            maxima = np.max(energies, axis=0)
            cycle['maximum_symmetric_modal_energy'] = float(maxima[0])
            cycle['maximum_antisymmetric_modal_energy'] = float(maxima[1])
            cycle['modal_plane'] = ('symmetric' if maxima[1] < 1e-12*max(maxima[0], 1e-30)
                                    else 'antisymmetric' if maxima[0] < 1e-12*max(maxima[1], 1e-30)
                                    else 'off both modal planes')
            if cycle['maximum_transverse_multiplier'] < 1:
                output['interpretation'] = f"Saved endpoint converged to a numerically verified stable cycle ({cycle['modal_plane']}); no persistent chaos supported by this audit."
                output['cycle_verification_record'] = f'cpu-row{args.row}-cycle-check.json'
                output['maximum_leading_minus_flow_window_rate'] = float(np.max(abs(np.asarray(windows)[:, 0]-flow_windows)))
                save(ROOT/f'cpu-row{args.row}-spectrum', output)
        save(ROOT/f'cpu-row{args.row}-cycle-check',cycle)
        print({k:v for k,v in cycle.items() if not isinstance(v,np.ndarray)},flush=True)


if __name__=='__main__':
    main()
