"""Independent CPU full-spectrum and periodic-orbit diagnostics for NIFTE scouts.

Example: python research/check_nifte_spectrum.py --row 215
No GPU calls. Physical seconds; full analytic variational system with DOP853,
warm QR tangents, a separately propagated normalized flow control, and optional
Newton Poincare closure checked through dimensional mass balances with Radau.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from nifte import ROOT, COLUMNS, matrix, rhs_jac, mass_rhs, save


def cycle_check(p, endpoint, linear, loss, remaining):
    n = len(endpoint)
    keep = np.delete(np.arange(n), 1)

    def section(t, y):
        return y[1]

    section.direction = 1
    orbit = solve_ivp(lambda t, y: rhs_jac(y, p, loss, linear)[0], (0, 60), endpoint,
                      method="DOP853", rtol=1e-11, atol=1e-13, max_step=.03, events=section)
    if not orbit.success or len(orbit.t_events[0]) < 12:
        return {"status": "insufficient section returns", "section_count": len(orbit.t_events[0])}
    seed = orbit.y_events[0][-1]
    periods = np.diff(orbit.t_events[0])[-10:]
    section_spread = float(np.max(np.ptp(orbit.y_events[0][-10:], axis=0)))
    if section_spread > 1e-6 or np.ptp(periods) > 1e-5:
        return {"status": "no routine period-one closure evident", "last_periods": periods,
                "last_section_states": orbit.y_events[0][-10:], "section_spread": section_spread}
    section.terminal = True
    cache = {}

    def mapping(v):
        remaining()
        key = v.tobytes()
        if key in cache:
            return cache[key]
        initial = np.zeros(n)
        initial[keep] = v

        def augmented(t, z):
            f, j = rhs_jac(z[:n], p, loss, linear)
            return np.r_[f, (j @ z[n:].reshape(n, n)).ravel()]

        z = np.r_[initial, np.eye(n).ravel()]
        skip = solve_ivp(augmented, (0, 1e-6), z, method="DOP853", rtol=1e-12, atol=1e-14)
        sol = solve_ivp(augmented, (1e-6, max(10., 3 * float(periods.mean()))), skip.y[:, -1],
                        method="DOP853", rtol=1e-12, atol=1e-14, max_step=.015, events=section)
        if not sol.success or len(sol.t_events[0]) != 1:
            raise RuntimeError("No section return in Newton mapping")
        end = sol.y_events[0][0]
        f = rhs_jac(end[:n], p, loss, linear)[0]
        phi = end[n:].reshape(n, n)
        derivative = (phi - np.outer(f, phi[1]) / f[1])[np.ix_(keep, keep)]
        value = end[keep], derivative, float(sol.t_events[0][0]), phi
        cache[key] = value
        return value

    fit = root(lambda v: mapping(v)[0] - v, seed[keep],
               jac=lambda v: mapping(v)[1] - np.eye(n - 1), tol=1e-10, options={"maxfev": 30})
    end, derivative, period, phi = mapping(fit.x)
    initial = np.zeros(n)
    initial[keep] = fit.x
    closure = float(np.max(abs(end - fit.x)))
    independent = solve_ivp(lambda t, y: mass_rhs(y, p, loss), (0, period), initial,
                            method="Radau", rtol=1e-12, atol=1e-14, max_step=.002, dense_output=True)
    independent_error = float(np.max(abs(independent.y[:, -1] - initial)))
    multipliers = np.linalg.eigvals(derivative)
    times = np.linspace(0, period, 4001)
    return {"status": "closed periodic orbit" if closure < 1e-8 and independent_error < 1e-8 else "closure not verified",
            "solver_message": str(fit.message), "section_state": initial, "period_seconds": period,
            "section_closure": closure, "independent_mass_radau_closure": independent_error,
            "section_multipliers": multipliers, "full_multipliers": np.linalg.eigvals(phi),
            "maximum_transverse_multiplier": float(max(abs(multipliers))),
            "floquet_exponents_per_s": np.sort(np.r_[0., np.log(abs(multipliers)) / period])[::-1],
            "times": times, "states": independent.sol(times).T, "section_spread_before_newton": section_spread,
            "direct_monodromy_determinant": float(np.linalg.det(phi)),
            "liouville_determinant": float(np.exp(np.trace(linear[0]) * period)),
            "determinant_caveat": "Direct determinant is unreliable for unresolved strongly contracting modes, especially 9D.",
            "scope": "Numerical orbit closure and linear transverse stability, not a rigorous existence proof."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "random-5d.npz")
    parser.add_argument("--row", type=int, default=215)
    parser.add_argument("--warm", type=float, default=1000.)
    parser.add_argument("--observe", type=float, default=4000.)
    parser.add_argument("--qr", type=float, default=.5)
    parser.add_argument("--wall-limit", type=float, default=540.)
    parser.add_argument("--no-cycle", action="store_true")
    args = parser.parse_args()
    if args.warm < 0 or args.observe <= 0 or args.qr <= 0 or args.wall_limit <= 0:
        raise ValueError("Require nonnegative warming and positive observation, QR interval, and wall limit")
    window_steps = args.observe / (4 * args.qr)
    if window_steps < 1 or abs(window_steps - round(window_steps)) > 1e-8:
        raise ValueError("Each of four observation windows must have an integer number of QR intervals")
    started = time.perf_counter()
    record = np.load(args.input)
    resolved_input = args.input.resolve()
    try:
        input_label = resolved_input.relative_to(ROOT.parent.parent).as_posix()
    except ValueError:
        input_label = str(resolved_input)
    input_hash = hashlib.sha256(resolved_input.read_bytes()).hexdigest()
    p = record["parameters"][args.row].copy()
    initial = record["final_states"][args.row].copy()
    y = initial.copy()
    n = len(y)
    loss = n == 9
    linear = matrix(p, loss)
    trace = float(np.trace(linear[0]))
    q = np.eye(n)
    samples, diagonal_logs, flow_logs, flow_residuals, timestamps = [], [], [], [], []
    t = 0.
    total = args.warm + args.observe
    warm_steps = int(round(args.warm / args.qr))
    total_steps = int(round(total / args.qr))
    if abs(warm_steps * args.qr - args.warm) > 1e-8 or abs(total_steps * args.qr - total) > 1e-8:
        raise ValueError("Warm and total durations must be integer multiples of QR interval")

    def remaining():
        if time.perf_counter() - started > args.wall_limit:
            raise TimeoutError("Bounded CPU wall-clock limit reached")

    def augmented(tau, z):
        f, j = rhs_jac(z[:n], p, loss, linear)
        return np.r_[f, (j @ z[n:n + n*n].reshape(n, n)).ravel(), j @ z[n + n*n:]]

    termination = "completed"
    try:
        for step in range(total_steps):
            remaining()
            flow_start = rhs_jac(y, p, loss, linear)[0]
            speed_start = np.linalg.norm(flow_start)
            if speed_start < 1e-14:
                raise RuntimeError("Flow-neutral control undefined near equilibrium")
            z0 = np.r_[y, q.ravel(), flow_start / speed_start]
            sol = solve_ivp(augmented, (0, args.qr), z0, method="DOP853", rtol=1e-10,
                            atol=1e-12, max_step=.05)
            if not sol.success:
                raise RuntimeError(sol.message)
            y = sol.y[:n, -1]
            raw = sol.y[n:n + n*n, -1].reshape(n, n)
            q, r = np.linalg.qr(raw)
            signs = np.where(np.diag(r) < 0, -1., 1.)
            q *= signs
            logs = np.log(abs(np.diag(r)))
            flow_end = rhs_jac(y, p, loss, linear)[0]
            propagated = sol.y[n + n*n:, -1]
            flow_log = float(np.log(np.linalg.norm(propagated)))
            expected = flow_end / speed_start
            residual = float(np.linalg.norm(propagated - expected) / (1 + np.linalg.norm(expected)))
            t += args.qr
            if step >= warm_steps:
                timestamps.append(t - args.warm)
                samples.append(y.copy())
                diagonal_logs.append(logs)
                flow_logs.append(flow_log)
                flow_residuals.append(residual)
            if (step + 1) % max(1, int(round(1000 / args.qr))) == 0:
                print({"time_seconds": t, "wall_seconds": time.perf_counter() - started,
                       "spectrum_so_far": (np.sum(diagonal_logs, axis=0) / (t - args.warm)).tolist() if t > args.warm else None}, flush=True)
    except (TimeoutError, RuntimeError) as exc:
        termination = str(exc)

    logs = np.asarray(diagonal_logs).reshape(-1, n)
    observed = len(logs) * args.qr
    spectrum = logs.sum(axis=0) / observed if observed else np.full(n, np.nan)
    count_per_window = int(round(args.observe / (4 * args.qr)))
    windows = []
    flow_windows = []
    for i in range(4):
        chunk = logs[i * count_per_window:(i + 1) * count_per_window]
        if len(chunk) != count_per_window:
            break
        windows.append(chunk.sum(axis=0) / (len(chunk) * args.qr))
        flow_windows.append(float(np.sum(flow_logs[i * count_per_window:(i + 1) * count_per_window]) / (len(chunk) * args.qr)))
    output = dict(parameters=p, initial_saved_endpoint=initial, endpoint=y, parameter_columns=COLUMNS,
                  input_file=input_label, input_sha256=input_hash, row=args.row, model_dimension=n,
                  warm_seconds=args.warm, actual_warm_seconds=min(t, args.warm),
                  requested_observation_seconds=args.observe, actual_observation_seconds=observed,
                  qr_interval_seconds=args.qr, method="DOP853 analytic full variational system",
                  rtol=1e-10, atol=1e-12, max_step=.05, time_unit="physical seconds",
                  spectrum_per_s=spectrum, spectrum_windows_per_s=np.asarray(windows),
                  flow_neutral_window_rates_per_s=np.asarray(flow_windows),
                  spectrum_sum=float(spectrum.sum()) if observed else None, exact_divergence_per_s=trace,
                  spectrum_sum_minus_trace=float(spectrum.sum() - trace) if observed else None,
                  normalized_flow_propagation_max_scaled_error=float(max(flow_residuals, default=0.)),
                  flow_control_scope="Each QR interval starts with normalized actual flow; growth telescopes to log speed ratio. Its finite-window rate need not be exactly zero.",
                  qr_times=np.asarray(timestamps), qr_log_diagonals=logs,
                  qr_flow_log_norms=np.asarray(flow_logs), qr_sampled_states=np.asarray(samples),
                  termination=termination, wall_seconds=time.perf_counter() - started,
                  physical_scope="Native mathematical ODE; no vessel/interface operating limits imposed.",
                  interpretation="Finite-time full spectrum; compare neutral flow control and periodic-orbit Floquet spectrum before any chaos claim.")
    prefix = ROOT / f"cpu-row{args.row}-spectrum"
    save(prefix, output)
    print({"termination": termination, "spectrum": spectrum.tolist(), "sum_error": output["spectrum_sum_minus_trace"],
           "flow_windows": flow_windows, "wall_seconds": output["wall_seconds"]}, flush=True)
    if termination == "completed" and not args.no_cycle:
        try:
            remaining()
            cycle = cycle_check(p, y, linear, loss, remaining)
            cycle.update(parameters=p, input_file=input_label, input_sha256=input_hash,
                         row=args.row, wall_seconds=time.perf_counter() - started)
            save(ROOT / f"cpu-row{args.row}-cycle", cycle)
            print({k: v for k, v in cycle.items() if not isinstance(v, np.ndarray)}, flush=True)
            if "floquet_exponents_per_s" in cycle:
                print("Floquet exponents", cycle["floquet_exponents_per_s"].tolist(), flush=True)
            if cycle["status"] == "closed periodic orbit" and cycle["maximum_transverse_multiplier"] < 1:
                output["interpretation"] = "Verified stable periodic orbit; finite-window leading rates follow the neutral flow direction. No chaos supported for this saved endpoint."
                output["cycle_verification_record"] = f"cpu-row{args.row}-cycle.json"
                output["maximum_leading_minus_flow_window_rate"] = float(np.max(abs(np.asarray(windows)[:, 0] - flow_windows)))
                save(prefix, output)
        except (TimeoutError, RuntimeError) as exc:
            save(ROOT / f"cpu-row{args.row}-cycle", {"status": str(exc), "row": args.row})
    if termination != "completed":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
