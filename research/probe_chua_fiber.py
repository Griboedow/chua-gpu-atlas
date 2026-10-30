"""A bounded, independent frozen-scalar-fiber diagnostic for scout row 2198.

This probes five already saved amplitudes, not a parameter survey.  It does not
classify the finite-alpha trajectory.  All dynamics here set epsilon=0.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad, simpson, solve_ivp
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/chua-scaling/scout.npz"
OUT = ROOT / "results/chua-scaling/fiber-probe"


def main():
    started = time.monotonic()
    data = np.load(SOURCE)
    row = 2198
    alpha, kappa, m0, m1 = map(float, data["parameters"][row])
    a = 1 + m1
    xe = (m1 - m0) / a
    omega = np.sqrt(kappa)
    period = 2 * np.pi / omega
    boundaries = data["boundary_states"][row]
    amplitudes = np.sqrt(boundaries[:, 1] ** 2 + boundaries[:, 2] ** 2 / kappa)
    deadline = started + 110
    arrays = {
        "source_boundary_states": boundaries,
        "source_boundary_original_times": data["boundary_original_times"][row],
        "amplitudes": amplitudes,
        "phase_time": np.linspace(0, period, 1001),
    }
    report = {
        "source": "results/chua-scaling/scout.npz",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_row": row,
        "source_parameters_alpha_kappa_m0_m1": [alpha, kappa, m0, m1],
        "scope": "epsilon=0 scalar fibers; no finite-alpha classification",
        "amplitude_selection": "all five saved observation-window boundaries, without interpolation",
        "forcing": "y=R*sin(sqrt(kappa)*s), Z=R*sqrt(kappa)*cos(sqrt(kappa)*s)",
        "section": "y=0 with positive y derivative",
        "period_fast_time": float(period),
        "outer_equilibrium": xe,
        "outer_affine_grazing_amplitude": float((xe - 1) * np.sqrt(a * a + kappa)),
        "inner_residence_contraction_threshold": float(a / (m1 - m0)),
        "settings": {"max_period_iterations": 70, "iteration_rtol": 1e-11,
                     "iteration_atol": 1e-12, "max_step_fraction_period": 1/200,
                     "final_rtol": 2e-12, "final_atol": 2e-13,
                     "stop_tolerance_relative": 2e-10, "wall_budget_seconds": 110},
        "probes": [],
        "shooting_controls": [],
    }

    for j, R in enumerate(amplitudes):
        R = float(R)

        def rhs(t, state):
            x = state[0]
            return [R * np.sin(omega * t) - a * x - (m0 - m1) * np.clip(x, -1, 1)]

        def lower(t, state):
            return state[0] + 1

        def upper(t, state):
            return state[0] - 1

        def averages(sol):
            switches = sorted(float(t) for events in sol.t_events for t in events)

            def integrand(t):
                y = R*np.sin(omega*t)
                return y*(float(sol.sol(t)[0])-y)

            integral, error = quad(integrand, 0, period, points=switches,
                                   epsabs=1e-10, epsrel=1e-10, limit=100)
            mesh = np.linspace(0, period, 4001)
            y = R*np.sin(omega*mesh)
            fine = simpson(y*(sol.sol(mesh)[0]-y), x=mesh)/period
            H = integral/period
            return {"average_H_y_times_x_minus_y": float(H),
                    "average_action_drift_original_time": float(kappa*H),
                    "average_amplitude_drift_original_time": float(H/R),
                    "average_H_quad_reported_error": float(error/period),
                    "average_H_quad_vs_4001_point_simpson": float(abs(H-fine)),
                    "average_method": "adaptive quadrature of dense DOP853 output split at every x=+/-1 event; checked against 4001-point Simpson"}

        for sign in [-1, 1]:
            name = f"boundary{j}_{'minus' if sign < 0 else 'plus'}"
            record = {"name": name, "boundary_index": j, "R": R,
                      "initial_x": sign * xe, "status": "unresolved"}
            report["probes"].append(record)
            x0 = sign * xe
            history = [x0]
            converged = False
            stable_count = 0
            for iteration in range(70):
                if time.monotonic() > deadline:
                    record["status"] = "wall_budget_before_convergence"
                    break
                sol = solve_ivp(rhs, (0, period), [x0], method="DOP853",
                                rtol=1e-11, atol=1e-12, max_step=period/200)
                if not sol.success:
                    record["status"] = "integration_failed"
                    record["message"] = sol.message
                    break
                new = float(sol.y[0, -1])
                residual = abs(new - x0)
                history.append(new)
                stable_count = stable_count + 1 if residual < 2e-10 * max(1, abs(new)) else 0
                x0 = new
                if stable_count >= 2:
                    converged = True
                    break
            arrays[name + "_iterates"] = np.asarray(history)
            record["iterations"] = len(history) - 1
            record["iteration_converged"] = converged
            if not converged:
                record["last_iterate"] = x0
                continue

            sol = solve_ivp(rhs, (0, period), [x0], method="DOP853", rtol=2e-12,
                            atol=2e-13, max_step=period/500, dense_output=True,
                            events=[lower, upper])
            other = solve_ivp(rhs, (0, period), [x0], method="Radau", rtol=2e-11,
                              atol=2e-12, max_step=period/100)
            if not sol.success or not other.success:
                record["status"] = "verification_failed"
                continue
            crossings = sorted([0.0, period] + [float(t) for events in sol.t_events for t in events])
            inner_time = sum(right-left for left, right in zip(crossings[:-1], crossings[1:])
                             if abs(float(sol.sol((left+right)/2)[0])) < 1)
            fraction = inner_time / period
            log_multiplier = -period * (a + (m0-m1)*fraction)
            trace_t = arrays["phase_time"]
            trace_x = sol.sol(trace_t)[0]
            shifted = sol.sol((trace_t + period/2) % period)[0]
            arrays[name + "_trajectory"] = trace_x
            arrays[name + "_switch_times"] = np.asarray(crossings[1:-1])
            record.update({
                "status": "verified_periodic_response",
                "section_x": x0,
                "dop853_stroboscopic_residual": abs(float(sol.y[0, -1])-x0),
                "radau_stroboscopic_residual": abs(float(other.y[0, -1])-x0),
                "independent_endpoint_difference": abs(float(other.y[0, -1]-sol.y[0, -1])),
                "x_min_sampled": float(trace_x.min()),
                "x_max_sampled": float(trace_x.max()),
                "inner_residence_time": float(inner_time),
                "inner_residence_fraction": float(fraction),
                "event_based_log_multiplier": float(log_multiplier),
                "event_based_multiplier": float(np.exp(log_multiplier)),
                "fiber_exponent_fast_time": float(log_multiplier/period),
                "breakpoint_crossings": len(crossings)-2,
                "half_wave_antisymmetry_error": float(np.max(np.abs(trace_x+shifted))),
            })
            record.update(averages(sol))
            print(name, json.dumps(record), flush=True)

        # The half-wave response satisfies flow(T/2, x0)=-x0.  This scalar
        # equation is strictly increasing, so its bracketed root is unique,
        # although the resulting full-period orbit need not be attracting.
        def half_residual(initial):
            half = solve_ivp(rhs, (0, period/2), [initial], method="DOP853",
                             rtol=2e-12, atol=2e-13, max_step=period/500)
            if not half.success:
                raise RuntimeError(half.message)
            return float(half.y[0, -1] + initial)

        if time.monotonic() < deadline - 5:
            center = brentq(half_residual, -20, 20, xtol=2e-13)
            controls = [(f"boundary{j}_half_wave_shoot", center, "root of x(T/2)+x(0)")]
            plus_name = f"boundary{j}_plus"
            if plus_name + "_trajectory" in arrays:
                partner = -float(arrays[plus_name + "_trajectory"][500])
                controls.append((f"boundary{j}_symmetry_partner", partner,
                                 "negative half-period shift of verified plus-start orbit"))
            for name, x0, construction in controls:
                sol = solve_ivp(rhs, (0, period), [x0], method="DOP853", rtol=2e-12,
                                atol=2e-13, max_step=period/500, dense_output=True,
                                events=[lower, upper])
                other = solve_ivp(rhs, (0, period), [x0], method="Radau", rtol=2e-11,
                                  atol=2e-12, max_step=period/100)
                crossings = sorted([0.0, period] + [float(t) for events in sol.t_events for t in events])
                inner_time = sum(right-left for left, right in zip(crossings[:-1], crossings[1:])
                                 if abs(float(sol.sol((left+right)/2)[0])) < 1)
                fraction = inner_time / period
                log_multiplier = -period * (a+(m0-m1)*fraction)
                trace_x = sol.sol(arrays["phase_time"])[0]
                arrays[name + "_trajectory"] = trace_x
                arrays[name + "_switch_times"] = np.asarray(crossings[1:-1])
                residual = abs(float(sol.y[0, -1])-x0)
                residual_other = abs(float(other.y[0, -1])-x0)
                control = {"name": name, "boundary_index": j, "R": R,
                           "construction": construction, "section_x": float(x0),
                           "status": "verified_periodic_response" if sol.success and other.success and max(residual,residual_other) < 1e-7 else "unresolved",
                           "dop853_stroboscopic_residual": residual,
                           "radau_stroboscopic_residual": residual_other,
                           "independent_endpoint_difference": abs(float(other.y[0,-1]-sol.y[0,-1])),
                           "inner_residence_fraction": float(fraction),
                           "event_based_log_multiplier": float(log_multiplier),
                           "event_based_multiplier": float(np.exp(log_multiplier)),
                           "fiber_exponent_fast_time": float(log_multiplier/period),
                           "x_min_sampled": float(trace_x.min()),
                           "x_max_sampled": float(trace_x.max()),
                           "breakpoint_crossings": len(crossings)-2,
                           "half_wave_antisymmetry_error": float(np.max(np.abs(trace_x+sol.sol((arrays["phase_time"]+period/2)%period)[0])))}
                report["shooting_controls"].append(control)
                control.update(averages(sol))
                print(name, json.dumps(control), flush=True)

    for j in range(len(amplitudes)):
        minus, plus = report["probes"][2*j:2*j+2]
        if minus["status"] == plus["status"] == "verified_periodic_response":
            minus["opposite_start_section_separation"] = plus["opposite_start_section_separation"] = abs(minus["section_x"]-plus["section_x"])
            plus_trace = arrays[plus["name"]+"_trajectory"]
            minus_trace = arrays[minus["name"]+"_trajectory"]
            # Uniform 1000-interval mesh has an exact half-period shift.
            symmetry = float(np.max(np.abs(minus_trace[:-1]+np.roll(plus_trace[:-1], 500))))
            minus["exchange_half_wave_error"] = plus["exchange_half_wave_error"] = symmetry
    report["elapsed_wall_seconds"] = time.monotonic()-started
    report["limitations"] = [
        "Five instantaneous amplitudes do not track the actual finite-alpha orbit's phase or basin.",
        "Two prescribed well starts plus symmetry/shooting controls are not a complete census of scalar periodic solutions.",
        "A contraction multiplier is only for the frozen scalar direction, not the finite-alpha full spectrum.",
        "No saddle-node location or global finite-alpha chaos exclusion is established.",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    np.savez_compressed(OUT.with_suffix(".npz"), **arrays)
    print("saved", OUT, "wall", report["elapsed_wall_seconds"], flush=True)


if __name__ == "__main__":
    main()
