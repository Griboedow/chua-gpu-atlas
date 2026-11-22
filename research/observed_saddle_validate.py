"""Measure the predeclared voltage-start forecasts and compare a full-field IVP."""
import csv
import json
import argparse
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov
from scipy.optimize import brentq, root

from observed_saddle_flights import NOMINAL, Section
from observed_saddle_residence import OUT


def direct(start, duration, rtol=2e-12, max_step=.025, dense=False):
    """Unsplit ODE. Switching times and affine-flight matrices are not inputs."""
    a, b, g = NOMINAL["alpha"], NOMINAL["beta"], NOMINAL["gamma"]
    m0, m1 = NOMINAL["m0"], NOMINAL["m1"]
    def field(t, xyz):
        x, y, z = xyz
        diode = m1*x + (m0-m1)*(abs(x+1)-abs(x-1))/2
        return [a*(y-x-diode), x-y+z, -b*y-g*z]
    def section(t, xyz):
        return xyz[0]+1
    section.direction = 1
    sol = solve_ivp(field, [0., duration], start, method="DOP853", rtol=rtol,
                    atol=rtol/25, max_step=max_step, dense_output=dense, events=section)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def exit_measurement(times, section_points, point, left, radius):
    sep = np.asarray(section_points)-point
    u = sep @ left
    n = int(np.flatnonzero(abs(u)>=radius)[0])
    if n == 0:
        raise ValueError("Launch missed the declared entry neighborhood")
    ratio = np.log(radius/abs(u[n-1])) / np.log(abs(u[n]/u[n-1]))
    return dict(first_exit=n, fractional=n-1+float(ratio),
                interpolated_time=float(times[n-1]+ratio*(times[n]-times[n-1])),
                entry_time=float(times[0]), entry_u=float(u[0]), exit_u=float(u[n]),
                section_time=list(times[:n+1]), section_yz=np.asarray(section_points[:n+1]).tolist(),
                section_u=u[:n+1].tolist())


def flight_measurement(section, x, point, left, radius):
    arcs = []
    t, initial, _ = section.launch(x, arcs)
    times, yz = [t], [initial[1:].copy()]
    for _ in range(60):
        dt, new_yz, _ = section.full_return(yz[-1], arcs)
        times.append(times[-1]+dt)
        yz.append(new_yz)
        if abs(left@(new_yz-point)) >= radius:
            result = exit_measurement(times, yz, point, left, radius)
            result["peak_abs_x_through_exit"] = section.absolute_x_peak(arcs)
            return result
    raise ArithmeticError("Exit from saddle neighborhood not detected")


def ivp_measurement(x, duration, point, left, radius):
    sol = direct([x, 0., 0.], duration)
    return exit_measurement(sol.t_events[0], sol.y_events[0][:, 1:], point, left, radius)


def independent_returns(point, left):
    def one_return(yz):
        sol = direct(np.r_[-1., yz], 2.6, rtol=3e-13, max_step=.015)
        k = np.flatnonzero(sol.t_events[0]>1e-6)[0]
        return sol.t_events[0][k], sol.y_events[0][k, 1:]
    t, yz = one_return(point)
    multipliers = {}
    for h in (0.0001, 0.00003):
        jac = np.column_stack([(one_return(point+e*h)[1]-one_return(point-e*h)[1])/(2*h)
                               for e in np.eye(2)])
        multipliers[str(h)] = sorted(float(m.real) for m in np.linalg.eigvals(jac))
    def cross(x):
        sol = direct([x, 0., 0.], 6.5, rtol=4e-13, max_step=.015)
        return float(left@(sol.y_events[0][1, 1:]-point))
    x_ivp = brentq(cross, 7.14, 7.16, xtol=1e-11)
    return dict(period=t, saddle_closure_norm=float(np.linalg.norm(yz-point)),
                centered_finite_difference_multipliers=multipliers,
                two_hit_threshold_x=x_ivp, two_hit_residual=cross(x_ivp))


def two_destinations(section, center):
    """Two post-exit examples in the three-segment ODE, outside the forecast test."""
    def close_outer(yz):
        _, end, dp = section.full_return(yz)
        return end-yz, dp-np.eye(2)
    outer = root(close_outer, [7.2, 55.], jac=True, tol=1e-12)
    arcs = []
    period, end, dp = section.full_return(outer.x, arcs)
    outer_gap = float(np.linalg.norm(end-outer.x))
    if outer_gap > 1e-8:
        raise ArithmeticError("Outer-cycle shooting did not converge")
    upper = direct([center+.001, 0., 0.], 750., rtol=2e-10, max_step=.1)
    lower = direct([center-.001, 0., 0.], 250., rtol=2e-10, max_step=.1)
    a, focus, _, _, _ = section.matrices[1]
    p = solve_continuous_lyapunov(a.T, -np.eye(3))
    # Largest level ellipsoid about the equilibrium contained in the outer affine cell.
    level = (focus[0]-1)**2 / np.linalg.inv(p)[0,0]
    delta = lower.y[:, -1]-focus
    return dict(offset_x=.001, lower_duration=250., lower_final=lower.y[:, -1].tolist(),
                positive_focus=focus.tolist(), focus_capture_level_ratio=(delta@p@delta)/level,
                focus_P_min_eigenvalue=np.linalg.eigvalsh(p)[0],
                focus_lyapunov_residual=float(np.linalg.norm(a.T@p+p@a+np.eye(3))),
                upper_duration=750., outer_section=outer.x.tolist(),
                upper_last_five_section_distances=np.linalg.norm(
                    upper.y_events[0][-5:,1:]-outer.x, axis=1).tolist(),
                outer_period=period, outer_multipliers=sorted(np.linalg.eigvals(dp).real),
                outer_closure=outer_gap, outer_peak_abs_x=section.absolute_x_peak(arcs),
                outer_diode_breakpoint_ratio=14.408/1.1606)


def measure(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    settings = json.loads((OUT / "forecast.json").read_text())
    point = np.array(settings["point"])
    left = np.array(settings["left_unstable_unit"])
    radius = settings["projected_exit_radius"]
    lam = settings["multipliers"][-1]
    gain = settings["entry_gain_per_x"]
    section = Section()
    second_grid = Section(step=.0075)
    records, direct_checks = [], []
    for pred in settings["predictions"]:
        offset, x = pred["offset_x"], pred["launch_x"]
        row = flight_measurement(section, x, point, left, radius)
        row.update(pred)
        row["fractional_error"] = row["fractional"]-pred["predicted_fractional_return"]
        row["integer_error"] = row["first_exit"]-pred["predicted_first_exit"]
        row["control_fractional"] = float(np.log(radius/(abs(offset)*gain))/np.log(1.5))
        fine = flight_measurement(second_grid, x, point, left, radius)
        row["halved_step_fractional_difference"] = fine["fractional"]-row["fractional"]
        row["halved_step_integer_difference"] = fine["first_exit"]-row["first_exit"]
        records.append(row)
        if abs(offset) in (.008, .001, .0001, .00001):
            control = ivp_measurement(x, row["section_time"][-1]+1., point, left, radius)
            shared = min(len(control["section_yz"]), len(row["section_yz"]))
            direct_checks.append(dict(offset_x=offset, first_exit=control["first_exit"],
                                      fractional=control["fractional"],
                                      integer_difference=control["first_exit"]-row["first_exit"],
                                      fractional_difference=control["fractional"]-row["fractional"],
                                      max_section_distance=float(np.max(np.linalg.norm(
                                          np.array(control["section_yz"][:shared])-row["section_yz"][:shared],axis=1)))))
    result = dict(records=records, independent_full_field=direct_checks,
                  independent_saddle=independent_returns(point, left),
                  post_exit_examples=two_destinations(section, settings["crossing_x"]))
    slopes = []
    for side in (-1, 1):
        for large, small in ((.001, .0001), (.0001, .00001)):
            lo, hi = [next(r for r in records if r["offset_x"]==side*v) for v in (large, small)]
            change = hi["fractional"]-lo["fractional"]
            expected = np.log(10)/np.log(lam)
            slopes.append(dict(side=side, larger=large, smaller=small, observed_change=change,
                               expected_change=expected, difference=change-expected))
    result["decade_tests"] = slopes
    result["summary"] = dict(
        exact_count=sum(r["integer_error"]==0 for r in records),
        within_one_return=sum(abs(r["integer_error"])<=1 for r in records),
        max_abs_fractional_error=max(abs(r["fractional_error"]) for r in records),
        rms_fractional_error=float(np.sqrt(np.mean([r["fractional_error"]**2 for r in records]))),
        control_rms_fractional_error=float(np.sqrt(np.mean(
            [(r["control_fractional"]-r["fractional"])**2 for r in records]))),
        direct_agree=sum(r["integer_difference"]==0 for r in direct_checks),
        max_decade_error=max(abs(r["difference"]) for r in slopes),
        fine_grid_agree=sum(r["halved_step_integer_difference"]==0 for r in records))
    scores = result["summary"]
    result["predeclared_tests_passed"] = {
        "main": scores["within_one_return"] == 24 and scores["max_abs_fractional_error"] < .35,
        "independent": scores["direct_agree"] == 8,
        "slope_control": bool(scores["max_decade_error"] < .20),
    }
    (output_dir/"measurements.json").write_text(
        json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    fields = ["offset_x", "launch_x", "predicted_first_exit", "first_exit",
              "predicted_fractional_return", "fractional", "fractional_error"]
    with (output_dir/"scores.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps(result["summary"], indent=2))
    if not all(result["predeclared_tests_passed"].values()):
        raise SystemExit("A declared numerical test failed; inspect measurements.json.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=OUT.parents[1]/"outputs/chua-saddle-residence",
                        help="Directory for regenerated measurements and figure")
    parser.add_argument("--figure", action="store_true", help="Render the figure after scoring")
    args = parser.parse_args()
    measure(args.output_dir)
    if args.figure:
        from observed_saddle_plot import main as plot
        plot(args.output_dir)
