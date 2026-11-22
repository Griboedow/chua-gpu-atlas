"""Finite-time voltage-preparation screen for the 2023 two-op-amp Chua circuit."""

import argparse
import json
from pathlib import Path

import numpy as np
from numba import njit, prange
from scipy.integrate import solve_ivp
from scipy.stats import qmc

OUT = Path(__file__).resolve().parents[1] / "outputs" / "chua-precharge"
COMPONENTS = ("R", "R1", "R2", "R3", "R4", "R5", "R6", "C1", "C2", "L", "RL")
NOMINAL = np.array(
    [
        788.0,
        8520.0,
        123.0,
        12020.0,
        1008.0,
        1003.0,
        24410.0,
        37e-9,
        311e-9,
        0.01579,
        2.8,
    ]
)
PARAMETERS = ("alpha", "beta", "gamma", "m0", "m1", "B_volts", "outer_break", "m2")
TABLE = np.array(
    [
        [8.3886, 8.4222],
        [12.1812, 12.2790],
        [0.04328, 0.04363],
        [-1.3689, -1.3634],
        [-0.1715, -0.1700],
    ]
)
PHYSICAL = np.array([8.4054, 12.2301, 0.0435, -1.3661, -0.1708])
MODEL = np.array([8.41, 12.23, 0.0435, -1.366, -0.17])
SEED = 20261122


def circuit(values):
    r, r1, r2, r3, r4, r5, r6, c1, c2, l, rl = values
    k1, k2 = 1 + r3 / r4, 1 + r5 / r6
    return np.array(
        [
            c2 / c1,
            c2 * r * r / l,
            c2 * r * rl / l,
            r * (-r3 / (r4 * r1) - r5 / (r6 * r2)),
            r * (1 / r1 - r5 / (r6 * r2)),
            15 / k1,
            k1 / k2,
            r * (1 / r1 + 1 / r2),
        ]
    )


@njit(inline="always")
def field(x, y, z, p):
    a, b, g, m0, m1, _, outer, m2 = p
    h = m1 * x + (m0 - m1) * min(1.0, max(-1.0, x))
    h += (m2 - m1) * (max(0.0, x - outer) + min(0.0, x + outer))
    return a * (y - x - h), x - y + z, -b * y - g * z


@njit(inline="always")
def cubic(f, v0, v1, d0, d1, dt):
    return (1 - f) ** 2 * ((1 + 2 * f) * v0 + f * dt * d0) + f * f * (
        (3 - 2 * f) * v1 + (f - 1) * dt * d1
    )


@njit(cache=True)
def integrate(u1, u2, p, dt, end):
    x, y, z = u1 / p[5], u2 / p[5], 0.0
    n = round(end / dt)
    count = 0
    mean = 0.0
    lower, upper = np.inf, -np.inf
    xmax = abs(x)
    # 24 upward zero crossings in physical x suffice to check an external cycle.
    crossings = np.zeros((24, 3))
    ncross = 0
    for i in range(n):
        k = field(x, y, z, p)
        l = field(x + 0.5 * dt * k[0], y + 0.5 * dt * k[1], z + 0.5 * dt * k[2], p)
        m = field(x + 0.5 * dt * l[0], y + 0.5 * dt * l[1], z + 0.5 * dt * l[2], p)
        j = field(x + dt * m[0], y + dt * m[1], z + dt * m[2], p)
        xn = x + dt * (k[0] + 2 * l[0] + 2 * m[0] + j[0]) / 6
        yn = y + dt * (k[1] + 2 * l[1] + 2 * m[1] + j[1]) / 6
        zn = z + dt * (k[2] + 2 * l[2] + 2 * m[2] + j[2]) / 6
        xmax = max(xmax, abs(xn))
        if i >= 3 * n // 4:
            count += 1
            mean += xn
            lower, upper = min(lower, xn), max(upper, xn)
            if x < 0 <= xn:
                final = field(xn, yn, zn, p)
                left, right = 0.0, 1.0
                for _ in range(18):
                    f = 0.5 * (left + right)
                    if cubic(f, x, xn, k[0], final[0], dt) < 0:
                        left = f
                    else:
                        right = f
                f = 0.5 * (left + right)
                crossings[ncross % 24] = [
                    (i + f) * dt,
                    cubic(f, y, yn, k[1], final[1], dt),
                    cubic(f, z, zn, k[2], final[2], dt),
                ]
                ncross += 1
        x, y, z = xn, yn, zn
    if ncross >= 24:
        crossings = np.roll(crossings.flatten(), -3 * (ncross % 24)).reshape(24, 3)
        period = (crossings[-1, 0] - crossings[0, 0]) / 23
        spread = max(
            crossings[:, 1].max() - crossings[:, 1].min(),
            crossings[:, 2].max() - crossings[:, 2].min(),
        )
    else:
        period, spread = 0.0, -1.0
    return np.array([mean / count, lower, upper, xmax, period, spread, x, y, z])


@njit(parallel=True)
def batch(grid, settings, dt, end):
    values = np.empty((len(settings), len(grid), 9))
    for i in prange(len(settings) * len(grid)):
        j, k = i // len(grid), i % len(grid)
        values[j, k] = integrate(grid[k, 0], grid[k, 1], settings[j], dt, end)
    return values


def classify(row, p):
    mean, low, high, _xmax, period, spread, x, y, z = row
    if max(abs(x), abs(y), abs(z)) < 1.0e-10 and high - low < 1.0e-9:
        return "origin"
    # Distance to the actual signed equilibrium guards against slowly decaying oscillations.
    eq = -(p[3] - p[4]) / (1 + p[4] - p[2] / (p[1] + p[2]))
    s = np.sign(mean)
    xyz = np.array([eq, eq * p[2] / (p[1] + p[2]), -eq * p[1] / (p[1] + p[2])]) * s
    if high - low < 0.003 and max(abs(np.array([x, y, z]) - xyz)) < 0.002:
        return "focus+" if mean > 0 else "focus-"
    if low < -3 and high > 3 and period > 0 and 0 <= spread < 0.01:
        return "outer-cycle"
    if high - low > 0.2 and abs(mean) > 0.3 and max(abs(low), abs(high)) < 3:
        return "inner+" if mean > 0 else "inner-"
    return "unresolved"


def rows(grid, settings, dt=0.02, end=600.0):
    vals = batch(np.ascontiguousarray(grid), np.ascontiguousarray(settings), dt, end)
    return np.array(
        [[classify(v, p) for v in group] for p, group in zip(settings, vals)]
    ), vals


@njit(parallel=True)
def refine_batch(grid, settings, dt, end):
    values = np.empty((len(grid), 9))
    for i in prange(len(grid)):
        values[i] = integrate(grid[i, 0], grid[i, 1], settings[i], dt, end)
    return values


def refine_outer(grid, settings, labels, values, end):
    # The second pair of breakpoints reduces the order of a fixed-step RK4
    # crossing. Keep the recurrence criterion fixed and recompute its failures.
    i, k = np.nonzero(
        (labels == "unresolved") & (values[:, :, 1] < -3) & (values[:, :, 2] > 3)
    )
    if len(i):
        refined = refine_batch(grid[k], settings[i], 0.005, end)
        for j, state in enumerate(refined):
            labels[i[j], k[j]] = classify(state, settings[i[j]])
            values[i[j], k[j]] = state
    return {
        "reintegrated": len(i),
        "still_unresolved": int(np.sum(labels[i, k] == "unresolved")),
        "step": 0.005,
    }


def nominal_settings():
    p = circuit(NOMINAL)
    settings = [p.copy(), p.copy(), p.copy()]
    settings[1][:5] = PHYSICAL
    settings[2][:5] = MODEL
    return np.array(settings)


def designs():
    # Fractions from section 4. These are bounded design coordinates, not a
    # probability model of the error of any individual measured component.
    fractions = np.full(11, 0.01)
    fractions[7:10] = 0.001
    candidates = qmc.Sobol(11, scramble=True, seed=SEED).random_base2(16)
    accepted = []
    for index, u in enumerate(candidates):
        components = NOMINAL * (1 + fractions * (2 * u - 1))
        p = circuit(components)
        if np.all((p[:5] >= TABLE[:, 0]) & (p[:5] <= TABLE[:, 1])):
            accepted.append((index, components, p))
            if len(accepted) == 64:
                break
    if len(accepted) != 64:
        raise RuntimeError("Insufficient compatible component draws")
    # Only the first five coordinates are reported in Table 2. An independent
    # envelope in these five coordinates is a second, deliberately different test.
    box = []
    for u in qmc.Sobol(5, scramble=True, seed=SEED + 1).random_base2(5):
        p = circuit(NOMINAL)
        p[:5] = TABLE[:, 0] + u * np.diff(TABLE, axis=1)[:, 0]
        box.append(p)
    return accepted, np.array(box), fractions


def first_entry(p, step=0.015, end=600.0):
    lo, hi = 2.93, 3.46

    def outcome(u):
        return classify(integrate(u, 0.0, p, step, end), p)

    if outcome(lo) != "focus+" or outcome(hi) != "inner+":
        return None
    for _ in range(17):
        mid = (lo + hi) / 2
        target = outcome(mid)
        if target == "focus+":
            lo = mid
        elif target == "inner+":
            hi = mid
        else:
            return None
    return [lo, hi]


def writes(data, name):
    OUT.mkdir(exist_ok=True, parents=True)
    (OUT / name).write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def survey():
    import csv

    accepted, box, _fractions = designs()
    nominal = nominal_settings()
    settings = np.vstack([nominal, [r[2] for r in accepted], box])
    u1 = np.round(np.arange(0.1, 10.01, 0.1), 2)
    grid = np.array([(u, v) for v in [-0.05, 0.0, 0.05] for u in u1])
    OUT.mkdir(exist_ok=True, parents=True)
    with (OUT / "settings.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["index", "design", "sobol_index", *COMPONENTS, *PARAMETERS])
        for i, p in enumerate(settings):
            if i < 3:
                w.writerow(
                    [
                        i,
                        ["components", "physical-rounded", "mathematical"][i],
                        "",
                        *(NOMINAL if i == 0 else [""] * 11),
                        *p,
                    ]
                )
            elif i < 67:
                idx, comps, _ = accepted[i - 3]
                w.writerow([i, "coupled", idx, *comps, *p])
            else:
                w.writerow([i, "table-box", i - 67, *([""] * 11), *p])
    print(
        "Compatible component candidates",
        len(accepted),
        "in",
        accepted[-1][0] + 1,
        flush=True,
    )
    labels, values = rows(grid, settings)
    retry = refine_outer(grid, settings, labels, values, 600.0)
    # Compact text records contain every launch, including all intermediate outcomes.
    with (OUT / "outcomes.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["U1_V", "U2_V", *[f"case_{i}" for i in range(len(settings))]])
        for k, (u, v) in enumerate(grid):
            w.writerow([u, v, *labels[:, k]])
    maps = {}
    for i in [0, 1, 2]:
        segment = labels[i, len(u1) : 2 * len(u1)]
        starts = np.r_[0, np.where(segment[1:] != segment[:-1])[0] + 1]
        maps[str(i)] = [
            [float(u1[a]), float(u1[b - 1]), str(segment[a])]
            for a, b in zip(starts, np.r_[starts[1:], len(u1)])
        ]
    # For each tick all 3 transverse preparation errors and all 99 settings agree.
    all_agree = []
    for k, u in enumerate(u1):
        unique = np.unique(labels[:, [k, k + len(u1), k + 2 * len(u1)]])
        all_agree.append(unique[0] if len(unique) == 1 else "mixed")
    all_agree = np.array(all_agree)
    starts = np.r_[0, np.where(all_agree[1:] != all_agree[:-1])[0] + 1]
    summary = {
        "nominal_bands": maps,
        "component_candidates_tested": accepted[-1][0] + 1,
        "settings": len(settings),
        "launches_per_setting": len(grid),
        "dt": 0.02,
        "horizon": 600,
        "recurrence_refinement": retry,
        "common_grid_runs": [
            [float(u1[a]), float(u1[b - 1]), str(all_agree[a])]
            for a, b in zip(starts, np.r_[starts[1:], len(u1)])
        ],
        "largest_tail_abs_x": float(np.abs(values[:, :, [1, 2]]).max()),
        "largest_abs_x_over_second_break": float(
            (values[:, :, 3] / settings[:, 6, None]).max()
        ),
    }
    writes(summary, "survey.json")
    print("Voltage runs", summary["common_grid_runs"], flush=True)


def sensitivity():
    accepted, box, _fractions = designs()
    groups = {
        "nominal": list(nominal_settings()),
        "coupled": [v[2] for v in accepted],
        "table-box": box,
    }
    edges = {}
    for name, params in groups.items():
        vals = [first_entry(p) for p in params]
        edges[name] = vals
        mid = np.array([np.mean(row) for row in vals if row is not None])
        print(name, len(mid), np.min(mid), np.median(mid), np.max(mid), flush=True)
    elasticity = {}
    for j, key in enumerate(COMPONENTS):
        xs = []
        for rel in [-0.0002, 0.0002]:
            val = NOMINAL.copy()
            val[j] *= 1 + rel
            xs.append(first_entry(circuit(val)))
        elasticity[key] = {
            "brackets": xs,
            "mV_per_plus_0.1_percent": float((np.mean(xs[1]) - np.mean(xs[0])) * 2500),
        }
    p = circuit(NOMINAL)
    controls = {str(step): first_entry(p, step, 900.0) for step in [0.02, 0.01, 0.005]}
    print("Component derivatives", elasticity, flush=True)
    writes(
        {
            "transition_brackets_V": edges,
            "component_derivatives": elasticity,
            "nominal_step_controls": controls,
            "threshold_start_bracket": [2.93, 3.46],
            "bisections": 17,
            "dt": 0.015,
            "horizon": 600,
        },
        "sensitivity.json",
    )


def independent(u, v, p, full=True):
    a, b, g, m0, m1, B, outer, m2 = p

    def rhs(t, state):
        x, y, z = state
        ax = abs(x)
        h = m0 * min(ax, 1.0) + m1 * max(0.0, ax - 1.0)
        if full and ax > outer:
            h += (m2 - m1) * (ax - outer)
        return a * (y - x - np.copysign(1.0, x) * h), x - y + z, -b * y - g * z

    def crossing(t, state):
        return state[0]

    crossing.direction = 1
    times = np.linspace(450.0, 600.0, 6001)
    sol = solve_ivp(
        rhs,
        (0.0, 600.0),
        (u / B, v / B, 0.0),
        method="DOP853",
        rtol=2e-9,
        atol=2e-11,
        max_step=0.15,
        t_eval=times,
        events=crossing,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    xs = sol.y[0]
    ts = sol.t_events[0]
    events = sol.y_events[0]
    mask = ts >= 450
    ts = ts[mask]
    events = events[mask] if len(events) else np.empty((0, 3))
    periods = np.diff(ts[-24:])
    spread = (
        float(np.ptp(events[-24:, 1:], axis=0).max()) if len(periods) == 23 else -1.0
    )
    period = float(np.mean(periods)) if len(periods) == 23 else 0.0
    metrics = [
        float(xs.mean()),
        float(xs.min()),
        float(xs.max()),
        float(max(np.abs(xs))),
        period,
        spread,
        *sol.y[:, -1].tolist(),
    ]
    return {
        "U1_V": u,
        "U2_V": v,
        "diode_segments": 5 if full else 3,
        "outcome": classify(metrics, p),
        "x_range": [metrics[1], metrics[2]],
        "period": period,
        "section_spread": spread,
        "period_spread": float(np.ptp(periods)) if len(periods) == 23 else None,
    }


def controls():
    nominal = circuit(NOMINAL)
    accepted, _, _ = designs()
    checks = []
    for u, v in [
        (2.0, 0.05),
        (3.1902, 0.0),
        (3.1912, 0.0),
        (3.5, -0.05),
        (3.6, 0.0),
        (3.7, 0.05),
        (4.8, 0.05),
        (5.8, 0.0),
        (8.5, 0.0),
        (-3.6, 0.0),
        (9.0, 0.0),
    ]:
        check = independent(u, v, nominal)
        check["rk4"] = classify(integrate(u, v, nominal, 0.01, 900.0), nominal)
        checks.append(check)
    three = independent(9.0, 0.0, nominal, False)
    # Search the interior of the two-voltage rectangles and include every edge
    # grid line. Dimensions: 9 longitudinal × 5 transverse launches per rectangle.
    groups = [
        ("focus+", 0.3, 2.8),
        ("inner+", 3.5, 3.7),
        ("focus-", 4.4, 5.3),
        ("outer-cycle", 9.0, 10.0),
    ]
    grid = np.array(
        [
            (u, v)
            for _, lo, hi in groups
            for u in np.linspace(lo, hi, 9)
            for v in np.linspace(-0.05, 0.05, 5)
        ]
    )
    settings = np.vstack([nominal_settings(), [v[2] for v in accepted], designs()[1]])
    labels, values = rows(grid, settings, 0.0125, 800.0)
    retry = refine_outer(grid, settings, labels, values, 800.0)
    details = {}
    for k, (expected, lo, hi) in enumerate(groups):
        local = labels[:, 45 * k : 45 * (k + 1)]
        keys, counts = np.unique(local, return_counts=True)
        details[expected] = {
            "U1_range_V": [lo, hi],
            "U2_range_V": [-0.05, 0.05],
            "trials": local.size,
            "counts": dict(zip(keys.tolist(), counts.tolist())),
        }
    writes(
        {
            "independent": checks,
            "three_segment_control": three,
            "rectangle_checks": details,
            "recurrence_refinement": retry,
            "rectangle_grid": [9, 5],
            "dt": 0.0125,
            "horizon": 800,
        },
        "controls.json",
    )
    print("Controls", details, flush=True)


def comparison():
    import csv

    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    accepted, box, _ = designs()
    pt = circuit(NOMINAL)
    saved = json.loads((OUT / "sensitivity.json").read_text())["transition_brackets_V"]
    ideal = np.array([np.mean(v) for v in saved["table-box"]])
    physical = np.array([np.mean(v) for v in saved["coupled"]])
    center = float(np.mean(saved["nominal"][0]) / pt[5])
    slopes = np.linalg.lstsq(box[:, :5] - pt[:5], ideal / pt[5] - center, rcond=None)[0]
    coupled = np.array([item[2] for item in accepted])
    forecast = coupled[:, 5] * (center + (coupled[:, :5] - pt[:5]) @ slopes)
    residual = 1000 * (forecast - physical)
    constant = 1000 * (center * pt[5] - physical)
    compare = {
        "reference_parameters": pt[:5].tolist(),
        "normalized_reference_threshold": center,
        "derivative_order": list(PARAMETERS[:5]),
        "dx_threshold_dp": slopes.tolist(),
        "training_cases": len(box),
        "separate_component_cases": len(physical),
        "affine_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
        "constant_RMSE_mV": float(np.sqrt(np.mean(constant**2))),
        "affine_max_error_mV": float(max(abs(residual))),
        "worst_case_component_grid_index": int(np.argmax(abs(residual))) + 3,
        "forecasts_V": forecast.tolist(),
    }
    finer = np.array([np.mean(first_entry(p, 0.005, 900.0)) for p in coupled])
    residual_fine = 1000 * (forecast - finer)
    compare["separate_numerical_controls"] = {
        "dt": 0.005,
        "horizon": 900,
        "reference_mV_shift_max": float(max(abs(finer - physical)) * 1000),
        "affine_RMSE_mV": float(np.sqrt(np.mean(residual_fine**2))),
        "affine_max_error_mV": float(max(abs(residual_fine))),
        "constant_RMSE_mV": float(
            np.sqrt(np.mean((1000 * (pt[5] * center - finer)) ** 2))
        ),
        "fine_thresholds_V": finer.tolist(),
        "dop853_bracket_checks": [],
    }
    selected = np.unique(
        [np.argmin(physical), np.argmax(physical), np.argmax(abs(residual))]
    )
    for j in selected:
        points = [
            independent(float(forecast[j] + d), 0.0, coupled[j])
            for d in [-0.001, 0.001]
        ]
        compare["separate_numerical_controls"]["dop853_bracket_checks"].append(
            {
                "settings_index": int(3 + j),
                "predicted_V": float(forecast[j]),
                "checks": points,
            }
        )
    p = nominal_settings()[2]
    compare["mathematical_column_near_outer_boundary"] = [
        independent(p[5] * (7.147684638806898 + 0.001), 0.0, p, full)
        for full in [True, False]
    ]
    writes(compare, "transfer.json")
    print(
        "Separate-design forecast (mV)",
        compare["affine_RMSE_mV"],
        compare["affine_max_error_mV"],
        "vs constant",
        compare["constant_RMSE_mV"],
        flush=True,
    )
    with (OUT / "outcomes.csv").open(encoding="utf-8") as f:
        rows0 = [r for r in csv.reader(f) if r[1] == "0.0"]
    names = ["focus+", "inner+", "focus-", "inner-", "outer-cycle", "unresolved"]
    colors = ["#a9cce5", "#1d70a2", "#eec38b", "#a33725", "#6c527e", "#eeeeee"]
    pixels = np.array([[names.index(v) for v in row[5:69]] for row in rows0]).T
    fig, (ax, bx) = plt.subplots(
        2,
        1,
        figsize=(8, 5.8),
        layout="constrained",
        gridspec_kw={"height_ratios": [1.3, 1.0]},
    )
    ax.imshow(
        pixels,
        origin="lower",
        extent=(0.05, 10.05, 0.5, 64.5),
        interpolation="nearest",
        cmap=ListedColormap(colors),
        vmin=-0.5,
        vmax=5.5,
        aspect="auto",
    )
    ax.set(
        xlabel="Initial capacitor voltage U1 (V), U2 = 0",
        ylabel="Component realization",
        title="Outcomes for 64 component realizations compatible with Table 2",
    )
    handles = [
        Patch(color=c, label=n)
        for c, n in zip(
            colors,
            ["Focus +", "Inner +", "Focus −", "Inner −", "Outer cycle", "Unresolved"],
        )
        if n != "Unresolved" or (pixels == 5).any()
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, 1.03),
        ncol=len(handles),
        fontsize=8,
        framealpha=0.96,
        columnspacing=0.7,
    )
    limits = (min(physical) - 0.005, max(physical) + 0.005)
    bx.plot(limits, limits, color="#5e5e5e", lw=1)
    bx.scatter(physical, forecast, s=16, color="#226c75")
    bx.set(
        xlim=limits,
        ylim=limits,
        xlabel="Integrated first transition (V)",
        ylabel="Predicted transition (V)",
        title="Threshold prediction on separate component settings",
    )
    bx.grid(alpha=0.15)
    fig.savefig(OUT / "voltage-preparation.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=OUT,
        help="Directory for regenerated records and comparison inputs (default: %(default)s).",
    )
    parser.add_argument(
        "--part",
        choices=["pilot", "survey", "sensitivity", "controls", "comparison"],
        default="survey",
    )
    args = parser.parse_args()
    OUT = args.output_dir.resolve()
    if args.part == "survey":
        survey()
        raise SystemExit
    if args.part == "sensitivity":
        sensitivity()
        raise SystemExit
    if args.part == "controls":
        controls()
        raise SystemExit
    if args.part == "comparison":
        comparison()
        raise SystemExit
    grid = np.column_stack([np.arange(0.025, 10.001, 0.025), np.zeros(400)])
    p = nominal_settings()
    labels, values = rows(grid, p)
    for j, name in enumerate(["components", "table physical", "table model"]):
        print(name, p[j], flush=True)
        for cut in np.r_[0, np.where(labels[j, 1:] != labels[j, :-1])[0] + 1]:
            print(grid[cut, 0], labels[j, cut], values[j, cut, :6])
