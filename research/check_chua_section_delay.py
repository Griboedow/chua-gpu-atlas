"""Score the corrected section-marker forecast at the five held-out scale parameters."""

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from check_chua_fold_delay import integrate
from chua_fiber_switching import FrozenFiber
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/chua-dynamic-fold"
DIGEST = "04bab0d992eb8b56344fecf7fd5d2f93a668c8d561aa14069726d26b23a333b1"


def check_frozen_derivatives(coeff):
    # Fix forcing phase y=R*sin(w*t) and solve both crossing times for
    # arbitrary initial x and R. This differentiates the map off its fixed branch.
    mp.mp.dps = 45
    f = FrozenFiber(mp.mpf(".1425"), mp.mpf(-2), -mp.mpf(5) / 7)
    T, w = f.period, f.omega
    d = f.a * T / (f.a - f.b)
    C, S = f.forcing(d, T, 1)
    offset = ((mp.atan2(S, C) - mp.pi / 2) % (2 * mp.pi)) / w
    r, x = mp.mpf(coeff["R_fold"]), mp.mpf(coeff["x_fold"])

    def P(h, j):
        amplitude = r + j

        def flight(k, start, end, initial, eq):
            return (
                eq
                + (initial - eq) * mp.exp(-k * (end - start))
                + amplitude * f.forcing_weights(k, start, end)[1]
            )

        t1 = mp.findroot(lambda t: flight(f.b, 0, t, x + h, 0) - 1, d - offset)
        t2 = mp.findroot(lambda t: flight(f.a, t1, t, 1, f.xe) - 1, T - offset)
        return flight(f.b, t2, T, 1, 0)

    measured = {
        name: float(mp.diff(P, (0, 0), order) / factor)
        for name, order, factor in [
            ("A", (0, 1), 1),
            ("B", (2, 0), 2),
            ("C", (1, 1), 1),
            ("D", (3, 0), 6),
        ]
    }
    for key, value in measured.items():
        if abs(value / coeff[key] - 1) > 1e-10:
            raise ArithmeticError(f"Frozen map derivative mismatch: {key} = {value}")
    return measured


def check_feedback_derivatives(k):
    """Independent central differences of full returns and the integrated radial drift."""
    r, x = k["R_fold"], k["x_fold"]
    w = np.sqrt(0.1425)
    T = 2 * np.pi / w

    def end(h, eps):
        def field(_, state):
            X, y, z, _ = state
            g = 2 / 7 * X - 9 / 7 * np.clip(X, -1, 1)
            return [y - g, z + eps * (X - y), -w * w * y, y * (X - y) / r]

        def rising(_, state):
            return state[1]

        rising.direction = 1
        s = solve_ivp(
            field,
            (0, T + 0.02),
            [x + h, 0, w * r, 0],
            method="DOP853",
            rtol=2e-13,
            atol=3e-15,
            max_step=0.07,
            events=rising,
        )
        if not s.success or len(s.t_events[0]) != 2:
            raise ArithmeticError("Feedback derivative control missed the section")
        return s.y_events[0][-1]

    eps, dx = 2e-6, 0.0005
    E = (end(0, eps)[0] - end(0, -eps)[0]) / (2 * eps)
    H = (end(dx, 0)[3] - end(-dx, 0)[3]) / (2 * dx)
    if abs(E / k["E"] - 1) > 2e-5 or abs(H / k["radial_drift_x"] - 1) > 2e-5:
        raise ArithmeticError(f"Feedback derivative disagreement: {E}, {H}")
    return {"E": float(E), "radial_drift_x": float(H)}


def plot(first, held, coefficients, target):
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.5, 3.2), layout="constrained")
    K, L = coefficients["leading_coefficient"], coefficients["correction_coefficient"]
    R = coefficients["R_fold"]
    grid = np.geomspace(2048, 4194304, 160)
    left.plot(
        grid, K * grid ** (-2 / 3), color="#af5733", ls="--", label="Leading order"
    )
    left.plot(
        grid, K * grid ** (-2 / 3) + L / grid, color="#24586c", label="With correction"
    )
    right.axhline(L, color="#24586c", label="Predicted coefficient")
    for rows, marker, name in [(first, "o", "First test"), (held, "s", "Held out")]:
        alpha = np.array([r["alpha"] for r in rows])
        excess = np.array([r["section_R"] - R for r in rows])
        for ax, y in [
            (left, excess),
            (right, alpha * (excess - K * alpha ** (-2 / 3))),
        ]:
            ax.scatter(
                alpha,
                y,
                marker=marker,
                facecolor="white",
                edgecolor="#263340",
                s=22,
                zorder=3,
                label=name,
            )
    left.axhline(0, lw=0.6, color=".65")
    left.set(yscale="symlog", ylabel="Section amplitude minus frozen fold")
    left.set_yscale("symlog", linthresh=1e-5)
    right.set(ylabel="α × error of leading prediction")
    for ax in (left, right):
        ax.set(xscale="log", xlabel="Scale parameter α")
        ax.legend(frameon=False, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(target, dpi=175)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "outputs/chua-dynamic-fold"
    )
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    raw = (SOURCE / "section-forecast.json").read_bytes().replace(b"\r\n", b"\n")
    if hashlib.sha256(raw).hexdigest() != DIGEST:
        raise ValueError("The held-out forecast changed")
    design = json.loads(raw)
    k = design["coefficients"]
    primary = json.loads((SOURCE / "forecast.json").read_text())["fold"]
    old = json.loads((SOURCE / "measurements.json").read_text())["runs"]
    rows = design["runs"]
    if {r["alpha"] for r in rows} & {r["alpha"] for r in old}:
        raise ValueError("Held-out settings overlap the original test")
    checks = check_frozen_derivatives(k)
    tasks = []
    for r in rows:
        # The common integrator reports two markers. Only its interpolated
        # section_R is scored for this second forecast.
        row = dict(
            r,
            section_prediction=k["R_fold"] + r["leading_excess"],
            two_well_prediction=k["R_fold"]
            + 2.33810741046 * r["leading_excess"] / k["Airy_marker"],
        )
        tasks.append((row, primary, "DOP853", 2e-11))
    repeats = [(row, fold, "LSODA", 2e-12) for row, fold, _, _ in tasks[-3:]]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        runs = list(pool.map(integrate, tasks + repeats))
    by_alpha = {r["alpha"]: r for r in rows}
    scores = []
    for run in runs:
        r = by_alpha[run["alpha"]]
        excess = run["section_R"] - k["R_fold"]
        base_error = excess - r["leading_excess"]
        error = excess - r["corrected_excess"]
        scores.append(
            {
                "alpha": r["alpha"],
                "method": run["method"],
                "rtol": run["rtol"],
                "section_R": run["section_R"],
                "measured_excess": excess,
                "leading_error": base_error,
                "corrected_error": error,
                "corrected_error_over_leading_excess": error / r["leading_excess"],
                "within_gate": abs(error) < 0.1 * r["leading_excess"]
                and abs(error) < abs(base_error),
            }
        )
    result = {
        "forecast_sha256": DIGEST,
        "independent_frozen_map_derivatives": checks,
        "full_field_feedback_derivatives": check_feedback_derivatives(k),
        "runs": scores[: len(rows)],
        "independent_controls": scores[len(rows) :],
        "all_three_held_out_checks_pass": all(
            r["within_gate"] for r in scores[len(rows) - 3 :]
        ),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "section-measurements.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    plot(old, result["runs"], k, args.output_dir / "section-delay.png")
    for r in scores:
        print(r["alpha"], r["method"], r["corrected_error_over_leading_excess"])


if __name__ == "__main__":
    main()
