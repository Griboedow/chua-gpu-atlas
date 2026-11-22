"""Integrate the autonomous Chua field and score the frozen dynamic-fold predictions."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parents[1]
FORECAST = ROOT / "results/chua-dynamic-fold/forecast.json"
DIGEST = "9ab3dd23bcfd0e8f5eb4e3fbd4aaedb41f7f7695ea0926abc1ec6cf2bf9ac0b1"


def integrate(args):
    row, fold, method, rtol = args
    alpha, x0 = row["alpha"], row["x0"]
    kappa, a, b = .1425, 2/7, -1
    omega = np.sqrt(kappa)

    def dx(state):
        x, y, _ = state
        return y-a*x-(b-a)*max(-1., min(1., x))

    def field(_, state):
        x, y, z = state
        return [dx(state), z+(x-y)/alpha, -kappa*y]

    def section(_, state):
        return state[1]

    section.direction = 1

    def extremum(_, state):
        return dx(state)

    speed = fold["radial_drift_per_turn_times_alpha"]/alpha
    turns = 17/(fold["normal_form_product"]*speed)**(1/3) + 8
    solution = solve_ivp(field, (0., turns*2*np.pi/omega),
                         [x0, 0., row["Z0"]], method=method,
                         rtol=rtol, atol=rtol/50, max_step=.35,
                         events=(section, extremum))
    if not solution.success:
        raise ArithmeticError(solution.message)
    times, peaks = solution.t_events
    section_states, extremes = solution.y_events
    x = section_states[:, 0]
    radius = np.hypot(section_states[:, 1], section_states[:, 2]/omega)
    x_fold, r_fold = fold["x_at_y0_rising"], fold["R"]
    crossings = np.flatnonzero((x[:-1] > x_fold) & (x[1:] <= x_fold))
    if not len(crossings):
        raise ArithmeticError("The test ended before the section marker")
    i = int(crossings[0])
    fraction = (x_fold-x[i])/(x[i+1]-x[i])
    marker = radius[i] + fraction*(radius[i+1]-radius[i])
    for turn, (start, end) in enumerate(itertools.pairwise(times)):
        segment = np.r_[x[turn:turn+2], extremes[(peaks>start) & (peaks<end), 0]]
        if segment.min() < -1 and segment.max() > 1:
            return {"alpha": alpha, "method": method, "rtol": rtol,
                    "section_turn_interpolated": float(i+fraction), "two_well_turn": turn,
                    "section_R": float(marker), "two_well_R": float(radius[turn]),
                    "section_ratio": float((marker-r_fold)/(row["section_prediction"]-r_fold)),
                    "two_well_ratio": float((radius[turn]-r_fold)/(row["two_well_prediction"]-r_fold)),
                    "turn_min_x": float(segment.min()), "turn_max_x": float(segment.max())}
    raise ArithmeticError("The test ended before a turn crossed both breakpoints")


def plot(rows, out):
    fig, ax = plt.subplots(figsize=(6.4, 3), layout="constrained")
    for marker, label, shape in [("section", "Fold coordinate on the section", "o"),
                                  ("two_well", "First turn spanning both wells", "s")]:
        ax.plot([r["alpha"] for r in rows], [r[marker+"_ratio"] for r in rows],
                marker=shape, label=label)
    ax.axhline(1, color=".4", linestyle=":", label="Airy prediction")
    ax.set(xscale="log", xlabel="Scale parameter α", ylabel="Measured / predicted excess amplitude")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"outputs/chua-dynamic-fold")
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    raw = FORECAST.read_bytes().replace(b"\r\n", b"\n")
    if hashlib.sha256(raw).hexdigest() != DIGEST:
        raise ValueError("The archived forecast changed")
    predictions = json.loads(raw)
    rows = predictions["runs"]
    fold = predictions["fold"]
    tasks = [(row, fold, "DOP853", 2e-11) for row in rows]
    tasks += [(row, fold, "LSODA", 2e-12) for row in (rows[0], rows[-2], rows[-1])]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        all_runs = list(pool.map(integrate, tasks))
    main_runs, controls = all_runs[:len(rows)], all_runs[len(rows):]
    results = {"forecast_sha256": DIGEST, "runs": main_runs, "independent_controls": controls,
               "checks": {
                   "section_within_20_percent_at_two_largest": all(
                       abs(r["section_ratio"]-1) < .2 for r in main_runs[-2:]),
                   "two_well_within_25_percent_at_largest": abs(main_runs[-1]["two_well_ratio"]-1) < .25}}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"measurements.json").write_text(json.dumps(results, indent=2)+"\n",
                                                     encoding="utf-8", newline="\n")
    plot(main_runs, args.output_dir/"fold-delay.png")
    for r in all_runs:
        print(r["alpha"], r["method"], r["section_ratio"], r["two_well_ratio"])


if __name__ == "__main__":
    main()
