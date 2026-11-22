"""Measure local passage outside three symmetric folds of the classical Chua circuit."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from chua_block_return import BlockReturn
from chua_fast_flights import FastFlights
from estimate_chua_homoclinic_tangency import HalfMap
from scipy.integrate import solve_ivp
from scipy.optimize import root
from trace_chua_window import Chua

ROOT = Path(__file__).resolve().parents[1]
GATE = 0.2


def inputs():
    directory = ROOT / "results/chua-layer-census"
    edges = json.loads((directory / "family-edges.json").read_text())
    second = json.loads((directory / "second-boundaries.json").read_text())
    registry = json.loads(
        (ROOT / "results/chua-layer-atlas/family-registry.json").read_text()
    )
    families = {r["id"]: r for r in registry["families"]}
    return [
        (name, families[name], row)
        for name, row in (
            ("S01", next(r for r in edges if r["beta"] == 14.259875)),
            ("S08", next(r for r in edges if r["beta"] == 14.32)),
            ("S11", second),
        )
    ]


def centered_fold(hmap, family, row):
    fold = row["fold"]
    flights = FastFlights(Chua(), family["regions"], family["targets"], True)
    beta, q = flights.boundary(fold["beta"], np.array(fold["q"]), 1)
    point = q[:2]
    end, jac, aux = hmap.apply(point, beta)
    vals, vectors = np.linalg.eig(jac)
    idx = int(np.argmin(abs(vals - 1)))
    if np.max(np.abs(end - point)) > 1e-10 or abs(vals[idx] - 1) > 1e-8:
        raise ArithmeticError("Fold closure or neutral multiplier mismatch")
    right = np.real(vectors[:, idx])
    left = np.real(np.linalg.inv(vectors)[idx])

    def derivatives(step):
        # Differentiating the analytic event Jacobian once avoids a second
        # difference of nearly equal section positions in the quadratic term.
        beta_plus = hmap.apply(point, beta + step)[0]
        beta_minus = hmap.apply(point, beta - step)[0]
        a = left @ (beta_plus - beta_minus) / (2 * step)
        jac_plus = hmap.apply(point + step * right, beta)[1]
        jac_minus = hmap.apply(point - step * right, beta)[1]
        b = left @ (jac_plus - jac_minus) @ right / (4 * step)
        return float(a), float(b)

    raw = [derivatives(h) for h in (3e-7, 1e-7, 3e-8)]
    a, b = raw[1]
    if b < 0:
        a, b, left, right = -a, -b, -left, -right
    if a <= 0 or b <= 0 or any(abs(aa * bb / (a * b) - 1) > 1e-4 for aa, bb in raw):
        raise ArithmeticError("Unfolding direction or derivative step check failed")
    check = BlockReturn(
        {"regions": family["regions"], "targets": family["targets"], "q": q}, digits=48
    )
    h = mp.mpf("1e-12")
    p, v, w = (mp.matrix(x.tolist()) for x in (point, right, left))
    beta_mp = mp.mpf(float(beta))
    slope = (
        check.evaluate(beta_mp + h, p)["point"]
        - check.evaluate(beta_mp - h, p)["point"]
    ) / (2 * h)
    curvature = (
        (
            check.evaluate(beta_mp, p + h * v)["matrix"]
            - check.evaluate(beta_mp, p - h * v)["matrix"]
        )
        * v
        / (4 * h)
    )
    ab_check = float(mp.fdot(w, slope) * mp.fdot(w, curvature))
    if abs(ab_check / (a * b) - 1) > 1e-5:
        raise ArithmeticError("Independent high-precision coefficient mismatch")
    return {
        "beta": beta,
        "width": beta - row["symmetry"]["beta"],
        "point": point.tolist(),
        "right": right.tolist(),
        "left": left.tolist(),
        "a": a,
        "b": b,
        "ab_step_check": [aa * bb for aa, bb in raw],
        "ab_48_digit_check": ab_check,
        "half_period": aux["period"],
        "itinerary": aux["itinerary"],
        "fold_defect": float(max(abs(point - end))),
        "transverse_multiplier": float(vals[1 - idx].real),
    }


def passage(hmap, fold, delta, gate=GATE):
    a, b, beta = fold["a"], fold["b"], fold["beta"]
    origin, right, left = (np.array(fold[k]) for k in ("point", "right", "left"))
    point = origin - right * (gate / b)
    section = [point]
    times = [0.0]
    for _ in range(10000):
        end, _, audit = hmap.apply(point, beta + delta)
        if audit["itinerary"] != fold["itinerary"]:
            raise ArithmeticError("Diode itinerary changed inside the gate")
        section.append(end)
        times.append(times[-1] + audit["period"])
        if left @ (end - origin) >= gate / b:
            break
        point = end
    else:
        raise ArithmeticError("Passage did not exit the gate")

    def residual(p):
        end, jac, _ = hmap.apply(p, beta - delta)
        return end - p, jac - np.eye(2)

    inside = root(
        residual, origin - right * np.sqrt(a * delta / b), jac=True, tol=1e-11
    )
    end, jac, audit = hmap.apply(inside.x, beta - delta)
    lam = np.max(np.linalg.eigvals(jac).real)
    if (
        max(abs(end - inside.x)) > 1e-10
        or not 0 < lam < 1
        or audit["itinerary"] != fold["itinerary"]
    ):
        raise ArithmeticError("Inside attracting orbit or shooting closure mismatch")
    frequency = np.sqrt(a * b * delta)
    count = len(times) - 1
    row = {
        "delta_beta": delta,
        "delta_over_width": delta / fold["width"],
        "gate": gate,
        "returns": count,
        "elapsed": times[-1],
        "inside_multiplier": lam,
        "leading_returns": np.pi / frequency,
        "finite_gate_returns": 2 * np.arctan(gate / frequency) / frequency,
        "multiplier_product": count * (1 - lam),
        "last_coordinate": b * (left @ (section[-1] - origin)),
    }
    return row, np.array(section), np.array(times)


def direct_check(fold, delta, section, times):
    model = Chua()

    def crossing(_, state):
        return state[0]

    sol = solve_ivp(
        lambda _, s: model.field(fold["beta"] + delta, s),
        (0.0, times[-1] + 0.1 * fold["half_period"]),
        [0.0, *section[0]],
        method="DOP853",
        rtol=3e-13,
        atol=3e-15,
        max_step=0.06,
        events=crossing,
    )
    if not sol.success or len(sol.t_events[0]) != len(times):
        raise ArithmeticError("Direct full-field crossing count mismatch")
    points = sol.y_events[0][:, 1:] * (-1.0) ** np.arange(len(times))[:, None]
    gate = fold["b"] * ((points - fold["point"]) @ fold["left"])
    if np.max(gate[:-1]) >= GATE or gate[-1] < GATE:
        raise ArithmeticError("Direct exit count disagrees")
    return {
        "delta_beta": delta,
        "half_returns": len(times) - 1,
        "max_section_error": float(np.max(np.abs(points - section))),
        "max_time_error": float(np.max(abs(sol.t_events[0] - times))),
    }


def survey_entries(fold, settings):
    """Observe the original survey's starts at an off-grid fold detuning."""
    model = Chua()
    beta = fold["beta"] + 1e-4 * fold["width"]
    finish = settings["transient"] + settings["duration"]

    def crossing(_, state):
        return state[0]

    records = []
    for seed in settings["seeds"]:
        for rtol in (1e-9, 1e-11):
            sol = solve_ivp(
                lambda _, s: model.field(beta, s),
                (0.0, finish),
                seed,
                method="DOP853",
                rtol=rtol,
                atol=rtol / 100,
                max_step=0.1,
                events=crossing,
            )
            if not sol.success:
                raise ArithmeticError(sol.message)
            times, states = sol.t_events[0], sol.y_events[0]
            delta = states[:, 1:] * np.sign(states[:, 1:2]) - fold["point"]
            longitudinal = delta @ fold["left"]
            transverse = delta - np.outer(longitudinal, fold["right"])
            inside = (fold["b"] * np.abs(longitudinal) <= GATE) & (
                fold["b"] * np.linalg.norm(transverse, axis=1) <= GATE
            )
            # Only intervals bracketed by consecutive qualifying crossings count.
            intervals = inside[:-1] & inside[1:]
            lengths = np.maximum(
                0.0, times[1:] - np.maximum(times[:-1], settings["transient"])
            )
            marked_time = float(lengths[intervals].sum())
            change = np.diff(np.r_[False, intervals, False].astype(int))
            spans = [
                (float(times[i]), float(times[j]))
                for i, j in zip(
                    np.flatnonzero(change == 1), np.flatnonzero(change == -1)
                )
            ]
            longest = max(spans, key=lambda pair: pair[1] - pair[0], default=None)
            records.append(
                {
                    "seed": seed,
                    "rtol": rtol,
                    "crossings": len(times),
                    "passages": len(spans),
                    "longest_span": longest,
                    "last_crossing_inside": bool(inside[-1]),
                    "observation_time_in_gate": marked_time,
                }
            )
    return records


def figure(records, target):
    fig, (axis, product) = plt.subplots(1, 2, figsize=(8.2, 3.3), layout="constrained")
    for family, data in records.items():
        rows = data["width_scaled_runs"]
        x = [r["delta_over_width"] for r in rows]
        (line,) = axis.loglog(x, [r["elapsed"] for r in rows], "o-", label=family)
        product.semilogx(
            x,
            [r["multiplier_product"] for r in rows],
            "o-",
            label=family,
            color=line.get_color(),
        )
    axis.axhline(1100, linestyle=":", color=".4")
    axis.set(
        xlabel="Distance beyond fold / window width",
        ylabel="Passage time (model units)",
    )
    axis.annotate(
        "Survey horizon (1,100)",
        (0.00015, 1100),
        xytext=(0, 4),
        textcoords="offset points",
        fontsize=8,
    )
    axis.legend(frameon=False)
    product.axhline(2 * np.pi, color=".4", linestyle=":")
    product.set(
        xlabel="Distance beyond fold / window width",
        ylabel=r"Half-returns × $(1-\lambda_{\rm inside})$",
    )
    product.annotate(
        r"$2\pi$",
        (0.003, 2 * np.pi),
        xytext=(0, -13),
        textcoords="offset points",
        fontsize=9,
    )
    for ax in (axis, product):
        ax.invert_xaxis()
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(target, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "outputs/chua-fold-passage"
    )
    output = parser.parse_args().output_dir
    output.mkdir(parents=True, exist_ok=True)
    records = {}
    survey = json.loads((ROOT / "results/original/classic_scan.json").read_text())[
        "settings"
    ]
    for name, family, source in inputs():
        hmap = HalfMap()
        fold = centered_fold(hmap, family, source)
        rows = []
        for ratio in [1e-2, 1e-3, 1e-4, 1e-5]:
            delta = ratio * fold["width"]
            row, points, times = passage(hmap, fold, delta)
            rows.append(row)
            if ratio == 1e-3:
                direct = direct_check(fold, delta, points, times)
        common, _, _ = passage(hmap, fold, 1e-6)
        gate_checks = [
            passage(hmap, fold, 1e-4 * fold["width"], gate)[0] for gate in (0.1, 0.3)
        ]
        records[name] = {
            "fold": fold,
            "width_scaled_runs": rows,
            "common_detuning": common,
            "gate_checks": gate_checks,
            "full_field_check": direct,
            "survey_starts": survey_entries(fold, survey),
        }
        print(name, [(r["returns"], round(r["elapsed"], 1)) for r in rows], flush=True)
    file = output / "passages.json"
    file.write_text(
        json.dumps({"gate": GATE, "families": records}, indent=2, allow_nan=False)
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    figure(records, output / "passage-time.png")


if __name__ == "__main__":
    main()
