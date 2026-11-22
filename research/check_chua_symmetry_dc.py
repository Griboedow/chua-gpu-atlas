"""Measure both asymmetric cycles at the preselected symmetry-breaking detunings."""

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from chua_block_return import BlockReturn
from chua_fold_passage import inputs
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from trace_chua_window import Chua

ROOT = Path(__file__).resolve().parents[1]
FORECAST = ROOT / "results/chua-symmetry-dc/forecast.json"
SHA256 = "d4dc9f58a4d030856f19846dfece7c81d8b780cbcd000ff63b987d7223896810"


def fixed_point(block, beta, seed):
    p = seed.copy()
    for _ in range(15):
        u = block.evaluate(beta, p)
        residual = u["point"] - p
        if mp.norm(residual) < mp.mpf("1e-33"):
            return p, u
        p -= mp.lu_solve(u["matrix"] - mp.eye(2), residual)
    raise ArithmeticError("Parent did not close")


def period_two(block, beta, parent, right, left, distance):
    plus, minus = parent + right * distance / 2, parent - right * distance / 2
    for _ in range(25):
        u, v = block.evaluate(beta, plus), block.evaluate(beta, minus)
        residual = mp.matrix([*(u["point"] - minus), *(v["point"] - plus)])
        if mp.norm(residual) < mp.mpf("1e-33"):
            if mp.fdot(left, plus - minus) <= distance / 2:
                raise ArithmeticError(
                    "Newton converged to the parent or exchanged the labels"
                )
            return plus, minus, u, v
        jac = mp.zeros(4)
        jac[:2, :2], jac[2:, 2:] = u["matrix"], v["matrix"]
        jac[2:, :2] = jac[:2, 2:] = -mp.eye(2)
        correction = mp.lu_solve(jac, residual)
        plus -= correction[:2, :]
        minus -= correction[2:, :]
    raise ArithmeticError("Asymmetric orbit did not close")


def integrated_half_voltage(block, beta, point, durations):
    # Charge balance on each segment gives a second way to integrate x.
    state = mp.matrix([0, *point])
    total = mp.mpf(0)
    for region, dt in zip(block.regions, durations):
        matrix, center, *_ = block.matrices(beta)[region]
        end = center + mp.expm(matrix * dt) * (state - center)
        slope = -mp.mpf(8) / 7 if region == 0 else -mp.mpf(5) / 7
        offset = -mp.mpf(3) * region / 7
        total += (
            (state[2] - end[2]) / beta + (state[0] - end[0]) / 9 - offset * dt
        ) / (1 + slope)
        state = end
    return total


def full_field(beta, point, period):
    """Independent unsegmented integration of state, running integral and first variations."""
    model = Chua()

    def flow(t, v):
        u = v[:3]
        tangent = v[4:].reshape(3, 2)
        jacobian = model.jacobian(beta, abs(u[0]) < 1)
        return np.r_[model.field(beta, u), u[0], (jacobian @ tangent).ravel()]

    def section(t, v):
        return v[0]

    section.direction = 1
    initial = np.array([0.0, *[float(v) for v in point]])
    sol = solve_ivp(
        flow,
        (0, period + 0.02),
        np.r_[initial, 0.0, np.eye(3)[:, 1:].ravel()],
        method="DOP853",
        rtol=1e-12,
        atol=5e-15,
        max_step=0.04,
        events=section,
    )
    times = sol.t_events[0]
    if not sol.success or abs(times[-1] - period) > 1e-5:
        raise ArithmeticError("Failed to find the full return")
    v = sol.y_events[0][-1]
    field = model.field(beta, v[:3])
    tangent = v[4:].reshape(3, 2)
    derivative = (tangent - np.outer(field, tangent[0]) / field[0])[1:]
    lambdas = np.linalg.eigvals(derivative)
    if np.iscomplex(lambdas).any():
        raise ArithmeticError("Unexpected full-return multipliers")
    return {
        "period": float(times[-1]),
        "mean_x": float(v[3] / times[-1]),
        "multiplier": float(max(lambdas.real)),
        "closure": float(np.max(abs(v[:3] - initial))),
    }


def measure(block, forecast, row):
    beta = mp.mpf(forecast["beta_c"]) - mp.mpf(row["delta"])
    p0 = mp.matrix([mp.mpf(v) for v in forecast["parent_section_and_durations"][:2]])
    parent, center = fixed_point(block, beta, p0)
    right, left = mp.matrix(forecast["right"]), mp.matrix(forecast["left"])
    plus, minus, u, v = period_two(block, beta, parent, right, left, row["section_gap"])
    period = u["time"] + v["time"]
    parent_period = 2 * center["time"]
    mean = (
        integrated_half_voltage(block, beta, plus, u["durations"])
        - integrated_half_voltage(block, beta, minus, v["durations"])
    ) / period
    multiplier = max(mp.eig(v["matrix"] * u["matrix"], left=False, right=False))
    values = {
        "section_gap": float(mp.norm(plus - minus)),
        "positive_section_mean_x": float(mean),
        "daughter_minus_parent_period": float(period - parent_period),
        "one_minus_daughter_multiplier": float(1 - multiplier),
    }
    checks = {
        name: full_field(float(beta), point, float(T))
        for name, point, T in (
            ("plus", plus, period),
            ("minus", minus, period),
            ("parent", parent, parent_period),
        )
    }
    full_bias = checks["plus"]["mean_x"]
    errors = {k: abs(float(values[k] / row[k] - 1)) for k in values}
    agreement = {
        "mean_x": abs(full_bias - float(mean)),
        "opposite_dc": abs(checks["minus"]["mean_x"] + float(mean)),
        "period_difference": abs(
            checks["plus"]["period"]
            - checks["parent"]["period"]
            - float(period - parent_period)
        ),
        "daughter_multiplier": abs(checks["plus"]["multiplier"] - float(multiplier)),
    }
    if (
        max(check["closure"] for check in checks.values()) > 2e-6
        or max(agreement.values()) > 2e-6
    ):
        raise ArithmeticError(f"Independent integration disagrees {agreement}")
    return {
        "offset_fraction": row["offset_fraction"],
        "beta": float(beta),
        "measured": values,
        "relative_errors": errors,
        "plus": [float(v) for v in plus],
        "minus": [float(v) for v in minus],
        "full_field": checks,
        "full_field_absolute_differences": agreement,
    }


def plot(frozen, results, path):
    fig, (dc, error) = plt.subplots(1, 2, figsize=(8.5, 3.1), layout="constrained")
    for name, color in (("S01", "#aa4a30"), ("S08", "#155c72"), ("S11", "#567a44")):
        predicted = frozen["families"][name]["predictions"]
        actual = results[name]
        offset = [float(x["offset_fraction"]) for x in predicted]
        leading = np.array([x["positive_section_mean_x"] for x in predicted])
        dc.loglog(offset, abs(leading), "-", color=color)
        sign = "+" if leading[0] > 0 else "−"
        dc.plot(
            offset,
            [abs(x["full_field"]["plus"]["mean_x"]) for x in actual],
            "o",
            markerfacecolor="none",
            color=color,
            label=f"{name} ({sign} DC)",
        )
        error.loglog(
            offset,
            [100 * max(x["relative_errors"].values()) for x in actual],
            "o-",
            color=color,
        )
    dc.set(xlabel=r"$(\beta_c-\beta)/W$", ylabel=r"Absolute mean voltage $|\bar{x}|$")
    dc.legend(frameon=False)
    error.set(xlabel=r"$(\beta_c-\beta)/W$", ylabel="Largest prediction error (%)")
    for ax in (dc, error):
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/chua-symmetry-dc"
    )
    args = parser.parse_args()
    raw = FORECAST.read_bytes().replace(b"\r\n", b"\n")
    if hashlib.sha256(raw).hexdigest() != SHA256:
        raise ValueError("Forecast checksum changed")
    frozen = json.loads(raw)
    records = {}
    for name, family, _ in inputs():
        f = frozen["families"][name]
        block = BlockReturn({**family, "q": f["parent_section_and_durations"]}, 55)
        records[name] = [measure(block, f, row) for row in f["predictions"]]
        print(
            name,
            [max(x["relative_errors"].values()) for x in records[name]],
            flush=True,
        )
    worst_smallest = max(
        max(rows[-1]["relative_errors"].values()) for rows in records.values()
    )
    result = {
        "forecast_sha256": SHA256,
        "records": records,
        "worst_relative_error": max(
            max(row["relative_errors"].values())
            for rows in records.values()
            for row in rows
        ),
        "worst_at_smallest_offset": worst_smallest,
        "smallest_offset_target": 0.01,
        "all_smallest_offset_targets_met": worst_smallest < 0.01,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "checks.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    plot(frozen, records, args.output / "branch-voltage.png")


if __name__ == "__main__":
    main()
