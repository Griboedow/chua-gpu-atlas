"""Predict voltage bias, period and stability at three symmetry-breaking cycles."""

import argparse
import json
from pathlib import Path

import mpmath as mp

from chua_block_return import BlockReturn
from chua_fold_passage import inputs
from chua_precision_flights import PrecisionFlights

ROOT = Path(__file__).resolve().parents[1]
DIGITS = 55
OFFSETS = ("0.01", "0.001", "0.0001")


def half_area(block, beta, point, durations):
    """Integral of x along successive affine flights from the upward section."""
    state = mp.matrix([0, *point])
    result = mp.mpf(0)
    mats = block.matrices(beta)
    for region, dt in zip(block.regions, durations):
        matrix, center, eig, v, w = mats[region]
        end = center + (v * mp.diag([mp.exp(e * dt) for e in eig]) * w * (state-center)).apply(mp.re)
        result += center[0] * dt + mp.lu_solve(matrix, end-state)[0]
        state = end
    return result


def coefficients(block, beta, p, h):
    base = block.evaluate(beta, p)
    a = base["matrix"]
    eigen, vectors = mp.eig(a)
    j = min(range(2), key=lambda i: abs(eigen[i]+1))
    if abs(eigen[j]+1) > mp.mpf("1e-28") or abs(eigen[1-j]) >= 1:
        raise ArithmeticError("A simple neutral flip and a stable second direction are required")
    r = vectors[:, j].apply(mp.re)
    r /= mp.norm(r)
    if r[0] < 0:
        r = -r
    left = (vectors**-1)[j, :].apply(mp.re)
    left /= (left*r)[0]
    plus = block.evaluate(beta, p+h*r)
    minus = block.evaluate(beta, p-h*r)
    d2 = (plus["matrix"] - minus["matrix"])*r / (2*h)
    d3 = (plus["matrix"] - 2*a + minus["matrix"])*r / h**2
    midpoint_shift = mp.lu_solve(mp.eye(2)-a, d2/2)
    mixed = ((block.evaluate(beta, p+h*midpoint_shift)["matrix"]
              - block.evaluate(beta, p-h*midpoint_shift)["matrix"])*r / (2*h))
    cubic = (left*(mixed + d3/6))[0]
    h_beta = (block.evaluate(beta+h, p)["point"]-block.evaluate(beta-h, p)["point"]) / (2*h)
    p_beta = mp.lu_solve(mp.eye(2)-a, h_beta)
    a_beta = (block.evaluate(beta+h, p+h*p_beta)["matrix"]
              - block.evaluate(beta-h, p-h*p_beta)["matrix"])/(2*h)
    slope = (left*a_beta*r)[0]
    gain = mp.sqrt(slope/cubic)
    d_area = (half_area(block, beta, p+h*r, plus["durations"])
              - half_area(block, beta, p-h*r, minus["durations"]))/(2*h)
    d2_time = (plus["time"]-2*base["time"]+minus["time"])/h**2
    d_time_midpoint = (block.evaluate(beta, p+h*midpoint_shift)["time"]
                       -block.evaluate(beta, p-h*midpoint_shift)["time"])/(2*h)
    return {
        "right": list(r), "left": list(left), "base": list(p), "cubic": cubic,
        "lambda_beta": slope, "gap_gain": 2*gain, "dc_gain": gain*d_area/base["time"],
        "period_difference_gain": gain**2*(d2_time+2*d_time_midpoint),
        "stability_deficit_gain": 4*slope,
    }


def forecast():
    records = {}
    for name, family, source in inputs():
        boundary = source["symmetry"]
        solver = PrecisionFlights(family["regions"], family["targets"], True, DIGITS)
        refined = solver.boundary(boundary["beta"], boundary["q"], -1)
        beta = mp.mpf(refined["beta"])
        p = mp.matrix([mp.mpf(x) for x in refined["q"][:2]])
        block = BlockReturn({**family, "q": refined["q"]}, DIGITS)
        # Small absolute steps retain high precision even in the very narrow S01 window.
        original = coefficients(block, beta, p, mp.mpf("1e-10"))
        coeff = coefficients(block, beta, p, mp.mpf("1e-12"))
        measured = ["gap_gain", "dc_gain", "period_difference_gain", "stability_deficit_gain"]
        convergence = max(abs(original[k]/coeff[k]-1) for k in measured)
        if convergence > mp.mpf("1e-6"):
            raise ArithmeticError(f"Derivative step change in {name}: {convergence}")
        width = mp.mpf(str(source["fold"]["beta"]))-beta
        rows = []
        for ratio in OFFSETS:
            delta = width*mp.mpf(ratio)
            rows.append({
                "offset_fraction": ratio, "delta": float(delta),
                "section_gap": float(coeff["gap_gain"]*mp.sqrt(delta)),
                "positive_section_mean_x": float(coeff["dc_gain"]*mp.sqrt(delta)),
                "daughter_minus_parent_period": float(coeff["period_difference_gain"]*delta),
                "one_minus_daughter_multiplier": float(coeff["stability_deficit_gain"]*delta),
            })
        records[name] = {
            "beta_c": refined["beta"], "parent_section_and_durations": refined["q"],
            "width": float(width), "derivative_step_relative_change": float(convergence),
            "right": [float(v) for v in coeff["right"]], "left": [float(v) for v in coeff["left"]],
            **{k: float(coeff[k]) for k in ["cubic", "lambda_beta", *measured]},
            "predictions": rows,
        }
    return {
        "coordinate": "upward x=0 section (y,z); on the plus daughter the signed half-gap has positive projection on the left eigenvector",
        "parameters": {"alpha":9, "m0":"-8/7", "m1":"-5/7", "gamma":0},
        "detuning": "delta = beta_c - beta, in units of the symmetric fold-to-symmetry-breaking width",
        "score": "Absolute relative errors for the four signed predictions. Target below 1 percent at the smallest detuning in each family.",
        "digits": DIGITS, "families": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"outputs/chua-symmetry-dc/forecast.json")
    args = parser.parse_args()
    data = forecast()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")
    for name, record in data["families"].items():
        print(name, record["dc_gain"], record["predictions"][-1], flush=True)


if __name__ == "__main__":
    main()
