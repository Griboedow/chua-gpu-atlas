"""Predict voltage harmonics from diode crossings, without reading a voltage trace."""

import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from chua_precision_flights import PrecisionFlights

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "results/chua-knot-atlas"
BANDS = (4, 16, 64)
FAMILIES = ("S08", "E06", "A04")


def crossing_data(family, peaks):
    path = ARCHIVE / f"{family}-d70-n192.json"
    row = json.loads(path.read_text())
    flow = PrecisionFlights(row["regions"], row["targets"], row["symmetric"], 70)
    beta = mp.mpf(row["beta"])
    q = [mp.mpf(v) for v in row["q"]]
    arcs = flow.evaluate(beta, q, gradients=False)["arcs"]
    period = mp.mpf(row["period"])
    time = mp.mpf(0)
    hits = []
    for arc, following in zip(arcs[:-1], arcs[1:]):
        time += arc["duration"]
        x, y, z = arc["end"]
        wall = int(mp.nint(x))
        assert abs(wall) == 1 and abs(x - wall) < mp.mpf("1e-35")
        outgoing = int(arc["region"] == 0) - int(following["region"] == 0)
        assert abs(outgoing) == 1
        velocity = [
            float(flow.alpha * (y - (1 + flow.m0) * wall)),
            float(x - y + z),
            float(-beta * y),
        ]
        hits.append(
            {
                "t": float(time),
                "wall": wall,
                "direction": outgoing,
                "velocity": velocity,
            }
        )
    if row["symmetric"]:
        # The cut at x=0 has the same affine slope on both sides.
        hits.extend(
            [
                {
                    "t": h["t"] + float(period) / 2,
                    "wall": -h["wall"],
                    "direction": h["direction"],
                    "velocity": [-v for v in h["velocity"]],
                }
                for h in list(hits)
            ]
        )
    jumps = [
        -float(flow.alpha * (flow.m1 - flow.m0)) * h["direction"] * h["velocity"][0]
        for h in hits
    ]
    for h, jump in zip(hits, jumps):
        h["curvature_jump"] = jump
    return {
        "family": family,
        "peaks": peaks,
        "beta": float(beta),
        "alpha": float(flow.alpha),
        "m0": float(flow.m0),
        "m1": float(flow.m1),
        "gamma": 0.0,
        "period": float(period),
        "symmetric": row["symmetric"],
        "crossings": hits,
        "scaled_band_power_limit": sum(v * v for v in jumps),
    }


def coefficients(row, harmonics):
    """Exact finite-frequency expression and leading high-frequency term for X_k."""
    harmonics = np.asarray(harmonics)
    if np.any(harmonics <= 0):
        raise ValueError("Expected positive harmonic numbers.")
    s = 2j * np.pi * harmonics / row["period"]
    resolvents = []
    for slope in (row["m0"], row["m1"]):
        A = np.array(
            [
                [-row["alpha"] * (1 + slope), row["alpha"], 0.0],
                [1.0, -1.0, 1.0],
                [0.0, -row["beta"], -row["gamma"]],
            ]
        )
        resolvents.append(np.linalg.inv(s[:, None, None] * np.eye(3) - A)[:, 0, :])
    delta = resolvents[1] - resolvents[0]
    hits = row["crossings"]
    phase = np.exp(-s[:, None] * [h["t"] for h in hits])
    velocity = np.array([h["velocity"] for h in hits])
    weights = [h["direction"] for h in hits]
    full = np.sum(phase * (delta @ velocity.T) * weights, axis=1) / (row["period"] * s)
    leading = (phase @ np.array([h["curvature_jump"] for h in hits])) / (
        row["period"] * s**3
    )
    return full, leading


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs/chua-diode-spectrum/forecast.json",
    )
    args = parser.parse_args()
    cards = json.loads((ROOT / "results/chua-layer-cards/cards.json").read_text())[
        "families"
    ]
    peaks = {r["id"]: r["peaks"] for r in cards}
    records = [crossing_data(f, peaks[f]) for f in FAMILIES]
    for r in records:
        # Retain concrete predictions even without replaying the matrix formula.
        k = 2 * r["peaks"] + 1
        full, _ = coefficients(r, [k])
        r["example_harmonic"] = {
            "k": k,
            "real": float(full[0].real),
            "imag": float(full[0].imag),
        }
    info = {
        "records": records,
        "bands_per_peak": BANDS,
        "checks": {
            "low_band": "1 through 6p inclusive; relative complex L2 error below 1e-5 for the finite-frequency formula",
            "tail_band": "64p through 128p-1; relative complex L2 error below 5% for the leading term",
            "tail_power": "same upper band; scaled power within 5% of the infinite-band limit",
        },
        "validation": "Uniform Fourier transform of x alone. No spectral samples or fitted coefficients in this forecast.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(info, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    main()
