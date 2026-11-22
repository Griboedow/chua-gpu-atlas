"""Recover y and z from one period of evenly sampled x, or replay the atlas benchmark."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import expm
from scipy.signal import resample

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "results/chua-knot-atlas"
RATES = (*range(2, 13), 16, 24, 32)
PHASES = tuple(i / 16 for i in range(16))


def reconstruct(x, period, beta, gamma=0.0, count=None):
    """Return y,z for the periodic trigonometric interpolant of x, endpoint excluded."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("Supply at least three finite samples of one voltage channel.")
    if (
        not np.isfinite([period, beta, gamma]).all()
        or period <= 0
        or beta <= 0
        or gamma < 0
    ):
        raise ValueError("Expected positive period and beta, and nonnegative gamma.")
    count = len(x) if count is None else count
    if count < len(x):
        raise ValueError("Output grid must have at least as many points as the input.")
    # Upsample the voltage before filtering to retain both quadratures of its
    # Nyquist cosine when evaluating between the original sampling instants.
    signal = resample(x, count) if count != len(x) else x
    s = 2j * np.pi * np.fft.rfftfreq(count, period / count)
    denominator = (s + 1) * (s + gamma) + beta
    drive = np.fft.rfft(signal) / denominator
    return np.column_stack(
        (
            np.fft.irfft((s + gamma) * drive, n=count),
            np.fft.irfft(-beta * drive, n=count),
        )
    )


def time_domain_check():
    period, beta = 3.75, 14.3
    omega = 2 * np.pi / period

    def voltage(t):
        return 0.23 + np.cos(omega * t) - 0.4 * np.sin(2 * omega * t)

    errors = {}
    for gamma in (0.0, 0.08):
        matrix = np.array([[-1.0, 1.0], [-beta, -gamma]])

        def field(t, yz, matrix=matrix):
            return matrix @ yz + [voltage(t), 0.0]

        options = {"method": "DOP853", "rtol": 2e-13, "atol": 2e-14}
        free_start = solve_ivp(field, (0.0, period), [0.0, 0.0], **options)
        start = np.linalg.solve(np.eye(2) - expm(matrix * period), free_start.y[:, -1])
        t = np.arange(601) * period / 601
        check = solve_ivp(field, (0.0, period), start, t_eval=t, **options)
        if not (free_start.success and check.success):
            raise ArithmeticError("Time-domain control failed")
        yz = reconstruct(
            voltage(np.arange(19) * period / 19), period, beta, gamma, len(t)
        )
        errors[str(gamma)] = float(np.abs(yz - check.y.T).max())
    if max(errors.values()) > 1e-9:
        raise ArithmeticError(
            "Periodic reconstruction disagrees with time-domain integration"
        )
    return errors


def reference_orbit(card, precision):
    stem = ARCHIVE / f"{card['id']}-{precision}"
    with np.load(stem.with_suffix(".npz")) as data:
        t, xyz = data["t"], data["xyz"]
    orbit = CubicSpline(t, xyz, axis=0, bc_type="periodic", extrapolate="periodic")
    meta = json.loads(stem.with_suffix(".json").read_text())
    return orbit, t[-1], float(meta["beta"])


def benchmark_family(card, precision):
    orbit, period, beta = reference_orbit(card, precision)
    peaks = card["peaks"]
    rows = []
    designs = [(rate, 1.0) for rate in RATES] + [(32, 1.01)]
    for rate, beta_ratio in designs:
        n = rate * peaks
        fine = max(1024, 128 * peaks)
        errors = []
        for phase in PHASES:
            offset = phase * period / n
            x = orbit(np.arange(n) * period / n + offset)[:, 0]
            reference = orbit(np.arange(fine) * period / fine + offset)[:, 1:]
            estimate = reconstruct(x, period, beta_ratio * beta, count=fine)
            rms = np.sqrt(np.mean((estimate - reference) ** 2, axis=0))
            errors.append(rms / np.std(reference, axis=0))
        worst = np.max(errors, axis=0)
        rows.append(
            {
                "family": card["id"],
                "peaks": peaks,
                "samples": n,
                "samples_per_peak": rate,
                "beta_ratio": beta_ratio,
                "y_error": worst[0],
                "z_error": worst[1],
            }
        )
    return rows


def run_atlas(output):
    cards = json.loads((ROOT / "results/chua-layer-cards/cards.json").read_text())[
        "families"
    ]
    if len({c["id"] for c in cards}) != 34:
        raise ValueError("Expected the 34 archived reference families.")
    rows, controls = [], []
    for card in sorted(cards, key=lambda c: c["id"]):
        rows.extend(benchmark_family(card, "d70-n192"))
        controls.extend(benchmark_family(card, "d60-n96"))
    report = {
        "phases": PHASES,
        "time_domain_max_error": time_domain_check(),
        "coarse_reference_max_change": max(
            abs(a[k] - b[k])
            for a, b in zip(rows, controls)
            for k in ("y_error", "z_error")
        ),
        "rates": {},
    }
    for rate in RATES:
        subset = [
            r for r in rows if r["samples_per_peak"] == rate and r["beta_ratio"] == 1
        ]
        worst = max(subset, key=lambda r: max(r["y_error"], r["z_error"]))
        report["rates"][rate] = {
            "worst_family": worst["family"],
            "worst_error": max(worst["y_error"], worst["z_error"]),
            "families_below_one_percent": sum(
                1 for r in subset if max(r["y_error"], r["z_error"]) < 0.01
            ),
        }
    biased = [max(r["y_error"], r["z_error"]) for r in rows if r["beta_ratio"] > 1]
    report["one_percent_beta_error"] = {"minimum": min(biased), "maximum": max(biased)}
    output.mkdir(parents=True, exist_ok=True)
    example = next(c for c in cards if c["id"] == "S11")
    orbit, period, _ = reference_orbit(example, "d70-n192")
    n = 9 * example["peaks"]
    np.savetxt(output / "S11-x.csv", orbit(np.arange(n) * period / n)[:, 0])
    with (output / "sampling.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    fig, ax = plt.subplots(figsize=(6.7, 3.2), layout="constrained")
    for family in sorted({r["family"] for r in rows}):
        sample = [r for r in rows if r["family"] == family and r["beta_ratio"] == 1]
        ax.plot(
            RATES,
            [100 * max(r["y_error"], r["z_error"]) for r in sample],
            color="#7599ac",
            alpha=0.4,
            linewidth=0.8,
        )
    ax.plot(
        RATES,
        [100 * r["worst_error"] for r in report["rates"].values()],
        "o-",
        color="#145369",
        label="Worst of 34 families",
        markersize=4,
    )
    ax.axhline(1, color="#ad5430", linestyle="--", label="1% error")
    ax.set(
        yscale="log",
        xticks=[2, 4, 6, 8, 10, 12, 16, 24, 32],
        xlabel="Voltage samples per peak (period average)",
        ylabel="Larger of y and z RMS errors (%)",
    )
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "sampling.png", dpi=170)
    plt.close(fig)
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    study = sub.add_parser("atlas", help="replay both archived orbit resolutions")
    study.add_argument(
        "--output", type=Path, default=ROOT / "outputs/chua-voltage-reconstruction"
    )
    trace = sub.add_parser("trace", help="process a one-column CSV of x with no header")
    trace.add_argument("input", type=Path)
    trace.add_argument("output", type=Path)
    trace.add_argument(
        "--period", type=float, required=True, help="full period in model time"
    )
    trace.add_argument("--beta", type=float, required=True)
    trace.add_argument("--gamma", type=float, default=0.0)
    args = parser.parse_args()
    if args.mode == "atlas":
        run_atlas(args.output)
        return
    x = np.loadtxt(args.input, delimiter=",", ndmin=2)
    if x.shape[1] != 1:
        parser.error("Expected one column containing x.")
    if args.input.resolve() == args.output.resolve():
        parser.error("Choose a separate output file.")
    x = x[:, 0]
    yz = reconstruct(x, args.period, args.beta, args.gamma)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    t = np.arange(len(x)) * args.period / len(x)
    np.savetxt(
        args.output,
        np.column_stack((t, x, yz)),
        delimiter=",",
        header="t,x,y,z",
        comments="",
    )


if __name__ == "__main__":
    main()
