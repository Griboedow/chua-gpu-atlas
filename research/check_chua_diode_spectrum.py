"""Compare the frozen crossing predictions with Fourier transforms of full-field trajectories."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter
from predict_chua_diode_spectrum import ARCHIVE, ROOT, coefficients
from scipy.integrate import solve_ivp
from trace_chua_window import Chua

PREDICTION = ROOT / "results/chua-diode-spectrum/forecast.json"


def waveform(row, count, method, step):
    # Separate orbit reconstruction at 60 digits supplies only the starting
    # state, beta and period. The time integrator receives no crossing times.
    info = json.loads((ARCHIVE / f"{row['family']}-d60-n96.json").read_text())
    initial = np.array([0.0, *map(float, info["q"][:2])])
    period, beta = float(info["period"]), float(info["beta"])
    circuit = Chua()
    t = np.linspace(0.0, period, count + 1)
    result = solve_ivp(
        lambda _, u: circuit.field(beta, u),
        (0.0, period),
        initial,
        t_eval=t,
        method=method,
        rtol=4e-14,
        atol=2e-15,
        max_step=step,
    )
    if not result.success:
        raise ArithmeticError(result.message)
    return np.fft.rfft(result.y[0, :-1]) / count, float(
        np.max(abs(result.y[:, -1] - initial))
    )


def relative_error(prediction, measurement):
    return float(np.linalg.norm(prediction - measurement) / np.linalg.norm(measurement))


def evaluate(row, bands):
    p = row["peaks"]
    n = 2 ** int(np.ceil(np.log2(32 * 2 * max(bands) * p)))
    measured, gap = waveform(row, n, "DOP853", 0.015)
    control, control_gap = waveform(row, n // 2, "Radau", 0.025)
    k = np.arange(1, 2 * max(bands) * p)
    full, leading = coefficients(row, k)
    low = k <= 6 * p
    checks = []
    for rate in bands:
        select = (k >= rate * p) & (k < 2 * rate * p)
        modes = k[select]
        s = 2j * np.pi * modes / row["period"]
        observed = measured[modes]
        band_ratio = (
            np.mean(abs(row["period"] * s**3 * observed) ** 2)
            / row["scaled_band_power_limit"]
        )
        checks.append(
            {
                "first_mode": int(modes[0]),
                "last_mode": int(modes[-1]),
                "per_peak": rate,
                "full_error": relative_error(full[select], observed),
                "leading_error": relative_error(leading[select], observed),
                "independent_solver_difference": relative_error(
                    control[modes], observed
                ),
                "band_power_ratio": float(band_ratio),
            }
        )
    low_error = relative_error(full[low], measured[k[low]])
    return {
        "family": row["family"],
        "waveform_samples": n,
        "dop853_closure": gap,
        "radau_closure": control_gap,
        "low_band_error": low_error,
        "bands": checks,
        "tests": {
            "low_band": low_error < 1e-5,
            "tail_band": checks[-1]["leading_error"] < 0.05,
            "tail_power": abs(checks[-1]["band_power_ratio"] - 1) < 0.05,
        },
    }


def plot(measurements, target):
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.7, 3.1), layout="constrained")
    for r in measurements:
        x = [b["per_peak"] for b in r["bands"]]
        left.plot(
            x, [100 * b["leading_error"] for b in r["bands"]], "o-", label=r["family"]
        )
        right.plot(x, [b["band_power_ratio"] for b in r["bands"]], "o-")
    left.set(ylabel="Complex voltage error (%)", yscale="log")
    right.set(ylabel="Measured power / limiting value")
    right.axhline(1, color=".6", ls="--")
    for ax in (left, right):
        ax.set(
            xscale="log",
            xticks=x,
            xticklabels=x,
            xlabel="First harmonic / peaks per period",
        )
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.spines[["top", "right"]].set_visible(False)
    left.legend(frameon=False)
    fig.savefig(target, dpi=175)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--forecast", type=Path, default=PREDICTION)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/chua-diode-spectrum"
    )
    args = parser.parse_args()
    content = args.forecast.read_bytes().replace(b"\r\n", b"\n")
    forecast = json.loads(content)
    records = [evaluate(r, forecast["bands_per_peak"]) for r in forecast["records"]]
    args.output.mkdir(parents=True, exist_ok=True)
    result = {
        "forecast_sha256": hashlib.sha256(content).hexdigest(),
        "records": records,
    }
    (args.output / "measurements.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    plot(records, args.output / "spectrum.png")
    for r in records:
        print(r["family"], r["tests"], r["bands"][-1], flush=True)
    if not all(all(r["tests"].values()) for r in records):
        raise SystemExit("At least one spectral prediction missed its threshold.")


if __name__ == "__main__":
    main()
