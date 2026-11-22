"""First-exit counts with preparation jitter or one random kick per return."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtr

from observed_saddle_flights import Section


ROOT = Path(__file__).resolve().parents[1]
RADIUS = 0.02  # same projected section threshold as the voltage-start forecast
SIGMA = 1e-4  # rms initial displacement as a fraction of RADIUS
LIMIT = 64


def linear_counts(multiplier, offset, initial_sigma, kick_sigma, trials, rng):
    q = offset + initial_sigma * rng.standard_normal(trials)
    counts = np.zeros(trials, dtype=int)
    remaining = np.arange(trials)[abs(q) < 1]
    for n in range(1, LIMIT + 1):
        if not remaining.size:
            break
        q[remaining] = multiplier * q[remaining] + kick_sigma * rng.standard_normal(remaining.size)
        escaped = abs(q[remaining]) >= 1
        counts[remaining[escaped]] = n
        remaining = remaining[~escaped]
    if remaining.size:
        raise ArithmeticError(f"{remaining.size} censored exits")
    return counts


def compare_full_map(section, point, multiplier, left, right, kick_sigma, rng):
    # Identical random numbers in the scalar map and in the two-dimensional
    # return. The kick is along the unstable eigenvector after each full turn.
    differences = []
    for offset in [0., SIGMA, 10 * SIGMA]:
        for _ in range(96):
            q = offset
            state = point + RADIUS * q * right
            linear_exit = full_exit = None
            for n in range(1, LIMIT + 1):
                kick = kick_sigma * rng.standard_normal()
                if linear_exit is None:
                    q = multiplier * q + kick
                    if abs(q) >= 1:
                        linear_exit = n
                if full_exit is None:
                    _, state, _ = section.full_return(state)
                    state = state + RADIUS * kick * right
                    if abs(left @ (state - point)) >= RADIUS:
                        full_exit = n
                if full_exit is not None and linear_exit is not None:
                    break
            if full_exit is None or linear_exit is None:
                raise ArithmeticError("Incomplete paired exit")
            differences.append(full_exit - linear_exit)
    return {str(k): differences.count(k) for k in sorted(set(differences))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/chua-residence-noise")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    section = Section()
    point, _, values, left = section.saddle()
    multiplier = max(values)
    _, _, derivative = section.full_return(point)
    vals, vectors = np.linalg.eig(derivative)
    right = vectors[:, np.argmax(abs(vals))].real
    right /= left @ right
    kick_sigma = SIGMA * np.sqrt(multiplier**2 - 1)
    rng = np.random.default_rng(11222026)
    trials = 65536
    offsets = np.r_[0., np.logspace(-7, -1, 13)]
    rows = []
    preparation_cdf_error = 0.
    for name, prep, kick in [("preparation", SIGMA, 0.), ("per_return", 0., kick_sigma)]:
        for offset in offsets:
            counts = linear_counts(multiplier, offset, prep, kick, trials, rng)
            if kick == 0.:
                turns = np.arange(LIMIT + 1)
                boundary = multiplier ** (-turns)
                # With no subsequent kicks the magnitude grows monotonically.
                exact = ndtr((-boundary - offset) / prep) + ndtr((offset - boundary) / prep)
                empirical = np.cumsum(np.bincount(counts, minlength=LIMIT + 1)) / trials
                preparation_cdf_error = max(preparation_cdf_error, float(max(abs(exact - empirical))))
            p10, p50, p90 = np.quantile(counts, [0.1, 0.5, 0.9], method="inverted_cdf")
            rows.append(dict(source=name, offset=offset, p10=int(p10), median=int(p50),
                             p90=int(p90), mean=float(np.mean(counts))))
    comparison = compare_full_map(section, point, multiplier, left, right, kick_sigma, rng)
    if preparation_cdf_error > 0.015:
        raise ArithmeticError(f"Preparation CDF discrepancy: {preparation_cdf_error}")
    metadata = dict(seed=11222026, trials_per_row=trials, saddle_section_yz=point.tolist(),
                    full_return_multiplier=multiplier, left_unit=left.tolist(),
                    right_normalized=right.tolist(), section_exit_radius=RADIUS,
                    rms_initial_over_radius=SIGMA, rms_kick_over_radius=kick_sigma,
                    max_preparation_cdf_discrepancy=preparation_cdf_error,
                    paired_full_minus_linear_exit_counts=comparison,
                    paired_offsets=[0., SIGMA, 10 * SIGMA], pairs_per_offset=96,
                    cap_returns=LIMIT)
    (args.output / "settings.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    with (args.output / "counts.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(6.6, 3.6), layout="constrained")
    for mode, label, color, marker in [("preparation", "Preparation jitter", "#205988", "o"),
                                       ("per_return", "Kick after each return", "#b34727", "x")]:
        data = [r for r in rows if r["source"] == mode and r["offset"]]
        xs = [r["offset"] for r in data]
        ax.fill_between(xs, [r["p10"] for r in data], [r["p90"] for r in data],
                        color=color, alpha=0.10, linewidth=0)
        ax.plot(xs, [r["median"] for r in data], color=color, marker=marker,
                markersize=5, label=label, lw=1.5, linestyle="-" if mode=="preparation" else "--")
    ax.axvline(SIGMA, color="0.55", lw=1, ls=":")
    ax.set(xscale="log", xlabel="Imposed offset / exit radius",
           ylabel="Full returns to first exit", ylim=(0, 18))
    ax.legend(frameon=False, loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(args.output / "exit-counts.png", dpi=185)
    plt.close(fig)
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
