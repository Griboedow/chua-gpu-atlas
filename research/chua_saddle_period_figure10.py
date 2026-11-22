"""Measure the change in peak spacing within each published Figure 10 trace."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ("https://media.springernature.com/lw1200/springer-static/image/"
          "art%3A10.1007%2Fs11071-022-08078-y/MediaObjects/11071_2022_8078_Fig10_HTML.png")
DIGEST = "a3c79cd7dda92071297d069293a23c5dfe9173b930f2b1d6d2eb44ca6d5246b0"


def peak_centers(pixels, panel, cutoff):
    # Fixed interior crops remove labels and both plot frames. Only the upper
    # part of each positive lobe is needed for a horizontal peak position.
    top, bottom = (17, 160) if panel == "a" else (397, 527)
    left, right = 65, 706
    crop = pixels[top:bottom, left:right]
    ink = (crop.max(axis=2) < cutoff) & (np.ptp(crop, axis=2) < 25)
    if panel == "b":
        ink[:15, :30] = False  # the printed U1 label
    envelope = np.where(ink.any(axis=0), top + ink.argmax(axis=0), bottom)
    peaks, _ = find_peaks(bottom - envelope, prominence=20, distance=22)
    centers = []
    for peak in peaks:
        dx = np.arange(-3, 4)
        a, b, _ = np.polyfit(dx, envelope[peak + dx], 2)
        vertex = np.clip(-b / (2 * a), -3, 3) if a > 0 else 0.
        centers.append(float(left + peak + vertex))
    return centers


def plot_spacings(rows, target):
    fig, ax = plt.subplots(figsize=(7.0, 2.9), layout="constrained")
    for panel, label, style, color in [
        ("b", "Hardware (b)", "o-", "#245e73"),
        ("a", "Simulation (a)", "s--", "#b05c38"),
    ]:
        positions = [x for source_panel, _, x in rows if source_panel == panel]
        spacings = np.diff(positions)
        relative = spacings / np.mean(spacings[:4])
        ax.plot(np.arange(1, 1 + len(spacings)), relative, style, color=color,
                markersize=4, linewidth=1.5, label=label)
    ax.axhline(1, color="0.65", linewidth=0.8, zorder=0)
    ax.set(xlabel="Successive intervals in each published trace",
           ylabel="Spacing / early mean", xticks=[1, 4, 8, 12, 16, 19])
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.savefig(target, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image", type=Path, help="Published 1200-pixel PNG: " + SOURCE)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/chua-saddle-period")
    args = parser.parse_args()
    if hashlib.sha256(args.image.read_bytes()).hexdigest() != DIGEST:
        raise ValueError("The pixel coordinates require the cited 1200-pixel source image.")
    with Image.open(args.image) as image:
        pixels = np.asarray(image.convert("RGB"))
    rows, checks = [], []
    for panel, expected in [("a", 17), ("b", 20)]:
        for cutoff in (60, 110, 160):
            centers = peak_centers(pixels, panel, cutoff)
            if len(centers) != expected:
                raise ValueError(f"Unexpected peak count for panel {panel} at cutoff {cutoff}")
            if cutoff == 110:
                rows.extend((panel, i, x) for i, x in enumerate(centers, 1))
            for count in (3, 4, 5):
                first = (centers[count] - centers[0]) / count
                last = (centers[-1] - centers[-count - 1]) / count
                checks.append(dict(panel=panel, cutoff=cutoff, intervals=count,
                                   first_spacing_px=first, last_spacing_px=last,
                                   ratio=last / first))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "figure10-peaks.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["panel", "interior_peak", "horizontal_pixel"])
        writer.writerows(rows)
    results = dict(source=SOURCE, source_sha256=DIGEST, primary_cutoff=110,
                   primary_intervals=4, checks=checks)
    (args.output_dir / "figure10-ratios.json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8")
    plot_spacings(rows, args.output_dir / "figure10-spacings.png")
    for row in checks:
        if row["cutoff"] == 110 and row["intervals"] == 4:
            print(f"Panel {row['panel']}: late/early spacing = {row['ratio']:.4f}")


if __name__ == "__main__":
    main()
