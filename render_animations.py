"""Render measured data with a slow, level camera sweep.

The animation changes the view, never the parameters or the integration time.
Small illustration datasets are included so this works without a GPU.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parent
BG = "#091321"


def camera(phase, sweep):
    """Ease through a horizontal arc; keep the tilt and screen orientation fixed."""
    return 23., -55. + sweep/2*np.sin(2*np.pi*phase), 0.


def configure(ax, labels):
    ax.set_facecolor(BG)
    ax.set_box_aspect((1, 1, 1))
    ax.set_proj_type("ortho")
    ax.set(xlabel=labels[0], ylabel=labels[1], zlabel=labels[2])
    ax.tick_params(colors="#b7c9dc", labelsize=8, pad=1)
    for a in [ax.xaxis, ax.yaxis, ax.zaxis]:
        a.set_pane_color((.06, .11, .18, 1))
        a._axinfo["grid"]["color"] = (.25, .34, .45, .45)
        a.label.set_color("#dce7f3")


def render(kind, data_dir, output, frame_count, fps, sweep):
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": "#e1edf8", "font.size": 11})
    fig = plt.figure(figsize=(9.6, 6.4), dpi=100, facecolor=BG)
    ax = fig.add_axes([.08, .15, .78, .7], projection="3d")
    labels = ("x", "y", "z") if kind == "trajectory" else ("α", "β", "m₀")
    configure(ax, labels)
    if kind in ["extent", "expanded"]:
        expanded = kind == "expanded"
        p = np.load(data_dir / ("expanded_candidates.npz" if expanded else "extent_candidates.npz"))["parameters"]
        if expanded:
            previous = np.load(data_dir / "extent_candidates.npz")["parameters"]
            ax.scatter(*previous[:, :3].T, c="#9caac0", s=3, alpha=.25, linewidths=0)
            vertices = np.array([[a, a*r, m] for a in [2, 256] for r in [.25, 64] for m in [-8, -1.003]])
            for i in range(8):
                for bit in [1, 2, 4]:
                    if i < (j := i ^ bit):
                        ax.plot(*vertices[[i, j]].T, c="#f6be57", lw=1)
        dots = ax.scatter(*p[:, :3].T, c=p[:, 2], cmap="viridis", s=7, alpha=.85,
                          linewidths=0, depthshade=False)
        ax.set(xlim=(0, 1024 if expanded else 256),
               ylim=(0, np.ceil(p[:, 1].max()/25000)*25000 if expanded else 12500),
               zlim=(-29 if expanded else -8, -1))
        cb = fig.colorbar(dots, cax=fig.add_axes([.9, .29, .015, .35]))
        cb.set_label("Inner slope m₀", color="#dce7f3")
        cb.ax.tick_params(colors="#b7c9dc")
        title = "Four times farther in parameter space" if expanded else "How far does the chaotic region extend?"
        subtitle = f"{len(p):,} repeated GPU candidates · fixed m₁ = −5/7 · linear coordinates"
        footer = ("Gold: previous search limits. Finite-time candidates; no proof of global unboundedness."
                  if expanded else "Sparse finite-time evidence; neither a closed solid nor a proof of unboundedness.")
    elif kind == "window":
        d = np.load(data_dir / "window_regular_mesh.npz")
        ax.add_collection3d(Poly3DCollection(d["vertices"][d["faces"]], facecolors="#ee8caa",
                                            linewidths=0, alpha=.9, shade=True))
        ax.set(xlim=(8.94, 9.06), ylim=(14.24, 14.40), zlim=(-8/7-.006, -8/7+.006))
        ax.ticklabel_format(useOffset=False)
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            axis.set_major_locator(MaxNLocator(4))
        ax.plot([9, 9], [14.24, 14.40], [-8/7]*2, color="#ffffff", lw=1.5)
        ax.scatter([9], [14.32], [-8/7], c="#f6be57", s=40, depthshade=False)
        title = "Thin regular windows inside chaos"
        subtitle = "41³ local grid · near-zero / nonpositive LLE for all three starts"
        footer = "Unsmoothed mesh, clipped at the scan box. Gold point: the five-peak example."
    else:
        d = np.load(data_dir / "classic_trajectory.npz")
        xyz = d["xyz"]
        ax.plot(*xyz.T, color="#57d4c6", lw=.35, alpha=.8)
        ax.set_box_aspect((1, .65, 1))
        title = "The Chua double scroll"
        subtitle = "Phase space (x, y, z) · α = 9, β = 100/7, m₀ = −8/7, m₁ = −5/7"
        footer = "One saved trajectory; only the camera moves. This is not parameter space."
    fig.text(.06, .95, title, fontsize=21, weight="bold")
    fig.text(.06, .905, subtitle, fontsize=10, color="#b7c9dc")
    fig.text(.06, .09, f"Gentle {sweep:g}° camera sweep · {frame_count/fps:g} s loop · fixed tilt",
             fontsize=12, color="#f6be57")
    fig.text(.06, .045, footer, fontsize=9, color="#b7c9dc")
    frames = []
    palette = None
    for i in range(frame_count):
        ax.view_init(*camera(i/frame_count, sweep))
        fig.canvas.draw()
        rgb = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
        if palette is None:
            palette = Image.fromarray(rgb).quantize(colors=128, method=Image.Quantize.MEDIANCUT,
                                                    dither=Image.Dither.NONE)
        frame = Image.fromarray(rgb).quantize(palette=palette, dither=Image.Dither.NONE)
        frames.append(frame)
        if i == 0:
            Image.fromarray(rgb).save(output / f"{kind}-poster.png")
        if (i+1) % max(1, frame_count//4) == 0:
            print(f"{kind}: rendered {i+1}/{frame_count} frames", flush=True)
    # GIF durations are multiples of 10 ms; distribute rounding error over frames.
    durations = [10*(round((i+1)*100/fps)-round(i*100/fps)) for i in range(len(frames))]
    target = output / f"{kind}-rotation.gif"
    frames[0].save(target, save_all=True, append_images=frames[1:], loop=0,
                   duration=durations, disposal=1, optimize=True)
    plt.close(fig)
    with Image.open(target) as check:
        encoded_frames = check.n_frames
        encoded_duration = 0
        for i in range(encoded_frames):
            check.seek(i)
            encoded_duration += check.info["duration"]
        # Pillow may merge identical frames at the smooth turning points.
        if not 1 < encoded_frames <= len(frames) or encoded_duration != sum(durations):
            raise RuntimeError("Unexpected GIF frame count or playback duration")
    return {"file": target.name, "frames": encoded_frames, "rendered_frames": len(frames), "fps": fps,
            "seconds": sum(durations)/1000, "bytes": target.stat().st_size,
            "camera_motion": "sinusoidal azimuth sweep", "sweep_degrees": sweep,
            "azimuth_center_degrees": -55, "elevation_degrees": 23, "roll_degrees": 0,
            "max_angular_speed_degrees_per_second": np.pi*sweep/(sum(durations)/1000),
            "changes_data": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", choices=["all", "extent", "expanded", "window", "trajectory"], default="all")
    parser.add_argument("--frames", type=int, default=240)
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--sweep", type=float, default=60, help="Total horizontal arc in degrees; camera tilt stays fixed")
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "illustrations")
    parser.add_argument("--output", type=Path, default=ROOT / "docs" / "assets")
    args = parser.parse_args()
    if args.frames < 12 or not 1 <= args.fps <= 50 or not 0 < args.sweep <= 360:
        parser.error("Use at least 12 frames, 1–50 fps, and a sweep between 0 and 360 degrees")
    args.output.mkdir(parents=True, exist_ok=True)
    kinds = ["extent", "expanded", "window", "trajectory"] if args.kind == "all" else [args.kind]
    records = [render(k, args.data, args.output, args.frames, args.fps, args.sweep) for k in kinds]
    manifest = args.output / "animations.json"
    previous = json.loads(manifest.read_text(encoding="utf-8")) if manifest.exists() else []
    current = {r["file"]: r for r in previous + records}
    manifest.write_text(json.dumps(list(current.values()), indent=2)+"\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
