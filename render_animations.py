"""Render measured data while rotating the camera about three coordinate axes.

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


def camera(axis, angle):
    """Rotate the viewing direction and up vector about a fixed world axis."""
    elev, azim = np.deg2rad([23, -55])
    direction = np.array([np.cos(elev)*np.cos(azim), np.cos(elev)*np.sin(azim), np.sin(elev)])
    up = np.array([-np.sin(elev)*np.cos(azim), -np.sin(elev)*np.sin(azim), np.cos(elev)])
    a = np.eye(3)[axis]
    theta = np.deg2rad(angle)
    def rotate(v):
        return v*np.cos(theta)+np.cross(a, v)*np.sin(theta)+a*np.dot(a, v)*(1-np.cos(theta))
    direction, up = rotate(direction), rotate(up)
    el = np.arcsin(np.clip(direction[2], -1, 1))
    az = np.arctan2(direction[1], direction[0])
    default_up = np.array([-np.sin(el)*np.cos(az), -np.sin(el)*np.sin(az), np.cos(el)])
    roll = np.arctan2(np.dot(direction, np.cross(default_up, up)), np.dot(default_up, up))
    # Matplotlib rotates its screen basis by minus the supplied roll angle.
    return np.rad2deg([el, az, -roll])


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


def render(kind, data_dir, output, frames_per_axis, fps):
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": "#e1edf8", "font.size": 11})
    fig = plt.figure(figsize=(9.6, 6.4), dpi=100, facecolor=BG)
    ax = fig.add_axes([.08, .15, .78, .7], projection="3d")
    labels = ("x", "y", "z") if kind == "trajectory" else ("α", "β", "m₀")
    configure(ax, labels)
    if kind == "extent":
        p = np.load(data_dir / "extent_candidates.npz")["parameters"]
        dots = ax.scatter(*p[:, :3].T, c=p[:, 2], cmap="viridis", s=7, alpha=.85,
                          linewidths=0, depthshade=False)
        ax.set(xlim=(0, 256), ylim=(0, 12500), zlim=(-8, -1))
        cb = fig.colorbar(dots, cax=fig.add_axes([.9, .29, .015, .35]))
        cb.set_label("Inner slope m₀", color="#dce7f3")
        cb.ax.tick_params(colors="#b7c9dc")
        title = "How far does the chaotic region extend?"
        subtitle = f"{len(p):,} repeated GPU candidates · fixed m₁ = −5/7 · linear coordinates"
        footer = "Sparse finite-time evidence; neither a closed solid nor a proof of unboundedness."
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
    caption = fig.text(.06, .09, "", fontsize=12, color="#f6be57")
    fig.text(.06, .045, footer, fontsize=9, color="#b7c9dc")
    frames = []
    for axis in [2, 0, 1]:
        for i in range(frames_per_axis):
            ax.view_init(*camera(axis, 360*i/frames_per_axis))
            caption.set_text(f"Camera rotation about {labels[axis]} axis  ·  {360*i/frames_per_axis:03.0f}°")
            fig.canvas.draw()
            rgb = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
            frame = Image.fromarray(rgb).quantize(colors=128, method=Image.Quantize.MEDIANCUT,
                                                  dither=Image.Dither.NONE)
            frames.append(frame)
            if i == 0 and axis == 2:
                Image.fromarray(rgb).save(output / f"{kind}-poster.png")
        print(f"{kind}: completed rotation about axis {labels[axis]}", flush=True)
    # GIF durations are multiples of 10 ms; distribute rounding error over frames.
    durations = [10*(round((i+1)*100/fps)-round(i*100/fps)) for i in range(len(frames))]
    target = output / f"{kind}-rotation.gif"
    frames[0].save(target, save_all=True, append_images=frames[1:], loop=0,
                   duration=durations, disposal=2, optimize=False)
    plt.close(fig)
    with Image.open(target) as check:
        if check.n_frames != len(frames):
            raise RuntimeError(f"Unexpected GIF frame count: {check.n_frames}")
    return {"file": target.name, "frames": len(frames), "fps": fps,
            "seconds": sum(durations)/1000, "bytes": target.stat().st_size,
            "camera_axes": list(labels), "changes_data": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", choices=["all", "extent", "window", "trajectory"], default="all")
    parser.add_argument("--frames-per-axis", type=int, default=48)
    parser.add_argument("--fps", type=int, default=16)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "illustrations")
    parser.add_argument("--output", type=Path, default=ROOT / "docs" / "assets")
    args = parser.parse_args()
    if args.frames_per_axis < 4 or not 1 <= args.fps <= 50:
        parser.error("Use at least 4 frames per axis and 1–50 fps")
    args.output.mkdir(parents=True, exist_ok=True)
    kinds = ["extent", "window", "trajectory"] if args.kind == "all" else [args.kind]
    records = [render(k, args.data, args.output, args.frames_per_axis, args.fps) for k in kinds]
    manifest = args.output / "animations.json"
    previous = json.loads(manifest.read_text(encoding="utf-8")) if manifest.exists() else []
    current = {r["file"]: r for r in previous + records}
    manifest.write_text(json.dumps(list(current.values()), indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
