"""Figures from the numerical window data and the interval certificate."""

import json
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator, LogFormatterMathtext

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/chua-layer-infinity"
ASSETS = ROOT / "docs/assets"
INK, MUTED = "#193547", "#657985"
BLUE, ORANGE, TEAL = "#2879b5", "#cc6337", "#148476"
PALE, GRID = "#f0f5f7", "#dce5e9"


def read(name):
    return json.loads((OUT / name).read_text())


def interval(record):
    # Plotting conversion only. Decisions in the certificate use Arb, not floats.
    mid = float(record["midpoint"])
    radius = float(record["radius"].strip("[]").split("+/-")[0])
    return mid, radius


def save(fig, name):
    fig.savefig(
        ASSETS / (name + ".png"), dpi=180, facecolor="white", bbox_inches="tight"
    )
    fig.savefig(
        OUT / (name + ".pdf"),
        facecolor="white",
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


def overview(cert):
    data = json.loads(
        (ROOT / "results/chua-layer-prediction/analysis.json").read_text()
    )
    rows = data["rows"]
    n = np.array([r["n"] for r in rows])
    width = np.array([r["width"] for r in rows])
    beta = float(cert["parameter_interval"]["midpoint"])
    distance = np.array([beta - r["fold_beta"] for r in rows])
    fig = plt.figure(figsize=(13.4, 7.9))
    grid = fig.add_gridspec(
        2,
        2,
        left=0.08,
        right=0.97,
        bottom=0.13,
        top=0.80,
        width_ratios=[1.12, 1],
        height_ratios=[1, 1],
        hspace=0.42,
        wspace=0.34,
    )
    fig.text(
        0.055,
        0.955,
        "An infinite cascade in the Chua parameter atlas",
        fontsize=22,
        weight="bold",
    )
    fig.text(
        0.055,
        0.909,
        r"$\alpha=9,\quad m_0=-8/7,\quad m_1=-5/7$"
        "     |     Measured windows and a certified organizing tangency",
        fontsize=12,
        color=MUTED,
    )
    ax = fig.add_subplot(grid[:, 0])
    for vals, color, label, marker in [
        (distance, BLUE, r"Distance to $\beta_*$", "o"),
        (width, ORANGE, "Symmetric-cycle stability width", "s"),
    ]:
        ax.semilogy(n, vals, color=color, lw=1.8, marker=marker, ms=6, label=label)
        ax.scatter(
            n[-2:],
            vals[-2:],
            s=110,
            facecolors="white",
            edgecolors=color,
            linewidths=2,
            zorder=5,
        )
    ax.set(
        xlabel="Peaks per complete periodic motion",
        ylabel=r"Parameter distance / width in $\beta$",
        xticks=n,
        xticklabels=[str(r["peaks"]) for r in rows],
        ylim=(1.3e-10, 8e-3),
    )
    ax.set_title("A   Eight measured members", loc="left", fontsize=14, pad=17)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=9))
    ax.yaxis.set_major_formatter(LogFormatterMathtext())
    ax.grid(axis="y", which="major")
    ax.legend(loc="upper right", frameon=False, fontsize=10)
    ax.annotate(
        "75 / 85 peaks\nfound after forecasts",
        xy=(7.2, 5e-9),
        xytext=(3.1, 2.5e-8),
        color=ORANGE,
        fontsize=11,
        arrowprops=dict(arrowstyle="->", color=ORANGE),
        bbox=dict(boxstyle="round,pad=.35", fc="white", ec="none"),
    )
    ax.text(
        0.02,
        0.02,
        "85-peak stability width: 6.02 × 10⁻¹⁰\nThe eight points alone do not prove infinity.",
        transform=ax.transAxes,
        fontsize=10,
        color=MUTED,
        va="bottom",
    )
    ax = fig.add_subplot(grid[0, 1])
    ax.set_title("B   The interval certificate", loc="left", fontsize=14, pad=17)
    ax.axis("off")
    lines = [
        (
            "Exact invariant manifolds",
            "Graph bounds include curvature and parameter motion.",
        ),
        ("Quadratic contact", r"$g_{aa}\in[2.0956,\,2.1447]\times10^{13}$"),
        ("Transverse unfolding", r"$g_\beta\in[3.3728,\,3.3731]$"),
        ("Dissipative saddle", r"$|\det DP|<2.788\times10^{-9}$"),
    ]
    for i, (label, value) in enumerate(lines):
        y = 0.96 - i * 0.245
        ax.text(0, y, label, fontsize=12, weight="bold", va="top", color=TEAL)
        ax.text(0, y - 0.095, value, fontsize=10.5, va="top", color=INK)
    ax = fig.add_subplot(grid[1, 1])
    ax.axis("off")
    ax.add_patch(
        Rectangle(
            (-0.03, -0.04),
            1.06,
            1.05,
            transform=ax.transAxes,
            facecolor=PALE,
            edgecolor="none",
            clip_on=False,
        )
    )
    ax.text(
        0.02,
        0.91,
        "C   Consequence of the sink-cascade theorem",
        weight="bold",
        fontsize=12,
        va="top",
    )
    ax.text(
        0.02,
        0.68,
        "Infinitely many attracting windows\nwith increasing return periods",
        fontsize=17,
        weight="bold",
        va="top",
    )
    ax.text(
        0.02,
        0.32,
        "They accumulate at a certified contact near",
        fontsize=11,
        color=MUTED,
    )
    ax.text(0.02, 0.16, r"$\beta_*\approx14.30635909996146$", fontsize=17, color=TEAL)
    fig.text(
        0.08,
        0.038,
        "Numerical data: 15–85 peak sequence.  Certificate: outward-rounded interval inequalities.\n"
        "Infinity concerns windows at different parameters; the total number of connected 3D surfaces is not determined.",
        fontsize=10.5,
        color=MUTED,
    )
    save(fig, "chua-infinity-overview")


def signs(cert):
    fig = plt.figure(figsize=(13.4, 7.8))
    grid = fig.add_gridspec(
        2,
        2,
        left=0.075,
        right=0.96,
        bottom=0.17,
        top=0.79,
        width_ratios=[1.1, 1],
        hspace=0.8,
        wspace=0.4,
    )
    fig.text(
        0.055, 0.95, "A tangency forced by interval signs", fontsize=22, weight="bold"
    )
    fig.text(
        0.055,
        0.90,
        "The bounds include every allowed exact stable and unstable graph.",
        fontsize=12,
        color=MUTED,
    )
    ax = fig.add_subplot(grid[:, 0])
    ax.add_patch(Rectangle((-1, -1), 2, 2, facecolor=PALE, edgecolor="none"))
    ax.plot([-1, 1], [-1, -1], color=BLUE, lw=6, solid_capstyle="butt")
    ax.plot([-1, 1], [1, 1], color=ORANGE, lw=6, solid_capstyle="butt")
    ax.plot([-1, -1], [-1, 1], color=BLUE, lw=4, solid_capstyle="butt")
    ax.plot([1, 1], [-1, 1], color=ORANGE, lw=4, solid_capstyle="butt")
    ax.text(0, 1.16, r"$g>0$ for every $a$", ha="center", color=ORANGE, fontsize=12)
    ax.text(0, -1.22, r"$g<0$ for every $a$", ha="center", color=BLUE, fontsize=12)
    ax.text(
        -1.12,
        0,
        r"$g_a<0$",
        ha="right",
        va="center",
        rotation=90,
        color=BLUE,
        fontsize=12,
    )
    ax.text(
        1.12,
        0,
        r"$g_a>0$",
        ha="left",
        va="center",
        rotation=90,
        color=ORANGE,
        fontsize=12,
    )
    ax.text(0, 0.49, "Throughout this box", ha="center", fontsize=12, weight="bold")
    ax.text(
        0,
        0.1,
        r"$g_{aa}>0\quad\mathrm{and}\quad g_\beta>0$",
        ha="center",
        fontsize=16,
        color=TEAL,
    )
    ax.text(
        0,
        -0.48,
        "A unique convex minimum\ncrosses zero inside.",
        ha="center",
        fontsize=12,
        linespacing=1.55,
    )
    ax.set(
        xlim=(-1.5, 1.5),
        ylim=(-1.5, 1.5),
        xticks=[-1, 0, 1],
        yticks=[-1, 0, 1],
        xlabel=r"$(a-a_0)/10^{-14}$",
        ylabel=r"$(\beta-\beta_0)/(5\times10^{-14})$",
    )
    ax.set_title(
        "A   Existence argument on the tested box", loc="left", fontsize=13, pad=19
    )
    for i, (names, scale, ylabel, title) in enumerate(
        [
            (
                ["beta_lower", "beta_upper"],
                1e13,
                r"Gap $g$ in units of $10^{-13}$",
                "B   Gap on the two parameter faces",
            ),
            (
                ["a_lower", "a_upper"],
                1,
                r"Derivative $g_a$",
                "C   Derivative on the two state faces",
            ),
        ]
    ):
        ax = fig.add_subplot(grid[i, 1])
        for k, (name, color) in enumerate(zip(names, [BLUE, ORANGE])):
            m, e = interval(cert["face_signs"][name])
            ax.errorbar(
                m * scale,
                k,
                xerr=e * scale,
                fmt="o",
                color=color,
                capsize=6,
                lw=4,
                ms=5,
            )
        ax.axvline(0, color=INK, lw=1, ls="--")
        ax.set(
            yticks=[0, 1],
            yticklabels=["Lower face", "Upper face"],
            ylim=(-0.65, 1.65),
            xlabel=ylabel,
        )
        limit = 2.15 if i == 0 else 0.29
        ax.set_xlim(-limit, limit)
        ax.set_title(title, loc="left", fontsize=13, pad=12)
        ax.grid(axis="x")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    fig.text(
        0.075,
        0.056,
        r"$a_0\approx-7.4583788929913551\times10^{-8}\qquad"
        r"\beta_0\approx14.3063590999614598$",
        fontsize=13,
    )
    fig.text(
        0.075,
        0.018,
        "Bars enclose the entire corresponding face, not statistical errors. "
        "The box diagram does not assign an uncomputed point location.",
        fontsize=10.5,
        color=MUTED,
    )
    save(fig, "chua-infinity-signs")


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": "#b9c7ce",
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "grid.color": GRID,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    cert = read("tangency-certificate-100.json")
    if cert["status"] != "passed":
        raise ValueError("Refusing to plot a successful certificate from failed data")
    overview(cert)
    signs(cert)
    print("Rendered two certificate figures as PNG and PDF.")


if __name__ == "__main__":
    main()
