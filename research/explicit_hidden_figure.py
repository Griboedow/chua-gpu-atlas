"""Draw the voltage projections of all three nominal cycles at two scales."""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from scipy.integrate import solve_ivp

from explicit_hidden_cycles import OUT, PARAMETERS


def main():
    data = json.loads((OUT / "certificate-70.json").read_text())
    a, b, m0, m1, g = map(float, PARAMETERS)

    def rhs(t, state):
        x, y, z = state
        return [
            a * (y - x - m1 * x - (m0 - m1) * np.clip(x, -1, 1)),
            x - y + z,
            -b * y - g * z,
        ]

    paths = {}
    for name, record in data["cycles"].items():
        q = record["shooting_center"]
        start = [1, *[float(v["midpoint"]) for v in q[:2]]]
        T = sum(float(v["midpoint"]) for v in q[2::3])
        sol = solve_ivp(
            rhs,
            [0, T],
            start,
            t_eval=np.linspace(0, T, 1801),
            method="DOP853",
            max_step=0.01,
            rtol=1e-11,
            atol=1e-13,
        )
        assert sol.success
        paths[name] = sol.y[:2]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.25), layout="constrained")
    colors = ["#116A70", "#AC4E11", "#4D5184"]
    xeq = b * (m1 - m0) / (b * (1 + m1) + g * m1) + g * (m1 - m0) / (
        b * (1 + m1) + g * m1
    )
    yeq = g * xeq / (g + b)
    for i, ax in enumerate(axes):
        ax.plot(*paths["outer"], c=colors[2], lw=1.7, label="Outer cycle")
        ax.plot(*paths["right_lobe"], c=colors[0], lw=2, label="Right lobe")
        ax.plot(*(-paths["right_lobe"]), c=colors[1], lw=2, label="Left lobe")
        ax.plot([xeq, -xeq], [yeq, -yeq], ".", c="black", ms=5)
        ax.plot([0], [0], "x", c="black", ms=4)
        ax.set_xlabel("First capacitor voltage, x")
        ax.set_ylabel("Second capacitor voltage, y")
        ax.grid(alpha=0.16)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].add_patch(
        Rectangle(
            (-3, -0.75), 6, 1.5, facecolor="none", edgecolor=".45", lw=1, linestyle="--"
        )
    )
    axes[0].set_title("All three periodic orbits", loc="left", fontweight="bold")
    axes[0].legend(loc="upper left", frameon=False, fontsize=8)
    axes[1].set(xlim=(-3, 3), ylim=(-0.75, 0.75))
    axes[1].set_title("Two reflected lobes, enlarged", loc="left", fontweight="bold")
    fig.savefig(OUT / "cycles.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
