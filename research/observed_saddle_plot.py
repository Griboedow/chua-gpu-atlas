"""Render the voltage launch, phase-space return, and out-of-sample comparison."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from observed_saddle_residence import OUT
from observed_saddle_validate import direct


def main(output_dir=OUT):
    forecast = json.loads((OUT/"forecast.json").read_text())
    data = json.loads((output_dir/"measurements.json").read_text())
    records = data["records"]
    blue, orange, grey = "#006895", "#a74714", "#596673"
    plt.rcParams.update({"font.size":10, "axes.spines.right":False, "axes.spines.top":False,
                         "axes.labelcolor":"#263744", "axes.titleweight":"bold"})
    fig, axes = plt.subplots(2, 2, figsize=(11.3, 7.35), layout="constrained")
    (ax1, ax2), (ax3, ax4) = axes
    point = np.r_[-1., forecast["point"]]
    p = direct(point, forecast["period"], dense=True)
    curve = p.sol(np.linspace(0,forecast["period"],901))
    x0 = forecast["crossing_x"]
    entry = direct([x0,0.,0.], 6., dense=True)
    path = entry.sol(np.linspace(0,entry.t_events[0][0],601))
    ax1.plot(path[0],path[1],color=orange,label="Launch to first section hit",lw=1.9)
    ax1.plot(curve[0],curve[1],color=blue,label="Saddle cycle",lw=1.9)
    ax1.scatter(x0,0.,color=orange,s=30,zorder=3)
    ax1.scatter(point[0],point[1],facecolor="white",edgecolor=blue,s=40,zorder=4)
    ax1.set(xlabel="First capacitor voltage, x",ylabel="Second capacitor voltage, y",
            title="A   Approach from a voltage-only start")
    ax1.legend(frameon=False,loc="lower right",fontsize=9)
    pairs = [(r, color) for r, color in ((next(r for r in records if r["offset_x"]==-.0001),orange),
                                         (next(r for r in records if r["offset_x"]==.0001),blue))]
    for r, color in pairs:
        u = np.abs(r["section_u"])
        ax2.plot(np.arange(len(u)),u,marker="s" if r["offset_x"]<0 else "o",
                 color=color,ms=5,mfc="white" if r["offset_x"]<0 else color,
                 ls="--" if r["offset_x"]<0 else "-",alpha=.85,
                 label=f"δx = {r['offset_x']:+.4f}")
    ax2.axhline(forecast["projected_exit_radius"],color=grey,ls=":",label="Exit at |u| = 0.02")
    ax2.set(yscale="log",xlabel="Full return index, first hit = 0",ylabel="Absolute unstable coordinate |u|",
            title="B   A rounding boundary at the eighth return")
    ax2.legend(frameon=False,loc="lower right",fontsize=9)
    volts = forecast["breakpoint_V"]*1000
    all_offsets = np.array(sorted(set(abs(r["offset_x"]) for r in records)))
    pred_x = np.geomspace(all_offsets.min(),all_offsets.max(),200)
    predicted = np.log(forecast["projected_exit_radius"]/(pred_x*forecast["entry_gain_per_x"])) / np.log(forecast["multipliers"][-1])
    ax3.plot(pred_x*volts,predicted,color=grey,lw=1.7,label="Linearized fractional forecast")
    for side, color, marker in [(-1,orange,"s"),(1,blue,"o")]:
        rows = sorted((r for r in records if np.sign(r["offset_x"])==side),key=lambda r:abs(r["offset_x"]))
        x = [abs(r["offset_x"])*volts for r in rows]
        style = dict(marker=marker, s=40 if side<0 else 20,facecolors="none" if side<0 else color,
                     edgecolors=color, label="Below threshold" if side<0 else "Above threshold")
        ax3.scatter(x,[r["first_exit"] for r in rows],**style)
        ax4.scatter(x,[r["fractional_error"] for r in rows],**style)
    for axis in (ax3,ax4):
        axis.set_xscale("log")
        axis.invert_xaxis()
        axis.set_xlabel("Distance from threshold |δU₁| (mV)")
        axis.grid(axis="y",alpha=.17)
    ax3.set(ylabel="Exit return",title="C   One decade closer adds 2.90 returns")
    ax4.axhline(0,color=grey,lw=1.)
    ax4.set(ylabel="Observed minus predicted (returns)",ylim=(-.021,.021),
            title="D   Fractional exit errors on 24 held-out starts")
    ax3.legend(frameon=False,fontsize=9)
    ax4.legend(frameon=False,fontsize=9,loc="center right")
    fig.savefig(output_dir/"voltage-start-residence.png",dpi=165)
    plt.close(fig)


if __name__ == "__main__":
    main()
