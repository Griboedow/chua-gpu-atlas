"""Draw the checked implications and the explicit boundary of the Chua proof pilot."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]


def main():
    evidence = json.loads((ROOT/"results/chua-proof-audit/numerical-audit.json").read_text())
    formal = json.loads((ROOT/"results/chua-proof-audit/formal-check.json").read_text())
    assert evidence["status"] == "local_numerical_hypotheses_reproduced"
    assert formal["status"] == "passed"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                         "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(14, 8.3))
    fig.patch.set_facecolor("#f5f3ee")
    ax.set(xlim=(0, 14), ylim=(0, 8.3))
    ax.axis("off")
    ink, blue, green, amber = "#20313a", "#266a94", "#20735b", "#9b5b17"
    ax.text(.35, 7.87, "One argument, three kinds of evidence", fontsize=23, color=ink, weight="bold")
    ax.text(.35, 7.45, "Chua tangency certificate · coexistence audit · first formal proof", color="#52616a", fontsize=13)

    def box(x, y, w, h, title, body, color, fill):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.025,rounding_size=.12",
                                   linewidth=1.5, edgecolor=color, facecolor=fill))
        ax.text(x+.2, y+h-.37, title, color=color, fontsize=14, weight="bold", va="top")
        ax.text(x+.2, y+h-.87, body, color=ink, fontsize=12, va="top", linespacing=1.5)

    box(.35, 4.6, 4.15, 2.25, "1  NUMERICAL CERTIFICATE",
        "Actual Chua return map\n80 and 100 digits: all 15 gates pass\nTwo shifted boxes: both rejected\n135 transverse flights per connection", blue, "#e8f1f7")
    box(4.9, 4.6, 4.15, 2.25, "2  LEAN-CHECKED IMPLICATION",
        "Arbitrary continuous gap g(a, b)\nFour face signs + actual derivatives\n⇒ an interior g = 0 and gₐ = 0\nAlso: ten exact rational guards", green, "#e6f2eb")
    box(9.45, 4.6, 4.15, 2.25, "3  CLASSICAL THEOREM",
        "Newhouse–Robinson\nApplies within a fixed parameter family\nCan give infinitely many sinks\nat a single parameter value", blue, "#e8f1f7")
    ax.add_patch(FancyArrowPatch((4.52, 5.7), (4.87, 5.7), arrowstyle="-|>", mutation_scale=16,
                                color=amber, linestyle="--"))
    ax.add_patch(FancyArrowPatch((9.07, 5.7), (9.42, 5.7), arrowstyle="-|>", mutation_scale=16,
                                color=amber, linestyle="--"))
    box(.35, 1.7, 6.35, 2.12, "FORMALIZATION STILL TO CONNECT",
        "The interval evaluator and exact invariant manifolds\nare not encoded in Lean. The general theorem checks\nthe logical implication once its hypotheses hold.\nNo theorem about Chua is assumed inside the proof.", amber, "#fcf0dc")
    box(7.05, 1.7, 6.55, 2.12, "COEXISTENCE STEP STILL TO CLOSE",
        "Show that the theorem-generated sinks remain inside\nthe valid Chua return domain. A compact extension alone\nis insufficient without that confinement argument.\nInfinite coexistence remains conditional here.", amber, "#fcf0dc")
    ax.text(.35, .99, "New formal observation", color=green, fontsize=13, weight="bold")
    ax.text(3.1, .99, "Convexity is not needed for existence of the stationary zero.", color=ink, fontsize=13)
    ax.text(.35, .43, "This is an argument map, not a plot of infinitely many computed attractors.", color="#52616a", fontsize=11)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    fig.savefig(ROOT/"docs/assets/chua-proof-audit.png", dpi=160, facecolor=fig.get_facecolor())
    fig.savefig(ROOT/"results/chua-proof-audit/chua-proof-audit.pdf", facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    main()
