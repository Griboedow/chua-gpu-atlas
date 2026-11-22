"""Compare the published saddle period with periods at measured component settings."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from flint import arb
from scipy.integrate import solve_ivp
from scipy.optimize import root

from observed_saddle_flights import NOMINAL, Section


ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / "results/chua-precharge/settings.csv"
CERTIFICATE = ROOT / "results/chua-saddle-period/interval.json"
SAMPLING_US = 2.5
MEASURED_SAMPLES = 218
MEASURED_US = SAMPLING_US * MEASURED_SAMPLES


def follow_saddle(parameters, start):
    """Follow the four-flight orbit from the nominal model to this row."""
    point = start.copy()
    for fraction in (0.5, 1.0):
        p = {k: v + fraction * (parameters[k] - v) for k, v in NOMINAL.items()}
        section = Section(p)

        def residual(yz):
            _, end, derivative = section.full_return(yz)
            return end - yz, derivative - np.eye(2)

        sol = root(residual, point, jac=True, tol=1e-12)
        point = sol.x
        arcs = []
        period, end, jac = section.full_return(point, arcs)
        closure = float(np.linalg.norm(end - point, ord=np.inf))
        multipliers = np.linalg.eigvals(jac)
        if not np.isreal(multipliers).all():
            raise ValueError("Saddle multipliers must be real")
        stable, unstable = sorted(multipliers.real)
        if (closure > 1e-10 or [a[0] for a in arcs] != [0, 1, 0, -1]
                or not (abs(stable) < 1 < unstable)):
            raise ValueError("Continuation left the four-flight saddle branch")
    return point, period, closure, stable, unstable, section.absolute_x_peak(arcs)


def independent_return(parameters, row):
    """Integrate the continuous field directly to its next section crossing."""
    a, b, g, m0, m1 = [parameters[k] for k in ("alpha", "beta", "gamma", "m0", "m1")]

    def field(t, state):
        x, y, z = state
        diode = m1*x + (m0-m1)*(abs(x+1)-abs(x-1))/2
        return a*(y-x-diode), x-y+z, -b*y-g*z

    def positive_section(t, state):
        return state[0] + 1

    positive_section.direction = 1
    start = [-1., row["section_y"], row["section_z"]]
    sol = solve_ivp(field, (0, 2.7), start, method="DOP853", rtol=2e-13,
                    atol=1e-14, max_step=0.0125, events=positive_section)
    times = sol.t_events[0]
    if not sol.success or len(times) != 2 or abs(times[0]) > 1e-10:
        raise RuntimeError("Expected the initial and one subsequent upward crossing")
    error = float(abs(times[-1] - row["period_dimensionless"]))
    closure = float(np.linalg.norm(sol.y_events[0][-1] - start, ord=np.inf))
    if error > 1e-8 or closure > 1e-7:
        raise ArithmeticError("Full-field integration disagrees with shooting")
    return dict(settings_index=row["settings_index"], ivp_period_dimensionless=times[-1],
                absolute_period_difference=error, section_closure_inf=closure)


def render(rows, output, certificate):
    coupled = [row for row in rows if row["design"] == "coupled"]
    periods = [row["period_us"] for row in coupled]
    clocks = [row["clock_us"] for row in coupled]
    with certificate.open(encoding="utf-8") as handle:
        lower, upper = json.load(handle)["coupled_period_microseconds"]
    # Keep the enclosure outward rounded when converting Arb balls to plot coordinates.
    bounds = [np.nextafter(float(arb(lower).lower()), -np.inf),
              np.nextafter(float(arb(upper).upper()), np.inf)]
    fig, (ax, zoom) = plt.subplots(1, 2, figsize=(8.8, 3.1), layout="constrained",
                                   gridspec_kw={"width_ratios": [1.25, 1]})
    ax.plot(bounds, [3, 3], color="#26705c", marker="|", markersize=12,
            markeredgewidth=1.8, linewidth=2.5)
    ax.scatter(periods, np.full(len(coupled), 2), marker="|", color="#27577a", s=110)
    ax.scatter(MEASURED_US, 0, marker="D", s=38, color="#a44630")
    nominal = next(row for row in rows if row["design"] == "components")
    ax.scatter(nominal["period_us"], 1, s=34, color="#172f46")
    ax.set(yticks=[0, 1, 2, 3], yticklabels=["Published peak spacing", "Nominal circuit",
                                             "64 component settings", "Certified range"],
           ylim=(-0.5, 3.5), xlabel="Time (µs)")
    ax.margins(x=0.08)
    ax.annotate("218 samples", (MEASURED_US, 0), xytext=(7, 2),
                textcoords="offset points", fontsize=9)
    zoom.scatter(clocks, periods, s=15, color="#27577a", alpha=0.78)
    zoom.set(xlabel="Clock factor R C₂ (µs)", ylabel="Period (µs)",
             title="Component settings, enlarged")
    for panel in (ax, zoom):
        panel.grid(axis="x", color="0.9")
        panel.set_axisbelow(True)
        panel.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "period-comparison.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/chua-saddle-period")
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE,
                        help="Interval certificate used for the figure (default: archived result)")
    args = parser.parse_args()
    with SETTINGS.open(encoding="utf-8", newline="") as handle:
        settings = [r for r in csv.DictReader(handle) if r["design"] != "table-box"]
    nominal_components = next(r for r in settings if r["design"] == "components")
    nominal_point = Section().saddle()[0]
    rows = []
    parameters = {}
    for source in settings:
        index = int(source["index"])
        p = {k: float(source[k]) for k in NOMINAL}
        point, period, closure, stable, unstable, peak = follow_saddle(p, nominal_point)
        if peak >= float(source["outer_break"]):
            raise ArithmeticError("Orbit reaches the second diode breakpoint")
        # The rounded parameter rows have no separately measured clock.
        clock_source = source if source["R"] else nominal_components
        clock = float(clock_source["R"]) * float(clock_source["C2"]) * 1e6
        row = dict(settings_index=index, design=source["design"],
                   clock_settings_index=int(clock_source["index"]), clock_us=clock,
                   period_dimensionless=period, period_us=period*clock,
                   period_in_sample_spacings=period*clock/SAMPLING_US,
                   section_y=point[0], section_z=point[1], section_closure_inf=closure,
                   lambda_stable=stable, lambda_unstable=unstable,
                   turns_per_decade=np.log(10)/np.log(unstable), max_abs_x=peak,
                   outer_break=float(source["outer_break"]))
        rows.append(row)
        parameters[index] = p
    coupled = [row for row in rows if row["design"] == "coupled"]
    if len(coupled) != 64 or len(rows) != 67:
        raise ValueError("Expected three nominal conventions and 64 coupled settings")
    lower = min(coupled, key=lambda r: r["period_us"])
    upper = max(coupled, key=lambda r: r["period_us"])
    nominal = next(row for row in rows if row["design"] == "components")
    checks = [independent_return(parameters[row["settings_index"]], row)
              for row in [nominal, lower, upper]]
    ranges = {key: [min(r[key] for r in coupled), max(r[key] for r in coupled)]
              for key in ["clock_us", "period_dimensionless", "period_us", "lambda_unstable",
                          "turns_per_decade", "max_abs_x", "outer_break"]}
    summary = dict(source="https://doi.org/10.1007/s11071-022-08078-y", figure="10b",
                   published_samples=MEASURED_SAMPLES, sample_interval_us=SAMPLING_US,
                   sample_derived_period_us=MEASURED_US, rounded_published_period_us=550,
                   coupled_settings=len(coupled), coupled_sample_ranges=ranges,
                   min_and_max_period_settings=[lower["settings_index"], upper["settings_index"]],
                   nearest_gap_us=lower["period_us"]-MEASURED_US,
                   nearest_gap_in_sample_spacings=(lower["period_us"]-MEASURED_US)/SAMPLING_US,
                   gap_from_rounded_period_us=lower["period_us"]-550,
                   independent_DOP853_returns=checks)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    with (out / "component-periods.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (out / "period-summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)
                                             + "\n", encoding="utf-8")
    render(rows, out, args.certificate)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
