"""Frozen response thresholds and observed switching on the archived large-alpha trajectory."""

import argparse
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/chua-scaling"


class FrozenFiber:
    def __init__(self, kappa, m0, m1):
        self.omega = mp.sqrt(kappa)
        self.a, self.b = 1 + m1, 1 + m0
        self.xe = (self.a - self.b) / self.a
        self.period = 2 * mp.pi / self.omega

    def forcing_weights(self, slope, start, end):
        # Integrate exp(-slope*(end-t)) times cos(omega*t) and sin(omega*t).
        w = self.omega
        z = (mp.exp(1j * w * end) - mp.exp(-slope * (end - start) + 1j * w * start)) / (
            slope + 1j * w
        )
        return [z.real, z.imag]

    def forcing(self, d, horizon, initial):
        weights = mp.matrix(
            [
                self.forcing_weights(self.b, 0, d),
                self.forcing_weights(self.a, d, horizon),
            ]
        )
        rhs = mp.matrix(
            [
                1 - initial * mp.exp(-self.b * d),
                (1 - self.xe) * (1 - mp.exp(-self.a * (horizon - d))),
            ]
        )
        return mp.lu_solve(weights, rhs)

    def amplitude(self, d, horizon, initial):
        c, s = self.forcing(d, horizon, initial)
        return mp.sqrt(c * c + s * s)

    def threshold(self, family):
        symmetric = family == "central"
        horizon = self.period / 2 if symmetric else self.period
        initial = -1 if symmetric else 1
        d = self.a * horizon / (self.a - self.b)
        c, s = self.forcing(d, horizon, initial)

        def amplitude(duration):
            return self.amplitude(duration, horizon, initial)

        record = {
            "amplitude": float(amplitude(d)),
            "inner_duration_per_segment": float(d),
            "forcing_cos": float(c),
            "forcing_sin": float(s),
            "dR_dd": float(mp.diff(amplitude, d)),
            "d2R_dd2": float(mp.diff(amplitude, d, 2)),
        }
        if not symmetric:
            # Fixed-phase amplitude derivative of the full return. The field
            # is continuous, so the two switching-time terms cancel.
            unit = [c / amplitude(d), s / amplitude(d)]
            inner = mp.fdot(self.forcing_weights(self.b, 0, d), unit)
            outer = mp.fdot(self.forcing_weights(self.a, d, horizon), unit)
            record["return_derivative_wrt_R"] = float(
                mp.exp(-self.a * (horizon - d)) * inner + outer
            )
        record["full_field_check"] = self.check_orbit(float(c), float(s), initial)
        return record

    def check_orbit(self, c, s, initial):
        a, b, w, T = map(float, (self.a, self.b, self.omega, self.period))

        def field(t, state):
            x = state[0]
            return [
                c * np.cos(w * t)
                + s * np.sin(w * t)
                - a * x
                - (b - a) * np.clip(x, -1, 1)
            ]

        def lower(t, state):
            return state[0] + 1

        def upper(t, state):
            return state[0] - 1

        sol = solve_ivp(
            field,
            (0, T),
            [initial],
            method="DOP853",
            dense_output=True,
            events=(lower, upper),
            rtol=2e-12,
            atol=2e-13,
            max_step=T / 160,
        )
        crossings = sorted(
            (float(t), wall)
            for times, wall in zip(sol.t_events, [-1, 1])
            for t in times
            if 1e-7 < t < T - 1e-7
        )
        walls = [wall for _, wall in crossings]
        if walls != ([1, 1, -1] if initial == -1 else [1]):
            raise ValueError(f"Wrong itinerary from {initial}: {walls}")
        times = [0] + [t for t, _ in crossings] + [T]
        inside = sum(
            right - left
            for left, right in itertools.pairwise(times)
            if abs(sol.sol((left + right) / 2)[0]) < 1
        )
        closure = abs(float(sol.y[0, -1]) - initial)
        log_q = -a * T + (a - b) * inside
        speeds = [abs(field(t, sol.sol(t))[0]) for t in times]
        if not sol.success or max(closure, abs(log_q)) > 1e-7 or min(speeds) < 0.1:
            raise ValueError("Independent frozen-field control failed")
        return {
            "closure": closure,
            "log_multiplier": log_q,
            "inner_time": inside,
            "interior_crossings": crossings,
            "minimum_crossing_speed": min(speeds),
        }


def itinerary(samples, tolerance=0.0):
    low, high = np.min(samples), np.max(samples)
    if low < -1 - tolerance and high > 1 + tolerance:
        return "both"
    if low < -1 - tolerance and high < 1 - tolerance:
        return "negative"
    if high > 1 + tolerance and low > -1 + tolerance:
        return "positive"
    return "ambiguous"


def observed_itineraries(trace, stride=1, tolerance=0.0):
    times = trace["time"][::stride]
    x = trace["states"][::stride, 0]
    section = trace["section_states"][:, 0]
    cuts = np.searchsorted(times, trace["section_time"])
    return np.asarray(
        [
            itinerary(np.r_[section[i], x[lo:hi], section[i + 1]], tolerance)
            for i, (lo, hi) in enumerate(itertools.pairwise(cuts))
        ]
    )


def replay_switches(trace, switches, labels):
    alpha, kappa, m0, m1 = trace["parameters"]
    section = trace["section_states"]
    times = trace["section_time"]
    max_error = 0.0

    def field(t, state):
        x, y, z = state
        return [
            y - (1 + m1) * x - (m0 - m1) * np.clip(x, -1, 1),
            z + (x - y) / alpha,
            -kappa * y,
        ]

    for i in switches:
        start, middle, end = (times[i - 1 : i + 2] - times[i - 1]) * alpha
        sol = solve_ivp(
            field,
            (start, end),
            section[i - 1],
            method="Radau",
            rtol=2e-10,
            atol=2e-12,
            max_step=0.25,
            dense_output=True,
        )
        for index, left, right in [(i - 1, start, middle), (i, middle, end)]:
            x = sol.sol(np.linspace(left, right, 301))[0]
            if not sol.success or itinerary(x, 0.05) != labels[index]:
                raise ValueError(f"Direct replay disagrees at turn {index}")
        error = np.max(np.abs(sol.y[:, -1] - section[i + 1]))
        max_error = max(max_error, float(error))
    if max_error > 1e-5:
        raise ValueError(f"Two-turn replay endpoint mismatch: {max_error}")
    return max_error


def survey(trace):
    kappa = trace["parameters"][1]
    times = trace["section_time"]
    section = trace["section_states"]
    if np.max(np.abs(section[:, 1])) > 1e-9 or np.any(section[:, 2] <= 0):
        raise ValueError("Expected rising y=0 crossings")
    amplitudes = section[:, 2] / np.sqrt(kappa)
    labels = observed_itineraries(trace)
    for stride in (1, 2, 4):
        if not np.array_equal(labels, observed_itineraries(trace, stride, 0.05)):
            raise ValueError(f"Breakpoint labels changed at stride {stride}")
    switches = np.flatnonzero(labels[:-1] != labels[1:]) + 1
    events = [
        {
            "turn": int(i),
            "time": float(times[i]),
            "previous": str(labels[i - 1]),
            "current": str(labels[i]),
            "amplitude": float(amplitudes[i]),
        }
        for i in switches
    ]
    directions = {}
    for pair in sorted({(e["previous"], e["current"]) for e in events}):
        values = [
            e["amplitude"] for e in events if (e["previous"], e["current"]) == pair
        ]
        directions["_to_".join(pair)] = {
            "count": len(values),
            "minimum": min(values),
            "maximum": max(values),
        }
    return (
        {
            "complete_turns": len(labels),
            "counts": {
                str(k): int(n) for k, n in zip(*np.unique(labels, return_counts=True))
            },
            "transitions": directions,
            "events": events,
        },
        switches,
        labels,
        amplitudes,
    )


def plot(output, trace, amplitudes, record):
    fig, (series, events) = plt.subplots(
        1,
        2,
        figsize=(9.8, 3.5),
        sharey=True,
        layout="constrained",
        gridspec_kw={"width_ratios": [1.35, 1]},
    )
    times = trace["section_time"]
    visible = times <= 4
    series.plot(times[visible], amplitudes[visible], color="#394857", lw=1.3)
    series.set(
        xlabel="Original time",
        ylabel="Oscillator amplitude R",
        title="First two excursions",
    )
    for label, color, mark, title in [
        ("both", "#00746e", "o", "First two-well turn"),
        ("negative", "#ad5533", "s", "First negative-well turn"),
    ]:
        group = [e for e in record["finite_alpha"]["events"] if e["current"] == label]
        selected = [e for e in group if e["time"] <= 4]
        series.scatter(
            [e["time"] for e in selected],
            [e["amplitude"] for e in selected],
            color=color,
            marker=mark,
            s=26,
            zorder=3,
        )
        events.scatter(
            np.arange(1, 1 + len(group)),
            [e["amplitude"] for e in group],
            color=color,
            marker=mark,
            s=24,
            label=title,
        )
    for ax in (series, events):
        for name, text, style in [
            ("well", "Well fold", "--"),
            ("central", "Central neutral", ":"),
        ]:
            level = record["frozen"][name]["amplitude"]
            ax.axhline(level, color="#7a8188", linestyle=style, lw=1, zorder=0)
            if ax is series:
                ax.annotate(
                    text,
                    (0.04, level),
                    (0, 5),
                    textcoords="offset points",
                    fontsize=8,
                    color="#545b62",
                )
        ax.spines[["right", "top"]].set_visible(False)
    events.set(
        xlabel="Excursion number", title="All 16 excursions", xticks=[1, 4, 8, 12, 16]
    )
    events.legend(loc="center", bbox_to_anchor=(0.5, 0.60), frameon=False, fontsize=8)
    fig.savefig(output, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "outputs/chua-fiber-switching"
    )
    out = parser.parse_args().output_dir
    source = "cpu-2198-aligned.npz"
    with np.load(SOURCE / source) as f:
        trace = dict(f)
    # The atlas defines the outer slope as the exact rational -5/7.
    mp.mp.dps = 60
    fiber = FrozenFiber(mp.mpf("0.1425"), mp.mpf(-2), -mp.mpf(5) / 7)
    if not np.allclose(
        trace["parameters"], [2048, 0.1425, -2, -5 / 7], atol=1e-15, rtol=0
    ):
        raise ValueError("Unexpected source parameters")
    frozen = {name: fiber.threshold(name) for name in ("central", "well")}
    survey_result, switches, labels, amplitudes = survey(trace)
    survey_result["max_two_turn_radau_endpoint_error"] = replay_switches(
        trace, switches, labels
    )
    old = json.loads((SOURCE / "fiber-probe.json").read_text())
    errors = []
    for row in old["shooting_controls"]:
        if row["name"].endswith("_half_wave_shoot"):
            h = fiber.period / 2
            d = (mp.log(row["event_based_multiplier"]) / 2 + fiber.a * h) / (
                fiber.a - fiber.b
            )
            errors.append(abs(float(fiber.amplitude(d, h, -1)) - row["R"]))
    if len(errors) != 5 or max(errors) > 1e-7:
        raise ValueError(
            "Archived central orbit multipliers disagree with the matched flights"
        )
    record = {
        "parameters": {"alpha": 2048, "kappa": 0.1425, "m0": -2, "m1": "-5/7"},
        "trajectory": "results/chua-scaling/" + source,
        "section": "y=0 increasing; amplitude at the start of the first turn with the new itinerary",
        "frozen": frozen,
        "max_amplitude_error_from_five_archived_multipliers": max(errors),
        "finite_alpha": survey_result,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "switching.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    plot(out / "switching.png", trace, amplitudes, record)
    print(survey_result["transitions"])
    print(
        "Maximum two-turn replay endpoint error:",
        survey_result["max_two_turn_radau_endpoint_error"],
    )


if __name__ == "__main__":
    main()
