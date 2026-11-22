"""Count periods on simulated U1 sampled at the Fig. 10 oscilloscope cadence."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

from observed_saddle_residence import BREAK_VOLTAGE, OUT as RESIDENCE
from observed_saddle_validate import direct


TAU_US = 788 * 311e-3  # R in ohms, C2 in nF, result in microseconds
SAMPLE_US = 2.5
STEP = SAMPLE_US / TAU_US
OFFSETS = (1e-3, 1e-4, 1e-5, 1e-6)
PHASES = np.arange(16) / 16
QUANTA_MV = (0, 1, 5)


def digitize(volts, quantum_mv):
    if quantum_mv:
        return np.rint(volts * 1000 / quantum_mv) * (quantum_mv / 1000)
    return volts


def crossing_indices(volts):
    """First sample at or above -B and linear threshold time, in sample units."""
    level = -BREAK_VOLTAGE
    lo = np.flatnonzero((volts[:-1] < level) & (volts[1:] >= level))
    interpolated = lo + (level - volts[lo]) / (volts[lo + 1] - volts[lo])
    return lo + 1, interpolated


def limits(values):
    return [float(np.min(values)), float(np.max(values))]


def evaluate(record, radius):
    # The archived affine events supply ground truth and the eligibility label.
    # Neither enters the voltage crossing detector.
    turns = min(8, record["first_exit"] - 1)
    times = np.array(record["section_time"][:turns + 1])
    if turns < 1 or np.max(np.abs(record["section_u"][:turns + 1])) >= radius:
        raise ValueError("Archived crossings lie outside the residence window")
    exact = np.diff(times) / STEP
    stop = times[-1] + 3 * STEP
    sol = direct([record["launch_x"], 0., 0.], stop, dense=True)
    if len(sol.t_events[0]) != len(times):
        raise ArithmeticError("Independent solver found a different crossing count")
    solver_error = float(np.max(np.abs(sol.t_events[0] - times)) / STEP)
    if solver_error > 0.001:
        raise ArithmeticError(f"Independent event times differ by {solver_error} samples")

    traces = []
    for phase in PHASES:
        grid = np.arange(phase * STEP, stop, STEP)
        traces.append(BREAK_VOLTAGE * sol.sol(grid)[0])

    rows = []
    for quantum in QUANTA_MV:
        raw, linear = [], []
        for volts in traces:
            index, interpolated = crossing_indices(digitize(volts, quantum))
            if len(index) != len(times):
                raise ArithmeticError("Sampled voltage yielded a different crossing count")
            raw.append(np.diff(index))
            linear.append(np.diff(interpolated))
        raw, linear = np.array(raw), np.array(linear)
        row = dict(offset_x=record["offset_x"], launch_x=record["launch_x"],
                   first_exit=record["first_exit"], turns=turns, quantum_mV=quantum,
                   single_exact_min=exact.min(), single_exact_max=exact.max(),
                   single_index_min=int(raw.min()), single_index_max=int(raw.max()),
                   single_linear_min=linear.min(), single_linear_max=linear.max(),
                   max_abs_index_error=np.max(np.abs(raw - exact)),
                   max_abs_linear_error=np.max(np.abs(linear - exact)),
                   eight_exact=None, eight_index_min=None, eight_index_max=None,
                   eight_linear_min=None, eight_linear_max=None,
                   eight_max_abs_index_error=None, eight_max_abs_linear_error=None)
        if turns == 8:
            exact_mean = np.mean(exact)
            index_mean, linear_mean = np.mean(raw, axis=1), np.mean(linear, axis=1)
            row.update(eight_exact=exact_mean, eight_index_min=index_mean.min(),
                       eight_index_max=index_mean.max(), eight_linear_min=linear_mean.min(),
                       eight_linear_max=linear_mean.max(),
                       eight_max_abs_index_error=np.max(np.abs(index_mean - exact_mean)),
                       eight_max_abs_linear_error=np.max(np.abs(linear_mean - exact_mean)))
        rows.append(row)
    return rows, solver_error, sol, times


def rounded_seed_peaks(forecast):
    # Begin at the printed point of Fig. 10(a). The first period starts at the
    # first complete local maximum, so the initial partial turn is excluded.
    stop = 9.5 * forecast["period"]
    sol = direct(forecast["published_seed"], stop, dense=True)
    grid = np.arange(0, stop, 0.04)
    coarse_peaks = find_peaks(sol.sol(grid)[0])[0]
    peaks = np.array([minimize_scalar(lambda t: -sol.sol(t)[0], method="bounded",
                                     bounds=(grid[i-1], grid[i+1]),
                                     options={"xatol": 1e-12}).x for i in coarse_peaks])
    if len(peaks) != 10:
        raise ArithmeticError("Expected ten complete voltage maxima")
    voltages = BREAK_VOLTAGE * sol.sol(peaks)[0]
    exact = np.diff(peaks) / STEP
    phase_traces = []
    for phase in PHASES:
        t = np.arange(phase * STEP, stop, STEP)
        phase_traces.append(BREAK_VOLTAGE * sol.sol(t)[0])
    rows, first_eight_ranges = [], []
    for quantum in QUANTA_MV:
        # find_peaks selects the middle sample of a flat maximum, rounding down.
        indices = [find_peaks(digitize(v, quantum))[0] for v in phase_traces]
        if any(len(i) != len(peaks) for i in indices):
            raise ArithmeticError("Peak count changed after digitization")
        spans = np.diff(indices, axis=1)
        first_eight_ranges.append(dict(quantum_mV=quantum,
                                       samples=limits(spans[:, :8].mean(axis=1))))
        for k, span in enumerate(exact):
            rows.append(dict(quantum_mV=quantum, interval=k+1,
                             first_peak_t_us=peaks[k] * TAU_US, last_peak_t_us=peaks[k+1] * TAU_US,
                             first_peak_V=voltages[k], last_peak_V=voltages[k+1],
                             continuous_samples=span, index_min=int(spans[:, k].min()),
                             index_max=int(spans[:, k].max())))
    u = (sol.y_events[0][:, 1:] - forecast["point"]) @ forecast["left_unstable_unit"]
    exit_index = int(np.flatnonzero(np.abs(u) >= forecast["projected_exit_radius"])[0])
    info = dict(start=forecast["published_seed"], first_exit_section=exit_index,
                first_exit_us=sol.t_events[0][exit_index] * TAU_US,
                first_positive_maximum_us=peaks[0] * TAU_US,
                first_rising_crossing_interval_samples=np.diff(sol.t_events[0])[0] / STEP,
                first_peak_interval_samples=exact[0], eighth_peak_interval_samples=exact[7],
                ninth_peak_interval_samples=exact[8], first_eight_mean_samples=exact[:8].mean(),
                first_eight_index_mean_ranges=first_eight_ranges)
    return info, rows


def plot(output, forecast, trace, times, rows, seed_peaks):
    saddle_samples = forecast["period"] / STEP
    fig, (wave, periods, late) = plt.subplots(1, 3, figsize=(13.4, 3.5), layout="constrained")
    wave_time = np.arange(times[0] - 12 * STEP, times[1] + 11 * STEP, STEP)
    elapsed = (wave_time - times[0]) / STEP
    wave.plot(elapsed, BREAK_VOLTAGE * trace.sol(wave_time)[0], color="#236868")
    wave.axhline(-BREAK_VOLTAGE, color="#909a9a", linestyle=":", linewidth=1)
    wave.axvline(218, color="#b24e22", linestyle="--")
    wave.scatter([0, (times[1] - times[0]) / STEP], [-BREAK_VOLTAGE] * 2,
                 color="#236868", zorder=3)
    wave.annotate("218 samples", (218, wave.get_ylim()[1]), xytext=(-3, -6),
                  ha="right", va="top", textcoords="offset points", color="#954019")
    wave.set(title="Simulated turn, δx = +10⁻⁵", xlabel="Samples from first rising crossing",
             ylabel="Capacitor voltage U₁ (V)")

    base = [r for r in rows if r["quantum_mV"] == 0]
    periods.hlines(range(len(base)), [r["single_exact_min"] for r in base],
                   [r["single_exact_max"] for r in base], color="#236868",
                   linewidth=3, label="Single turns")
    qualifying = [i for i, r in enumerate(base) if r["turns"] == 8]
    periods.plot([base[i]["eight_exact"] for i in qualifying], qualifying, "o",
                 markerfacecolor="white", markeredgecolor="#262a32", markersize=5,
                 label="Eight-turn mean")
    periods.set_yticks(range(len(base)), [f'{r["offset_x"]:+.0e}' for r in base])
    periods.invert_yaxis()
    periods.axvline(saddle_samples, color="#909a9a", linestyle=":", linewidth=1,
                    label="Saddle orbit")
    periods.set(title="Eight voltage precharges", ylabel="Precharge offset δx",
                xlabel="Continuous crossing interval (samples)")
    periods.ticklabel_format(axis="x", style="plain", useOffset=False)
    periods.legend(loc="lower left", fontsize=8)
    seed = [r for r in seed_peaks if r["quantum_mV"] == 0]
    late.plot([r["interval"] for r in seed], [r["continuous_samples"] for r in seed],
              ".-", color="#236868", label="Successive maxima")
    late.axhline(218, linestyle="--", color="#b24e22", label="Reported 218")
    late.set(title="Rounded Fig. 10(a) start", ylabel="Peak-to-peak period (samples)",
             xlabel="Interval between positive maxima", xticks=[1, 3, 5, 7, 9])
    late.legend(loc="lower left", fontsize=8)
    fig.savefig(output / "trace-periods.png", dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "outputs/chua-saddle-period")
    output = parser.parse_args().output_dir
    forecast = json.loads((RESIDENCE / "forecast.json").read_text())
    records = json.loads((RESIDENCE / "measurements.json").read_text())["records"]
    records = [r for r in records if abs(r["offset_x"]) in OFFSETS]
    if len(records) != 2 * len(OFFSETS):
        raise ValueError("Archived launch set is incomplete")
    rows, checks = [], []
    for record in records:
        measured, error, trace, times = evaluate(record, forecast["projected_exit_radius"])
        rows.extend(measured)
        checks.append(error)
        if record["offset_x"] == 1e-5:
            sample = trace, times

    seed_info, seed_rows = rounded_seed_peaks(forecast)
    summary = []
    for quantum in QUANTA_MV:
        singles = [r for r in rows if r["quantum_mV"] == quantum]
        eights = [r for r in singles if r["turns"] == 8]
        summary.append(dict(
            quantum_mV=quantum,
            single_index=limits([r[k] for r in singles for k in ("single_index_min", "single_index_max")]),
            single_linear=limits([r[k] for r in singles for k in ("single_linear_min", "single_linear_max")]),
            eight_index=limits([r[k] for r in eights for k in ("eight_index_min", "eight_index_max")]),
            eight_linear=limits([r[k] for r in eights for k in ("eight_linear_min", "eight_linear_max")]),
            single_max_abs_linear_error=max(r["max_abs_linear_error"] for r in singles),
            eight_max_abs_linear_error=max(r["eight_max_abs_linear_error"] for r in eights)))
    meta = dict(source="https://doi.org/10.1007/s11071-022-08078-y",
                scope="Synthetic U1 at nominal parameters for voltage precharges and the printed Fig. 10(a) start. Periods in samples per turn.",
                reference="Affine return times from results/chua-saddle-residence/measurements.json.",
                eligibility="Both endpoints and all intermediate crossings satisfy |u| < 0.02. At most eight turns from entry.",
                clock_phases_in_samples=PHASES.tolist(), hypothetical_quanta_mV=QUANTA_MV,
                R_ohm=788, C2_nF=311, tau_us=TAU_US, sample_us=SAMPLE_US, breakpoint_V=BREAK_VOLTAGE,
                published_samples_per_turn=218, nominal_orbit_samples=forecast["period"] / STEP,
                starts=len(records), eligible_single_turns=sum(r["turns"] for r in singles),
                eligible_eight_turn_starts=len(eights), max_affine_vs_ivp_event_error_samples=max(checks),
                continuous_single=limits([r[k] for r in singles for k in ("single_exact_min", "single_exact_max")]),
                continuous_eight=limits([r["eight_exact"] for r in eights]), by_quantum=summary,
                rounded_seed_peaks=seed_info)
    output.mkdir(parents=True, exist_ok=True)
    for name, data in [("trace-periods.csv", rows), ("trace-rounded-seed-peaks.csv", seed_rows)]:
        with (output / name).open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=data[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)
    (output / "trace-summary.json").write_text(json.dumps(meta, indent=2, allow_nan=False) + "\n",
                                              encoding="utf-8")
    plot(output, forecast, *sample, rows, seed_rows)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
