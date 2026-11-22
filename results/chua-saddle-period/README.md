# Saddle-cycle period records

The [research note](../../docs/chua-saddle-period.md) compares the marked hardware interval in Fig. 10 with the model orbit, voltage traces sampled at the oscilloscope cadence, and peak positions in the published image.

## Interval certificate

[interval.json](interval.json) gives the parameter ranges, grid, physical clock, global bounds and one row per cell. The `coupled_period_microseconds` field uses the identity $R C_2=\sqrt{\beta L C_2}$. The wider `physical_period_microseconds` field varies the clock independently of the five model parameters. Each printed Arb endpoint includes its rounding radius.

The [certificate program](../../research/chua_saddle_period_bound.py) matches two successive affine flights to the central inversion of their starting point. A Krawczyk inclusion and a contraction bound establish a unique root in each shooting box at each parameter value. Flight guards check the crossing order and transversality. The signed half-return derivative certifies saddle stability. On each shared face, both boxes contain the same enclosed root at its midpoint. Local uniqueness continues that agreement over the connected face.

## Numerical records

| Files | Contents |
|---|---|
| [Component periods](component-periods.csv) · [summary](period-summary.json) | Three nominal conventions and 64 compatible draws from the [component settings](../chua-precharge/settings.csv). Columns include dimensionless period, clock and period in µs, section closure and return multipliers. The summary includes independent DOP853 checks at the nominal components and sample extrema. |
| [Voltage returns](trace-periods.csv) · [summary](trace-summary.json) | Eight starts in four signed pairs. Periods are in sampling intervals. Each row summarizes 16 clock phases at one voltage resolution. Empty eight-turn fields mean the trajectory left the neighborhood sooner. |
| [Voltage peaks](trace-rounded-seed-peaks.csv) | Successive peak intervals from the rounded Fig. 10(a) starting state, peak voltages and integer readings across clock phases. |

## Published traces

[Peak positions](figure10-peaks.csv) give horizontal pixel coordinates for 17 interior peaks in panel (a) and 20 in (b). The [spacing ratios](figure10-ratios.json) compare the mean of the last and first complete spacings in each panel. The main comparison uses four spacings per end and a grayscale cutoff of 110 out of 255. The file repeats it for three to five spacings and cutoffs of 60, 110 and 160. The source image digest pins the pixel coordinates to the [publisher’s 1,200-pixel figure](https://media.springernature.com/lw1200/springer-static/image/art%3A10.1007%2Fs11071-022-08078-y/MediaObjects/11071_2022_8078_Fig10_HTML.png).

## Reproduce

Run from the repository root. The scripts write to `outputs/chua-saddle-period/`. Save the linked figure locally and supply its path for the final command.

```bash
python -m pip install -r requirements-certificate.txt
python research/chua_saddle_period_bound.py --workers 4
python research/chua_saddle_period.py --certificate outputs/chua-saddle-period/interval.json
python research/chua_saddle_period_trace.py
python research/chua_saddle_period_figure10.py path/to/11071_2022_8078_Fig10_HTML.png
```
