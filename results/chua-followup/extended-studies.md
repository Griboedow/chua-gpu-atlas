# Excitation, coexistence-box and extended-slope data

The [report](../../docs/chua-excitation-robustness-prediction.md) separates numerical results from uniform interval statements.

| Directory | Contents |
|---|---|
| [excitation-plane](excitation-plane/) | Polar basin grid at two time steps; directional energy thresholds; rectangular-pulse grid and refined crossings; local pulse-cost minimum; independent DOP853/Radau checks; demonstration trajectories |
| [coexistence-box](coexistence-box/) | Four-cycle uniform interval certificate, rejected wider enclosures, 31 × 31 numerical slope section, eight selection trajectories and a conservative interval component image |
| [extended-slopes](extended-slopes/) | Previously frozen forecasts, training data, subsequent endpoint solves, recalculated scores and twelve independent orbit checks |

The extension forecasts were published in [0c39117](https://github.com/Griboedow/chua-gpu-atlas/commit/0c39117), before target searches. Their calibration data and central bias correction must not be refitted when reproducing the scores. `initialization-attempts` retains superseded training guesses; these are not successful branch identifications.

## Reproduction

Run from the repository root, using the existing environment on drive D. Numba accelerates only the nonrigorous basin/pulse scan. The interval proof uses python-flint/Arb, as specified in `requirements-certificate.txt`.

```powershell
.\.venv\Scripts\python.exe -X utf8 -m pip install -r requirements-excitation.txt
.\.venv\Scripts\python.exe -X utf8 research/study_chua_excitation_plane.py plane
.\.venv\Scripts\python.exe -X utf8 research/study_chua_excitation_plane.py pulses
.\.venv\Scripts\python.exe -X utf8 research/study_chua_excitation_plane.py check
.\.venv\Scripts\python.exe -X utf8 research/derive_chua_pulse_limit.py
.\.venv\Scripts\python.exe -X utf8 research/certify_chua_coexistence_box.py --digits 100 --radius 1e-8 --state-radius auto --output replay-100.json
.\.venv\Scripts\python.exe -X utf8 research/check_chua_coexistence_box.py
.\.venv\Scripts\python.exe -X utf8 research/check_chua_extended_predictions.py
.\.venv\Scripts\python.exe -X utf8 research/review_chua_extended_studies.py
.\.venv\Scripts\python.exe -X utf8 research/plot_chua_extended_studies.py
```

The basin and pulse commands replace their numerical output files; use a separate checkout to compare against the publication. The certificate command above writes a new file and preserves `certificate.json`. The independent prediction checker resumes retained checks, verifies the frozen forecast against Git and recomputes all scores; remove its output in a disposable checkout to replay every selected integration.

The slope/component numerical screen can be repeated with `research/study_chua_coexistence_tolerances.py`. Its failed solves are retained as unresolved. The physical component image is an interval calculation, not a random tolerance sample.

The slope grid stores section points as well as continuation status. Counts merge attracting returns within 10⁻⁶: continuation references are not themselves distinct attractors. There are 124 nodes with four distinct attracting returns and 332 nodes where at least two reference continuations converge to the same return. These are numerical counts; the uniform certificate separately proves distinctness by interval separation.

To reproduce the prospective experiment from scratch in a separate checkout, run the `init`, `calibrate`, `freeze`, then `search` phases of `research/extend_chua_slope_predictions.py`. The existing frozen file is deliberately protected against overwrite. Do not run `search` before publishing the new frozen forecast when proposing another out-of-sample test.

## Meaning of the checks

The new record review checks Krawczyk inclusions, contractions, nominal Schur margins, determinant bounds, exclusion of ±1, separations and component containment. Its negative controls alter a contraction bound and an inclusion interval. It does not replace rerunning the interval program or reading the stability argument.

For inclusion tests, read the retained lower and upper endpoint enclosures. The short `ball` display strings can be much wider after decimal formatting; using those displays in place of the endpoints can lose an otherwise verified strict inclusion. The reviewer uses inward bounds for the containing box and outward bounds for its enclosed image.

The basin labels and bisection endpoints are numerical. Subdivision to a small scalar bracket does not account for integration error. Threshold-side checks with two adaptive integrators are retained separately. One preliminary destination check was unresolved at a section tolerance of 10⁻⁷: DOP853's final distance was 1.33 × 10⁻⁷ while Radau's was 2.63 × 10⁻⁹. The final declared tolerance is 10⁻⁶; this change affects the finite-time label, not the basin-grid or bisection calculation. No interval certificate is claimed for those new thresholds.

The uniform coexistence certificate is distinct: every parameter in its stated box is covered by interval enclosures, including flight admissibility. Its engineering tolerance image is extremely conservative and is not a practical component specification.
