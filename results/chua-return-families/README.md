# Return geometry in two further Chua sequences

This directory extends the S08-sequence study to long cycles near the S11 and E05 short-cycle branches. [Report](../../docs/chua-return-families.md).

The protocol was published in commit [bd9e698](https://github.com/Griboedow/chua-gpu-atlas/commit/bd9e698); calibration and numerical forecasts were frozen in [a2444a0](https://github.com/Griboedow/chua-gpu-atlas/commit/a2444a0). Six subsequent roots have the predicted peak counts. Corrected forecasts pass the prescribed gates for 6/6 periods, 4/6 gains and 0/6 widths. Passing numerical verification does not turn a failed forecast into a success.

## Records

| File pattern | Contents |
|---|---|
| `protocol.md`, `calibration.json`, `forecasts-frozen.json` | Prospective protocol, two known long members per sequence, short-saddle measurements and unaltered forecasts |
| `*-reference.json`, `*-edges.json` | 70-digit roots at half multiplier +1/2 and endpoint roots at +1 and −1 |
| `*-search-failure.json`, `*-search-retries.json` | Original failed Newton trial and subsequent predictor attempts |
| `*-audit.json` | 90-digit repeats, independent eigen-flow reconstruction, full-field counts, section derivatives and eight numerical gates |
| `*-gain-s8.*`, `*-gain-s128.*` | Coarse and fine transverse propagators, sampled maxima and coarse-seeded local refinement |
| `*-gain-control.json` | Independent full-field variational checks at two step limits |
| `*-companion.json`, `*-normal-form-d*.json` | Same-parameter saddle roots and actual nonlinear secant returns |
| `*-geometry.json`, `geometry-summary.json`, `*-parent.*` | Continued short saddles, section distances and plotting arrays |
| `*-directions.json` | One transported Floquet direction and the sampled maximizing input direction; not a new optimum at every elapsed time |
| `scores.json`, `review.json` | All forecast outcomes and numerical record audit |
| `*.pdf` | Vector versions of the four publication figures |

`S11-n3` through `S11-n5` have 35, 45 and 55 peaks. `E05-n3` through `E05-n5` have 49, 63 and 77. These identifiers distinguish the sequences even when their peak counts coincide with other atlas entries. The short reference branch label does not imply that its original reference orbit was a saddle.

## Reproduction

Run locally from the repository root using its Python environment. Keep the working copy, environment and temporary output on D:. No GitHub Actions are needed. Search and measurement stages retain completed outputs. For an independent replay, use a separate checkout and change the `OUT` constant in `test_chua_return_families.py` there. Copy only the original protocol, calibration and frozen forecasts into that alternative directory, then run from `search` onward. This preserves the predictions and recomputes their targets. Repeat calibration separately if desired; do not overwrite the historical calibration or regenerate its timestamped forecast file.

The original order was:

```powershell
.\.venv\Scripts\python.exe -X utf8 research/test_chua_return_families.py calibrate
.\.venv\Scripts\python.exe -X utf8 research/test_chua_return_families.py freeze
# Publish the forecast commit before running the search.
.\.venv\Scripts\python.exe -X utf8 research/test_chua_return_families.py search
.\.venv\Scripts\python.exe -X utf8 research/test_chua_return_families.py gains
.\.venv\Scripts\python.exe -X utf8 research/test_chua_return_families.py gains --subdivisions 128
.\.venv\Scripts\python.exe -X utf8 research/audit_chua_return_families.py orbits
.\.venv\Scripts\python.exe -X utf8 research/audit_chua_return_families.py companions
.\.venv\Scripts\python.exe -X utf8 research/audit_chua_return_families.py normal --ids S11-n5 E05-n5
.\.venv\Scripts\python.exe -X utf8 research/audit_chua_return_families.py normal --ids S11-n5 E05-n5 --digits 90
.\.venv\Scripts\python.exe -X utf8 research/audit_chua_return_families.py gain --ids S11-n5 E05-n5
.\.venv\Scripts\python.exe -X utf8 research/analyze_chua_return_families.py
.\.venv\Scripts\python.exe -X utf8 research/review_chua_return_families.py
.\.venv\Scripts\python.exe -X utf8 research/plot_chua_return_families.py
```

The initial direct high-precision S11-n3 solve failed. The current search first uses a double-precision predictor and retains its attempts. The original trial can be replayed with the search implementation in `a2444a0`. No numerical forecast was changed during this correction.

The audit command can be run directly on the published records. It compares the frozen file to its original Git commit, so the checkout must include that commit. It checks scientific outcomes separately from implementation validity: a correct audit is expected to retain all six width misses and two gain misses.

## Review targets

1. Confirm the parameter convention and distinguish signed half-return multipliers from full-period multipliers.
2. Reconstruct forecast coefficients exclusively from the two old members and the short saddle at the second member's β. Keep corrected root locations out of the prediction score.
3. Check region itineraries, primitive section points and independent peak counts. A converged algebraic shooting root alone is insufficient.
4. Inspect whole-period drift separately from independently restarted next-section errors. The latter are not a global orbit error bound.
5. Compare fixed-direction amplitude curves with operator-norm gain curves. The section norm and the three-dimensional transverse quotient norm are also different measurements.
6. Check the companion at the identical β and recompute event times for secant samples. Nine samples do not bound a global normal-form remainder.
7. Keep the 34 illustrated references, six additional roots, three candidate organizing mechanisms and the separately proved infinity result distinct.

These are numerical results. They add no interval certificate, no complete atlas classification and no claim of infinitely many layers for either new sequence.
