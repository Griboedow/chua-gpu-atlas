# Return geometry: records and reproduction

This directory supports the [return-geometry report](../../docs/chua-return-geometry.md). It contains numerical measurements of transient transverse amplification for 34 displayed reference cycles, a ten-member period-adding sequence with same-parameter saddle companions, direct nonlinear return samples, and prospective tests on new 95- and 105-peak cycles. No new interval or Lean certificate is claimed.

## Evidence and provenance

| Record | Contents |
|---|---|
| [`protocol.md`](protocol.md) | Questions, definitions, exploratory gates and the original section-separation conjecture, published in [fb72c48](https://github.com/Griboedow/chua-gpu-atlas/commit/fb72c48) before these measurements |
| [`forecasts-frozen.json`](forecasts-frozen.json) | P95/P105 width, gain, section-separation and multiplier forecasts, with source hashes, published in [8df8da8](https://github.com/Griboedow/chua-gpu-atlas/commit/8df8da8) before those searches |
| `ID-d70-s8.json`, `.npz` | Reference source/hash, reconstruction diagnostics, coarse phase maxima and full sampled propagator data |
| `ID-d85-s128.*`, selected `ID-d85-s32.*` | Finer phase grids for all ten sequence members; the digits label is working precision, not a certified error bound |
| `ID-optimized.json` | Bounded local refinement of start/end times around the coarse maximum; retains the coarse value and search bounds |
| `ID-companion-d70.json` | Same-β saddle root, signed half-return multiplier, section displacement, path samples and Newton residuals; P85 also has a 90-digit repeat |
| [`sequence-analysis.json`](sequence-analysis.json) | Matched section points, distances to the organizing saddle and between the two cycles, aligned at the largest excursion |
| `P55-normal-form-d75.json`, `P85-normal-form-d75.json`, `P85-normal-form-d90.json` | Nine actual secant-return samples per run, compared with (s + s²)/2; the 90-digit repeat also records sampled region checks |
| `P95-reference.json`, `P105-reference.json` | New symmetric representatives with signed half-return multiplier +1/2 |
| `P95-edges.json`, `P105-edges.json` | Fold/symmetry endpoints, widths and frozen-forecast hash |
| `ID-dDIG-sSUB-check.json` | Independent full-field variational integrations at two step limits, plus finite perturbation controls |
| `P95-independent-orbit.json`, `P105-independent-orbit.json` | Whole-period full-field integration and separately restarted next-section variation; their errors have different meanings |
| `P95-precision-check.json`, `P105-precision-check.json` | 95-digit root repeats and independent eigendecomposition of each affine flow, with peak and section counts |
| `P95-edge-precision-check.json`, `P105-edge-precision-check.json` | 95-digit repeats of both window endpoints |
| [`energy-metric-check.json`](energy-metric-check.json) | Coarse-mesh amplitude gains in Euclidean and stored-energy norms |
| [`review.json`](review.json) | Recomputed forecast scores, source checks and aggregate numerical controls |
| [`environment.json`](environment.json) | Python and numerical-library versions used for the local computations |
| `chua-return-*.pdf` | Four vector figures; corresponding PNGs are displayed in the report |

The NPZ arrays store physical time `t`, reference states `xyz`, velocity, orthonormal normal frames `basis`, and the local 2 × 2 quotient propagators `local`. `gain` is the best sampled gain for each initial phase; `maximizing_steps` and `maximizing_duration` identify its endpoint. `worst_history` contains operator norms along the period starting at the largest sampled phase. Its optimal initial direction can change with the displayed endpoint time.

The original reference roots and earlier widths remain in their original directories. The new cycles are additional numerical test objects; they have not been added to the 34-card display registry, assigned knot types or given three-dimensional family meshes.

## Quick audit of the published records

Run from the repository root with NumPy, SciPy, mpmath 1.3.0 and Matplotlib installed. The existing `requirements.txt` and `requirements-precision.txt` specify compatible versions. This study uses the CPU and does not require OpenCL, the GPU surveys or GitHub Actions.

```bash
python research/review_chua_return_geometry.py
python research/manifest_refined_chua_atlas.py --check
```

The first command recomputes scores and consistency gates from the stored measurements. The second verifies file integrity. Neither replaces replaying a numerical calculation or checking the mathematical interpretation.

## Replay without replacing published evidence

The programs honor `CHUA_RETURN_OUTPUT`. On the Windows workstation, run the following PowerShell setup from a checkout on drive D:. It creates a new local directory and preserves the original forecasts and their calibration records. It refuses to reuse an existing replay directory.

```powershell
$returnReplay = Join-Path (Get-Location) '.publication/return-geometry-replay'
if (Test-Path -LiteralPath $returnReplay) { throw 'Choose a new replay directory.' }
New-Item -ItemType Directory -Path $returnReplay | Out-Null
Copy-Item -LiteralPath 'results/chua-return-geometry/protocol.md' -Destination $returnReplay
Copy-Item -LiteralPath 'results/chua-return-geometry/forecasts-frozen.json' -Destination $returnReplay
$env:CHUA_RETURN_OUTPUT = $returnReplay
```

Use the checkout's Python environment. For the commands below, `python` can be replaced by `.\.venv\Scripts\python.exe -X utf8`. Keep the published directory present: the frozen source hashes identify those original calibration files, including run metadata. A replay is compared with those records; it does not rewrite their history.

First reconstruct the original references and measure the initial sequence:

```powershell
python research/study_chua_return_geometry.py
python research/study_chua_return_geometry.py --ids S04 S05 S06 S07 P55 P65 P75 P85 --subdivisions 128 --digits 85
python research/study_chua_return_geometry.py --ids S04 P65 P75 P85 --subdivisions 32 --digits 85
python research/refine_chua_transient_maxima.py
python research/study_chua_companion_saddles.py
python research/study_chua_companion_saddles.py --ids P85 --digits 90
python research/analyze_chua_return_geometry.py
python research/sample_chua_return_normal_form.py --ids P55 P85 --digits 75
python research/sample_chua_return_normal_form.py --ids P85 --digits 90
```

Then solve the two prospective targets using the retained forecast file. Do not run `--freeze` over the published forecasts.

```powershell
python research/extend_chua_return_tests.py
$returnRoots = @((Join-Path $returnReplay 'P95-reference.json'), (Join-Path $returnReplay 'P105-reference.json'))
python research/study_chua_return_geometry.py --root-files $returnRoots
python research/study_chua_return_geometry.py --root-files $returnRoots --subdivisions 128 --digits 85
python research/refine_chua_transient_maxima.py --ids P95 P105
python research/study_chua_companion_saddles.py --root-files $returnRoots
```

Run the independent controls and produce the audit and figures:

```powershell
python research/check_chua_return_geometry.py --labels S04-d85-s32 P85-d85-s32 P95-d70-s8 P105-d70-s8 --finite
python research/check_chua_gain_metric.py
python research/check_chua_new_return_orbits.py
python research/check_chua_new_return_orbits.py --edges
python research/check_chua_new_return_orbits.py --full-field
python research/review_chua_return_geometry.py
python research/plot_chua_return_geometry.py
Remove-Item Env:CHUA_RETURN_OUTPUT
```

With the output override, figures are written beside the replay data, preserving the published assets. Most expensive producers skip records that already exist; the analysis, audit, precision checks and plotting commands recompute their summaries. Use a new replay directory when changing the experiment. The fine-grid gain scan is quadratic in its number of phase samples, and the high-precision event solves can take several minutes.

## Retained failures and corrections

- The initial protocol proposed a uniform Λ⁻ⁿ companion-separation scale. Phase-aligned measurements instead distinguish an inherited-section scale Λ⁻²ⁿ from a maximum matched section separation of order Λ⁻ⁿ. The corrected rule was frozen before the P95/P105 tests.
- [`P105-initial-trial.json`](P105-initial-trial.json) retains the 18-iteration Newton failure. Only the fold-solve iteration budget was increased to 36; the residual gate was unchanged. The original 18-iteration default remains for earlier callers.
- [`P95-full-field.json`](P95-full-field.json) retains an initial whole-period integration with maximum step 0.005 and closure error about 1.22 × 10⁻⁶. The later file records the step-0.00125 repeat.
- P105's whole-period integration still drifts by about 4.91 × 10⁻⁶. This cannot resolve the 3.72 × 10⁻¹² companion separation; high-precision roots and separate section controls address different parts of the numerical problem.
- Finite perturbation errors are nonmonotone across the tested amplitudes. All amplitudes are retained rather than selecting only the best agreement.

The 75-digit return samples precede the addition of sampled path diagnostics; the later P85 90-digit record includes them. A replay using the final program adds those fields to its own 75-digit files. Timing, floating-point roundoff and such additional diagnostics can differ from the historical records; scientific quantities and stated gates are the comparison targets.

## Review targets

Check the quotient-flow construction and the chosen norm before interpreting the gain. Verify that every propagated flight belongs to the prescribed continuous Chua field, that a matched companion is a distinct root at identical β, and that section phase is specified before comparing separations. Score forecasts against the subsequently computed values, not against a rescaled copy of the forecasts.

The independent ODE controls validate selected sampled gains, not a global maximum over all start/end phases. The normal-form samples are a finite local test, not a uniform error bound. A saddle root alone does not identify a basin boundary. The report's relation W G² ≈ constant is a numerical sequence result; the previous infinity theorem and the global completeness of the atlas are separate questions.
