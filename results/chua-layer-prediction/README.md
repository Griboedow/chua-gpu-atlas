# Evidence for theory-first layer predictions

The [research report](../../docs/chua-layer-prediction.md) compares forecasts published in [0fc47ea](https://github.com/Griboedow/chua-gpu-atlas/commit/0fc47ea) with a later targeted search. The model is the classical piecewise-linear Chua system at alpha = 9, m0 = -8/7, m1 = -5/7. Two new attracting reference cycles have 75 and 85 maxima, extending the operational catalogue from 32 to 34. Their complete three-dimensional sheets were not continued in this update.

The 75-peak forecast overlaps the measured core. The 85-peak forecast misses it by several measured widths; orbit correction subsequently finds it. The width discrepancies are 1.89% and 2.70% of the forecasts. These numerical checks preceded the separate [interval certificate](../chua-layer-infinity/README.md), which now establishes a tangency and an infinite cascade. The frozen forecasts and recorded numerical outcomes here remain unchanged.

## Evidence files

| File | Meaning |
|---|---|
| [predictions-frozen.json](predictions-frozen.json) | Unchanged pre-search forecasts, input hashes, calibration and engineering gates. |
| [search-results.json](search-results.json) | Prescribed itineraries, rejected float64 trials, 55-digit roots at fold/symmetry/core and forecast errors. Decimal strings preserve the computed values. |
| [precision-implementation-checks.json](precision-implementation-checks.json) | Independent finite-difference Jacobian checks and comparisons with two previously known roots. |
| [implementation-trial-01.json](implementation-trial-01.json) | Failed initial implementation: negative indexing of an mpmath matrix accessed the wrong element. Preserved for provenance; not an absence-of-orbit result. |
| [independent-checks.json](independent-checks.json) | 75-digit repeats, original affine implementation, DOP853/Radau full-field integrations and variational controls. Includes rejected whole-period multiplier checks and the passing section-by-section checks. |
| [analysis.json](analysis.json) | Forecast scoring, scaling ratios, normalized core coordinates and measured section sequences. |
| [organizer-approximants.json](organizer-approximants.json) | Six finite-iterate tangent-line approximations, derivative controls and sampled return folds. Post-search exploratory evidence, in float64; their agreement is not a rigorous error bound. |
| [motion-descriptions.json](motion-descriptions.json), [CSV](motion-descriptions.csv) | Original 32-family route words, residence times, symmetry, periods and coarse-code groups. |
| [extended-family-registry.json](extended-family-registry.json) | Original 32 entries plus P75/P85. Compatible float64 core fields are accompanied by the high-precision core and evidence links. The original registry is unchanged. |
| [predicted-orbits.npz](predicted-orbits.npz) | Sampled trajectories for the two new portraits. These illustrate the checked solutions; the independent solvers count peaks and section returns. |

All five figures are available as vector PDFs: [forecast comparison](chua-layer-prediction-proof.pdf), [scaling](chua-layer-prediction-scales.pdf), [motion descriptions](chua-layer-motion-description.pdf), [new motions](chua-layer-predicted-motions.pdf), and [organizing geometry](chua-layer-prediction-organizer.pdf). Their PNG versions appear in the report and under `docs/assets/`. The [atlas manifest](../chua-layer-atlas/manifest.json) covers the study, code and input evidence.

## What passed, and what did not

All final orbit-audit gates pass. A 55-to-75-digit repeat changes beta by less than 7.4e-42 for these records; this is precision sensitivity within one implementation, not a certified accuracy statement. Direct full-field DOP853 and Radau integrations independently recover the peak counts, primitive section sequences and inversion symmetry. Their full-period state closure errors are below 2e-8.

For 85 peaks, three whole-period variational calculations miss the multiplier-error gate of 2e-5. The record retains errors of 7.66e-5, 2.22e-5 and 7.45e-5. Projecting out the neutral flow direction at each successive section crossing, then composing the two-dimensional return derivatives, gives errors of 3.21e-7 and 2.29e-7 at two step settings. The local next-section errors below 2.2e-14 must not be reported as whole-period closure errors. The independent full-period state tests remain separate.

Both forecasts pass the frozen location gate of 2e-8, but that tolerance is wider than either new stable core. The 85-peak point forecast is therefore still a **fine-scale location miss**, with zero interval overlap. No forecast was replaced after the search.

The organizer script solves a fold of a finitely propagated initial tangent line relative to a terminal tangent stable line. It checks multiple numbers of local/terminal iterates and samples nonzero curvature and parameter unfolding. It does not bound the difference between those lines and the exact invariant manifolds. The contraction argument in the report proves an implication under uniform return-map assumptions; it does not assert those assumptions have been proved here.

## Reproduction

Use Python with the repository's `requirements.txt` and the additional pinned `requirements-precision.txt`. The published run used mpmath 1.3.0. Commands below are run from the repository root. The separate output directory preserves the original frozen predictions and search; those two programs refuse to overwrite their existing result files.

```sh
python -m pip install -r requirements.txt -r requirements-precision.txt
python research/check_chua_precision_flights.py
python research/predict_chua_layers.py --output results/chua-layer-prediction-reproduction
python research/verify_chua_layer_predictions.py --output results/chua-layer-prediction-reproduction
python research/audit_chua_layer_predictions.py --output results/chua-layer-prediction-reproduction
```

The initial audit may exit with status 2 when a numerical gate fails, as it did for the published 85-peak whole-period variational control. Preserve that record, then run the explicit refinements separately:

```sh
python research/audit_chua_layer_predictions.py --output results/chua-layer-prediction-reproduction --refine-failed
python research/audit_chua_layer_predictions.py --output results/chua-layer-prediction-reproduction --section-controls
```

The first refinement also failed in the published run; the second changed the conditioning strategy without loosening the gate. Reference values and numerical errors can vary with the numerical libraries. The reproduction forecast has a fresh creation timestamp and therefore a different file hash; its numerical inputs and predictions should agree. The original publication commit establishes the chronology, not a newly generated timestamp.

The following commands regenerate descriptions, analysis and figures **from the published evidence directory**. They do not automatically substitute a separate reproduction run:

```sh
python research/describe_chua_layer_words.py
python research/analyze_chua_prediction_evidence.py
python research/trace_chua_layer_organizer.py
python research/plot_chua_layer_predictions.py
python research/manifest_refined_chua_atlas.py --check
```

The arbitrary-precision solver is implemented specifically for this fixed classical-parameter transect. Its exact affine-flow formula does not turn its floating-point evaluations into interval certification. Forecast calibration uses four known locations and a short-saddle multiplier; the two new cycles are a small targeted test, not an exhaustive classification or a claim of scientific priority over all prior Chua work.
