# Prediction, impulse and hidden-coexistence evidence

[Combined report](../../docs/chua-prediction-and-hidden-coexistence.md) · [Hidden-infinity proof](../../docs/chua-hidden-infinity-certificate.md) · [Computed summary](analysis.json)

## Slope-sheet predictions

- [`sheets/protocol.json`](sheets/protocol.json) and [`transport-protocol.json`](sheets/transport-protocol.json): training sets, two models, held-out definition and scoring rules.
- [`forecasts-frozen.json`](sheets/forecasts-frozen.json): predictions published in [112b01f](https://github.com/Griboedow/chua-gpu-atlas/commit/112b01f), before the 48 off-centre endpoint searches. Central target folds are calibration, not validation.
- `sheets/*-calibration.json`: the 18 local organizers, additional-iterate checks, training folds and coefficients.
- `sheets/*-search.json`: all 54 targets, including the six central calibration targets, with fold, symmetry and interior roots. The uncorrected 85-peak forecasts miss all sixteen held-out cores. Failed initial guesses remain recorded even when later continuation succeeds.
- [`sheets/independent-checks.json`](sheets/independent-checks.json): twelve 65-digit direct-exponential refinements, independent ODE/variational checks, peak counts and displaced-parameter rejection controls.

The corrected predictor passes 48/48 midpoint tests. The maximum corrected fold error is 0.022103 measured widths; the maximum width error is 6.6803%. These are local numerical validations, not interval certificates for the entire interpolated surfaces.

```bash
python research/predict_chua_slope_sheets.py search
python research/check_chua_followup.py sheets
```

The `init`, `calibrate`, `central` and `freeze` phases reconstruct the prediction stage. Run them in a separate checkout/output copy if preserving the published files. Re-fitting after inspecting targets would not constitute a fresh prospective test.

## Impulse certificates

- [`trigger/screen.json`](trigger/screen.json): 567 finite-amplitude tests over seven rays.
- [`thresholds.json`](trigger/thresholds.json): finite-iterate crossing approximants.
- [`local-certificate-80.json`](trigger/local-certificate-80.json), [`local-certificate-100.json`](trigger/local-certificate-100.json): seven exact stable-graph intersections, opposite basin destinations and full-flight checks. Each bracket has width 2 × 10⁻¹¹.
- [`below-certificate-80.json`](trigger/below-certificate-80.json): complete lower-amplitude coverage along the y ray; its `passed`, `partial` and `coverage` fields state completion separately from the local certificates. `below-prefix-80.json` and `below-part-80-*.json` retain the locally computed chunks.
- [`delay-forecast.json`](trigger/delay-forecast.json), [`delay.json`](trigger/delay.json): multiplier-based slope prediction, 100 measurements and four fitted slopes. This prediction was saved before its numerical measurements but was not separately published as a blind test.
- [`independent-checks.json`](trigger/independent-checks.json): twelve long DOP853/Radau destination controls.
- `local-certificate-80-four-returns-rejected.json`: the insufficient preliminary enclosure, preserved as a rejection.

```bash
python research/study_chua_trigger.py --phase thresholds
python research/study_chua_trigger.py --phase delay
python research/certify_chua_trigger.py local --digits 80
python research/certify_chua_trigger.py local --digits 100
python research/certify_chua_trigger.py below --digits 80
python research/check_chua_followup.py trigger
```

The all-below calculation passes with 662 adjoining intervals after the central segment [0,0.299]. It is slower than the local certificate because intervals meeting a diode contact require continuous-flow enclosures and subdivision. It was divided into independent local CPU processes after original bands 0–30. `finalize_chua_trigger.py` joins those retained chunks only after checking their completion and every shared endpoint. A fresh sequential `below` run covers the same domain without this merging step. The full lower-amplitude proof was run at 80 digits; the seven local crossing proofs were repeated at both precisions.

## Hidden-region search and certificate

- `hidden-search/*screen*.json` and `.npz`: stated search sectors, starts, settings, candidates and full finite-time classifications. Stable-origin surveys contain 18,432 parameter samples and no passing chaos candidates. The separate stable-outer-equilibrium survey contains 4,096 samples and six passing candidates. Neither finding proves chaotic attraction or absence.
- [`candidate-controls.json`](hidden-search/candidate-controls.json): longer runs at two smaller steps and starts near both unstable-origin exits.
- [`classical-continuation.json`](hidden-search/classical-continuation.json): the stopped numerical continuation of the earlier tangency. A failed corrector is not a proof of branch termination.
- [`lobe-saddles.json`](hidden-search/lobe-saddles.json), `connection-screen-818*.json`, `connection-newton-818*.json`: saddles, finite-iterate scans and Newton searches. These exploratory records make no manifold-existence claim.
- [`connection-center-14-2.json`](hidden-search/connection-center-14-2.json): high-precision numerical center, fixed decimal coordinates and parameters for the certificate. The failed lower-precision refinement is retained separately.
- [`hidden-tangency-140.json`](hidden-search/hidden-tangency-140.json), [`hidden-tangency-160.json`](hidden-search/hidden-tangency-160.json): all interval conditions pass. The return is F on x = 1; the inherited `full_return_*` output fields describe F². The recorded 36 F returns use 72 affine flights and 14,555 global time cells.
- [`hidden-tangency-140-failure.json`](hidden-search/hidden-tangency-140-failure.json): retained fixed-coordinate wrapping failure. The subsequent rational moving-frame calculation succeeds without dropping interval uncertainty.
- [`uniform-equilibrium-exclusion-80.json`](hidden-search/uniform-equilibrium-exclusion-80.json), [`uniform-equilibrium-exclusion-100.json`](hidden-search/uniform-equilibrium-exclusion-100.json): complete equilibrium neighborhoods excluded from every periodic-attractor basin on the surrounding β interval.
- `equilibrium-exclusion-80.json`, `equilibrium-exclusion-100.json`: the corresponding pointwise tests for all six original candidates. These alone do not prove existence of the apparent chaotic attractors.
- [`independent-tangency-checks.json`](hidden-search/independent-tangency-checks.json): separate 90-digit mpmath flights, two finite-difference scales, direct ODE controls and two rejected shifted β boxes. `classical-regression-80.json` checks the old classical certificate after shared implementation changes.

```bash
python research/refine_chua_hidden_connection.py --digits 140
python research/certify_chua_hidden_tangency.py --digits 140
python research/certify_chua_hidden_tangency.py --digits 160
python research/certify_chua_hidden_neighborhoods.py --uniform --digits 80
python research/certify_chua_hidden_neighborhoods.py --uniform --digits 100
python research/check_chua_hidden_tangency.py
```

Exploratory scans can be rerun with `search_chua_hidden_tangencies.py` phases `continuation`, `screen`, `broad`, `outer`, followed by `check_chua_hidden_candidates.py`. They require the local OpenCL environment. The interval and mpmath checks run on the CPU. Numerical seed reconstruction may produce different last digits while the stored center remains inside the certified box; it is a separate operation from replaying the published certificate.

## Figures and integrity

The [browser check](browser-check.json) records inspection of the rendered GitHub reports and README, including the displayed equations and figures. The local screenshots accompany the working copy; their hashes are retained in the record.

- [Slope-sheet comparison, PDF](chua-followup-sheets.pdf)
- [Impulse and delay comparison, PDF](chua-followup-trigger.pdf)
- [Hidden-connection certificate, PDF](chua-followup-hidden.pdf)

```bash
python research/analyze_chua_followup.py
python research/audit_chua_followup.py
python research/plot_chua_followup.py
python research/manifest_refined_chua_atlas.py --check
```

`pilot-fold.json` and `pilot-spectral.json` are preliminary solver comparisons, not held-out forecast evidence. All published results and programs enter the repository manifest. Checksums establish provenance, not mathematical truth. The [proof report](../../docs/chua-hidden-infinity-certificate.md) identifies the analytic and implementation assumptions; no complete Lean proof or historical-priority claim is made.
