# Recorded evidence

The directories below contain the settings, data, certificates and numerical controls supporting the reports. Recent layer studies are included directly in this repository; the historical release supplies the larger survey arrays. The [research guide](../docs/research-guide.md) connects the studies and explains their different evidence standards.

- [`chua-layer-infinity/`](chua-layer-infinity/README.md): interval certificates at 80 and 100 digits, two rejected parameter boxes, independent local/global derivative controls, a reconstructed numerical center and two figures. The [proof report](../docs/chua-infinity-certificate.md) derives an infinite cascade from the certified quadratic tangency.

- [`chua-layer-prediction/`](chua-layer-prediction/README.md): forecasts published before the 75/85-peak search, 55/75-digit orbit checks, independent full-field controls, preserved forecast and numerical failures, the preliminary organizing geometry, a 34-entry reference registry and five scientific figures. See [the prediction report](../docs/chua-layer-prediction.md).

- [`chua-layer-atlas/`](../docs/chua-refined-atlas.md): the 32-family mesh snapshot, registry, ancestry, continuation checks, dense README view and [file manifest](chua-layer-atlas/manifest.json).

- [`chua-slope-plane/`](../docs/chua-slope-plane.md), [`chua-slope-volumes/`](../docs/chua-slope-volumes.md), [`chua-slope-families/`](../docs/chua-slope-families.md): sampled slope sections and volumes, continued families and basin controls.

- [`chua-overlap-geometry/`](chua-overlap-geometry/README.md): symmetry-resolved P5/P7 coexistence, 18 contact curves, projection extrema, checked two-/three-/four-cycle examples, predictive tests, six figures and a self-contained visual guide. See [the research note](../docs/chua-overlap-geometry.md).

- `original/`: two 51³ scans, reference trajectories and eight numerical validation checks.
- `detailed/`: two 81³ grids, 41³ local zoom, beta transects, example regimes, connectivity summaries and independent validation.
- `expanded/`: adaptive solver controls, expanded-scout settings, coordinate extrema and corrected full-field DOP853 checks.
- `extent/`: wide-scout settings, repeated-candidate interpretation, per-parameter seed mapping and historical v1.0 affine-branch event checks (superseded; see correction below).

The large basin metadata, complete NPZ arrays, all offline interactive reports and browser-check records are in [release v1.1.0](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.1.0). Extract the archive into the repository root. Its internal `RESULTS_MANIFEST.json` provides per-file SHA-256 hashes; `SHA256SUMS.txt` verifies the ZIP itself.

These are finite-time numerical experiments. The old affine-branch CPU reference in `extent/event_resolved_checks.json` has been superseded after a missed-crossing failure at a distant point. Those files are historical evidence, not definitive reference values; use the new full-field checks under `expanded/`. Quantitative convergence over the entire cloud is not established. See the main README and methods for interpretation.

- `chua-scaling/`: fixed-slope alpha/kappa sections, preserved clock failures and compensated-clock resumptions, CPU spectra/sections, affine-transient counterexample and frozen-fiber controls. See [the numerical report](../docs/chua-scaling-research.md) and [the structural note](../docs/chua-scaling-structure.md).

## Exploratory searches beyond Chua

These datasets and their independent Python/OpenCL implementations now live in [Chaos Attractor Lab](https://github.com/Griboedow/chaos-attractor-lab). See its [results index](https://github.com/Griboedow/chaos-attractor-lab/blob/main/results/README.md) and [findings overview](https://github.com/Griboedow/chaos-attractor-lab/blob/main/docs/research-findings.md). Chua datasets and all parameter-scaling work remain in this repository. Original historical commits are unchanged.
