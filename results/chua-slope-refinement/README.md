# Refined slope-volume evidence

The [report](../../docs/chua-slope-refinement.md) distinguishes periodic-cycle birth, symmetry breaking and an internal contact with the diode breakpoint. Two 41 × 41 meshes per family contain 6,713 accepted cross-sections. Eleven P5 nodes remain unresolved; both P7 meshes are complete within the requested box.

## Reproduction

Run from the repository root using the existing NumPy, SciPy and Matplotlib environment. All computations run locally on the CPU. The recorded environment used Python 3.14.3, NumPy 2.5.4 and SciPy 1.18.1. Set `OPENBLAS_NUM_THREADS=1` for small-matrix continuation.

```bash
python research/refine_chua_slope_atlases.py --family P5 --kind alpha
python research/refine_chua_slope_atlases.py --family P5 --kind beta
python research/refine_chua_slope_atlases.py --family P7 --kind alpha
python research/refine_chua_slope_atlases.py --family P7 --kind beta
python research/check_chua_slope_refinement.py
python research/study_chua_slope_grazing.py
python research/check_chua_slope_branches.py
python research/plot_chua_slope_refinement.py
python research/manifest_refined_chua_atlas.py --check
```

The four mesh jobs are independent. They reuse and recheck the previous 21 × 21 accepted roots from `../chua-slope-families/`. Rerunning a study changes its elapsed-time metadata; the manifest checks the published snapshot and will then require an intentional refresh. The original grazing forecast is preserved and checked for consistency on rerun.

## Files

| File | Meaning |
|---|---|
| `P5-alpha.json`, `P5-beta.json`, `P7-alpha.json`, `P7-beta.json` | Corrected fold, symmetry edge and attracting interior cycle at every accepted node, each with its own itinerary; failed attempts and completeness |
| `audit.json` | Whole-mesh residual/path checks, 84 independent ODE controls, sampled fold coefficients and old-cell interpolation discrepancies |
| `heldout-interpolation.json` | 32 additional corrected cell centres, coarse/fine interpolation errors and full witnesses |
| `bifurcation-anchors.json`, `bifurcation-controls.json` | Twelve branch-switch locations, stable/saddle pairs, asymmetric inversion pairs, unstable parents and unsuccessful initial guesses |
| `grazing-prediction.json` | Preserved initial constant-multiplier protocol; its proposed asymptotic coefficient was not validated on that constraint |
| `grazing-constraint-control.json` | Actual outcomes of the unsuccessful first protocol, including a failed smallest point |
| `grazing-linear-prediction.json` | Fixed smooth path, Taylor coefficients and their step checks; formulated after preliminary exploration, not a blind prediction |
| `grazing-results.json` | Three central contact points, six sampled contact curves, nine cycle witnesses, DOP853/Radau controls and acceptance gates |
| `chua-slope-*.pdf` | Vector versions of the three report figures |

`q` stores the two initial section coordinates followed by affine-flight durations. Mesh records are symmetric half-period solutions; the period is twice the duration sum. `regions` and `targets` belong to that particular root. They must not be replaced by a single family-wide itinerary. Branch-control records explicitly state whether they are symmetric or full-period.

Independent segmented controls integrate the actual piecewise field and variational equation with no state resets. The affine flight times only partition the time axis so a new short central excursion cannot be stepped over. Root closure, sampled extrema and independent integration are numerical evidence, not interval proofs. Neither mesh completeness nor correct file hashes establish an exhaustive attractor count.
