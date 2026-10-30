# Refined Chua layer atlas: evidence and reproduction

Read the [research report](../../docs/chua-refined-atlas.md) for interpretation. The catalogue contains 32 numerically distinguished family identities, not a proven exhaustive number of periodic windows. All calculations use the continuous piecewise-linear model in this repository with m1 = -5/7.

## Main files

| Files | Contents |
|---|---|
| `atlas.html` | Offline interactive atlas; download and open in a browser |
| `family-registry.json` | 32 normalized roots, identities, source candidates, and ancestry |
| `family-catalog.csv`, `family-behavior.json`, `family-orbits.npz` | Periods, maxima, symmetry, state statistics, contraction, sampled orbits |
| `*-sheet.json` | Continued representatives and available fold/symmetry boundaries; missing samples stay explicit |
| `display-samples.npz`, `plot-summary.json` | Arrays and counts used by the figures; a missing value does not mean no family exists |
| `atlas-audit.json` | Independent matrix-exponential and full-field variational checks |
| `fast-flight-checks.json`, `event-map-checks.json` | Comparisons of numerical implementations |
| `interpolation-checks.json` | Fresh midpoint checks for S01–S11 and P55 only |
| `extended-adding.json` | New 55- and 65-maximum cycles, boundaries, independent checks, failed and refined Radau trials |
| `offcenter-raw.npz`, `offcenter-summary.json` | Compact four-transect GPU evidence; the large peak matrix is omitted |
| `offcenter-candidates.json`, `offcenter-cycles.json`, `offcenter-ancestry.json` | Candidate selection, accepted and unresolved CPU checks, identity resolution |
| `off21-*.json` | Local parent recovery, small-step transport, and inversion-symmetry checks for the difficult E09 source |
| `four-cycles.json`, `four-cycles.npz` | Four attracting cycles at exactly the same parameters |
| `overlap-basins*`, `basin-capture*`, `basin-perturbations*`, `basin-transients*` | Initial-condition outcomes, step comparisons, capture and sensitivity diagnostics |
| `*-initial-continuation.json`, `*-daughter.json`, `P65-boundary-attempt.json`, `fast-flight-initial-checks.json` | Retained preliminary/rejected attempts; these are not final displayed sheets |
| `raw-provenance.json`, `manifest.json` | Original raw hashes, compacting rules, versions, and final published checksums |

Arrays in `offcenter-raw.npz` retain the name used by the scanner so the verification scripts can read the compact file directly. It is not byte-identical to the original raw output. `raw-provenance.json` records the original hash and retained arrays. Sampled arc geometry is removed from published JSON; flight times, region itineraries, parameters and checks are retained. Regenerate an arc with `Flights.record(beta, q, True)` using the stored model and itinerary. A node-specific `regions`/`targets` pair takes precedence over the sheet-level template.

E09's initial continuation is excluded from the displayed atlas because local branch continuity was not established. Its rebuilt samples use the automatic event map and bounded small-step affine corrections, with five paths cross-checked between the two methods. E13 remains partially sampled. Missing points, including any stopped by a per-node continuation time budget, must not be interpreted as dynamical boundaries. P65 has a complete representative grid but no claimed complete boundary grid. `P65-boundary-attempt.json` records that unsuccessful boundary attempt.

## Recheck the published evidence

Use the repository Python environment with `requirements.txt` installed. Copy this directory to a scratch location before rerunning scripts that replace their outputs. The examples below use PowerShell and an arbitrary writable output directory:

```powershell
python research/manifest_refined_chua_atlas.py --check
Copy-Item -Recurse results/chua-layer-atlas D:/tmp/chua-atlas-recheck
python research/select_chua_offcenter.py --output D:/tmp/chua-atlas-recheck --check
python research/check_fast_flights.py --output D:/tmp/chua-atlas-recheck
python research/check_chua_event_map.py --output D:/tmp/chua-atlas-recheck
python research/audit_refined_chua_atlas.py --output D:/tmp/chua-atlas-recheck
python research/plot_refined_chua_atlas.py --output D:/tmp/chua-atlas-recheck --assets D:/tmp/chua-atlas-figures
```

The audit checks all saved sheet points using the original 4 by 4 matrix-exponential affine flow. It separately integrates the full state and variational equations for two in-box representatives per family. Gates are closure below 2e-8 for the original affine implementation; multiplier differences below 2e-5; sampled region violations below 1e-7; independent full-field closure below 2e-7; and independent attraction. Region sampling and ordinary floating-point solvers do not provide interval enclosures. Passing these gates does not establish exhaustive enumeration or global family identity.

## Regenerate discovery and continuation

The preceding central-census evidence in `results/chua-layer-census` and the first-family continuation in `results/chua-layers` are inputs. Their own READMEs explain regeneration. For a fresh refinement directory:

```powershell
python research/refine_chua_sheets.py --output D:/tmp/chua-atlas-new
python research/extend_chua_adding.py --output D:/tmp/chua-atlas-new
python research/refine_chua_adding_check.py --output D:/tmp/chua-atlas-new
python research/refine_chua_satellites.py --output D:/tmp/chua-atlas-new
python research/refine_chua_satellites.py --output D:/tmp/chua-atlas-new --adding
python research/chua_window_ancestry.py --output D:/tmp/chua-atlas-new

python research/scan_chua_offcenter.py --output D:/tmp/chua-atlas-new
python research/select_chua_offcenter.py --output D:/tmp/chua-atlas-new
python research/verify_chua_offcenter.py --output D:/tmp/chua-atlas-new
python research/chua_window_ancestry.py --output D:/tmp/chua-atlas-new --local-offcenter-index 21
python research/bridge_chua_itineraries.py --output D:/tmp/chua-atlas-new --recover-off21
python research/chua_window_ancestry.py --output D:/tmp/chua-atlas-new --offcenter
python research/build_chua_family_registry.py --output D:/tmp/chua-atlas-new

Copy-Item D:/tmp/chua-atlas-new/A03-sheet.json D:/tmp/chua-atlas-new/A03-daughter.json
Copy-Item D:/tmp/chua-atlas-new/A05-sheet.json D:/tmp/chua-atlas-new/A05-daughter.json
python research/refine_chua_satellites.py --output D:/tmp/chua-atlas-new --registry D:/tmp/chua-atlas-new/parent-templates.json
python research/refine_chua_satellites.py --output D:/tmp/chua-atlas-new --registry D:/tmp/chua-atlas-new/additional-templates.json
python research/refine_chua_sheets.py --output D:/tmp/chua-atlas-new --adding --families P55
Copy-Item D:/tmp/chua-atlas-new/E09-sheet.json D:/tmp/chua-atlas-new/E09-initial-continuation.json
Copy-Item D:/tmp/chua-atlas-new/E13-sheet.json D:/tmp/chua-atlas-new/E13-initial-continuation.json
python research/bridge_chua_itineraries.py --output D:/tmp/chua-atlas-new --rebuild --families E09 --seconds-per-node 25
python research/bridge_chua_itineraries.py --output D:/tmp/chua-atlas-new --families E13 --seconds-per-node 25
python research/bridge_chua_itineraries.py --output D:/tmp/chua-atlas-new --families E09 E13 --fast-local --seconds-per-node 25
python research/check_chua_event_map.py --output D:/tmp/chua-atlas-new

python research/check_sheet_interpolation.py --output D:/tmp/chua-atlas-new --families S01 S02 S03 S04 S05 S06 S07 S08 S09 S10 S11 P55
python research/audit_refined_chua_atlas.py --output D:/tmp/chua-atlas-new
```

The off-center scan requires an FP64 OpenCL device. Candidate index 21 refers to the saved selection order; if hardware or numerical changes alter that order, identify the candidate near m0 = -8/7 + 0.006, beta = 14.35328125 before invoking the targeted recovery. The current adding constructor stops when a candidate fails its checks; the recorded run reached 65 and needed the separate tightened Radau replay. It does not promise to verify every requested longer cycle. Per-node time limits can change how many unresolved mesh points remain on different machines. Increase that budget to attempt more points, retaining the earlier result for comparison.

The original E09 rebuild was interrupted after 18 accepted nodes and resumed with a 25-second per-node continuation budget. The remaining 19 holes were then filled by bounded local affine continuation; previous timeout records remain in `prior_bridge_failures`. Each saved node is independently checked by the final audit. The original P65 boundary attempt is retained as evidence of a numerical limitation but is not required to reproduce its verified representative grid. The initial strict absolute Jacobian comparison is also retained; the final fast-flow check reports both absolute discrepancies and scale-aware relative errors for long-flight sensitivities.

To rerun ancestry after adding E-family sheets, pass `--families A01 A02 A03 A04 A05 A06 P55 P65 S01 S02 S03 S04 S05 S06 S07 S08 S09 S10 S11` when rebuilding the central ancestry file; otherwise the new sheets would also be treated as central inputs. For the original daughter ancestry, use the saved daughter sheets or rerun the initial satellite continuation before replacing A03 and A05 with their parent sheets.

## Regenerate the overlap experiment

```powershell
python research/verify_four_chua_cycles.py --output D:/tmp/chua-atlas-new
python research/chua_overlap_basins.py --output D:/tmp/chua-atlas-new
python research/chua_overlap_basins.py --output D:/tmp/chua-atlas-new --tight
python research/chua_basin_capture.py --output D:/tmp/chua-atlas-new
python research/chua_basin_capture.py --output D:/tmp/chua-atlas-new --full
python research/chua_basin_capture.py --output D:/tmp/chua-atlas-new --perturb
python research/chua_basin_capture.py --output D:/tmp/chua-atlas-new --transients
```

The basin and transient scripts require the GPU. Step refinement does not converge individual labels on this grid; retain both the images and their disagreement rates. The eight-return capture criterion includes a confirmation delay. It is not a direct measurement of chaotic-transient lifetimes.

## Export and render

```powershell
python research/export_refined_chua_atlas.py --source D:/tmp/chua-atlas-new --output D:/tmp/chua-atlas-publication
python research/plot_refined_chua_atlas.py --output D:/tmp/chua-atlas-publication --assets D:/tmp/chua-atlas-figures
```

The plotting script rebuilds the CSV, display arrays, three PNG/PDF scientific figures, and self-contained interactive HTML from the compact evidence. It uses the earlier original-atlas regular mask for comparison. Numeric results and wall-clock-limited coverage can differ slightly by platform; compare tolerances and branch identities, not binary file equality.
