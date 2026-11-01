# Continued slope families and coexistence evidence

Read the [report](../../docs/chua-slope-families.md). Local labels P5 and P7 follow the specific closed five- and seven-peak reference cycles from `results/chua-slope-plane/closed-cycles.json`; they are not a global numbering of all periodic families.

| Files | Contents |
|---|---|
| `anchors.json` | Two central symmetric periodic-orbit representatives and their m0 fold/symmetry edges |
| `P5-alpha.json`, `P5-beta.json`, `P7-alpha.json`, `P7-beta.json` | Continued boundary and interior-orbit records, axes, successful node indices, attempted failures and scope |
| `independent-checks.json` | 60 selected full-field state-and-variational controls |
| `overlaps.json` | 13 common-parameter witnesses, each with two independently checked attracting cycles |
| `coexistence-witness.json`, `coexisting-traces.npz` | Highlighted alpha = 8.88 witness and one full period of each trajectory |
| `bifurcations.json` | Fold nondegeneracy controls, stable/saddle witnesses and symmetry-breaking pairs |
| `itinerary-change.json` | Actual-field continuation of P7 past a fixed-itinerary failure; 8 to 12 breakpoint crossings |
| `basin.npz`, `basin.json` | 121-square finite-time initial-state capture slice and three additional standard starts |
| `basin-refinement.npz`, `basin-time-only.npz`, `basin-time-only.json` | Matched-start step/time sensitivity and a same-step longer-run control |
| `audit.json`, `animation.json` | Consistency checks and calm camera settings |

Each sheet node is a `(control, m1)` cross-section. Its m0 fold and symmetry values are solved by affine-flight continuation, and `representative` is a closed attracting cycle at their midpoint. The control is alpha or beta according to the filename. Missing nodes must remain missing in plots. Affine arcs are checked at sampled points and at numerically located x extrema; derivative roots are bracketed with maximum spacing 0.025. This is not interval-certified root enumeration. The surfaces interpolate stability-interval midpoints only, without thickening thin stability intervals.

The basin arrays have flat node order `(z, y)` with y fastest, followed by the three standard starts. All initial states are retained. Labels 5 and 7 require a repeated ordered peak signature matching the independently closed reference; 0 is unresolved or another outcome. Full peak streams are compacted to the last 28 maxima, with original counts and statuses retained. Refinement files store indices into the baseline file. The detailed basin map is sensitive to integration step, as quantified in the report.

## Reproduce

Run from the repository root with the existing Python dependencies. Only the capture-map stages require an FP64 OpenCL GPU.

```powershell
python research/continue_chua_slope_families.py --output D:/tmp/chua-slope-families --stage anchors
python research/continue_chua_slope_families.py --output D:/tmp/chua-slope-families --stage sheets
python research/check_chua_slope_families.py --output D:/tmp/chua-slope-families --stage overlaps
python research/check_chua_slope_families.py --output D:/tmp/chua-slope-families --stage boundaries
python research/check_chua_slope_families.py --output D:/tmp/chua-slope-families --stage sampled
python research/check_chua_slope_families.py --output D:/tmp/chua-slope-families --stage itinerary
python research/map_chua_slope_coexistence.py --output D:/tmp/chua-slope-families
python research/map_chua_slope_coexistence.py --output D:/tmp/chua-slope-families --stage time-control
python research/check_chua_slope_families.py --output D:/tmp/chua-slope-families --stage audit
python research/plot_chua_slope_families.py --source D:/tmp/chua-slope-families --assets D:/tmp/chua-slope-family-figures
```

The default sheet meshes contain 21 values per axis. `--nodes` changes that experimental mesh; the published audit and highlighted witness use the default indices. Successful complete sheets are skipped on a repeat; incomplete sheets are recalculated. Plotting accepts `--stage static` or `--stage animation`. The GIF has a 24-second cycle, a 40-degree horizontal sweep, fixed 24-degree elevation and zero roll.

The [combined atlas manifest](../chua-layer-atlas/manifest.json) hashes these data, scripts, report and assets. Verify with `python research/manifest_refined_chua_atlas.py --check`.
