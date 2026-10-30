# Two-dimensional slope-plane evidence

Read the [report](../../docs/chua-slope-plane.md). This experiment fixes alpha = 9 and beta = 100/7 while varying m0 and m1 independently. It complements the earlier atlas, whose m1 was fixed; no correspondence to all 32 earlier family identities is claimed.

| Files | Contents |
|---|---|
| `wide.npz`, `wide.json` | 241 by 241 broad grid, all three starts, Lyapunov windows, states, statistics, classifications, GPU details and settings |
| `zoom.npz`, `zoom.json` | 241 by 241 grid within 0.006 of each classical slope; spacing 5e-5 |
| `*-refinement.npz`, `refinement.json` | 512 stratified nodes per panel rerun at dt/2 and doubled transient/observation, including all original and refined labels |
| `zoom-peaks.npz`, `zoom-peaks.json` | Peak buffers and shortest detected repeating patterns for seed zero at the 9,158 uniformly regular zoom nodes |
| `cpu-controls.json` | Six full-field DOP853 checks started from saved late GPU states |
| `closed-cycles.json` | Independently closed and checked five- and seven-peak cycles |
| `plot-summary.json` | Displayed domains and category counts |

Parameter rows are `[alpha, beta, m0, m1]`. Flattening follows an `xy` mesh: m0 varies fastest, and maps reshape to `(len(m1), len(m0))`. Per-start class codes 0–5 follow `chua_gpu.CLASS_NAMES`. The displayed `labels` additionally use 6 for disagreements across starts. Uniform displayed classes require agreement of all three starts; any numerical failure would override them. There were no numerical-failure codes in these runs.

The broad grid uses dt = 0.005, transient = 600 and observation = 1,600. The zoom uses dt = 0.005, transient = 1,000 and observation = 2,000. Three original starts are `(0.7,0,0)`, `(-0.2,0.1,0)` and `(1.5,0.02,-1.5)`. Escape means a coordinate exceeds 100; it is not a proof of unboundedness. The parent [methods guide](../../docs/methods.md#numerical-method-and-interpretation) defines the finite-time Lyapunov classes.

The refinement sample selects up to 64 nodes from each displayed category, then fills uniformly to 512 using RNG seed 20261030. It deliberately overrepresents uncertain and mixed points. Do not interpret its agreement percentages as whole-grid accuracy. The published main images retain the original complete grids; refined samples are separate evidence.

The peak diagnostic continues seed zero for 1,200 additional time units at dt = 0.0025, retains up to 512 peaks, and examines the last 384 for at least four repeats. Its maximum lag is 96, with relative height tolerance 5e-4. Zero means unresolved. These are recurrence diagnostics; rare large counts are not automatically new closed cycles or new family identities.

## Reproduce

Use the repository environment and an FP64 OpenCL GPU. Run from the repository root into a fresh directory:

```powershell
python research/scan_chua_slope_plane.py --output D:/tmp/chua-slope-plane
python research/check_chua_slope_plane.py --output D:/tmp/chua-slope-plane --stage refine
python research/check_chua_slope_plane.py --output D:/tmp/chua-slope-plane --stage peaks
python research/check_chua_slope_plane.py --output D:/tmp/chua-slope-plane --stage cpu
python research/check_chua_slope_plane.py --output D:/tmp/chua-slope-plane --stage cycles
python research/plot_chua_slope_plane.py --source D:/tmp/chua-slope-plane --assets D:/tmp/chua-slope-figures
```

The CPU and cycle stages do not require a GPU calculation. Plotting only reads saved arrays. `--kind wide` or `--kind zoom` runs one grid. The combined [atlas manifest](../chua-layer-atlas/manifest.json) includes these data, scripts, report and images; verify it with `python research/manifest_refined_chua_atlas.py --check`.
