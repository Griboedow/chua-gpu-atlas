# Three-dimensional slope-volume evidence

Read the [report](../../docs/chua-slope-volumes.md). The two volumes vary m0 and m1 together with alpha or beta, keeping the other parameter fixed at its classical value. Each has shape **(41, 61, 61)** in `(control, m1, m0)` order; m0 varies fastest.

| Files | Contents |
|---|---|
| `alpha.npz`, `alpha.json` | Alpha = 7 to 11, fixed beta = 100/7; per-start diagnostics, final states, class counts and settings |
| `beta.npz`, `beta.json` | Beta = 100/7 +/- 4, fixed alpha = 9; the same fields |
| `*-refinement.npz`, `refinement.json` | 512 stratified nodes per volume rerun at dt/2 and doubled transient/observation |
| `cpu-controls.json`, `low-alpha-recheck.json` | Independent full-field controls, including the unresolved weak low-alpha chaos candidate and its longer matched-start checks |
| `render-summary.json`, `animation.json` | Mesh sizes, exact crop bounds, plotted context counts and slow-camera settings |
| `volumes.html` | Self-contained interactive viewer; download and open in a browser |
| `viewer-summary.json`, `viewer-checks.json` | Trace mapping and controller checks; browser rendering was not available for verification |

The NPZ files retain axes, aggregate labels, three per-start classes/statuses, mean and four-window Lyapunov estimates, outer-region visitation flags and final states. Raw accumulation statistics are omitted. Reconstruct parameter tuples with `parameters(m0, m1, control, kind)` from `research/scan_chua_slope_volumes.py`; the order is `[alpha, beta, m0, m1]`.

Displayed classes follow the earlier slope-plane experiment. Uniform regular or chaotic labels require agreement of all three starts. Class 6 denotes different classifications across starts, not automatically distinct attractors. Escape uses the operational coordinate cutoff 100. The central control index is 20. Its 61-square classes and finite exponents match the earlier wide 2D scan sampled with stride four.

The baseline uses FP64, dt = 0.005, transient = 600 and observation = 1,600. The three starts are `(0.7,0,0)`, `(-0.2,0.1,0)` and `(1.5,0.02,-1.5)`. The stratified refinement chooses up to 48 nodes per displayed category, then fills to 512 using RNG seed 20261030. Do not interpret this selection's agreement percentages as whole-volume accuracy.

## Reproduce

Run from the repository root using the existing Python dependencies and an FP64 OpenCL GPU. Scan into a fresh folder:

```powershell
python research/scan_chua_slope_volumes.py --output D:/tmp/chua-slope-volumes
python research/check_chua_slope_volumes.py --output D:/tmp/chua-slope-volumes --stage refine
python research/check_chua_slope_volumes.py --output D:/tmp/chua-slope-volumes --stage cpu
python research/check_chua_slope_volumes.py --output D:/tmp/chua-slope-volumes --stage low-alpha
python research/plot_chua_slope_volumes.py --source D:/tmp/chua-slope-volumes --assets D:/tmp/chua-slope-volume-figures
python research/build_chua_slope_volume_viewer.py --source D:/tmp/chua-slope-volumes
node research/check_chua_slope_volume_viewer.cjs D:/tmp/chua-slope-volumes
```

The scan reads the earlier `results/chua-slope-plane/wide.npz` to use exactly matching slope coordinates. `--kind alpha` or `--kind beta` selects one volume; `--resume` skips existing completed volumes. Rendering can use `--stage static` or `--stage animation`. No GPU integration is needed for plotting or HTML export. Node.js is needed only for the optional controller check, which uses state stubs and does not test browser/WebGL rendering.

Marching cubes uses every node of the binary mask at level 0.5. The interior view selects m0 <= -1.06 and -0.90 <= m1 <= -0.20; its exact retained ranges are in `render-summary.json`. There is no smoothing or artificial closure at cropped boundaries. Finite resolution can create missing or fragmented layers. The previous 32-family catalogue is not a count for either new volume.

The combined [atlas manifest](../chua-layer-atlas/manifest.json) covers these data, scripts, report and assets. Verify with `python research/manifest_refined_chua_atlas.py --check`.
