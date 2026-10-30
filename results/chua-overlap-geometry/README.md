# Symmetry-resolved coexistence geometry

Read the [research note](../../docs/chua-overlap-geometry.md) or download the self-contained [visual guide](guide.html). Figures are available as PNG files in `docs/assets/` and as PDF files here, with vector lines/text and rasterized dense color fields. The four selectable SVG sections are bundled in the guide; no network connection is required after download.

This dataset follows the specific P5 and P7 families from `results/chua-slope-families/anchors.json`. It is not an exhaustive attractor count. Primary branches are included only up to the first period doubling of their asymmetric pair.

## Files

| File | Contents |
|---|---|
| `alpha-XX.json`, `beta-XX.json` | Eighteen contact curves, each with 66 corrected records, full cycle/edge data, tangent vectors and residuals |
| `turning-points.json` | Projection extrema of FF, FS, SS and SF, both sections; curvature controls and independent orbit checks |
| `turning-extended.json` | Projection extrema involving first period-doubling edges: FD, DF, SD, DS and DD |
| `predictions.json` | 32 predictions from tip curvature and separately corrected pairs of horizontal-cut intersections |
| `slice-8.84.json`, `slice-8.85.json`, `slice-8.86.json`, `slice-8.865.json` | Six boundary values at each of 101 m1 positions per section; four stable interior checks per position; corrected contact roots |
| `witnesses.json`, `witness-traces.npz` | Three common-parameter examples supporting at least 2, 3 and 4 attracting cycles; section coordinates, full-field checks and trajectories |
| `independent-checks.json` | 180 full-field cycle controls; 90 surface-transversality and Jacobian-rank controls |
| `branch-controls.json` | 2,340 interior samples on 36 selected stable branch segments, signed-multiplier monotonicity and 36 independent midpoint checks |
| `audit.json` | Executable evidence gates and numerical summary |
| `chua-overlap-*.pdf` | Six scientific figures in PDF format |
| `slice-*.svg`, `guide.html` | Four selectable sections and the standalone visual explanation |

The two-letter contact code selects one edge of P5 followed by one edge of P7. F denotes the symmetric-cycle fold, S the symmetry-breaking edge, and D the first period-doubling edge of an asymmetric cycle. For example, DF means `D5 = F7` at the same m0 and alpha/beta, on two distinct trajectories. Tangents store derivatives of `w` with respect to m1. The vector `w` concatenates P5 flight unknowns, P7 flight unknowns, common m0, and the free control (alpha or beta). Each flight vector consists of initial y, initial z and positive region-flight durations.

P5 and P7 symmetric itineraries each have five half-orbit arcs. Their full asymmetric itineraries have nine arcs after joining the two central flights. F and S use half-return targets +1 and -1; D uses a full-return target -1. Physical peak counts account for the difference between half and full flights. The signed half-return target must not be confused with a full-return period-doubling multiplier.

The displayed quantity `1000(m0-C5)` subtracts the local midpoint of the P5 symmetric interval. It does not multiply physical widths by an undisclosed exaggeration factor. The PDF/PNG color maps interpolate the 101 corrected boundary samples to a finer display mesh; those interpolated cells are not independent solves.

## Reproduce

Use Python with the repository's NumPy, SciPy and Matplotlib dependencies. No GPU is required for this continuation study. From the repository root:

```powershell
python research/resolve_chua_overlap_geometry.py --stage contacts
python research/resolve_chua_overlap_geometry.py --stage contacts --extended
python research/resolve_chua_overlap_geometry.py --stage turning
python research/resolve_chua_overlap_geometry.py --stage turning --extended
python research/resolve_chua_overlap_geometry.py --stage predictions
python research/check_chua_overlap_geometry.py --stage slices
python research/check_chua_overlap_geometry.py --stage witnesses
python research/check_chua_overlap_geometry.py --stage controls
python research/check_chua_overlap_geometry.py --stage branches
python research/check_chua_overlap_geometry.py --stage audit
python research/plot_chua_overlap_geometry.py
python research/manifest_refined_chua_atlas.py --check
```

Contact and slice stages skip existing result files. For a fresh calculation, pass the same new `--output` directory to every numerical stage, followed by `--source` and a separate `--assets` directory for plotting. Input anchors remain the published slope-family records. The default contact mesh contains 65 uniformly spaced m1 values plus the classical slope -5/7. The audit expects the published four detailed slices and 18 contact curves.

The offline guide has a separate script/data check:

```powershell
node research/check_chua_overlap_guide.cjs
```

The [combined manifest](../chua-layer-atlas/manifest.json) includes this dataset and its source dependencies. Numerical thresholds and hashes support reproducibility; they are not interval certificates or a proof of novelty.
