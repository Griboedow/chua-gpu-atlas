# Denser layers for the README

This add-on makes the existing 32-family catalogue easier to see. It does not discover another 32 families, establish an exhaustive count, or replace the earlier 9 by 9 atlas grids. Read the [main report](../../../docs/chua-refined-atlas.md#sampling-and-visualization-of-individual-families) for interpretation.

## What the figures show

The [animation](../../../docs/assets/chua-readme-layers.gif) highlights all 32 families once, with an overview before and after. The camera, coordinates and alpha=9 slice stay fixed. The [static overview](../../../docs/assets/chua-readme-layers-poster.png) and [width zoom](../../../docs/assets/chua-readme-layer-zoom.png) are available separately; [the zoom PDF](readme-layer-zoom.pdf) retains vector lines and text.

Dots are numerically closed attracting cycles. Colored facets locate a representative surface at a selected transverse multiplier; they do not measure the full window's physical thickness. Triangles connect adjacent accepted mesh vertices and are geometrically clipped at beta = 14.24 and 14.40. Neither the facets nor straight slice segments certify attraction between their vertices.

The width zoom uses a categorical row for each of six period-adding families, then physical beta offsets for P55 and P65. Its open circles mark the earlier boundary estimates. The 25 newly solved interior points per core run from 2% to 98% of that interval; endpoints are excluded. No vertical separation of categorical rows represents a parameter distance.

## Data and numerical checks

| File | Meaning |
|---|---|
| `A01.json` through `S11.json` (32 family files) | Selected 17 by 17 mesh locations, closed representatives, itineraries and independent affine checks |
| `audit.json` | Final surface counts, coverage and numerical gates; all selected points recovered |
| `thin-cores.json` | 25 fresh interior cycles each for P55 and P65, plus full-field midpoint controls |
| `P65-positive-half-multiplier-attempt.json` | Retained preliminary display attempt, with 35 unresolved locations; excluded from the figures and final totals |
| `render-summary.json` | Frame order, face counts and plotted sample totals |
| `readme-layer-zoom.pdf` | Vector version of the multiscale figure |

There are **5,447 checked surface samples**, including **3,975 new parameter locations** relative to the old mesh; **4,126** fall inside the displayed beta box. The total includes 1,321 halo samples used to clip adjacent facets. This set overlaps the original 2,587 representatives and must not be added to it as a distinct-sample count. There are a further **50 thin-core profile samples** at alpha = 9 and m0 = -8/7.

A linear interpolant of the old beta samples selects candidate cells intersecting beta = [14.238, 14.402]; it does not generate the accepted orbit data. All four cell vertices are selected. Each candidate is closed by bounded local continuation from a nearby earlier representative. Missing source coverage is not extrapolated. Selection near the visible box is not an exhaustive search of the complete 17 by 17 mesh or of all possible periodic orbits.

The representative target is half-return multiplier +0.5 for symmetric families and full-return multiplier -0.5 for asymmetric families, except P65. P65 retains its existing interior half-return multiplier **-0.48332857085597425**. The initial +0.5 display attempt closed 169 of 204 selected locations, leaving 35 numerical corrector failures. Keeping the existing interior target recovers all 204. P65's catalogue identity remains the normalized +0.5 root; changing its display representative does not create or merge a family.

Every accepted point passes the original 4 by 4 matrix-exponential implementation: closure below 2e-8, return-multiplier magnitude discrepancy below 2e-5, and attracting return magnitude below 0.999. Region sampling is gated at 1e-7. The largest observed surface closure error is **2.09e-10**, and the largest multiplier discrepancy is **2.17e-7**. These checks compare two affine implementations, not independent physical models.

Two new in-box points each for A05, S08, E09, E13 and P65 also receive independent full-field state and variational integration. All ten controls pass; maximum closure is **1.69e-10** and maximum multiplier discrepancy **6.13e-7**. Each thin-core midpoint receives the same independent check; their largest closure is **5.08e-10** and multiplier discrepancy **8.16e-7**. The full-field gates are closure below 2e-7, multiplier error below 2e-5 and attraction. These are numerical checks, not interval certificates.

The P55 and P65 profile spacings are approximately **8.9653e-9** and **1.2704e-9** in beta. Their widths, approximately 2.2413e-7 and 3.1759e-8, come from the previously published fold and symmetry-breaking estimates in `../extended-adding.json`. Descendant branches can extend beyond the symmetric core.

## Reproduce

Run from the repository root with its Python dependencies installed. These scripts read the earlier published family sheets and registry. Use a fresh scratch directory so existing evidence remains available:

```powershell
python research/refine_readme_layers.py --output D:/tmp/chua-readme-view
python research/refine_readme_thin_cores.py --output D:/tmp/chua-readme-view
python research/plot_readme_layers.py --source D:/tmp/chua-readme-view --assets D:/tmp/chua-readme-figures
```

To resume completed families, add `--resume` to the first command. To recompute individual families, use `--families E09 P65`, then run the first command with `--resume` and without `--families` to rebuild the aggregate audit. Resume trusts existing family summaries; use a fresh directory for an independent recomputation. Replaying the historical P65 attempt is optional and writes a separate result:

```powershell
python research/refine_readme_layers.py --output D:/tmp/chua-readme-p65-trial --families P65 --p65-half-multiplier 0.5
```

Counts in a replay can differ with numerical-library versions, especially in the retained preliminary attempt. Unresolved locations remain explicit. The parent directory's [manifest](../manifest.json) hashes this entire add-on, the new scripts, README/report text, and rendered assets. Run `python research/manifest_refined_chua_atlas.py --check` to verify the published files.
