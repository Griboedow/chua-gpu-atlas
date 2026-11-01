# Numerical knot atlas: data and reproduction

The [report](../../docs/chua-knot-atlas.md) classifies 34 reference cycles into five numerical knot types. These files do not certify the ambient isotopy between an exact ODE orbit and its sampled polygon.

| Files | Contents |
|---|---|
| `ID-d60-n96.json/.npz` | 60-digit shooting root; 96 subdivisions per section-to-section arc; xyz and t arrays |
| `ID-d70-n192.json/.npz` | Independent refinement at 70 digits and twice the sampling density, for all 34 entries |
| `ID-d80-n384.json/.npz` | Third resolution and precision for A04, A05 and P85 |
| `diagrams/` | 108 projections: raw and simplified planar-diagram codes, two Alexander calculations, recognition paths, margins and source hashes |
| `review.json` | Family classification, numerical checks, counts and scope |
| `controls.json`, `standard-knot-curves.npz` | 22 analytic-curve controls and standard embeddings used for illustrations |
| `sampling-counterexample.json`, `A05-uniform-*.npz` | Discrepant coarse sampling and the finer uniform-time comparison at the same root |
| `knot-atlas.html` | Self-contained interactive comparison of every reference orbit with a standard embedding |

The NPZ polygons include a repeated endpoint for section-aligned data; the uniform-time controls use implicit closure. The diagram program removes the repeated endpoint before indexing segments. Crossing heights are measured after positive coordinate scaling and a proper rotation; they are diagnostics of a projection, not lower bounds on all nonadjacent orbit separations. `simplified_crossings` is the size of the obtained diagram, not a general crossing-number minimization algorithm.

Polynomial coefficients are stored from lowest to highest degree, after removing powers of t and choosing a positive leading coefficient. `knot_type` preserves handedness relative to the standard Spherogram torus convention. The identifier `unresolved` in the coarse uniform-time counterexample means the limited torus recognizer did not identify that polygon; its degree-12 polynomial already distinguishes it from the converged A05 result.

## Reproduction

Run locally in a separate checkout when preserving the published records matters. Keep the environment and package cache on D: on the author's Windows machine. Python 3.14, Spherogram 2.4.1, SymPy 1.14.0 and Numba 0.68.1 were used. Install the repository requirements, `requirements-precision.txt`, and `requirements-knots.txt`; no SageMath or GitHub Actions is required.

```bash
python research/reconstruct_chua_knot_orbits.py --digits 60 --subdivisions 96
python research/reconstruct_chua_knot_orbits.py --digits 70 --subdivisions 192
python research/reconstruct_chua_knot_orbits.py --digits 80 --subdivisions 384 --ids A04 A05 P85
python research/analyze_chua_knots.py --pattern '*-d60-n96.npz' --seeds 0
python research/analyze_chua_knots.py --pattern '*-d70-n192.npz' --seeds 0 2
python research/analyze_chua_knots.py --pattern '*-d80-n384.npz' --seeds 0 2
python research/check_chua_knot_sampling.py
python research/review_chua_knots.py --replay-diagrams
python research/plot_chua_knots.py
```

The reconstruction and classification programs resume by skipping completed records. To recompute every value, use an empty `results/chua-knot-atlas` directory in the separate checkout. The review command does not skip the diagram replay. Omit `--replay-diagrams` for a faster algebra, recognition-path and record-consistency audit; this does not repeat the large-diagram reductions.

## Review targets

- Confirm that the refined β values and flight itineraries agree with the referenced catalogue entries. β changes by less than 2.1 × 10⁻¹³.
- Compare the two complete sampling resolutions and both fine-grid projections. Check the third resolution for A04, A05 and P85.
- Replay raw-diagram reductions and the saved type-III paths to the standard summands. A matching polynomial alone is insufficient.
- Verify the analytic-curve controls, including reversed traversal and a figure-eight knot.
- Repeat the A05 uniform-time counterexample. Its coarse mismatch must remain visible.
- Treat the missing certified polygon-to-orbit isotopy as an open proof obligation. The previous interval and Lean results do not fill it.

## Provenance

The first 32 starting roots come from [`family-registry.json`](../chua-layer-atlas/family-registry.json). P75 and P85 use the high-precision roots in [`search-results.json`](../chua-layer-prediction/search-results.json). Earlier records and display arrays are unchanged. The [repository manifest](../chua-layer-atlas/manifest.json) includes these new data, scripts, figures and documentation.
