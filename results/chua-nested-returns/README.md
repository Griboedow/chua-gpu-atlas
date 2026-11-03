# Nested return experiment

Six new primitive attracting cycle representatives, with 105, 175 and 245 peaks near each of two 35-peak parent families. The [report](../../docs/chua-nested-returns.md) distinguishes one additional numerically supported level from an arbitrary-depth conjecture.

The protocol was published in `bd7429b`; calibration and forecasts were frozen in `eaa2eca`. Preserve those records. The forecast digest is `d737a5cd38a20e90b47aa50ab22192fbeb0cfb1b211510678976407035b81f21`. Location forecasts pass 2/6, widths 6/6 and periods 6/6; numerical orbit acceptance is separate from forecast success.

## Records

| Files | Contents |
|---|---|
| `protocol.md`, `calibration.json`, `forecasts-frozen.json` | Pre-search selection, old parent roots and six unmodified forecasts |
| `*-reference.json`, `*-edges.json` | 60-digit multiple-shooting representatives and fold/symmetry endpoints |
| `*-precision.json` | 85-digit repetitions of all 18 roots |
| `*-exponential.json`, `*-orbit.npz` | Independent matrix-exponential controls and plot trajectories |
| `*-audit.json` | Full-field peak counts, whole-period drifts, independent section variational equations and numerical acceptance gates |
| `*-geometry.json`, `geometry-summary.json` | Same-parameter long and short parents, seven-point block comparisons and coarse words |
| `negative-controls.json` | Ordinary 105-peak comparison and intentionally repeated nonprimitive parents |
| `parent-branch-continuation.json` | Attracting-to-saddle connection through the parent fold |
| `parent-block-map.json` | Fixed-coordinate return samples and independent exponential checks |
| `moving-coordinate-map.json` | Explicitly post-search secant normalization; not a replacement forecast |
| `scores.json`, `review.json` | All forecast errors and recomputed record checks |
| `*.pdf` | Vector versions of the four report figures |

The `block-map-control.json` record checks the new return implementation's Jacobian against a central difference at the old S06 root. The independent orbit audits provide broader controls. Region membership is sampled, not enclosed by interval arithmetic.

## Reproduction

Use Python with NumPy, SciPy, mpmath and Matplotlib, as installed for the preceding precision studies. Run locally; no GPU or GitHub Actions are required. Keep the checkout, virtual environment and Matplotlib cache on the selected work drive.

The record review does not rerun integrations:

```bash
python research/review_chua_nested_returns.py
python research/manifest_refined_chua_atlas.py --check
```

To replay expensive stages, use a separate checkout and move generated target files out of this result directory while retaining `protocol.md`, `calibration.json` and `forecasts-frozen.json`. Search, audit and geometry stages retain existing results. The forecast generator deliberately refuses to overwrite a frozen forecast.

```bash
python research/test_chua_nested_returns.py search --digits 60
python research/audit_chua_nested_returns.py
python research/analyze_chua_nested_returns.py
python research/analyze_chua_nested_returns.py --branches
python research/sample_chua_nested_returns.py
python research/sample_chua_nested_returns.py --moving
python research/review_chua_nested_returns.py
python research/plot_chua_nested_returns.py
```

`--ids S06-k3 S06-k5 S06-k7` selects a subset for search or audit. The independent audit integrates entire orbits and individual section transitions; the latter errors must not be reported as whole-period drift. The plotting program places its cache in `.publication/matplotlib` and writes PNGs under `docs/assets` and PDFs here.

## Review targets

1. Compare the frozen records and source digests with the pre-search commit. Recompute location errors in **measured** widths, including all four misses. No searched child may enter the forecast calibration.
2. Reconstruct representative and endpoint roots independently. Check positive flight times, diode-region membership, closure, section crossings and multipliers. A repeated itinerary alone is insufficient.
3. Reject multiple traversal: actual normalized section points must be distinct. The intentionally repeated-parent control demonstrates why an apparently correct peak count and multiplier do not establish primitiveness.
4. Verify the long saddle at the child's identical β, distinguish it from the other continued long branch, and compare all seven points in each block. The binary word is only a coarse diagnostic.
5. Recompute period increments, increasing close residence and the positive long-parent multipliers. Check the fold connection rather than assuming an old attracting reference is already a saddle.
6. Preserve the poor fixed-coordinate quadratic approximation, off-secant components and whole-period integration drift. The moving-coordinate comparison is calibrated after the searches.
7. Do not infer arbitrary recursive depth, a new infinity theorem, global layer completeness, basin size or knot type from these numerical records.

The review script checks record consistency and the numerical controls; a `passed` field is not a proof object. Decimal reproduction and independent implementation agreement do not replace interval existence and uniqueness bounds.
