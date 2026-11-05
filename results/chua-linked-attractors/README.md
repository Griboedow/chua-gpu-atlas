# Joint topology of four coexisting attractors

[Report](../../docs/chua-linked-attractors.md) · [Interactive 3D viewer](linked-attractors.html) · [Review](review.json)

The four certified cycles F1a, F1b, F2a and F2b are reconstructed at the same exact parameter point. Every one of their 15 nonempty sublinks reduces numerically to an unlink at all tested resolutions and projections. No Borromean or other Brunnian triple was found in this finite configuration.

## Evidence

- `protocol.md`: prospective protocol, published in commit `39a84c8` before the joint Chua classifications.
- `F*-d*-n*.json`: fixed-parameter shooting roots, periods, multipliers, residuals, itinerary and source digest. The suffix gives working digits and subdivisions per affine section-to-section arc.
- Matching `.npz`: physical coordinates `xyz` and sampling times `t`. The repeated closing endpoint is omitted when counting polygon vertices.
- `diagrams/`: 105 complete extracted planar diagrams, component bookkeeping, linking matrices, simplification stages, reduced diagrams and polynomial values.
- `analysis-d60-n96.json`, `analysis-d80-n192.json`, `analysis-d90-n384.json`: all 15 sublinks at each level and six independently calculated solid-angle Gauss integrals.
- `controls.json`, `standard-link-curves.npz`: 34 analytic geometric controls and Hopf orientation/reflection checks.
- `control-development.json`: retained control failures and their corrections.
- `review.json`: source and protocol integrity, root-box containment, refinement agreement, eight independent ODE integrations, all diagram replays, control polynomial checks and minimum polygon segment distances.

The interval certificate is an earlier input in `../chua-accessibility/`. It certifies the periodic cycles, not an isotopy between those cycles and the new polygons. All geometry and reduction claims in this directory retain that distinction.

## Reproduction

Run locally from the repository root, using the Python environment on D:. Dependencies are listed in `requirements.txt`, `requirements-precision.txt` and `requirements-knots.txt`.

```powershell
.\.venv\Scripts\python.exe -X utf8 research/study_chua_linked_attractors.py controls
.\.venv\Scripts\python.exe -X utf8 research/study_chua_linked_attractors.py reconstruct --digits 60 --subdivisions 96
.\.venv\Scripts\python.exe -X utf8 research/study_chua_linked_attractors.py reconstruct --digits 80 --subdivisions 192
.\.venv\Scripts\python.exe -X utf8 research/study_chua_linked_attractors.py reconstruct --digits 90 --subdivisions 384
.\.venv\Scripts\python.exe -X utf8 research/study_chua_linked_attractors.py analyze --digits 60 --subdivisions 96
.\.venv\Scripts\python.exe -X utf8 research/study_chua_linked_attractors.py analyze --digits 80 --subdivisions 192
.\.venv\Scripts\python.exe -X utf8 research/study_chua_linked_attractors.py analyze --digits 90 --subdivisions 384 --seeds 0 2 5
.\.venv\Scripts\python.exe -X utf8 research/review_chua_linked_attractors.py
.\.venv\Scripts\python.exe -X utf8 research/plot_chua_linked_attractors.py
```

Existing reconstruction files and individual diagram records are reused. To recompute them independently, use a separate checkout and a separate output directory, preserving the published originals. The review replays all reductions from their unsimplified integer diagrams; it does not replace those inputs. It writes a new review summary, while controls and plotting regenerate their respective outputs.

## Review targets

1. Confirm identical exact parameters across all components. A multiplier-constrained reconstruction that moves β would answer a different question.
2. Check that each refined root corresponds to its existing certified cycle and that independent full-field integration follows the reconstructed curve.
3. Preserve all component coordinates under one common spatial transformation. Independent centering or scaling can change a link.
4. Track components with no crossings during extraction and reduction. An empty crossing list with a lost component is not the intended unlink.
5. Require crossing-free reduction for the unlink conclusion. Neither zero linking number nor a trivial polynomial is sufficient by itself.
6. Reproduce the positive Borromean control and its unlinked proper sublinks. A polynomial match to a named link does not generally identify an unknown link.
7. Inspect the retained resolution-sensitive raw crossing counts, including 138/98/50 for F2a/F2b in projection 2. Do not interpret agreement of final diagrams as an interval isotopy certificate.
8. Keep the finite scope: these are four preselected coexisting cycles, not an exhaustive attractor search or a proof that the entire Chua system has no Brunnian attractors.
