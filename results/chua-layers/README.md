# Classical periodic-layer continuation

See [the research note](../../docs/chua-layer-bifurcations.md) for interpretation,
parameter conventions, reproduction commands, and limitations.

- `initial.json`: exploratory full-field cycle used to initialize affine shooting.
- `affine-central.json`: symmetric branch and its fold/symmetry thresholds.
- `extended.json`: asymmetric branch, fold partners, fold nondegeneracy controls,
  and 49 points on each of three local bifurcation surfaces.
- `cascade.json`: doubled cycle and its next multiplier crossing.
- `verification.json`: all independent state replays and an asymmetric partner.
- `verified-orbits.npz`: time and state arrays from selected Radau replays.
- `audit.json`: executable gates, independent variational checks, sensitivities,
  and source hashes. Read the measured errors, not only the final status.
- `plot-data.json`: additional saddle continuation and symmetry-partner sections.
- `bifurcations.pdf`, `surfaces.pdf`: standalone figures; PNG versions are in
  `docs/assets`.
- `manifest.json`: provenance, dependency versions, and artifact hashes.

All values are floating-point numerical results. This is a local study of one
window, not an interval proof, exhaustive basin survey, or a complete account
of the atlas's layers. The second period-doubling threshold is located, but its
daughter branch and the final loss of periodic attraction are not continued.
