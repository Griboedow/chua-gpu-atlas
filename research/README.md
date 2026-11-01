# Research programs

Run programs from the repository root. The [research guide](../docs/research-guide.md) connects their scientific questions to reports and datasets.

## Cross-parameter prediction, impulses and hidden coexistence

`predict_chua_slope_sheets.py` separates calibration, frozen forecasts and endpoint searches. `study_chua_trigger.py` and `certify_chua_trigger.py` compute impulse crossings and verify their relation to the exact basin separator. `chua_interval_lobe.py` provides the second return geometry; `certify_chua_hidden_tangency.py` and `certify_chua_hidden_neighborhoods.py` verify the hypotheses for infinite hidden coexistence. `check_chua_followup.py` and `check_chua_hidden_tangency.py` use separate numerical implementations. The shared composition supports rational moving frames to reduce interval wrapping. [Reports and exact reproduction commands](../results/chua-followup/README.md).

## Local basin and coexistence results

`certify_chua_basin_separator.py` covers both saddle exit strips, with fixed-time flow enclosures at diode contacts. `check_chua_basin_proof.py` checks event derivatives, a rejected root and independent flow integrations. `audit_chua_local_results.py` checks coverage and the finite premises of both theorem arguments; `plot_chua_local_results.py` draws the figures. [Report](../docs/chua-basin-and-coexistence.md) · [Reproduction](../results/chua-basin-proof/README.md).

## Periodic-attractor accessibility

The follow-up [global-stability and basin study](../docs/chua-hidden-boundary.md) uses `derive_chua_global_sector.py`, `study_chua_hidden_boundary.py`, `certify_chua_basin_saddle.py`, `check_chua_hidden_boundary.py`, `check_chua_storage_formal.py` and `plot_chua_hidden_boundary.py`. These derive the storage identities, follow both cycles, certify saddle/sink roots, audit numerical evidence, check three Lean implications and draw the figures. [Commands and records](../results/chua-hidden-boundary/README.md).

- `certify_chua_coexisting_cycles.py` validates four first-return sinks and distinguishes their orbits.
- `certify_chua_hidden_cycle.py` encloses a periodic attractor and proves exclusion neighborhoods at all three equilibria.
- `study_chua_stability.py` validates another hidden cycle, a common Lyapunov sector and local-stability predictions.
- `check_chua_accessibility.py` independently checks all six cycles and repeated-precision enclosures.
- `prepare_chua_accessibility.py`, `probe_chua_accessibility.py` and `check_chua_stability_prediction.py` reproduce the exploratory and finite-time numerical experiments.
- `study_chua_circuit_constraints.py` constructs positive component values and tests nominal and perturbed lossy models.
- `plot_chua_accessibility.py` produces the three figures from these records.

[Report](../docs/chua-accessibility.md) · [Evidence and commands](../results/chua-accessibility/README.md)

## Interval proof

| Program | Purpose |
|---|---|
| `certify_chua_homoclinic_tangency.py` | Main interval certificate and shifted-box controls |
| `chua_interval_map.py` | Affine event maps, implicit first/second derivatives, all-time path checks |
| `chua_interval_connection.py` | Enclosed composition in coordinates following the saddle |
| `certify_chua_saddle.py` | Local saddle, parameter derivatives and coordinate bounds |
| `locate_chua_homoclinic_tangency.py` | Reconstruct the numerical center used by the certificate |
| `check_chua_interval_derivatives.py` | Independent local derivative checks |
| `check_chua_global_connection.py` | Independent 27-fold composition checks |
| `plot_chua_infinity.py` | Figures from the recorded intervals and measured windows |

[Proof and derivation](../docs/chua-infinity-certificate.md) · [Commands and evidence](../results/chua-layer-infinity/README.md)

## Coexistence audit and formal proof

- `audit_chua_coexistence.py` recomputes the local tangency hypotheses at two precisions and rejects two shifted boxes.
- `check_chua_formal.py` builds the [Lean project](../formal/README.md), checks the five declarations' axiom dependencies and runs rejection controls.
- `check_chua_formal_local.ps1` runs that check on the Windows research workstation using tools, caches and temporary files on D:.
- `plot_chua_proof_audit.py` draws the checked implications and remaining obligations from those records.

[Audit report](../docs/chua-coexistence-audit.md) · [Evidence and reproduction](../results/chua-proof-audit/README.md)

## Numerical studies

- `refine_chua_slope_atlases.py`, `check_chua_slope_refinement.py`, `study_chua_slope_grazing.py`, `check_chua_slope_branches.py` and `plot_chua_slope_refinement.py` reproduce the [refined slope volumes and contact study](../results/chua-slope-refinement/README.md).

- `predict_chua_layers.py`, `verify_chua_layer_predictions.py` and the associated precision, audit and plotting programs implement the [forecast protocol](../results/chua-layer-prediction/README.md).
- The `refine_chua_*`, registry, gallery and atlas programs implement the [family refinement](../docs/chua-refined-atlas.md).
- The `scan_chua_slope_*`, `continue_chua_slope_*` and associated checks resolve [slope sections and volumes](../docs/chua-slope-volumes.md).
- `resolve_chua_overlap_geometry.py` and its checks implement the [symmetry and coexistence study](../results/chua-overlap-geometry/README.md).
- `manifest_refined_chua_atlas.py --check` verifies the published file hashes.

Older programs remain alongside their data so the earlier numerical stages can be reproduced. Exact per-study commands are given in the reports and results guides; the interval proof does not require rerunning the GPU surveys.

The field-name consistency check, `python research/check_chua_schema.py`, compares the current records with their preceding Git publication. It verifies exact numerical values, array contents, mathematical blocks and local links, and writes its report under `.publication/`.
