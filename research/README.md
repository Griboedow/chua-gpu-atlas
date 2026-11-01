# Research programs

Run programs from the repository root. The [research guide](../docs/research-guide.md) connects their scientific questions to reports and datasets.

## Interval proof

| Program | Purpose |
|---|---|
| `certify_chua_homoclinic_tangency.py` | Main interval certificate and shifted-box controls |
| `chua_interval_map.py` | Affine event maps, implicit first/second derivatives, all-time path checks |
| `chua_interval_connection.py` | Enclosed composition in coordinates following the saddle |
| `certify_chua_layer_organizer.py` | Local saddle, parameter derivatives and coordinate bounds |
| `locate_chua_homoclinic_tangency.py` | Reconstruct the numerical center used by the certificate |
| `check_chua_interval_derivatives.py` | Independent local derivative checks |
| `check_chua_global_connection.py` | Independent 27-fold composition checks |
| `plot_chua_infinity.py` | Figures from the recorded intervals and measured windows |

[Proof and derivation](../docs/chua-infinity-certificate.md) · [Commands and evidence](../results/chua-layer-infinity/README.md)

## Coexistence audit and formal proof

- `audit_chua_coexistence.py` recomputes the local tangency hypotheses at two precisions and rejects two shifted boxes.
- `check_chua_formal.py` builds the [Lean project](../formal/README.md), checks the five declarations' axiom dependencies and runs rejection controls.
- `plot_chua_proof_audit.py` draws the checked implications and remaining obligations from those records.

[Audit report](../docs/chua-coexistence-audit.md) · [Evidence and reproduction](../results/chua-proof-audit/README.md)

## Numerical studies

- `predict_chua_layers.py`, `verify_chua_layer_predictions.py` and the associated precision, audit and plotting programs implement the [forecast protocol](../results/chua-layer-prediction/README.md).
- The `refine_chua_*`, registry, gallery and atlas programs implement the [family refinement](../docs/chua-refined-atlas.md).
- The `scan_chua_slope_*`, `continue_chua_slope_*` and associated checks resolve [slope sections and volumes](../docs/chua-slope-volumes.md).
- `resolve_chua_overlap_geometry.py` and its checks implement the [symmetry and coexistence study](../results/chua-overlap-geometry/README.md).
- `manifest_refined_chua_atlas.py --check` verifies the published file hashes.

Older programs remain alongside their data so the earlier numerical stages can be reproduced. Exact per-study commands are given in the reports and results guides; the interval proof does not require rerunning the GPU surveys.
