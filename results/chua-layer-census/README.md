# Layer census evidence

See the [research report](../../docs/chua-layer-census.md) for interpretation, limitations, and reproduction commands. These are numerical findings, not validated interval proofs or an exhaustive census of basins or globally connected layers.

| Files | Contents |
|---|---|
| `line.npz`, `tight.npz`, `adding.npz` | Compact GPU scans: parameter axes, classifications, Lyapunov observations, final states, peak counts/status, and candidate recurrence periods. Full peak matrices omitted. |
| `line-summary.json` | Original baseline report with the overly tight 2e−5 recurrence threshold; retained as a diagnostic trial. Do not use its 173 fragments as a layer count. |
| `line-recurrence.json/.npz` | Saved baseline reanalysis at 5e−4. The exported `line.npz` uses these revised periods. |
| `tight-summary.json`, `adding-summary.json` | Scan settings, device/timing information, candidate intervals and nested-resolution counts. The adding scan uses one initial state. |
| `cycles.json/.npz`, `cycles-extra.json/.npz`, `adding-cycles.json/.npz` | Full CPU verification records and sampled independently integrated orbits. Include two non-recurrent candidates, a duplicate nine-maximum example, and both Radau trials of the 45-maximum cycle. |
| `family-edges.json` | Symmetric fold and symmetry-breaking thresholds for selected verified families, including both representatives of the same nine-maximum window. |
| `second-boundaries.json`, `second-surfaces.json` | Three central thresholds and 49 × 3 continued boundary records for the second broad family. |
| `census-analysis.json` | Seventeen deduplicated window representatives, 16 finite-time gap witnesses, the four-member adding sequence, section-return distances, and original-grid comparison. |
| `atlas-comparison.npz` | Original 41³ regular mask and masks for the two interpolated five-maximum cores. |
| `coexistence.json/.npz` | Two attracting cycles at identical parameters, independent state/multiplier checks, perturbed-start checks, and sampled trajectories. |
| `audit.json` | Executable numerical acceptance gates, independent boundary integrations, fold controls, and the second fold's stable/saddle partners. |
| `chua-layer-*.pdf` | Exportable figures; PNG versions are in `docs/assets`. |
| `raw-provenance.json` | SHA-256 hashes and sizes of the complete raw GPU scans and details of compact export. |
| `manifest.json` | Hashes of published code and artifacts, environment versions, and source commit. |

The NPZ files can be read with `numpy.load`. For the survey files, `parameters` has columns `(alpha, beta, m0, m1)`. `classes`, `lle_windows`, and `final_states` retain a start axis (three starts in line/tight, one in adding). `period` is the candidate peak recurrence for the first start. The existing [methods guide](../../docs/methods.md) defines the finite-time class codes. The published compact scans can seed CPU verification, but recurrence reanalysis requires regenerating the omitted raw `peaks` arrays.

The current scripts include stricter solver-success gating and the corrected candidate tolerance. Original per-run hashes in the reports identify earlier script revisions; they are preserved as provenance, not asserted to equal the final script hashes in `manifest.json`.
