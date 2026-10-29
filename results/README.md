# Recorded evidence

These JSON files are copies of the computation settings and validation summaries supplied in the complete results release. They allow checking the headline claims without downloading the large arrays.

- `original/`: two 51³ scans, reference trajectories and eight numerical validation checks.
- `detailed/`: two 81³ grids, 41³ local zoom, beta transects, example regimes, connectivity summaries and independent validation.
- `extent/`: wide-scout settings, repeated-candidate interpretation, per-parameter seed mapping and independent event-resolved checks.

The large basin metadata, complete NPZ arrays, all offline interactive reports and browser-check records are in [release v1.0.0](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.0.0). Extract the archive into the repository root. Its internal `RESULTS_MANIFEST.json` provides per-file SHA-256 hashes; `SHA256SUMS.txt` verifies the ZIP itself.

These are finite-time numerical experiments. The far-field exponent discrepancy documented in `extent/event_resolved_checks.json` is unresolved, not a passed convergence test. See the main README and methods for interpretation.
