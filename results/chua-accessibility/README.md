# Coexistence, hiddenness and circuit accessibility: evidence

The [report](../../docs/chua-accessibility.md) gives the analytical arguments, parameter conventions, figures and limits. The new certificates are computer-assisted Python/FLINT results, not Lean theorems.

## Files

| Files | Contents |
|---|---|
| `four-cycle-certificate-{80,100}.json` | Four first-return fixed points, complete itineraries, Schur stability, disjointness and a rejected shifted box |
| `hidden-cycle-certificate-{80,100}.json` | One periodic sink, equilibrium exclusion neighborhoods, saddle exit tube and invariant ellipsoids |
| `stability-certificate-{80,100}.json` | A second hidden cycle with a unique stable equilibrium, a global Lyapunov sector, a rejection control and local-stability predictions |
| `independent-checks.json` | DOP853/Radau closure, variational multipliers and repeated-precision comparisons |
| `hidden-cycle-seed.json` | Numerical guess used for the first hidden-cycle enclosure |
| `hidden-inner-continuation.json` | Exploratory m₀ continuation used to select and initialize the prediction tests |
| `hidden-slope-continuation.json` | Supplementary numerical m₁ continuation; no complete parameter-domain certificate |
| `stability-basin-check.json` | Fixed-step checks at four m₀ values, including unresolved outcomes |
| `adaptive-prediction-check.json` | Independent adaptive checks of 48 unstable-side starts, including nine longer follow-ups |
| `adaptive-prediction-first-horizon.json` | Preserved time-1000 result: 39 of 48 pass the distance criterion |
| `equilibrium-screen.json`, `.npz` | 768-start four-cycle basin screen; 557 step-refinement disagreements |
| `*-linear.json` | Superseded section classifications using linear event interpolation |
| `circuit-constraints.json` | Component map, nominal lossy model, 128 random perturbations and three independent five-segment-diode checks |
| `*.pdf` | Exportable versions of the three report figures |

The `ball` representation of an Arb enclosure includes decimal serialization error. Separate `midpoint` and `radius` strings are diagnostic and must not be interpreted as a lossless interval pair. A floating-point plot or a `passed` field alone is not the proof; reproduce the executable inequalities and read the analytical reduction in the report.

## Reproduction

Run locally from the repository root using `requirements-certificate.txt`. The commands below require no GPU, hosted workflow or new dependency beyond that environment. On the research workstation, use the repository on D: and its existing Python environment.

The stored numerical seeds suffice for the certificates:

```bash
python research/certify_chua_coexisting_cycles.py --digits 80
python research/certify_chua_coexisting_cycles.py --digits 100
python research/certify_chua_hidden_cycle.py --digits 80
python research/certify_chua_hidden_cycle.py --digits 100
python research/study_chua_stability.py --digits 80
python research/study_chua_stability.py --digits 100
python research/check_chua_accessibility.py
```

To reconstruct exploratory seeds before repeating those checks:

```bash
python research/prepare_chua_accessibility.py
```

For the numerical prediction and physical-model experiments:

```bash
python research/study_chua_stability.py --basins
python research/check_chua_stability_prediction.py
python research/probe_chua_accessibility.py
python research/study_chua_circuit_constraints.py
python research/plot_chua_accessibility.py
python research/manifest_refined_chua_atlas.py --check
```

The last command checks identity with the recorded publication. Recomputed floating-point seeds, figures or metadata can change their file hashes without invalidating a mathematical enclosure. Run the interval and independent checks to assess numerical reproduction; do not regenerate the publication manifest merely to conceal a mismatch.

All six cycles pass independent DOP853 and Radau integration. The component sample is finite and does not establish a tolerance-box theorem. The four-cycle basin screen is explicitly inconclusive. The global Lyapunov sector and the separate hidden-cycle point do not locate an exact boundary between them.
