# Local basin separation and infinite periodic coexistence

The [report](../../docs/chua-basin-and-coexistence.md) contains two completed mathematical arguments, their domains and their limitations. This directory supplies the new finite computations and the application audit. Neither complete theorem is formalized in Lean.

| File | Role |
|---|---|
| `separator-80.json`, `separator-100.json` | Complete covers of two local saddle exit strips: 64 rectangles per side, all-time event checks, capture domains and continuous-flow enclosures where needed |
| `independent-checks.json` | mpmath value/Jacobian/Hessian checks at two difference scales; rejection of a shifted fixed-point box; five independent DOP853 checks across a diode contact |
| `audit.json` | Coverage, adjacency, destinations, capture inequalities, repeated-precision summaries, tangency premises and source hashes |
| `preliminary-failures.json` | Unsuccessful enclosure strategies and the reasons for replacing or refining them |
| `chua-basin-proof.pdf` | Vector version of the local basin and coverage figure |
| `chua-coexistence-localization.pdf` | Labelled schematic of the localized return and nested parameter construction |
| `browser-check.json` | Live browser review of the published report and main page |

From the repository root, using the existing CPU certificate dependencies:

```bash
python research/certify_chua_basin_separator.py --digits 80
python research/certify_chua_basin_separator.py --digits 100
python research/check_chua_basin_proof.py
python research/audit_chua_local_results.py
python research/plot_chua_local_results.py
python research/plot_chua_proof_audit.py
python research/manifest_refined_chua_atlas.py --check
```

`--only INDEX` runs one u-interval on each side and produces a `-pilot.json` file marked `partial`. Such a file is not a complete proof. The audit requires both full covers and rejects partial runs, missing cells, unresolved destinations or uncovered adjacent intervals. Continuous-flow subdivisions must certify both children. Rejected intermediate enclosures are retained inside accepted continuous-flow records.

The two basin certificates use exact rational parameters α = 42/5, β = 12, m₀ = −19/20, m₁ = −1/20, γ = 0. They enclose all points in a section square of coordinate radius 10⁻⁷. The curve drawn over a larger domain in the earlier basin study remains a numerical approximation.

Infinite coexistence is a separate result on α = 9, m₀ = −8/7, m₁ = −5/7, with β varying. Its interval premises are the unchanged [80/100-digit homoclinic certificates](../chua-layer-infinity/README.md). The domain/localization argument and the imported Newhouse–Robinson theorem are analytical dependencies, not assertions established by `audit.json`. The conclusion is existential: no finite decimal witness or exhaustive list of the infinitely many attracting cycles is supplied.

The independently integrated trajectories test the implementation. They do not replace the uniform enclosures, graph transform, Lyapunov argument or classical theorem. All scratch work and computations for this study run locally on D:; no hosted workflow is used.
