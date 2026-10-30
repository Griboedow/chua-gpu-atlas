# HHS single-cell source controls: time units and a figure mismatch

Audit: 30 October 2026 UTC. **Six independent integrations of the literal single-cell equations do not reproduce the constant traces in supplementary Fig. S1 under the fixed-reference interpretation used here.** This is an unresolved reproduction question, not evidence of a new attractor or a conclusion about the authors' unseen implementation.

## Source settings and scope

Both the original manuscript [1] and latest public version [2] use the five-state single-cell equation S2. We set its optional force to zero. The main text defines `E0` as the uncoupled fixed-point concentration; here it is named `E_reference` to separate it from the initial concentration. SI initialization also uses `E0` for sampled initial values, so the notation needs care.

The source sets `tau_R=1` as its time unit, not one second. Defaults are `a=c=1`, `b=3`, `D0=.1`, `tau_D=1.1`, `tau_L=1.2`. Chemical starts are sampled from the uncoupled limit cycle, with `D=.1`, `L0=L=1`. Fig. S1 displays times 975–1000, panels labeled `mu=.80` and `1.50`, and rounded `(M,E)` starts `(1.05,3.43)`, `(.44,2.71)`, `(.49,4.21)`. Its caption instead describes equal `alpha,beta` values. We interpret the panel label as `mu=alpha*beta`.

The source's `mu=.01…9` sweeps and durations `5e4–1e5 tau_R` concern the spatial model. That model adds cubic regularizers and conserves total length; its chaos claims and simulation durations are not single-cell controls.

## Independent literal-equation controls

The fixed reference is `E_reference=3.8626272360560576`, calculated from the uncoupled chemical equilibrium in the [equation audit](late-candidate-screening-3.md). It was held fixed for every initial condition. We used `alpha=beta=sqrt(mu)` and the original physical-coordinate state order `(M,E,D,L0,L)`, without calling our GPU model or changing the source equations.

SciPy Radau integrated all six cases from 0 to 1000 with relative tolerance `1e-8`, absolute tolerance `1e-10` and maximum step `.2`. Each reached the requested endpoint. There were 501 dense-output samples in the last 25 units and 10,001 samples over the full interval. The six runs took approximately 79 seconds in total.

| `mu` | Initial `(M,E)` | Late `E` range, rounded | Late peak-to-peak `E` |
| ---: | --- | --- | ---: |
| .8 | (1.05, 3.43) | 1.385645–4.374427 | 2.988782 |
| .8 | (.44, 2.71) | 1.385644–4.374428 | 2.988784 |
| .8 | (.49, 4.21) | 1.385645–4.374424 | 2.988779 |
| 1.5 | (1.05, 3.43) | 1.236618–4.452877 | 3.216258 |
| 1.5 | (.44, 2.71) | 1.236816–4.452926 | 3.216110 |
| 1.5 | (.49, 4.21) | 1.237065–4.452923 | 3.215858 |

All six retain substantial late oscillations. Fig. S1 visually shows three essentially constant traces at `.8`, and a constant middle trace at `1.5`. Rounded starts, the panel/caption parameter interpretation, and the overloaded reference notation limit an exact comparison. We did not replace the reference with each initial `E` or fit parameters to the picture. No Lyapunov exponent, periodic-orbit closure or asymptotic attractor classification is claimed from these six controls.

Sampled lengths stayed positive. Nevertheless, maxima reached approximately `L=2.20–2.67` and `L0=2.45–3.02`; these are large deformations relative to the unit reference length. Positivity alone therefore does not validate a small-strain physical interpretation. The minima and maxima are sampled numerical diagnostics, not rigorous continuous-time bounds. No negative length was clipped.

The complete settings, endpoints, ranges, software versions and exact independent source code are embedded in [source-audit.json](../results/hhs/source-audit.json). To reproduce the six controls from the repository root:

```bash
python -c "import json; from pathlib import Path; a=json.loads(Path('results/hhs/source-audit.json').read_text()); exec(compile(a['independent_controls']['implementation_source'], '<archived control>', 'exec'))"
```

This executes our archived short control, not downloaded author code. It writes its numerical output under ignored `outputs/candidate_research/`.

## Source-code provenance and priority limits

The public arXiv v3 source archive was downloaded through its ordinary endpoint and inspected without executing or extracting its contents. Its 43 members comprise one TeX manuscript, four bibliography files, 22 PDF figures, 15 PNG figures and one submission metadata JSON. We found no simulation implementation or raw trajectory file there. The TeX independently confirms the reference definition, time unit and Fig. S1 caption. The manuscript links to supplementary movies, but no author simulation repository was identified in the inspected material or targeted identifier/title searches. This is not proof that no public or private code exists.

The latest paper already claims spatial chaos. Related Brusselator chaos on curved surfaces is also published [3]; it is a spatial model with a different coupling, not the exact S2 system. This bounded follow-up found no additional explicit chaotic parameter set for the literal autonomous single-cell equations. That negative search result does not establish novelty. The [earlier screening](late-candidate-screening-3.md) retains the broader version and citation qualifications.

## References

1. P. Dewan, S. Mondal, and S. Sarkar. *Oscillation death by mechanochemical feedback*. [arXiv:2504.19655v1](https://arxiv.org/html/2504.19655v1). Original SI sections I–II inspected.
2. S. Mondal, P. Dewan, L. S. Kumar, and S. Sarkar. *Symmetry-protected phases in a 1D active solid with mechanochemical feedback*. [Public v3 and SI](https://arxiv.org/html/2504.19655v3); [author source archive](https://arxiv.org/src/2504.19655v3); [DOI](https://doi.org/10.1103/p1bw-v87g). SI Fig. S1 visually inspected. Source figures and manuscript files remain in the ignored local cache and are not relicensed as MIT project material.
3. *Oscillatory and chaotic pattern dynamics driven by surface curvature*. *Physical Review E* **111**, L022202 (2025). [Primary publisher text](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.111.L022202). Related spatial Brusselator model, not evidence of chaos in S2.
