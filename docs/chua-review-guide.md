# Reviewing the prediction, excitation and coexistence results

This guide identifies what must be checked to accept each conclusion of the [combined study](chua-prediction-and-hidden-coexistence.md). It separates arithmetic that can be replayed from mathematical implications that require reading. It also identifies failures that would invalidate a conclusion. A file containing `passed: true` is a record of an execution, not an independent proof object.

## 1. Questions, quantities and evidence

| Result | Quantity being determined | Evidence | What remains undetermined |
|---|---|---|---|
| R1: prediction of thin windows | Location and width of three specified periodic families on two small slope patches | Forecasts published before 48 held-out searches; refined roots and independent orbit controls | A global predictor, a rigorous remainder bound and completeness of the atlas |
| R2: excitation of a hidden cycle | A basin-boundary intersection on seven rays; the first excitation infimum on the positive y ray | Exact local stable graph, transverse intersections and complete coverage of smaller y amplitudes | An optimum over all directions or finite-duration pulses; the entire basin boundary |
| R3: infinite hidden periodic coexistence | Existence of individual parameter values supporting infinitely many simultaneous hidden cycles | Exact-manifold tangency certificate, uniform equilibrium exclusion and localized parameter theorem | A decimal witness, a count of chaotic attractors, component tolerance or global sheet count |

The original mesh samples 32 catalogue entries under the stated branch-matching convention; its sheets interpolate selected representatives. A periodic window is an attracting interval along a specified parameter transect. Coexistence counts distinct attracting motions at one parameter point. These three counts need not agree. Symmetry-related cycles are distinct motions, although a family convention can group them together.

The model and hidden-attractor definition are stated before the computations in the [combined report](chua-prediction-and-hidden-coexistence.md#model-questions-and-interpretation). All three results use the ideal three-segment model with γ = 0, but different parameter sets. The coincidence of some parameter names does not make these the same experiment.

## 2. Predictions and excitation

### R1. Reconstruct the denominator and the score

Start with the [frozen forecasts](../results/chua-followup/sheets/forecasts-frozen.json) and [transport protocol](../results/chua-followup/sheets/transport-protocol.json), published in [112b01f](https://github.com/Griboedow/chua-gpu-atlas/commit/112b01f). There are two sections, nine slope nodes per section and three target families. The two central nodes provide six calibration targets. Only the remaining 48 targets enter the held-out score.

For each target, read `fold.coordinate`, `symmetry.coordinate` and `core.coordinate` in its `*-search.json` record. With f and s denoting the two boundaries, the measured stability interval of the symmetric cycle is the open interval between them and its width is |f − s|. Check the forecast midpoint against those endpoints. Recompute the corrected fold error divided by |f − s|, and the absolute relative width error. Do not score a forecast against its own predicted interval or count the central calibration targets as validation.

The [record checker](../research/review_chua_followup.py) recomputes these quantities, yielding 32/48 raw midpoint hits, 48/48 corrected hits, maximum corrected fold error 0.022103 widths and maximum width error 6.6803%. It also checks that the stored interior root lies between the endpoints and satisfies the numerical two-dimensional Schur stability inequalities. These are checks on computed roots, not interval existence proofs. The [independent orbit controls](../results/chua-followup/sheets/independent-checks.json) additionally test twelve corner cycles with a separate exponential calculation, ODE/variational integration and peak counting.

The displayed surfaces join nine solved folds per family. Their plotted thickness is not the stability width. In the lower figure panels the numerical organizing surface has been subtracted; the vertical coordinate is a residual parameter coordinate, not the original α or β.

### R2. Distinguish a local crossing from the first crossing

At α = 8.4, β = 12, m₀ = −0.95 and m₁ = −0.05, inspect the [local certificate](../results/chua-followup/trigger/local-certificate-100.json). For every direction it retains the image of the entire amplitude interval, its amplitude derivative and both endpoint images. The proof uses |χ′| ≤ 0.02 for the previously certified stable graph. The image must remain in its local square; the gap derivative must be positive; and the two endpoint gap signs must be opposite.

The [basin-boundary proof](chua-basin-and-coexistence.md#2-proof-of-theorem-a) is a dependency: local transversality alone cannot identify the two destinations. Its complete exit strips reach the origin and hidden cycle. Six half returns then transfer that statement to each tested initial-state ray.

Only the y ray has the additional [complete lower-amplitude coverage](../results/chua-followup/trigger/below-certificate-80.json). Check the central segment [0,0.299], the 662 adjoining captured intervals, their ordering and their final endpoint. Each retained chunk must be complete. Checking a selection of trajectories would not exclude an earlier small basin interval. The local crossing interval overlaps this coverage and supplies the remaining segment up to A*.

These ingredients prove that A* is an infimum; at A* the trajectory approaches the saddle cycle. The seven local brackets are not seven globally minimal excitation amplitudes. For nonzero angle labels the certified direction is the stored decimal vector, not an exact trigonometric expression.

The logarithmic delay is a separate asymptotic prediction. Its slope uses the saddle multiplier and period; the plotted intercepts are fitted. The 100 delay measurements and twelve long integration controls support that prediction, without proving a finite-time error bound for every impulse.

## 3. Infinite hidden coexistence

Read the [theorem statement and notation](chua-hidden-infinity-certificate.md#statement-and-scope) first. The fixed parameters are exact decimals; β is the only varying parameter. The conclusion asserts the existence of one β∞ supporting the entire infinite collection, arbitrarily close to the tangency β*. It is not a list of parameters supporting one cycle each.

The proof is the implication **T + H + N ⇒ R3**, with the following independently inspectable obligations.

| Obligation | Stored evidence | Producing routine | Reason to reject the inference |
|---|---|---|---|
| T1. A saddle of the correct return branch exists uniformly in β. | `parameter_family_saddle.krawczyk`, `root`, `path_cells` | `enclose_saddle` in [certify_chua_layer_organizer.py](../research/certify_chua_layer_organizer.py) | A small numerical residual without a strict inclusion, or an unvalidated earlier crossing. |
| T2. Exact local invariant manifolds satisfy the stated derivative bounds. | `local_map`, `graph_constants`, `graph_transform_bounds` | `local_coordinates` and `certify` in [certify_chua_homoclinic_tangency.py](../research/certify_chua_homoclinic_tangency.py) | Bounds on a tangent line substituted for bounds on the invariant graph. |
| T3. The complete excursion is valid and nontrivial. | `itinerary_cells`, `departure_witness`, `return_value` | `connect` in [chua_interval_connection.py](../research/chua_interval_connection.py), [event map](../research/chua_interval_lobe.py) | A correct endpoint reached after an invalid itinerary, or departure/arrival outside a manifold domain. |
| T4. The exact gap has a quadratic, transverse zero. | Four `face_signs`, `normal_curvature`, `parameter_unfolding` | `certify` and the stationary-minimum argument in proof §3 | A face sign contains zero; curvature or parameter speed is unproved; manifold or moving-saddle derivatives are omitted. |
| T5. The return saddle is dissipative. | `full_return_multipliers`, `full_return_determinant` | `certify`, using interval derivatives | Using the wrong map iterate or replacing an enclosure by a floating eigenvalue. |
| H. All equilibrium neighborhoods avoid every periodic basin uniformly near the tangency. | `uniform_beta_interval`, `spectrum`, `modal_basis`, `outer_lyapunov`, `captured_energy`, path counts | [certify_chua_hidden_neighborhoods.py](../research/certify_chua_hidden_neighborhoods.py) | Testing only eigenvector starts, only one parameter, or failing to cover the entire exit rectangle. |
| N. The parameter theorem applies to these partial return maps and its sinks remain in their domain. | Proof §§1 and 5; the [localization argument](chua-basin-and-coexistence.md#32-localization-and-extension-argument) | Mathematical argument, using the cited theorem | Extending a map and then counting sinks that visit the arbitrary part of the extension. |

T1–T5 refer to the [160-digit tangency record](../results/chua-followup/hidden-search/hidden-tangency-160.json), repeated at 140 digits. H refers to the [100-digit uniform exclusion record](../results/chua-followup/hidden-search/uniform-equilibrium-exclusion-100.json), repeated at 80 digits. The bound on H covers a larger interval J, with the whole tangency box strictly inside it. Pointwise hiddenness only at β* would not establish hiddenness at nearby β∞.

The proof also checks global forward boundedness in this family: the complete field is a Hurwitz outer linear system plus a bounded saturation term. The verified Lyapunov matrix makes its quadratic energy decrease outside a bounded ellipsoid. This is distinct from the local determinant condition in T5 and does not require assuming global boundedness from a finite simulation.

### Map conventions and interval serialization

F is one ordinary positive-lobe first return with two flights; P = F² has positive saddle multipliers. The connection uses F³⁶ = P¹⁸. For compatibility, the tangency JSON retains the older names `half_eigenvalues` for F and `full_return_multipliers` for P. In the uniform exclusion JSON, `outer_lyapunov.P` is instead the rational Lyapunov matrix, denoted R in the proof. These are different mathematical objects.

Arb's short `ball` string is an outward enclosure, but its display rounding can erase a strict inequality. In the central-segment record, for example, the displayed ball is `[1.0 +/- 6.26e-3]`, while its separately enclosed upper endpoint is approximately 0.999994186601. The checker forms the union of `lower` and `upper`, each itself an outward-rounded ball. It never uses a rounded midpoint as an exact interval. For extremely tight Krawczyk inclusions even those printed endpoints need not retain sufficient precision; replay the producing calculation.

### What is independent, and what is shared

The mpmath control uses its own matrix exponentials and event-time solves, and the ODE control integrates the original vector field. They can expose errors in the return implementation. They still share the model, seeds and intended itinerary. Neither validates the exact manifold bounds or the imported theorem independently.

The record checker recomputes scores, graph inequalities, strict signs, capture-energy comparisons and Lyapunov principal minors. It rejects both shifted parameter boxes and two deliberately corrupted records. It does not reconstruct every time tube or prove the meaning of a capture label. Certificate replay recomputes those enclosures using the same interval implementation. External mathematical review and a complete independent implementation remain distinct tasks.

The Lean development checks limited algebraic and stationary-zero implications documented in [the formalization report](chua-coexistence-audit.md). It does not check the new Chua tangency, graph transforms, return-domain localization or infinite hidden coexistence.

## 4. Reproduction without modifying published evidence

Run from the repository root using an environment with [requirements-certificate.txt](../requirements-certificate.txt) and [requirements-precision.txt](../requirements-precision.txt). These checks use the CPU. The commands below store new output on D:; replace only that output location when working elsewhere. Keep the published inputs unchanged.

First recompute the record-level implications:

```powershell
python research/review_chua_followup.py --output D:/chua-review/record-check.json
python research/manifest_refined_chua_atlas.py --check
```

Then replay the central interval arguments and independent numerical controls:

```powershell
python research/certify_chua_hidden_tangency.py --digits 160 --output D:/chua-review/tangency-160.json
python research/certify_chua_hidden_neighborhoods.py --uniform --digits 100 --output D:/chua-review/equilibria-100.json
python research/certify_chua_trigger.py local --digits 100 --output D:/chua-review/thresholds-100.json
python research/check_chua_hidden_tangency.py --output-directory D:/chua-review/independent
```

The much slower complete lower-ray calculation can be reproduced sequentially:

```powershell
python research/certify_chua_trigger.py below --digits 80 --output D:/chua-review/below-80.json
```

An unfinished output has `partial: true` and must not be treated as a successful full-coverage result. A different environment can change numerical basis choices, subdivision counts and the last displayed digits. Accept strict inclusions and the specified domains, not byte equality of a new run. The [revision audit](../results/chua-followup/review/README.md) records which checks were actually repeated for this revision.

The [full evidence index](../results/chua-followup/README.md) also gives reconstruction commands for the exploratory searches and prediction stage. Refitting the frozen predictor after inspecting the target data would not reproduce a prospective prediction test.

## 5. References and scope of comparison

1. N. V. Stankevich, N. V. Kuznetsov, G. A. Leonov and L. O. Chua, *Scenario of the Birth of Hidden Attractors in the Chua Circuit*, International Journal of Bifurcation and Chaos 27(12), 1730038 (2017). [DOI](https://doi.org/10.1142/S0218127417300385), [author manuscript](https://arxiv.org/abs/1710.02677). System (1), the hidden-attractor definition and §5 locate the relevant model and question.
2. N. Kuznetsov, T. Mokaev, V. Ponomarenko, E. Seleznev, N. Stankevich and L. Chua, *Hidden attractors in Chua circuit: mathematical theory meets physical experiments*, Nonlinear Dynamics 111, 5859–5887 (2023). [Article](https://doi.org/10.1007/s11071-022-08078-y). Sections 2 and 6 separate model assumptions and the coexistence-count question.
3. C. Robinson, *Bifurcation to infinitely many sinks*, Communications in Mathematical Physics 90, 433–459 (1983). [DOI](https://doi.org/10.1007/BF01206892). The parameter theorem used in R3.
4. E. Catsigeras, M. Cerminara and H. Enrich, *Simultaneous continuation of infinitely many sinks at homoclinic bifurcations*, Discrete and Continuous Dynamical Systems 29, 693–736 (2011). [Published article](https://www.aimsciences.org/article/doi/10.3934/dcds.2011.29.693), [author manuscript, titled *Simultaneous Continuation of Infinitely Many Sinks Near a Quadratic Homoclinic Tangency*](https://www.fing.edu.uy/~eleonora/Recopilacion/Archivos/InfinitosPozos_Premat-2008-108.pdf). Theorem 2.6 in that manuscript and its local return construction give the formulation used here.

The comparison identifies a mathematical question and checks the applicability of an established mechanism. It is not a systematic priority review. The present ideal-model results do not settle the broader Hopf conjecture, all hidden-chaos scenarios, the global basin geometry or experimental observability.
