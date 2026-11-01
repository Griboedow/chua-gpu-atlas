# Research guide

The central question is how the periodic layers of the Chua atlas are organized: which motions they represent, how they appear and overlap, and whether a finite catalogue can exhaust them.

The investigation proceeds from sampled parameter maps to continued periodic orbits, forecasts tested against new searches, and an interval proof of one organizing homoclinic tangency. These stages support different kinds of claims.

## Main line of the study

| Stage | Question and result | Report | Evidence |
|---|---|---|---|
| Identify the layers | 32 sampled family sheets; orbit ancestry, symmetry and behaviour | [Refined atlas](chua-refined-atlas.md) | [Atlas data](../results/chua-layer-atlas/) |
| Compare their motions | One complete orbit, signal and peak sequence per reference | [Family gallery](chua-family-gallery.md) | [Portrait data](../results/chua-layer-atlas/) |
| Predict thinner windows | Published forecasts followed by 75/85-peak discoveries; 34 reference identities | [Prediction study](chua-layer-prediction.md) | [Frozen forecasts and checks](../results/chua-layer-prediction/README.md) |
| Establish infinity | An exact quadratic tangency plus the sink-cascade theorem | [Proof report](chua-infinity-certificate.md) | [Interval certificate](../results/chua-layer-infinity/README.md) |
| Audit infinite coexistence | Local hypotheses pass; localization remains open; Lean proves the tangency-box implication | [Coexistence audit](chua-coexistence-audit.md) | [Audit records](../results/chua-proof-audit/README.md), [formal source](../formal/README.md) |
| Explain overlaps | Symmetry breaking, coexistence of 2–4 cycles, 18 contact curves | [Overlap geometry](chua-overlap-geometry.md) | [Continuation and checks](../results/chua-overlap-geometry/README.md) |
| Vary both slopes | Two-dimensional sections and two three-dimensional volumes | [Slope plane](chua-slope-plane.md), [volumes](chua-slope-volumes.md) | [Slope-plane data](../results/chua-slope-plane/), [volume data](../results/chua-slope-volumes/) |
| Resolve slope-family structure | Continued five-/seven-peak families and missed basins | [Slope families](chua-slope-families.md) | [Family data](../results/chua-slope-families/) |

For a first reading, start with the [main page](../README.md), compare the [motion portraits](chua-family-gallery.md), then read the [prediction study](chua-layer-prediction.md) and the [certificate](chua-infinity-certificate.md). The latter contains the mathematical argument and separates its conclusion from the remaining classification questions.

## Counts and terminology

**17** was an early census on one transect. **32** is the subsequent catalogue with sampled three-dimensional sheets. **34** includes two later reference cycles predicted before discovery. These are operational counts under a periodic-root and ancestry convention.

**Infinitely many** refers to attracting periodic windows with unbounded return periods near a certified parameter on the fixed classical transect. It does not turn the finite mesh into an exhaustive atlas or prove that each visual sheet is globally disconnected from every other family.

A family can have multiple attracting descendants at the same parameters, particularly after symmetry breaking. A parameter projection can also merge distinct branches. Counts of families, individual attracting cycles, connected components and projected islands must therefore be kept separate.

## Strength of the evidence

| Type | What supports the claim | What it does not establish |
|---|---|---|
| Finite-time survey | Repeated integrations, multiple starts, step and duration controls | An asymptotic chaotic attractor or exhaustive basin coverage |
| Computed periodic cycle | Closed orbit, valid itinerary, attracting multipliers and independent controls | A global component count |
| Interval certificate | Outward-rounded inclusions and inequalities for complete parameter boxes and flight intervals | A classification outside the certified neighborhoods |
| Theorem consequence | Verified hypotheses plus the established sink-cascade theorem | Exact widths and peak labels for every uncomputed member |
| Partial formal proof | Lean checks the general boundary-sign implication and rational guards | The Chua numerical implementation, invariant-manifold construction or infinite coexistence |

The 85-peak forecast miss, failed whole-period multiplier checks and two rejected tangency boxes are retained. The independent checks and the interval proof are distinct forms of evidence.

## Broader surveys and scaling

The original GPU experiments explored local chaos candidates, the double scroll, regular windows and basins. Wider surveys then enlarged the search domain to α = 1024 and m₀ = −29. A separate fixed-slope study followed β = κα² through α = 2048.

| Topic | Report or data |
|---|---|
| Numerical model, solver methods, thresholds and complete reproduction | [Methods](methods.md) |
| Local dynamics, basins, wide and expanded surveys | [Experiment notes](experiment-notes.md), [recorded evidence](../results/README.md) |
| Early periodic branches and births | [Bifurcations](chua-layer-bifurcations.md) |
| Earlier central count and first adding sequence | [Layer census](chua-layer-census.md) |
| Large-α numerical sections and solver controls | [Scaling experiments](chua-scaling-research.md) |
| Exact rescaling, limiting equations and amplitude drift | [Scaling structure](chua-scaling-structure.md) |

The expanded survey retained 1,832 finite-time candidates from 16,384 prescribed Sobol samples. Its cloud reaches the imposed search limits; it does not establish global unboundedness. Independent full-field checks replace an older affine reference routine that missed a crossing at a distant point. The [results index](../results/README.md) identifies the superseded records.

The large-α problem is separate from the infinitely many windows proved near β ≈ 14.30636. An infinite hierarchy in a bounded parameter neighborhood says nothing by itself about chaos at arbitrarily large α.

## Files and provenance

- `docs/` contains the reports and displayed figures.
- `research/` contains continuation, analysis, interval and plotting programs.
- `formal/` contains the pinned Lean project and its explicit verification scope.
- `results/` contains settings, small datasets, certificates and checks. Each recent study has a reproduction guide.
- The root scanner and OpenCL kernels support the original broad surveys.
- [Release v1.1.0](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.1.0) supplies large historical arrays and offline HTML reports. It predates the recent layer studies.

The [manifest](../results/chua-layer-atlas/manifest.json) hashes the finite atlas, subsequent studies and relevant programs. The frozen forecasts retain their original contents and commit. The [transfer record](research-transfer.md) documents the move of non-Chua work to a separate repository.
