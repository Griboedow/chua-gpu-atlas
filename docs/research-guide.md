# Study outline and verification scope

The central question is how the periodic layers of the Chua atlas are organized: which motions they represent, how they appear and overlap, and whether a finite catalogue can exhaust them.

The investigation proceeds from sampled parameter maps to continued periodic orbits, forecasts tested against new searches, and an interval proof of one organizing homoclinic tangency. These stages support different kinds of claims.

For the current prediction, excitation and hidden-coexistence results, start with the [claim-by-claim review guide](chua-review-guide.md). It states how each result can be used, what would invalidate it and how to repeat the checks without modifying the published evidence.

The [terminology guide](terminology.md) distinguishes periodic orbits, their stability regions, plotted sheets and catalogue entries.

## Studies and principal results

| Stage | Question and result | Report | Evidence |
|---|---|---|---|
| Predict slope layers, quantify excitation and prove hidden coexistence | 48 held-out windows; seven certified impulse crossings; a second tangency and uniform hiddenness imply infinite hidden coexistence | [Combined study](chua-prediction-and-hidden-coexistence.md), [hidden-infinity proof](chua-hidden-infinity-certificate.md) | [Forecasts, certificates and controls](../results/chua-followup/README.md) |
| Identify the layers | 32 sampled family sheets; periodic-orbit branch relationships, symmetry and behaviour | [Refined atlas](chua-refined-atlas.md) | [Atlas data](../results/chua-layer-atlas/) |
| Compare their motions | One complete orbit, signal and peak sequence per reference | [Family gallery](chua-family-gallery.md) | [Portrait data](../results/chua-layer-atlas/) |
| Predict thinner windows | Published forecasts compared with 75/85-peak computations; 34 reference entries | [Prediction study](chua-layer-prediction.md) | [Frozen forecasts and checks](../results/chua-layer-prediction/README.md) |
| Periodic-window accumulation | An exact quadratic tangency plus the sink-cascade theorem | [Proof report](chua-infinity-certificate.md) | [Interval certificate](../results/chua-layer-infinity/README.md) |
| Audit infinite coexistence | Local hypotheses pass; localization is now proved in the companion study; Lean proves only the stated general implications | [Coexistence audit](chua-coexistence-audit.md) | [Audit records](../results/chua-proof-audit/README.md), [formal source](../formal/README.md) |
| Explain overlaps | Symmetry breaking, coexistence of 2–4 cycles, 18 contact curves | [Overlap geometry](chua-overlap-geometry.md) | [Continuation and checks](../results/chua-overlap-geometry/README.md) |
| Certify finite coexistence and hiddenness | Four distinct periodic sinks; two hidden-cycle examples with excluded equilibrium neighborhoods; a global-convergence sector and circuit constraints | [Accessibility](chua-accessibility.md) | [Certificates and controls](../results/chua-accessibility/README.md) |
| Vary both slopes | Two-dimensional sections and two three-dimensional volumes | [Slope plane](chua-slope-plane.md), [volumes](chua-slope-volumes.md) | [Slope-plane data](../results/chua-slope-plane/), [volume data](../results/chua-slope-volumes/) |
| Resolve slope-family structure | Continued five-/seven-peak families and missed basins | [Slope families](chua-slope-families.md) | [Family data](../results/chua-slope-families/) |
| Refine two-slope families | 6,713 cross-sections; recovered P7 nodes; diode-contact curves and a tested duration law | [Slope refinement](chua-slope-refinement.md) | [Refinement data and controls](../results/chua-slope-refinement/README.md) |
| Prove local separation and infinite coexistence | Uniform capture of both saddle exit strips; localized Newhouse–Robinson application | [Two local results](chua-basin-and-coexistence.md) | [Certificates and checks](../results/chua-basin-proof/README.md) |
| Exclude oscillations and locate a basin separator | Explicit global-convergence sector; certified saddle beside a hidden sink; tested two-thirds/one-third contact laws | [Global stability and basins](chua-hidden-boundary.md) | [Storage identities, certificates and controls](../results/chua-hidden-boundary/README.md) |

The [main page](../README.md) summarizes the study. The [orbit catalogue](chua-family-gallery.md) describes the reference motions, the [prediction study](chua-layer-prediction.md) compares forecasts with subsequent computations, and the [certificate](chua-infinity-certificate.md) presents the mathematical argument and its scope.

## Counts and terminology

**17** was an early census on one transect. **32** is the subsequent catalogue with sampled three-dimensional sheets. **34** includes two later reference cycles predicted before discovery. These are operational counts under a periodic-root and branch-matching convention.

**Infinitely many** refers to attracting periodic windows with unbounded return periods near a certified parameter on the fixed classical transect. The localized parameter theorem also gives individual nearby parameters with infinitely many simultaneous periodic attractors. Neither result makes the finite mesh exhaustive or proves that every motion belongs to a globally disconnected family surface.

The later hidden-region tangency supplies a separate fixed-α, fixed-slope family in which infinitely many simultaneous attracting cycles are all hidden. Its equilibrium-neighborhood certificate is uniform on an interval surrounding the tangency. The distinction between window count, simultaneous cycle count and connected family-sheet count remains necessary.

A family can have multiple attracting descendants at the same parameters, particularly after symmetry breaking. A parameter projection can also merge distinct branches. Counts of families, individual attracting cycles, connected components and projected islands must therefore be kept separate.

## Verification methods and scope

| Type | What supports the claim | What it does not establish |
|---|---|---|
| Finite-time survey | Repeated integrations, multiple starts, step and duration controls | An asymptotic chaotic attractor or exhaustive basin coverage |
| Computed periodic cycle | Closed orbit, valid itinerary, attracting multipliers and independent controls | A global component count |
| Interval certificate | Outward-rounded inclusions and inequalities for complete parameter boxes and flight intervals | A classification outside the certified neighborhoods |
| Theorem consequence | Verified hypotheses plus the established sink-cascade theorem | Exact widths and peak labels for every uncomputed member |
| Partial formal proof | Lean checks the general boundary-sign implication, rational guards and three storage-algebra implications | The Chua numerical implementation, invariant-manifold construction, full global-convergence argument or infinite coexistence |

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
