# Additional candidate screening

Status: 29 October 2026. This is a research queue, not a list of newly discovered chaotic systems. An abstract that discusses periodic oscillation without mentioning chaos does not establish a publication gap.

The current best numerical lead is the [four-oscillator synchronematic model](synchronematic-research.md). Its scope is narrow: the underlying single-particle physical model already has reported chaos, and the article and author code disagree on a coefficient sign. Both coefficient variants nevertheless produce independently checked collective-chaos candidates in our calculations.

A second independently checked lead is the [light-switching LCE slider](lce-slider-research.md), with positive examples beyond the paper's numerical stiffness/damping ranges. Its full equations are available and admit an exact event-map implementation. General relay and LCE chaos are already known; exact-model novelty remains unverified.

## Other leads and exclusions

| Specific target | What was checked | Current decision |
| --- | --- | --- |
| Boiardi et al., three-dimensional morphoelastic hydrogel filament, [DOI 10.1007/s00707-026-04851-x](https://doi.org/10.1007/s00707-026-04851-x) | Open 2026 article gives a continuum rod model and periodic 2D/3D motions after primary and secondary instabilities; targeted chaos searches found no matching result | A substantive future target, but physical 3D motion is not a three-state ODE. The COMSOL implementation is available only by author request. A faithful independent continuum solver and spatial convergence checks are needed; no reduced substitute has been presented as the same model. |
| Dhiman et al., low-excitation Kerr optomechanics, [DOI 10.1038/s41467-026-73259-x](https://doi.org/10.1038/s41467-026-73259-x) | Open 2026 paper, semiclassical model and public code; observes self-sustained motion at low excitation | Lower priority for first-chaos searching: neighboring optomechanical/Kerr models already have an extensive chaos literature. A distinct low-excitation question would need a careful validity and equivalence analysis. No simulation performed here. |
| Oshima et al., molecular nickel-dithiolene memristive oscillator, [DOI 10.1002/apxr.202300117](https://doi.org/10.1002/apxr.202300117) | Publisher article describes temperature/current-dependent resistance, apparent inductance and a parallel capacitor; exact-title/material/chaos searches found no matching result | Analyze state dimension before scanning. With one thermal state and one capacitor, a smooth single-valued circuit branch is planar and cannot have a classical autonomous strange attractor. The implicit current-voltage relation may fold, and additional physical relaxation may matter; these must be checked rather than counting the apparent inductance as an automatically independent state. No numerical chaos claim. |
| Wu et al., LCE spring in a linear light field, [DOI 10.1016/j.chaos.2024.114587](https://doi.org/10.1016/j.chaos.2024.114587) | Title/abstract and targeted searches describe equilibrium/self-oscillation and a bifurcation analysis | Exact model/full text not retrieved in this pass. The OpenAlex snapshot reports no open full-text location. Do not invent the equations from the abstract. Any claimed gap remains unverified. |
| Zhao, Yang and Dai, three-line LCE planar oscillator, [DOI 10.1016/j.tws.2026.115328](https://doi.org/10.1016/j.tws.2026.115328) | Publisher preview describes one active and two passive fibers, in-plane oscillation, geometry-dependent switching and an experimental periodic regime | A potentially interesting multi-state mechanical/thermal target. The preview does not provide the complete equations; no open full-text location was returned. Its own introduction cites established LCE chaos, so the entire material class is not a novelty target. |
| Modi et al., enzyme/coacervate feedback, [DOI 10.1039/D2SC03838B](https://doi.org/10.1039/D2SC03838B) | Open publisher paper describes fuel-triggered pH change and droplet dissolution | Low priority for a sustained-attractor search in the reported closed-fuel setup. Negative feedback alone does not supply a recurrent driven dynamical system. A continuously fed extension would be a changed model. |
| Petrzela, higher-order oscillatory equations, [DOI 10.1038/s41598-024-72034-6](https://doi.org/10.1038/s41598-024-72034-6) | Primary article explicitly studies chaotic higher-order circuits and raises a question about structurally stable hyperchaotic behavior | Not a first-chaos target. A second positive finite-time exponent would not settle structural stability or uniform hyperbolicity. Potential mathematical project, not yet simulated here. |
| Periodically driven binary Quincke colloids, [arXiv:2610.19593](https://arxiv.org/abs/2610.19593) | October 2026 preprint abstract describes size-dependent mobility and a dynamic square lattice | Different system from our fixed-center phase/orientation model. No exact-model novelty assessment or simulation completed. |

These are bounded screening decisions, not global impossibility claims. The two-dimensional exclusion applies only to a smooth autonomous regular branch with the stated state variables. A physical device can contain unmodeled memory, forcing, switching, delay or spatial degrees of freedom.

## Selection rule for further work

Favor an explicit published model with enough effective dynamical states, a meaningful bounded operating domain, and a reason to expect loss of regular motion. Check exact identifiers, author work, supplements, dissertations, forward citations and equivalent equations. Then separate three questions:

1. Does the model actually show persistent deterministic irregularity under independent numerical checks?
2. Does that regime survive the approximation choices needed to relate the model to a physical system?
3. Is the precise result absent from the literature checked so far?

The answers need not coincide. Our optical-feedback search supplied a useful mechanism and a negative chaos result. The synchronematic investigation supplied a stronger numerical candidate, a source inconsistency, and close prior art that narrows any novelty claim. All three are useful outcomes of the search.
