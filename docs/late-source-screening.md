# Further source screening: reported chaos and dimensional exclusions

These are literature and model-selection results, not new attractor discoveries.

## Modular opto-magnetic oscillator: chaos is already reported

Fu et al. [1] was initially retained from its periodic-oscillation abstract. Recovering the complete author-hosted PDF changes that decision. Page 4 explicitly calls the irregular high-intensity regime a **chaotic oscillation mode**, and page 7 makes another chaos claim under unfocused light. The first discussion concerns laser intensity above 39.4 mW/mm² and points to supplementary Figure S8; the second points to Figure S15.

This is enough to exclude the work from a search for a system with no prior chaos report. It is not enough to independently verify a deterministic strange attractor. Our inspection of those passages found irregular amplitude and timing, without a Lyapunov spectrum there. The supporting DOCX returned HTTP 403, so its contents and full dynamical closure remain unaudited. No fitted three-state surrogate was substituted. The PDF is cached locally for research and is not redistributed under this project's MIT license.

## Thermo-responsive gel: the published lumped model is planar

The primary model of Wang et al. [2], Eqs. 3–8, evolves polymer volume fraction and temperature:

```math
\dot\phi=\frac{\Pi_{\rm el}(\phi)-\pi_{\rm osm}(\phi,T)}{\tau_\phi},
\qquad \dot T=\gamma f(\phi)-k(T-T_0).
```

The pressure laws are explicit functions of these two variables; the heating gate does not have its own dynamics. Their stated Hill gate is continuously differentiable for exponent at least two. On a compact regular region with `0 < phi < 1` and positive absolute temperature, this is an autonomous planar flow. The Poincare–Bendixson restriction therefore rules out a classical chaotic strange attractor of this model. A fold, hysteresis or coexistence of a stable equilibrium and a periodic orbit does not evade that restriction.

This is our structural exclusion, not a claim about all thermo-responsive gels. Spatial deformation, thermal diffusion, forcing or another relaxation variable could change the state space. No GPU sweep is warranted for the literal two-state equations.

## A mathematical lead: minimal two-strain dynamics

Gavish [3] resolves an eigenvalue degeneracy that invalidated an earlier stability approximation. The minimal lifelong partial-cross-immunity model has seven independent compartments; the published oscillatory examples are unusually narrow in parameter space. The author repository was inspected as source material, without copying or executing its implementation.

The related preprint [4] supplies a five-state baseline with one direction of secondary infection suppressed. Its population simplex is invariant and bounded. It presents local and large cycles, calls some patterns possibly aperiodic, and mentions a planned companion analysis of a Shilnikov-like connection. This substantially narrows any possible priority claim. Its linked author repository returned HTTP 404 during this audit.

Before simulation, compare the exact equations with older two-strain chaos models, recover precise example parameters and handle slow demographic timescales. The first paper's Figure 3 author script uses a final time of `1e6`; the indexed paper describes `1e5`, so that reproducibility detail also needs reconciliation. General epidemic chaos is already established, for example in [5]. No biological prediction or exact-model novelty claim is made here.

**Subsequent work:** the [five-state numerical study](two-strain-research.md) now supplies GPU/CPU candidates, precision controls, a periodic orbit and an analytic boundary obstruction. The seven-state first paper was not surveyed. Exact-model priority and a biological interpretation remain unresolved.

## References

1. Fu, Y.-H., et al. (2026). Modular Opto-Magnetic Oscillators: Harnessing Light to Drive Versatile Materials and Functionalities. *Advanced Materials* **38**(16), e20348. [doi:10.1002/adma.202520348](https://doi.org/10.1002/adma.202520348). Author-hosted PDF, pp. 4 and 7; page 4 visually checked.
2. Wang, Y., Yuan, L., Ren, L., Liu, Z., and Gao, Q. (2026). Thermo-responsive self-oscillating gel: mathematical model and theoretical analysis. [arXiv:2606.21225v1](https://arxiv.org/html/2606.21225v1), Eqs. 3–8. Preprint.
3. Gavish, N. (2026). Hidden Oscillations: How Minimal Immunological Interactions Destabilize Two-Strain Epidemic Models. *SIAM Journal on Life Sciences* **1**(3), 420–443. [doi:10.1137/25M1826160](https://doi.org/10.1137/25M1826160). [Author repository](https://github.com/NGavish/HiddenOscillations); no root license was listed.
4. Gavish, N. (2026). An Analytically Tractable Framework for Multi-Strain Epidemics: Resolving Algebraic Complexity to Map Oscillatory Dynamics. [arXiv:2608.13995v1](https://arxiv.org/html/2608.13995v1), Eqs. 2 and 12, discussion. Preprint; local full HTML saved.
5. Aguiar, M., Kooi, B., and Stollenwerk, N. (2008). Epidemiology of dengue fever: A model with temporary cross-immunity and possible secondary infection shows bifurcations and chaotic behaviour in wide parameter regions. *Mathematical Modelling of Natural Phenomena* **3**(4), 48–70. [doi:10.1051/mmnp:2008070](https://doi.org/10.1051/mmnp:2008070). Prior-art warning; equation-by-equation comparison remains pending.
