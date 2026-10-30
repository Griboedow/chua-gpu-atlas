# Oscillating heat pipes: source and priority gate

Audit: 30 October 2026 UTC. **Do not launch a search of a reconstructed 2026 model.** Its complete equations and numerical control were not retrieved. The accessible 2024 predecessor already describes prospective chaotic motion, so even that older exact model would be a quantitative validation target rather than an unsuspected chaotic family. No simulation or author-code execution was performed here.

## The requested 2026 model

Chicone, Feng and Lombardo [1] describe a reduced model with three ODEs per liquid-slug/vapor-plug pair: two mechanical states and one phase-change state. It studies startup, rest, stagnation, stopovers, multiple attractors and thermal conductance. Their introduction distinguishes this finite-temperature-difference extension from an earlier isothermal model and explicitly acknowledges earlier single-branch chaos. A prescribed spatial temperature profile supplies the thermal driving. These are statements about the published model, not a reconstruction of its vector field.

The publisher exposed the abstract, introduction and section excerpts, but not the equations. A 31-page SSRN preprint [2] is indexed; ordinary requests to its landing page and public delivery URL returned HTTP 403. The final article is marked closed in OpenAlex. Neither a verified public code repository nor a complete final-model parameter tuple was recovered. These access observations do not imply that no public version exists.

Three states per pair do not establish the number of independent asymptotic degrees of freedom. A conservation/redundancy audit must precede a Lyapunov-spectrum interpretation. Likewise, slug mergers cannot silently become arbitrary numerical reset rules. Both issues remain gates for the final model.

## What was actually recovered: the 2024 predecessor

The author-upload manuscript [3], equations (5), (6), (12), (13), supplies the following **isothermal predecessor**, not a substitute for [1]. Let $\theta_i$ denote a slug meniscus, $v_i$ its velocity, $\gamma_i$ its fixed length, $\phi_i$ scaled vapor mass and $\ell$ total scaled tube length. Indices are cyclic.

```math
\dot\theta_i=v_i,\qquad
\dot v_i=\frac{p_i-p_{i+1}}{\gamma_i}-\nu v_i,\qquad
\dot\phi_i=e(\theta_{i-1}+\gamma_{i-1})+e(\theta_i).
```

```math
p_i=\frac{\phi_i}{\beta_i},\qquad
\beta_1=\theta_1+\ell-\theta_N-\gamma_N,\qquad
\beta_i=\theta_i-\theta_{i-1}-\gamma_{i-1}\quad(i\ge2).
```

```math
e(\theta)=-R\,\mathrm{clip}\!\left(\frac{Y(\theta)}{\epsilon},-1,1\right).
```

Here $Y$ is the vertical coordinate along the specified serpentine tube, and clip limits its argument to the stated interval. The source already derives $S=\sum_i\gamma_i\theta_i=c_1+c_2e^{-\nu t}$. Section 3.2 describes apparent torus-to-chaotic motion for three slugs and stronger irregularity for five; it explicitly leaves precise attractor characterization unresolved.

The figures use $r/L=0.1$, $h/L=0.03$, $\nu=0.5$, $\sigma=R/\epsilon$: three slugs at approximately $0.617$ and $1$, five at $0.65$ and approximately $0.834$. These captions alone are **not a complete reproducible tuple**: the global geometry, separate $R,\epsilon$, initial menisci and masses must also be specified. Published oscillation shapes are not a Lyapunov calculation.

## Low-cost mathematical consequences

These are direct checks of the displayed predecessor equations, using unwrapped ordered positions before mergers, not new theorems claimed for the unavailable extension. Summing the momentum equations gives

```math
P=\sum_i\gamma_i v_i,\qquad \dot P=-\nu P.
```

For positive damping, every compact recurrent set must lie in $P=0$: any nonzero value decays forward and grows without bound backward. On that set $S$ is constant. Thus the attractor analysis must respect an invariant leaf, rather than treating all $3N$ plotted coordinates as independent recurrent directions. This is consistent with the source's existing reduction, not a new discovery.

For one pair on the closed ring, pressure differences telescope to zero in its sole mechanical equation. The velocity decays; this one-pair limit cannot be identified with the **open-ended** single-branch oscillator used in older chaos studies. Equal nominal state counts do not establish equivalence.

Positive plug lengths and vapor masses are necessary physical conditions. The denominator $\beta_i$ makes collision/plug collapse a singular boundary. The displayed evaporation law does not automatically enforce $\phi_i\ge0$: its right-hand side can be negative at zero mass. Consequently, numerical boundedness in unconstrained coordinates would not establish a physically admissible attractor. A future implementation needs boundary detection and an explicit, sourced decision about termination or merger. No closure is supplied here.

The serpentine drawing is the **fixed apparatus geometry**, not an attractor. Any three-dimensional plot of a many-slug state would be a projection. Visual resemblance to a coiled pipe would therefore not be evidence of a new attractor topology.

## Priority and forward-citation check

| Source | Evidence actually inspected | Consequence |
|---|---|---|
| Chicone et al. 2024 [3] | Author-upload full text, especially §§2–3.2 | The direct predecessor already proposes chaotic multi-slug regimes; no priority claim based on an absent keyword in its title. |
| Kumar and Singh 2023 [4] | Primary publisher abstract | Explicit period doubling, Lyapunov exponent and correlation-dimension analysis in a single-branch model. Full equations were not retrieved, so no exact equivalence to [1] or [3] is asserted. |
| Kumar and Singh 2023 [5] | Primary abstract/introduction and section excerpts | Separate bifurcation study identifies period-doubling and Neimark–Sacker routes. It should not be conflated with the numerical chaos diagnosis of [4]. |
| Zhou et al. 2026 [6] | Primary abstract/highlights; OpenAlex forward link | Two-dimensional thermofluidic comparison of uniform/dual tube diameters. Retrieved material does not establish a chaos result for the exact reduced ODE model. Full text remains unchecked. |

The OpenAlex snapshot returned one indexed citing work for [1], namely [6]. Its online/publication metadata and December 2026 issue date differ; it was already searchable at this audit. Citation-index coverage is incomplete and this single result cannot clear scientific priority.

**Disposition:** stop this candidate at the source/priority gate. A later defensible project could quantify the predecessor's already-suggested irregular regimes, distinguish long transients and collisions from sustained chaos, and compare the finite-temperature extension after obtaining its actual equations. No new attractor, new mechanism, or absence of prior publication is established by this audit.

## References

1. C. Chicone, Z. C. Feng and S. J. Lombardo, “Modeling of oscillating heat pipes: Equilibria, stagnation, startup, bifurcation, oscillation, conductance, and statistics,” *Applied Thermal Engineering* **289**, 129761 (2026). [DOI](https://doi.org/10.1016/j.applthermaleng.2026.129761). [Publisher](https://www.sciencedirect.com/science/article/pii/S1359431126000694).
2. S. Lombardo, C. Chicone and Z. C. Feng, same-title preprint, SSRN 5732707, posted 10 November 2025. [Landing page](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5732707). [DOI](https://doi.org/10.2139/ssrn.5732707). Full PDF not retrieved.
3. C. Chicone, Z. C. Feng, S. J. Lombardo and D. G. Retzloff, “Start-up of oscillating heat pipes via Hopf bifurcation,” *International Journal of Mathematical Modelling and Numerical Optimisation* **14**, 32–52 (2024). [Final DOI](https://doi.org/10.1504/IJMMNO.2024.143238). [Author-upload manuscript actually read](https://www.researchgate.net/publication/377035225_Startup_of_oscillating_heat_pipes_via_Hopf_bifurcation), associated with early DOI 10.1504/IJMMNO.2024.10059141. Draft pagination/wording need not equal the final version.
4. A. Kumar and S. Singh, “On the interaction between liquid slug and vapor bubble in the chaotic operation of pulsating heat pipe,” *Physics of Fluids* **35**, 054113 (2023). [DOI](https://doi.org/10.1063/5.0147852).
5. A. Kumar and S. Singh, “Bifurcation analysis of a pulsating heat pipe,” *International Journal of Thermal Sciences* **192**, 108384 (2023). [Publisher](https://www.sciencedirect.com/science/article/abs/pii/S1290072923002454). [DOI](https://doi.org/10.1016/j.ijthermalsci.2023.108384).
6. G. Zhou, Y. Lv, Z. Wang and J. Qu, “Comparison of operational and thermofluidic characteristics between uniform-diameter and dual-diameter oscillating heat pipes: From a perspective of numerical simulation,” *International Journal of Heat and Mass Transfer* **271**, 129266 (2026). [Publisher](https://www.sciencedirect.com/science/article/abs/pii/S0017931026009427). [DOI](https://doi.org/10.1016/j.ijheatmasstransfer.2026.129266).

No publisher text, figures or code is relicensed by this repository. Machine-readable audit: [heat-pipe-source-audit.json](../results/candidate-audits/heat-pipe-source-audit.json).
