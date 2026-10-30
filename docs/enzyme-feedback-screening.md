# Enzyme feedback: structural screening before numerical searches

Audit date: 2026-10-29 UTC. This is a model audit, not a discovery claim. No parameter survey or new chaotic trajectory was produced in this block.

Later work: the [full urease survey](urease-research.md) and [structural audit](urease-structural-audit.md) implement the next steps below. The [esterase source audit](esterase-compartment-audit.md) identifies the remaining supplementary-model gap. Their results do not retroactively turn this initial screening into a numerical study.

The October 2026 artificial-cell model and the 2025 two-variable urease oscillator should not enter a search for an autonomous strange attractor. The original four-species urease network remains a concrete candidate for further investigation, subject to the exact dimensional constraint below. The three-state Milster microreactor also has a structural obstruction, even before its published planar reduction.

## Decisions

| Model | Exact restriction or unresolved issue | Decision |
|---|---|---|
| Göth–Dzubiella 2026 closed GOx/urease artificial-cell mixture | Two finite fuels decrease monotonically; pH and permeability have no independent dynamical state | Exclude sustained autonomous chaotic attractors in this published model |
| Straube et al. 2025 geometric urease oscillator | Smooth autonomous planar flow in the positive domain | Exclude a strange attractor; limit cycles are already established |
| Milster et al. 2023 single-enzyme microreactor, full three-state model | A sign change makes it competitive; compact limit dynamics embed in a plane | Exclude bounded chaotic limit attractors under the stated signs |
| Straube–Winkelmann–Höfling 2023 original four-species urease network | Exact attracting three-dimensional hyperplane when product and substrate outflow rates agree; the subsequent planar reduction is approximate | Eligible conditional candidate; simulate the full equations before claiming anything |
| Leathard–Beales–Taylor 2023 esterase compartment | Published seven-species model and known oscillations; equations and supplementary XPP code not yet audited here | Secondary source lead, not a screened candidate |

These conclusions concern the equations and parameter signs specified below. They do not exclude chaos in spatial, delayed, forced, coupled, or chemically extended versions.

## 1. Closed artificial cells: a monotone fuel budget

Göth and Dzubiella [1] describe glucose-oxidase and urease populations in a shared, well-mixed environment. Equation 1 consumes glucose and creates gluconic acid at a nonnegative rate. Urease consumes urea and creates two ammonia molecules plus carbon dioxide. Their rates include nonnegative activity, permeability, and saturation factors. The paper discusses switching and pH bistability, with external acid/base additions used to obtain repeated cycles.

Write the two reaction rates as `vG >= 0` and `vU >= 0`. With glucose `G`, urea `U`, gluconic acid `A`, ammonia `N`, and carbon dioxide `C`, the stoichiometry is

```text
G' = -vG       A' = vG
U' = -vU       N' = 2 vU       C' = vU.
```

Our direct consequence is

```math
L=G+U,\qquad \dot L=-v_G-v_U\leq 0.
```

For any compact recurrent set, this monotone coordinate must be constant. Then both nonnegative rates vanish, and all five derivatives vanish. Thus the closed published model has no nonstationary recurrent attractor. For one enzyme there is only one reaction extent after imposing the conserved totals.

The acid/base calculation cannot supply a hidden dynamical degree of freedom. With `H = [H+] > 0`, the paper's charge-balance equation has the form

```math
F(H)=H-\frac{K_w}{H}
-\sum_a\frac{c_aK_a}{K_a+H}
+\sum_b\frac{c_bK_b}{K_b+K_w/H}=p_0.
```

For nonnegative concentrations and positive dissociation constants,

```math
F'(H)=1+\frac{K_w}{H^2}
+\sum_a\frac{c_aK_a}{(K_a+H)^2}
+\sum_b\frac{c_bK_bK_w}{H^2(K_b+K_w/H)^2}>0.
```

Together with the limits at zero and infinity, this proves a unique positive algebraic `H` for each chemical state. A nearly constant pH during ongoing fuel consumption is not an equilibrium of the complete concentration state. The proof also survives a prescribed nonnegative startup-rate correction; such a correction does not create autonomous recurrence.

## 2. The 2025 geometric urease model is planar

The primary preprint [2], version 1, studies the two variables `s,h`. Its system is

```math
\dot s=-r(h)s+K_s,\qquad
\dot h=-q(s,h)+K_h(1-h),
```

```math
r(h)=\left[\frac{\beta\varepsilon_1}{h}+1+
\frac{\beta h}{\varepsilon_1}\right]^{-1},
\quad
q(s,h)=\frac{-v(h)+\sqrt{v(h)^2+4K\varepsilon_2^{-1}r(h)h^2s}}{2},
```

```math
v(h)=\alpha K\varepsilon_2^{-1}h^2-K_h(1-h).
```

The authors establish an attracting periodic orbit under their instability and trapping conditions using Poincaré–Bendixson. Multiple time scales and sharp pH changes do not increase the state dimension. This exclusion applies to compact regular dynamics in the positive planar domain. The reactor has an external reservoir, so the finite-fuel argument of Section 1 is not applicable here.

The label `nlin.CD` on an arXiv classification is not evidence that the article demonstrates chaos. Only version 1, dated 2025-08-10, was listed on arXiv at this audit.

## 3. Milster's full three-state model is competitive

Milster et al. [3], Eqs. 1–2, give

```math
\dot g=k_D^g g_0e^{-B\phi}-(k_D^g+k_R)g,
\qquad
\dot h=k_Rg-k_D^h(h-h_0),
```

```math
\tau_\phi\dot\phi=(\phi-\phi_c)
-a(\phi-\phi_c)^3+m(h-h_c).
```

The source's subsequent Eq. 4 is a planar approximation with FitzHugh–Nagumo-like dynamics. A stronger obstruction follows directly for the unreduced equations when `B,m,kR,kDg,g0,tau_phi > 0`.

Set `x = (g,-h,phi)`. The only nonzero off-diagonal derivatives become

```math
\frac{\partial\dot x_1}{\partial x_3}
=-Bk_D^g g_0e^{-Bx_3}<0,\qquad
\frac{\partial\dot x_2}{\partial x_1}=-k_R<0,\qquad
\frac{\partial\dot x_3}{\partial x_2}=-m/\tau_\phi<0.
```

All other off-diagonal entries are zero. This is a competitive smooth vector field on the convex domain `R^3`. Hirsch's Theorem A / 3.1 [4] embeds its limit sets, preserving the flow, into a Lipschitz flow in one lower dimension. Therefore any compact limit dynamics are planar, including sets containing equilibria. This rules out a bounded chaotic limit attractor for these equations and signs. It does not prove every orbit is bounded or physically meaningful: the Landau polynomial alone does not enforce `0 <= phi <= 1` for arbitrary parameters.

The related monotone cyclic feedback theorem [5] supplies the familiar equilibrium-free periodic-orbit conclusion. We use the stronger competitive embedding here so that an equilibrium-containing limit set is not silently omitted. This structural application is our audit observation; no claim of a new mathematical theorem is made.

The four-state two-enzyme Shinbo model [8] being investigated separately is **not** excluded by this calculation. Its additional substrate/product loop and pH-dependent activity cannot be replaced by Milster's three-variable loop.

## 4. Next candidate: original four-species urease network

The concrete source is Straube, Winkelmann, and Höfling [6], Section 2 of arXiv:2212.14503v1, published in J. Phys. Chem. B in 2023. State order is `(S,H,P,Q)`: urea, hydrogen ions, ammonia, and ammonium concentrations, all in mol/L. Time is in seconds. The source omits carbon dioxide chemistry in this minimal network.

Define

```math
v(S,H)=\frac{v_{\max}S}{K_M+S}
\frac{1}{1+H/K_{E1}+K_{E2}/H},\qquad
j=k_2PH-k_{2r}Q.
```

The exact four equations are

```math
\begin{aligned}
\dot S&=-v+k_S(S_{\rm ext}-S),\\
\dot H&=-j+k_H(H_{\rm ext}-H),\\
\dot P&=2v-j-kP,\\
\dot Q&=j-kQ.
\end{aligned}
```

The positive-domain extension at `H=0` is obtained without a singular division by writing the activity as `H / (H + H*H/KE1 + KE2)`.

| Parameter | Source reference value | Units |
|---|---:|---|
| `k2` | `4.3e10` | L mol^-1 s^-1 |
| `k2r` | `24` | s^-1 |
| `KM` | `3e-3` | mol/L |
| `KE1` | `5e-6` | mol/L |
| `KE2` | `2e-9` | mol/L |
| `kS` | `1.4e-3` | s^-1 |
| `kH` | `9e-3` | s^-1 |
| `k` | `kS` | s^-1 |
| `vmax` | `1.85e-4` | mol L^-1 s^-1 |
| `Sext` | `3.8e-4` | mol/L |
| `Hext` | `1.3e-4` | mol/L |
| `S(0)` | `5e-5` | mol/L |
| `H(0)` | `1e-5` | mol/L |

The inspected paragraph does not specify `P(0),Q(0)`. A later numerical control must state its own choices, unless the earlier 2021 source or supplementary material resolves them. Do not label an invented complete initial state as the source's state. The source shows sustained regular oscillation at the tabulated parameters; it is the first qualitative control to reproduce, not a chaos candidate by itself.

### Exact bounds and an exact dimension reduction

This is an open, continuously supplied reactor. Unlike the 2026 closed model, its substrate need not decrease. Nevertheless, it has useful exact identities.

For nonnegative initial concentrations, the vector field points into the nonnegative orthant at every boundary face. Write

```math
N=2S+P+Q.
```

Direct stoichiometric cancellation gives

```math
\dot N=2k_SS_{\rm ext}-2k_SS-k(P+Q)
\leq2k_SS_{\rm ext}-\min(k_S,k)N.
```

Thus `S,P,Q` are bounded for positive outflow rates. Also

```math
\dot H\leq k_{2r}N+k_HH_{\rm ext}-k_HH,
```

which bounds `H`. These statements establish boundedness of nonnegative solutions of this mathematical model, not experimental validity across arbitrary rate choices.

At the source setting `k=kS`,

```math
\dot N=k_S(2S_{\rm ext}-N),\qquad
N(t)=2S_{\rm ext}+[N(0)-2S_{\rm ext}]e^{-k_St}.
```

Every compact recurrent set therefore lies in the exact invariant hyperplane

```math
Q=2S_{\rm ext}-2S-P.
```

The long-time dynamics are exactly three-dimensional, on `S,H,P >= 0` and `2S+P <= 2Sext`. The scalar transverse coordinate has decay rate `-kS`; this is an exact normal contraction to check independently of any Lyapunov calculation. The full variational equations are block triangular in coordinates `(S,H,P,N)`, with this scalar normal block.

The divergence of the original four-state vector field is

```math
\nabla\cdot f=-v_S-k_S-k_H-k_{2r}-2k-k_2(P+H)<0.
```

Volume contraction does not settle whether its remaining three-dimensional dynamics can be chaotic.

### Why the published planar approximation is not an exclusion proof

Section 3.2 of [6] explicitly introduces `P(t) H(t) approximately constant` using a scale-separation argument. This is not an exact first integral of the three-state mass-action reduction. In that reduction,

```text
S' = -v + kS (Sext-S)
H' = -a P H + kH (Hext-H)
P' = 2v - a P H - k P,
a = k*k2/(k2r+k).
```

At the reference parameters and `(S,H,P)=(5e-5,1e-5,1e-8)`, direct substitution yields `(PH)' = 1.771719621009859e-11 mol^2 L^-2 s^-1`, not zero. This is an algebraic counterexample to a global first integral, not criticism of the reduction's accuracy near the source limit cycle. Therefore neither the 2023 approximate planar model nor the 2025 planar theorem excludes chaos in the original network for all parameters.

The activity derivative changes sign at `H=sqrt(KE1*KE2)=1e-7 mol/L`. The original network is not a single fixed-sign cyclic loop, so the Milster exclusion must not be transferred to it.

### Proposed next numerical block

1. Reproduce the reference regular oscillation with an implicit CPU solver and an analytic Jacobian. Supply and label the missing product initial concentrations explicitly.
2. Verify nonnegativity, the exact `N(t)` identity, and the divergence against variational traces. Compare four-state integration with exact three-state integration on the invariant hyperplane.
3. Vary experimentally motivated reservoir concentrations, enzyme-dependent `vmax`, and membrane rates. Keep protonation constants fixed in the primary survey. Changing them broadly would investigate a different chemical regime.
4. Only then consider GPU screening. Large `k2` makes an explicit fixed-step solver potentially impractical; an implicit or properly validated asymptotic method is needed. A convenient planar approximation would defeat the purpose of this search.
5. Any positive finite-time exponent still requires long continuation, several starts, independent stiff CPU solvers, convergence checks, a Poincaré section, and physical projections.

## 5. Literature audit and limits of novelty checking

Primary sources were read through official ACS/PMC XML and arXiv full text. The 2023 article's final ACS DOI returned HTTP 403 during this block; the inspected full equation version is explicitly the public preprint. The arXiv records listed a single version for both the 2023 and 2025 urease papers. The 2026 source references the earlier urease line, but supplies no deterministic-chaos result for the exact four-species network.

Searches included `urease deterministic chaos`, `urea-urease chaotic`, `urease strange attractor`, the exact reduced-model title plus `chaos`, and the Milster title plus `chaos`. No exact-model chaos demonstration was identified in this limited search. This is not evidence that none exists. A generic paper published in the journal *Chaos* is not a report of chaotic dynamics.

The Europe PMC citing-article endpoint returned one indexed citation for PMID 36976257: Mahato et al. 2025 [7], whose primary text concerns local inhibition and distal activation after a chemical injection, not a demonstration of chaos in the isolated vesicle model. This index is visibly incomplete: the relevant 2025 arXiv preprint is absent. The same endpoint returned three citations for Milster (PMID 37978612), with topics electron-beam nanoparticle patterning, tracer motion in polymers, and transient soft matter fueling. No exhaustive citation claim is made.

Leathard et al. [9] is a useful separate lead: its primary abstract and published text identify a seven-species esterase compartment with known transport-driven pH oscillations and supplementary XPP equations. Those equations have not been screened here. It must be audited before entering a numerical queue.

Machine-readable decisions, equations, source hashes, parameter values, and the small algebra checks are in [enzyme-feedback-screening.json](../results/candidate-audits/enzyme-feedback-screening.json). Source PDFs/HTML/XML are cached locally under ignored `outputs/`; they are not relicensed under this project's MIT license.

## References

1. N. Göth and J. Dzubiella, “Chemical Feedback and Collective pH-Dynamics in Enzyme-Powered Artificial Cells,” *J. Phys. Chem. B* **130**, 10269–10278 (2026). [DOI:10.1021/acs.jpcb.6c03926](https://doi.org/10.1021/acs.jpcb.6c03926). [Official full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13637022/fullTextXML).
2. A. V. Straube, G. Olicón Méndez, S. Winkelmann, F. Höfling, and M. Engel, “Unfolding the geometric structure and multiple timescales of the urea–urease pH oscillator” (2025). [arXiv:2508.07275v1](https://arxiv.org/html/2508.07275v1).
3. S. Milster, A. Darwish, N. Göth, and J. Dzubiella, “Synergistic chemomechanical dynamics of feedback-controlled microreactors,” *Phys. Rev. E* **108**, L042601 (2023). [DOI:10.1103/PhysRevE.108.L042601](https://doi.org/10.1103/PhysRevE.108.L042601). [Inspected preprint](https://arxiv.org/html/2306.05181v1).
4. M. W. Hirsch, “Systems of Differential Equations Which Are Competitive or Cooperative: I. Limit Sets,” *SIAM J. Math. Anal.* **13**, 167–179 (1982), Theorem A / 3.1. [DOI:10.1137/0513013](https://doi.org/10.1137/0513013). [Author-uploaded primary full text](https://www.researchgate.net/publication/247695745_Systems_of_Differential_Equations_Which_Are_Competitive_or_Cooperative_I_Limit_Sets).
5. J. Mallet-Paret and H. L. Smith, “The Poincaré–Bendixson theorem for monotone cyclic feedback systems,” *J. Dyn. Differ. Equ.* **2**, 367–421 (1990). [DOI:10.1007/BF01054041](https://doi.org/10.1007/BF01054041). The precise compact, equilibrium-free statement and convex-projection condition were also checked in L. A. Sánchez, [“Global asymptotic stability of the Goodwin system with repression,” Section 2](https://www.ugr.es/~ecuadif/files/cyclicneg.pdf).
6. A. V. Straube, S. Winkelmann, and F. Höfling, “Accurate reduced models for the pH oscillations in the urea–urease reaction confined to giant lipid vesicles,” *J. Phys. Chem. B* **127**, 2955–2967 (2023). [DOI:10.1021/acs.jpcb.2c09092](https://doi.org/10.1021/acs.jpcb.2c09092). [Inspected arXiv:2212.14503v1](https://arxiv.org/html/2212.14503v1).
7. R. R. Mahato, N. Adikessavane, N. K. Ettikkan, and S. Maiti, “A chemoenzymatic reaction-diffusion system with local inhibition and distal activation,” *Nat. Commun.* **16**, 9885 (2025). [DOI:10.1038/s41467-025-65605-2](https://www.nature.com/articles/s41467-025-65605-2).
8. Shinbo, “Theoretical Model of Chemomechanical Self-Oscillation in a Two-Enzyme Core–Shell Hydrogel With pH-Gated Microcapsules,” *ChemSystemsChem* (2026). [DOI:10.1002/syst.70050](https://doi.org/10.1002/syst.70050). Investigated separately in this project.
9. A. S. Leathard, P. A. Beales, and A. F. Taylor, “Design of oscillatory dynamics in numerical simulations of compartment-based enzyme systems,” *Chaos* **33**, 123128 (2023). [DOI:10.1063/5.0180256](https://doi.org/10.1063/5.0180256). [Institutional published version](https://eprints.whiterose.ac.uk/id/eprint/207599/).
