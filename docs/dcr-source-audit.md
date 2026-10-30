# DCR source audit: contact boundaries and optional callbacks

Audit: 30 October 2026 UTC. This follow-up to [screening VI](late-candidate-screening-6.md) reads the final paper, its two supplementary PDFs, Zenodo version 1.0.0, and pinned Julia source [1–3]. **The literal strict contact field has a conditional boundary-existence problem.** This does not invalidate the reported finite-step simulations, establish a global failure, or exclude chaos in a properly specified contact model. It means that event conventions must be stated before claiming a unique continuous-time attractor.

No author code was executed and no simulations were run in this audit. Independent numerical checks are recorded separately by the project.

## What the source actually implements

The pinned file is `DCR_Model_available.jl` at commit `796c5e39468b06fd6fc2c378437d618593948e98`, blob `85ccf19ab244815929045fb646814f54b56279e8`. Source line numbers below refer to this file. Both recovered code versions use the same contact and optional callback rules. The default `events=false` call uses IDA with a discontinuous DAE residual and **does not pass the contact callback**. The presence of a callback definition in the file therefore does not mean that default simulations locate contact events.

| Source element | Observed implementation |
|---|---|
| Angular contact window, lines 69–82 | Strictly between the two angular endpoints, after phase wrapping |
| Contact test, lines 117–120 | Ordinary wrapped-tooth distance strictly less than 0.02 |
| Contact phase velocity | Backbone velocity divided by the absolute phase cosine |
| Noncontact phase velocity | Exactly $2\pi$ |
| Default solver, lines 337–340 | IDA/LapackDense, maximum step $10^{-4}$, relative and absolute tolerances $10^{-7}$, BrownFullBasicInit |
| Optional callback | Only passed when `events=true` |

The main article discusses contact detection and numerical convergence, but neither the examined main text, the two SI PDFs, nor these code files specifies a Filippov convention, persistent contact mode, complementarity rule, position reset, grazing rule, or ordering rule for simultaneous contacts. This is a bounded source finding, not a statement about every version the authors may possess. No package manifest fixes the installed SciML/SUNDIALS versions.

## Contact captures a gap rather than reducing it

Use the **author-code** rotor position with the minus sine; the printed plus-sine discrepancy is documented in screening VI. For a particular tooth and rotor, away from tooth wrapping, write

```math
g=t_k-r_j,
\qquad r_j=j\ell_r-\sin\theta_j,
\qquad \dot t_k=v=\dot x.
```

Within the angular window, $\cos\theta_j\lt0$. Therefore

```math
\dot g=v+\cos\theta_j\dot\theta_j
=\begin{cases}
v-2\pi|\cos\theta_j|,&\text{free},\\
0,&\text{in contact}.
\end{cases}
```

The contact branch preserves whatever interior gap was captured. It does not dynamically relax that gap to zero. An existing contact can end at an angular boundary or a tooth-wrap discontinuity; it cannot leave the distance band by smooth gap evolution while the same tooth remains active. Changes in other contacts alter the common velocity but preserve this identity.

### A local obstruction at a strict distance boundary

Consider an isolated first entry at $g=d$, with $d=0.02$, away from angular boundaries, tooth wrapping, and other switching surfaces. Assume the free normal speed is strictly negative and remains bounded away from zero in a neighborhood. At equality the strict test selects **free motion**; just inside it selects the tangent contact branch. For $h=g-d$, the local rule has

```math
\dot h\le-\delta\quad\text{when }h\ge0,
\qquad \dot h=0\quad\text{when }h\lt0,
\qquad \delta\gt0.
```

There is no absolutely continuous solution starting at $h=0$ that satisfies this rule almost everywhere. Such a solution would be nonincreasing and could not become positive. It also cannot enter a negative interval: there its derivative is zero, making it constant, inconsistent with continuity from the first zero crossing. Remaining at zero contradicts the prescribed negative boundary derivative. The scalar example $\dot h=-1$ for $h\ge0$, and zero for $h\lt0$, has the same obstruction.

This is a **conditional local derivation**, not an assertion that every initial condition reaches this configuration. The opposite distance boundary has the analogous inward-entry issue. It assumes a locally regular algebraic velocity branch and excludes simultaneous and grazing events.

A finite step can move a point slightly inside the contact band and then preserve that penetrated gap. A boundary-tracking limit would need a closed-contact, hybrid-mode, or differential-inclusion convention. Choosing one can be reasonable mathematical work, but it is an explicit model interpretation. Reducing solver tolerances alone does not define that interpretation. A finite distance band does not remove the issue at its edges.

## The optional callback is a distinct numerical intervention

`condition1_safe`, lines 183–224, evaluates the nearest tooth distance minus **0.023**, only while the rotor is in its angular window and its current derivative exceeds $0.9(2\pi)$. Otherwise it returns 1. Its guard depends on the integrator's derivative and on discontinuous gates. Consequently it is not globally a continuous function of physical state, despite being supplied to `VectorContinuousCallback` with root finding enabled.

The detection shell lies outside the residual's 0.02 contact band. At a detected event, `affect_safe!`, lines 227–260, assigns

```math
\omega_j^{\mathrm{assigned}}
=\frac12\left|\omega_j^{\mathrm{old}}
-\frac{v^{\mathrm{old}}}{|\cos\theta_j|}\right|.
```

It changes only `integrator.du[j+1]`. It does not set the derivative to the contact velocity, reset a phase or displacement, set a contact flag, recompute the backbone force balance, or explicitly request consistent DAE reinitialization. For example, free speed $2\pi$ and zero backbone speed produce an assigned speed $\pi$, although a rotor at distance 0.023 remains free under the residual and requires $2\pi$. Thus the assigned derivative is **not generally a consistent DAE derivative**. What a particular library version subsequently accepts or overwrites has not been established here.

The assignment is nonnegative even though the stated contact law permits negative phase velocity during spring-driven backtracking. It should not be described as a derived restitution law or as an exact enforcement of the mechanical constraint. The source has no separately specified entry and exit callback laws. The file's comment describes its intended contact-detection purpose; the executable expressions above are the evidence for what it actually assigns.

## Consequences for tangent dynamics

After specifying a valid hybrid convention, isolated transverse events with identity state reset have the usual saltation update [4]

```math
\Xi=I+\frac{(f^+-f^-)n^{\mathsf T}}
{n^{\mathsf T}f^-},
\qquad n=\nabla g.
```

For capture into the gap-preserving branch, $n^{\mathsf T}f^+=0$, so $n^{\mathsf T}\Xi=0$: the update removes the normal variation to first order. Its singularity is expected for capture and is not itself a numerical defect or a proof of contraction of all directions. This formula requires a declared post-event mode; it does not repair the literal strict-field existence problem. At grazing, the denominator vanishes. Simultaneous events require separate sequence-sensitive treatment, not arbitrary multiplication of single-event updates. The within-branch Jacobian alone is insufficient for a Lyapunov spectrum.

For the sigmoid force, even within a fixed contact set one should check the scalar algebraic root. Its derivative can have a small positive part during backdrive. The independent [structural audit](dcr-structure.md) gives the derivative bound and uses the source spacing to sharpen the maximum contact count. That algebraic issue is separate from the distance-boundary convention.

## Priority and remaining gaps

The final paper already reports deterministic irregularity and strong initial-condition sensitivity. Its full main text and the recovered SI do not supply an event-aware Lyapunov spectrum. The 2024 discrete-sarcomere follow-up [5] remains a **full-text gap**: its metadata and abstract were recovered, while the public full-text and JATS endpoints were inaccessible in this bounded check. We cannot exclude an earlier chaos diagnosis hidden in that work. The author's thesis search record also points to chaos as future investigation, but its full PDF was not recovered and is not treated as a completed thesis audit.

The closest explicit-chaos prior remains Blanch-Mercader, Orlandi, and Casademunt's interacting quasideterministic ratchets [6]. That model has stochastic alternation of deterministic circle maps; no exact equivalence to this deterministic DCR network has been established. Generic novelty claims about motor-driven chaos would therefore be unjustified even if a well-defined DCR implementation later yields positive exponents.

**Decision:** keep short literal-code diagnostics as source-behavior checks. Before a mass search or attractor claim, record a defensible event convention, show that numerical contact handling converges to it, distinguish the optional callback variant, and validate event-aware sensitivity. No new-attractor claim follows from this audit.

## References

1. B. Warmington, J. Rossiter, H. Bloomfield-Gadêlha. *A discontinuously coupled network of phase oscillators replicates cooperation between motors in actomyosin systems*. J. R. Soc. Interface **23**, 20250438 (2026). [DOI](https://doi.org/10.1098/rsif.2025.0438), [public final PDF](https://aura.abdn.ac.uk/bitstreams/2e76dc81-a0f5-4fb5-aa9c-7ca1016f8730/download), [official supplements](https://doi.org/10.6084/m9.figshare.c.8319421).
2. *Discontinuously coupled rotor model*, version 1.0.0. [Zenodo](https://doi.org/10.5281/zenodo.18196853).
3. Author implementation, [pinned Julia file](https://github.com/polymaths-lab/actomyosin-kuramoto-rotor-model/blob/796c5e39468b06fd6fc2c378437d618593948e98/DCR_Model_available.jl).
4. N. J. Kong, J. J. Payne, J. Zhu, A. M. Johnson. *Saltation matrices: The essential tool for linearizing hybrid dynamical systems*. Proc. IEEE **112**, 585–608 (2024). [DOI](https://doi.org/10.1109/JPROC.2024.3440211), [author full text](https://arxiv.org/html/2306.06862v3), section III-A.
5. B. Warmington, J. Rossiter, H. Bloomfield-Gadêlha. *A discrete sarcomere model reveals internal unbinding dynamics and conditions for regularity of spontaneous oscillatory contractions*. [2024 preprint](https://doi.org/10.1101/2024.07.03.601948); metadata/abstract only.
6. C. Blanch-Mercader, J. G. Orlandi, J. Casademunt. *Unraveling the hidden complexity of quasideterministic ratchets: Random walks, graphs, and circle maps*. Phys. Rev. E **101**, 012203 (2020). [DOI](https://doi.org/10.1103/PhysRevE.101.012203).
