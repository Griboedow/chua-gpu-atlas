# Late screening VI: a discontinuously coupled motor network

Audit: 30 October 2026 UTC. **A complete executable author model has been recovered:** Warmington, Rossiter, and Bloomfield-Gadêlha's discontinuously coupled rotor (DCR) network [1–3]. It is a substantially better reproduction target than an abstract-only oscillator. However, its published controls include regular and irregular oscillations rather than a certified periodic orbit, and the printed geometry differs from the author implementation. No simulation or chaos diagnosis was performed here.

## Physical mechanism and state space

Rotors driven at a natural angular speed intermittently push teeth on a common translating backbone. Viscous resistance and a Hookean load provide mechanical feedback. Contact causes motor recruitment and loss; collective release allows the spring to pull the backbone back. This is an autonomous mechanical model with maintained motor power, not an added artificial feedback state or a finite-fuel chemical system.

The state is $(x,\theta_1,\ldots,\theta_N)$: **$N+1$ independent dynamic coordinates**, with a state-dependent contact set. The published physiological simulations use $N=100$, hence 101 coordinates. Phases should be understood modulo $2\pi$ when discussing compactness. A three-dimensional picture would only be a projection. Reducing to two or three motors is a parameter study of the network family, not reproduction of the 100-motor controls; source synchronization depends on network size and geometry.

## Authoritative sources recovered

The final 18-page 2026 article was retrieved from the University of Aberdeen's public repository, along with both official supplementary PDFs. The source code is a single Julia file on Zenodo, DOI 10.5281/zenodo.18196853. Its public GitHub counterpart was inspected at commit `796c5e39468b06fd6fc2c378437d618593948e98`. The GitHub repository declares CC0; no author code has been executed, copied into tracked implementation files, or relicensed as our work.

The final article calls the irregularity deterministic and sensitive to initial conditions in section 4.2. It does not provide a Lyapunov spectrum or label a computed strange attractor. Therefore, this would be an investigation of already reported irregular dynamics, not discovery of an unsuspected irregular regime.

## Printed equations versus executable variant

The final paper's equations 2.2–2.8 give the dimensionless force balance

```math
\alpha\dot x+\beta x+\gamma=
\sum_{j=1}^{N}F(\dot\theta_j)H(\theta_j)|\cos\theta_j|,
```

```math
F(\omega)=\begin{cases}1-\omega/(2\pi),&\omega\ge0,\\1,&\omega\lt0,\end{cases}
\qquad
\dot\theta_j=\begin{cases}\dot x/|\cos\theta_j|,&j\text{ in contact},\\2\pi,&\text{otherwise}.\end{cases}
```

Here $H=1$ for a wrapped angle strictly between $\pi-\arcsin c$ and $\pi+\arcsin c$, and zero otherwise. Contact additionally requires proximity to a tooth. With zero-based indices, the **author-code geometry** is

```math
r_j=j\ell_r-\sin\theta_j,\qquad
t_k=(k\ell_t+x)\bmod (N_t\ell_t),\qquad
N_t=\left\lceil\frac{N\ell_r}{\ell_t}\right\rceil,
```

with contact when $H(\theta_j)=1$ and at least one $|t_k-r_j|\lt d$. The code uses ordinary absolute separation after wrapping tooth positions; replacing it by shortest periodic distance would change the published implementation near boundaries.

**Important source distinction:** printed equation 2.7 has **plus** $\sin\theta_j$, whereas both recovered code versions use **minus**. This is visually confirmed on PDF page 4. With the source contact interval, $\cos\theta_j\lt0$, and its positive contact-speed rule, the code gives $\dot r_j=\dot x$, while the printed plus sign gives $\dot r_j=-\dot x$. The code geometry therefore preserves rotor–tooth contact under its own kinematic rule. A reproduction should name the chosen author-code variant and record this difference rather than silently calling both formulations identical. A common shift of indexing does not resolve this sign discrepancy.

The code also implements the following smooth saturation at `switch=0`:

```math
F_s(\omega)=1-\frac{\omega}{4\pi}\left[1+\tanh(40\omega)\right].
```

This is not exactly the piecewise function above, although the paper explicitly describes sigmoid smoothing in section 6.1. Its constant-load implementation is $\gamma[1+\tanh(2(x-1.3))]/2$, rather than a globally constant $\gamma$. The difference disappears for $\gamma=0$, the useful spring-oscillation case. Turning on the code's `switch` changes backdrivability and is explicitly unsuitable for its oscillatory mode. Optional event callbacks also alter derivative initialization; their behavior must be reported separately from the default `events=false` branch.

## Complete starting specification

Use the author-code variant with the following source-grounded control settings:

| Quantity | Setting |
|---|---:|
| Motor count $N$ | 100 |
| Contact angular parameter $c$ | 0.14 |
| Rotor spacing $\ell_r$ | 1.204 |
| Tooth spacing $\ell_t$ | 1.092 |
| Contact distance $d$ | 0.02 (author code) |
| Initial backbone $x(0)$ | 0 |
| Initial phases | Independent uniform phases on $[0,2\pi)$; save the actual vector |
| Constant load $\gamma$ | 0 |
| Force law | Author sigmoid with coefficient 40; `switch=0` |
| Contact events | Default code branch off; independently resolve switching in a later check |
| Source Fig. 5 parameter pairs $(\alpha,\beta)$ | $(1,1),(1,4),(4,1),(4,4)$ |

The geometry is from the Fig. 5 caption and matches `8.6*c` and `7.8*c` in the code. It differs slightly from rounded physiological values 1.16 and 1.05 discussed elsewhere in the paper; the figure-specific numbers should be used for reproduction. Source runs last 10,000 time units, except $(4,1)$, which lasts 30,000. Its largest $\alpha$, largest $\beta$ control has a narrow period distribution near 0.55 s, using the source dimensional frequency $f_0=30$ Hz. This is **not evidence of an exactly periodic full-state orbit**.

The supplied file's default values are instead $\alpha=2$, $\beta=0$, $\gamma=2$, duration 400 on Zenodo; GitHub uses $\gamma=3$, duration 100. Those defaults describe unloaded/sustained translation rather than the spring-oscillation controls. Initial phases are unseeded in the supplied file. Reproducible independent work must save phases and distinguish defaults from Fig. 5 settings.

The source uses IDA/SUNDIALS and a maximum step of $10^{-4}$ dimensionless time. Printed tolerances are relative $10^{-6}$, absolute $10^{-8}$; recovered code uses $10^{-7}$ for both. The authors warn that missed contacts create spurious disorder and that small drag/high load can defeat contact detection. A normal smooth RK4/Lyapunov pipeline cannot be applied without resolving these switching events.

## Cheap mathematical preparation

For the **sharp saturation law** in equation 2.3, a fixed contact set $C$ allows the DAE velocity to be eliminated analytically. Let $n=|C|$, $S=\sum_{j\in C}|\cos\theta_j|$, and $A=S-\beta x-\gamma$. Then

```math
\dot x=\begin{cases}
A/(\alpha+n/(2\pi)),&A\ge0,\\
A/\alpha,&A\lt0.
\end{cases}
```

Unconnected rotors have exactly zero force; connected rotors share the sign of $\dot x$. For $\alpha\gt0$, this gives a unique algebraic velocity on each contact branch. No extra state or closure is needed. This formula is **not** the exact solver for the sigmoid-smoothed code, whose scalar velocity equation must instead be solved consistently.

For the sharp law with $\beta\gt0$ and $\gamma=0$, $0\le S\le N$ implies that the interval $0\le x\le N/\beta$ is invariant: the velocity points inward at both boundaries. With wrapped phases, this provides a compact candidate region, conditional on a well-defined event solution. It does not prove chaos, attraction, or uniqueness through every grazing/simultaneous-contact event. The contact interval also implies $|\cos\theta_j|\ge\sqrt{1-c^2}$ during contact, so its velocity denominator is separated from zero for $0\lt c\lt1$.

The full hybrid network is not a planar oscillator, and contact switching cannot be discarded in a dimension reduction. Any Lyapunov calculation must include event-time sensitivity (or an independently converged finite-perturbation method); evolving only the within-branch smooth Jacobian can miss the dominant stretching. Before geometry claims, test contact-resolution convergence, persistence over long times, several starting phase vectors, and physical versus phase-quotient coordinates.

## Independent reproduction plan

1. Implement the pinned author-code contact geometry and smoothed force law independently; retain the supplied code only as a reference. Use $\gamma=0$, `switch=0`, the four Fig. 5 tuples, and explicit saved phase vectors.
2. First verify instantaneous contact kinematics and the scalar force residual on fixed contact branches. Compare a consistent implicit DAE solver with an independent event-resolved implementation. Use the sharp-law branch formula only as a separately labeled limit/control.
3. Reproduce the regular-versus-irregular source behavior, release cascades, contact fraction, and period distributions before any parameter survey. No claim of matching the exact source trajectory is possible without its unreported random phases.
4. Repeat with smaller maximum steps and tighter contact resolution. Check convergence of event sequences/statistics and distinguish the finite contact width from integration tolerance. Preserve the physical nonzero contact width until a separately labeled width-sensitivity study.
5. Only after these controls, estimate event-aware exponents and full-state recurrence. If plotting three coordinates, label the result a projection of the 101-state network. Do not call a narrow period histogram an attracting periodic orbit.

## Prior art and forward checks

The final paper's own closest complexity reference is Blanch-Mercader, Orlandi, and Casademunt (2020) [4], which explicitly studies piecewise chaos in interacting quasideterministic ratchets. Its model combines hard-core interaction with stochastic alternation of deterministic circle maps. That is relevant prior art, but no exact equivalence to the deterministic DCR contact network is established here. Familiar Kuramoto terminology likewise does not identify the DCR network with the usual smooth sinusoidal phase-coupling equations.

The OpenAlex forward query for the final paper returned two records: an author follow-up preprint on a discrete sarcomere model [5] and a news article. The preprint's metadata and abstract were checked, but its full equations were not recovered within this bounded block. The 2023 rotor preprint has at least two versions; the final 2026 paper, not its early title/abstract, was used for this audit. Exact-title searches with chaos/chaotic/Lyapunov did not return a verified exact-DCR chaos computation. Such a bounded search does not establish that none exists. The full 2024 follow-up, older motor models, and unpublished analyses remain priority gaps.

## Decision

The **author-code DCR network is ready for an independently implemented source-behavior check**, with a clear distinction between sharp printed equations and smoothed executable equations. It does not yet meet the stronger requirement of an independently verified periodic full-state control or a novelty-cleared chaos target. The narrow useful question is whether the reported deterministic irregularity yields robust positive exponents after contact events and numerical resolution are treated correctly. A distinctive physical mechanism is present; a new attractor has not been established.

## References

1. B. Warmington, J. Rossiter, H. Bloomfield-Gadêlha. *A discontinuously coupled network of phase oscillators replicates cooperation between motors in actomyosin systems*. J. R. Soc. Interface **23**, 20250438 (2026). [DOI](https://doi.org/10.1098/rsif.2025.0438), [public final PDF](https://aura.abdn.ac.uk/bitstreams/2e76dc81-a0f5-4fb5-aa9c-7ca1016f8730/download), [official supplements](https://doi.org/10.6084/m9.figshare.c.8319421).
2. B. Warmington. *Discontinuously coupled rotor model*, version 1.0.0. [Zenodo](https://doi.org/10.5281/zenodo.18196853).
3. Author repository, [pinned source](https://github.com/polymaths-lab/actomyosin-kuramoto-rotor-model/blob/796c5e39468b06fd6fc2c378437d618593948e98/DCR_Model_available.jl).
4. C. Blanch-Mercader, J. G. Orlandi, J. Casademunt. *Unraveling the hidden complexity of quasideterministic ratchets: Random walks, graphs, and circle maps*. Phys. Rev. E **101**, 012203 (2020). [Primary paper](https://doi.org/10.1103/PhysRevE.101.012203).
5. B. Warmington, J. Rossiter, H. Bloomfield-Gadêlha. *A discrete sarcomere model reveals internal unbinding dynamics and conditions for regularity of spontaneous oscillatory contractions*. [2024 preprint](https://doi.org/10.1101/2024.07.03.601948). Metadata-only follow-up check in this block.
