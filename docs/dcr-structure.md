# Structure of the discontinuously coupled rotor source model

This note audits the **author-code variant** of Warmington, Rossiter, and Bloomfield-Gadêlha's discontinuously coupled rotor (DCR) network [1–3]. It fixes `switch=0`, constant-load parameter `gamma=0`, and the Figure 5 geometry: 100 rotors, contact parameter $c=0.14$, contact distance $d=0.02$, rotor spacing $\ell_r=1.204$, and tooth spacing $\ell_t=1.092$. The source drag/load pairs are $(\alpha,\beta)=(1,1),(1,4),(4,1),(4,4)$. The pinned source is commit `796c5e39468b06fd6fc2c378437d618593948e98`; its file hash is recorded in the [audit results](../results/dcr/structure-audit.json).

The useful results are a unique algebraic velocity for this geometry and an inexpensive exact description within a fixed contact branch. **They do not establish a unique global event solution.** Two explicit boundary states below explain why capture and angular switching need a stated interpretation before event-based chaos diagnostics.

These are independent algebraic derivations and finite checks, not a novelty claim or a trajectory reproduction. The source paper already reports deterministic irregular oscillations. See [the source-screening note](late-candidate-screening-6.md) for the distinction between the printed and executable models, source behavior, and prior-art limits.

## Exact tooth geometry with ordinary distance

Use zero-based rotor and tooth indices. With $N_t=111$ and $L=N_t\ell_t=121.212$, the source positions are

```math
r_j=j\ell_r-\sin\theta_j,
\qquad
t_k=(k\ell_t+x)\bmod L.
```

The minus sign is the code's convention. It differs from the printed plus sign in equation 2.7 of [1]. Sorting the tooth set gives exactly

```math
\{t_k\}_{k=0}^{N_t-1}
=\{\delta+k\ell_t\}_{k=0}^{N_t-1},
\qquad \delta=x\bmod\ell_t.
```

Thus the nearest tooth in **ordinary absolute distance** has index

```math
k_* = \min\left(N_t-1,
\max\left(0,\mathrm{round}\frac{r_j-\delta}{\ell_t}\right)\right).
```

The tie convention is irrelevant to the contact predicate because $d\ll\ell_t/2$. Replacing ordinary distance by shortest periodic distance changes the source model: at $x=1.082$ and $\theta_0=\pi-\arcsin(0.01)$, the ordinary nearest distance is 1.092, while the periodic distance is zero. The source rotor is **not** in contact.

At $x$ equal to an integer multiple of $\ell_t$, the sorted indices change. Interior physical teeth remain continuous and must simply be relabeled; only a tooth crossing the finite domain boundary disappears or reappears at the opposite end. Automatically releasing every contact at a change of $\delta$ would be wrong.

## How many motors can contact simultaneously?

A contacting rotor has an angle strictly between $\pi-\arcsin c$ and $\pi+\arcsin c$, hence $|\sin\theta_j|\lt c$. Its center must lie within $c+d=0.16$ of a tooth. As a necessary condition, the center residue modulo $\ell_t$ lies in a circular interval of length $2(c+d)=0.32$.

The source ratio is exactly $\ell_r/\ell_t=43/39$. Its center residues occupy 39 equally spaced sites separated by $\ell_t/39=0.028$. An interval of length 0.32 contains at most 12 such sites, and each site occurs at most three times among 100 rotors. Therefore

```math
n=|C|\le36.
```

This simple bound needs no simulation. Counting the integer residues $4j\bmod39$, for $0\le j\lt100$, tightens the necessary-condition bound to **33**. The reproducible script records all site multiplicities and every cyclic 12-site count. Finite, nonperiodic tooth boundaries can only remove contacts from this residue-based eligibility test. These bounds rely on the specified geometry and motor count.

## The sigmoid force balance has a unique velocity here

For `switch=0`, the algebraic source force is

```math
F_s(\omega)=1-\frac{\omega}{4\pi}[1+\tanh(40\omega)].
```

Let $b_j=|\cos\theta_j|$ on a fixed contact set $C$. Then $b_j\ge\sqrt{1-c^2}>0$, and contacting angular velocities are $v/b_j$, where $v=\dot x$. Define the contribution of noncontacting rotors as $B$; its tiny analytic-versus-FP64 distinction is given below. The scalar residual is

```math
G(v)=\alpha v+\beta x-
\sum_{j\in C}b_jF_s(v/b_j)-B.
```

It is not correct to assume $F_s$ is everywhere decreasing. Writing $u=40\omega$ gives

```math
F_s'(\omega)=-\frac{1+\tanh u+u\,\mathrm{sech}^2u}{4\pi}.
```

Its largest positive derivative occurs at $u=-w$, where $w>0$ solves $w\tanh w=1$. At this point the expression in the numerator equals $1-w$. Hence

```math
M=\max_\omega F_s'(\omega)
=\frac{w-1}{4\pi}
\approx0.01588992131344,
\qquad w\approx1.199678640258.
```

The factors $b_j$ cancel when differentiating $G$ with respect to $v$, so

```math
G'(v)=\alpha-\sum_{j\in C}F_s'(v/b_j)
\ge\alpha-nM
\ge1-36M
\approx0.4279628327>0
```

for the source geometry and $\alpha\ge1$. Moreover, $G(v)\to-\infty$ as $v\to-\infty$, and $G(v)\to+\infty$ as $v\to+\infty$. Thus each fixed admissible contact branch has exactly one algebraic velocity. This does not prove monotonicity at arbitrary low drag, for arbitrary motor counts, or after changing the geometry.

For a convenient bracket, put $A=\sum_{j\in C}b_j+B-\beta x$. Since $F_s(0)=1$, the velocity has the sign of $A$. For positive arguments $F_s\le1$; for negative arguments $F_s\ge1$. Therefore the unique root lies between $0$ and $A/\alpha$, including the endpoints. A bracketed scalar solve avoids guessing an algebraic branch.

## Conditional backbone bounds and the free-force caveat

The actual Julia FP64 expression evaluates $F_s(2\pi)$ as zero because its hyperbolic tangent rounds to saturation. Analytically, however,

```math
\eta=F_s(2\pi)
=\frac{1}{1+\exp(160\pi)}>0,
\qquad
B=\eta\sum_{j\notin C}H(\theta_j)|\cos\theta_j|.
```

The script records both values. This exponentially small remainder has no numerical relevance at source precision, but it matters to the literal wording of an exact reduction. On a branch, $B$ does not depend on $v$ and therefore does not affect the uniqueness proof.

Because $G$ is increasing and $G(0)=\beta x-\sum_C b_j-B$, the sign of $v$ is the sign of $\sum_C b_j+B-\beta x$. For $\beta>0$, this gives inward velocities at the endpoints of

```math
0\le x\le\frac{N}{\beta}.
```

The analytic sharper bound is $[0,(36+64\eta)/\beta]$, while the FP64 source reduction gives $[0,36/\beta]$. They follow from $\sum_C b_j+B\le n+(N-n)\eta$ and $n\le36$. The sharper integer count can improve 36 to 33, with the corresponding free term changed to $67\eta$.

These are **conditional forward-invariance statements**: they assume a well-defined event continuation with no unmodeled jump in $x$ and with branch velocities obeying the force balance. They do not supply such a continuation at problematic boundaries. The source state has 101 coordinates; compactness of phases modulo $2\pi$ does not reduce it to a three-dimensional system.

## Scalar evolution within a fixed contact branch

While a rotor remains in the contact interval, $\cos\theta_j<0$, and

```math
\frac{d}{dt}\sin\theta_j
=\cos\theta_j\frac{v}{|\cos\theta_j|}=-v.
```

If it entered the current branch at $(x_a,\theta_{j,a})$, then

```math
\sin\theta_j(x)=\sin\theta_{j,a}-(x-x_a),
\qquad
\theta_j(x)=\pi-\arcsin[\sin\theta_{j,a}-(x-x_a)]
\pmod{2\pi}.
```

The branch is valid only while the reconstructed angle remains in the contact interval and its physical tooth has not disappeared. With an unwrapped tooth position locally equal to $t_{k,a}+x-x_a$, the separation obeys

```math
\frac{d}{dt}(r_j-t_k)=v-v=0.
```

Free phases are explicit: $\theta_j(t)=\theta_{j,a}+2\pi(t-t_a)$. For the FP64 source force, they exert no force until contact. Consequently a fixed contact branch gives the scalar autonomous equation $\dot x=V(x)$, with algebraic velocity solved as above. If the positive analytic remainder $\eta$ is retained, their explicitly known phases contribute $B(t)$, giving a scalar **nonautonomous** equation instead. Neither statement removes phase offsets or contact history from the full hybrid problem.

The computational opportunity is to solve one scalar backbone equation and reconstruct contacting phases, while locating the earliest distance, angular, and tooth-wrap events. It is not permission to freeze a contact set over a missed event or to replace the network by one global scalar oscillator.

## A strict-distance entry needs an interpretation

Consider the source-admissible state with $\alpha=\beta=1$, $\gamma=0$, $x=d=0.02$, rotor 0 at $\theta_0=\pi$, and every other rotor at angle zero. Tooth 0 is at $d$, the rotor tip is at zero, and its signed gap is $g=r_0-t_0=-d$. It is angularly eligible, but the source's strict distance test excludes contact.

On that free branch the FP64 backbone speed is $v=-0.02$, and

```math
\dot g=2\pi-v=6.303185307\ldots>0.
```

Thus the free field points into the contact band. The analytic free-force remainder changes $v$ by only $\eta$ and does not change this sign. On the limiting included-contact branch the algebraic speed is approximately $0.84544349$, and $\dot g=0$ exactly.

An event solver cannot simply stop at $g=-d$, apply the contact branch forever, and also claim to obey the source's strict predicate: the contact branch keeps $g=-d$, which remains excluded. Conversely, a finite numerical step can penetrate slightly into the open band and then keep that small penetration. Near this example, with no other event intervening, the normal-coordinate issue has the form “positive derivative outside/on the boundary, zero derivative inside.” It is not the usual transverse crossing into a smooth interior vector field.

This establishes a **local continuation ambiguity/regularity obstruction that needs explicit treatment**, not a claim that every source trajectory fails or a global nonexistence theorem. A closed-boundary capture rule, an identified limiting interpretation, or a specified contact regularization would resolve a modeling choice, but each must be declared rather than silently substituted for the author code. Small-step comparisons must test the resulting contact sequence and penetration dependence. No such substitute is implemented in this audit.

## An angular boundary can have opposing directions

There is a separate issue at the lower angular limit. Set $x=\ell_r-c=1.064$, rotor 1 at $\theta_1=\pi-\arcsin c$, all other angles zero, and again $\alpha=\beta=1$, $\gamma=0$. Its tip coincides with tooth 0, so the tooth-distance condition is strictly satisfied. The angular condition is not: the boundary itself has $H=0$.

The free branch has $\dot\theta_1=2\pi>0$, pointing into the angular contact interval. On the included-contact side,

```math
S=\sqrt{1-c^2}\approx0.9901515<\beta x,
\qquad
v\approx-0.07381839,
\qquad
\dot\theta_1\approx-0.07455262<0,
```

pointing back out. Thus a naive sequence of enter/exit events can chatter at the same angular boundary. Whether to use a sliding solution, a capture/release rule, a regularization, or a particular numerical limit is not fixed by the within-branch identities. Simultaneous contacts, grazing events, and tooth disappearance likewise require explicit handling. This note does not choose a physical contact law on the authors' behalf.

## Verification and a gate before chaos claims

Run `python research/audit_dcr_structure.py`. It performs only algebraic checks, in approximately 0.06 seconds after Python startup:

- 256 independent tooth-set and nearest-distance comparisons against explicit wrapping of all 111 teeth: largest error $1.8\times10^{-13}$.
- 128 source-admissible instantaneous states with deliberately recruited eligible motors: at most 33 contacts, largest force residual $1.1\times10^{-14}$, and no velocity-sign failures.
- Independent scalar-force derivative and finite-difference checks, exact residue counts, contact kinematic identities, and reconstruction checks.
- The two explicit switching-boundary states above, with both one-sided algebraic velocities retained.

The [JSON record](../results/dcr/structure-audit.json) includes source identity, assumptions, settings, and every result. No trajectories or exponents are computed.

The next scientific gate is event consistency: identify and document the source numerical limit or an explicitly separate hybrid interpretation, resolve capture/release and angular boundary behavior, and test event-sequence convergence before interpreting irregularity or positive exponents. A within-branch tangent calculation alone cannot capture event-time sensitivity. Uniqueness of the algebraic speed does not establish uniqueness of the discontinuous flow.

## References

1. B. Warmington, J. Rossiter, and H. Bloomfield-Gadêlha, *A discontinuously coupled network of phase oscillators replicates cooperation between motors in actomyosin systems*, **Journal of the Royal Society Interface 23**, 20250438 (2026). [DOI](https://doi.org/10.1098/rsif.2025.0438), [public final article](https://aura.abdn.ac.uk/bitstreams/2e76dc81-a0f5-4fb5-aa9c-7ca1016f8730/download).
2. B. Warmington, *Discontinuously coupled rotor model*, version 1.0.0. [Zenodo source release](https://doi.org/10.5281/zenodo.18196853).
3. Author repository, [pinned Julia implementation](https://github.com/polymaths-lab/actomyosin-kuramoto-rotor-model/blob/796c5e39468b06fd6fc2c378437d618593948e98/DCR_Model_available.jl). The locally cached pinned file was inspected; the independent audit does not execute or redistribute it.
