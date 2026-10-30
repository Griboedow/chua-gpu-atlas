# Chua along the quadratic parameter path: a regular singular limit with important caveats

The question is whether a chaotic family can persist as $\alpha\to\infty$ with $\beta=\kappa\alpha^2$ and both diode slopes fixed. An exact rescaling produces a useful limiting model: a harmonic oscillator drives a scalar periodic equation. That limiting equation has no recurrent chaotic attractor at fixed amplitude. **This does not prove that chaos terminates at a finite $\alpha$.** Long-time limits, weak amplitude drift, and losses of uniform attraction must be treated separately.

The derivations below are our explicit checks of the canonical equations, not a claim of a new Chua mechanism or a new singular-perturbation theory. Chua return maps and slow–fast analysis have substantial primary literature [1–4]. This note belongs to the original Chua project, rather than the separate search for other physical models.

## Exact scaling and an action identity

Start with the same normalization used by this repository:

```math
\dot x=\alpha[y-x-f(x)],\qquad
\dot y=x-y+z,\qquad
\dot z=-\beta y,
```

```math
f(x)=m_1x+\frac{m_0-m_1}{2}(|x+1|-|x-1|),\qquad g(x)=x+f(x).
```

Thus $m_0,m_1$ are the slopes of $f$; the slopes of $g$ are $1+m_0,1+m_1$. Some primary papers use slope symbols for $g$ instead, so numerical slopes cannot be transferred without checking conventions [1]. Assume $\kappa>0$, and put

```math
\varepsilon=\frac1\alpha,\qquad s=\alpha t,\qquad Z=\frac z\alpha,
\qquad \beta=\kappa\alpha^2.
```

Subscripts denote differentiation in the fast time $s$. The transformation is exact:

```math
x_s=y-g(x),\qquad
y_s=Z+\varepsilon(x-y),\qquad
Z_s=-\kappa y.
```

The oscillator action satisfies the exact balance

```math
I=\frac{\kappa y^2+Z^2}{2},\qquad
I_s=\varepsilon\kappa y(x-y).
```

For nonzero amplitude, write $y=R\cos\theta$, $Z=-\sqrt\kappa R\sin\theta$, so $I=\kappa R^2/2$. Then

```math
R_s=\varepsilon\cos\theta\,[x-R\cos\theta],\qquad
\theta_s=\sqrt\kappa-\varepsilon\frac{\sin\theta}{R}[x-R\cos\theta].
```

The polar description is singular at $R=0$; it is not a uniform coordinate system for starts whose oscillator amplitude vanishes with $\varepsilon$.

## What is rigorously simple at zero epsilon

At $\varepsilon=0$, the pair $(y,Z)$ is an independent harmonic oscillator. Fix its amplitude $R>0$ and phase. The remaining equation is

```math
x_s=R\cos(\sqrt\kappa s)-g(x),\qquad T=\frac{2\pi}{\sqrt\kappa}.
```

Because $g$ is globally Lipschitz, solutions are unique and two scalar solutions cannot cross. Consequently the time-$T$ map $P_R$ is continuous and strictly increasing. If $P_R(x_0)>x_0$, its iterates increase; if $P_R(x_0)<x_0$, they decrease. Every **bounded** sequence of these iterates therefore converges, and continuity makes its limit a fixed point of $P_R$. Continuous dependence on initial conditions over one period then gives convergence to a $T$-periodic response.

This excludes nontrivial subharmonic cycles and recurrent chaotic behavior on a regular fixed-amplitude fiber. It is an order argument, not a finite numerical search. For the usual positive outer slope $a=1+m_1>0$, write

```math
g(x)=ax+h(x),\qquad h(x)=(m_0-m_1)\mathrm{sat}(x),\qquad
|h(x)|\leq |m_0-m_1|,
```

where $\mathrm{sat}(x)=\max(-1,\min(1,x))$. Scalar comparison also bounds every forced solution by the linear envelope with forcing magnitude $R+|m_0-m_1|$. At $R=0$, the limiting problem is instead the scalar autonomous equation $x_s=-g(x)$.

An **unstable periodic response can nevertheless have a positive exponent**. On a regular periodic fiber its scalar multiplier is

```math
P_R'(x_*)=\exp\!\left[-\int_0^T g'(x_*(s))\,ds\right]>0.
```

This gives an exact residence-time diagnostic. If $d$ is the fraction of the period spent inside $|x|<1$, then

```math
\frac{\log P_R'(x_*)}{T}
=-[(1+m_1)+(m_0-m_1)d].
```

When $m_0<-1$ and $m_1>-1$, attraction requires $d<(1+m_1)/(m_1-m_0)$. Uniform attraction needs a strict uniform gap. At $m_1=-5/7$, the thresholds for $m_0=-2,-5,-10,-30$ are respectively $2/9$, $1/15$, $2/65$, and $2/205$. Thus very deep inner slopes can lose contraction after only a small fraction of the cycle in the inner region. A multiplier reaching one indicates lost linear attraction, not by itself a proved fold. This diagnostic belongs to the frozen scalar periodic response; it is not the full finite-$\varepsilon$ Lyapunov exponent.

For example, with $m_0=-8/7$, $m_1=-5/7$, $\kappa=10/49$, and $R=0.1$, the exact inner-branch periodic response stays inside $|x|<0.212$, has multiplier approximately 7.2930, and exponent $1/7$ per fast time. It is an unstable periodic orbit, not a chaotic attractor. The matching positive outer-branch response has multiplier approximately 0.0188012. Both are checked against their exact affine solutions in the accompanying script.

## The averaged drift is conditional on a selected periodic branch

Suppose a family of attracting scalar periodic responses $x_j(s;I)$ has been selected and remains uniformly attracting over a compact amplitude interval. The formal first averaged action equation, in the original model time $t=\varepsilon s$, is

```math
\frac{dI}{dt}\simeq\Gamma_j(I),\qquad
\Gamma_j(I)=\frac{\kappa}{T}\int_0^T y(s)[x_j(s;I)-y(s)]\,ds.
```

This is a scalar envelope equation on a chosen branch. If several scalar periodic solutions coexist, averaging cannot silently replace them with a single response. Folds, loss of uniform transverse attraction, switching between branches, and phase-coordinate degeneracy require separate analysis.

A uniformly attracting periodic family suggests an attracting two-dimensional amplitude–phase manifold for sufficiently small $\varepsilon$. **This remains a conditional route to a local exclusion, not a theorem established here.** One must verify an appropriate invariant neighborhood, uniform normal hyperbolicity, and enough regularity of the actual return map. The vector field is continuous but only piecewise differentiable. Smooth persistence theorems cannot simply be applied across every crossing or grazing of $x=\pm1$. Specialized nonsmooth averaging exists [4], but its precise hypotheses have not been verified for this proposed manifold.

## Outer-branch covariance and the exact stability line

Assume $m_0<-1$ and $-1<m_1<0$, and write $a=1+m_1\in(0,1)$. In the positive outer branch, $g(x)=a(x-x_e)$ with $x_e=(m_1-m_0)/a>1$. A periodic response that stays entirely in this branch is

```math
x(s)=x_e+\frac{R}{a^2+\kappa}
\left[a\cos(\sqrt\kappa s)+\sqrt\kappa\sin(\sqrt\kappa s)\right].
```

It follows exactly that

```math
\frac{\langle xy\rangle}{\langle y^2\rangle}
=\frac{a}{a^2+\kappa},\qquad
\langle I_s\rangle=\varepsilon I
\left[\frac{a}{a^2+\kappa}-1\right].
```

The constant offset averages out. The one-branch assumption holds, for example, if $R/\sqrt{a^2+\kappa}<x_e-1$. This exact formula must not be assigned to a switching double-scroll trajectory.

It also supplies a large-amplitude asymptote. For any bounded periodic response with $a>0$, subtract the periodic solution of $x_s+ax=y$. The difference solves a stable scalar linear equation with forcing $-h(x)$, hence has magnitude at most $|m_0-m_1|/a$. Therefore the covariance ratio differs from $a/(a^2+\kappa)$ by $O(1/R)$ as $R\to\infty$. This is an amplitude asymptote, separate from the $\alpha\to\infty$ limit.

The exact fast-time characteristic polynomial at either outer equilibrium is

```math
p(\lambda)=\lambda^3+(a+\varepsilon)\lambda^2
+[\kappa+\varepsilon(a-1)]\lambda+a\kappa.
```

The cubic stability criterion gives the exact linear boundary

```math
\kappa_c=(1-a)(a+\varepsilon),\qquad
\beta_c=a(1-a)\alpha^2+(1-a)\alpha.
```

For every finite $\alpha$ (so $\varepsilon>0$), the outer equilibria are linearly asymptotically stable above this boundary and unstable below it, under the stated slope assumptions. At $\varepsilon=0$, the oscillator pair is neutral for every $\kappa>0$. This is the same outer-equilibrium stability surface already used by the atlas, expressed in scaled variables. It is not a newly discovered chaos boundary.

There is also a qualification to the term *Hopf line*. At equality,

```math
p(\lambda)=(\lambda+a+\varepsilon)
[\lambda^2+a(1-a)].
```

The local field is exactly affine. Sufficiently small ellipses in its center plane stay entirely in the outer branch and are exact periodic orbits. They form a nonisolated family; generic nondegenerate Hopf conclusions do not follow. A full-field integration at $\alpha=64$ reproduces one such ellipse, with $1.45\leq x\leq1.55$, and closes to approximately $6.7\times10^{-16}$.

## Why the leading parabola hides very slow growth

For the standard outer slope $m_1=-5/7$, the limiting value and finite-$\alpha$ correction are

```math
\kappa_0=\frac{10}{49},\qquad
\kappa_c=\frac{10}{49}+\frac{5}{7\alpha}.
```

At fixed $\kappa=\kappa_0$, every finite $\alpha$ lies below the exact outer stability boundary. Expanding the oscillatory eigenvalues gives

```math
\mathrm{Re}\,\lambda_s
=\frac{1-a}{2a}\varepsilon^2+O(\varepsilon^3),\qquad
\mathrm{Re}\,\lambda_t
=\frac{1-a}{2a\alpha}+O(\alpha^{-2})
=\frac{1.25}{\alpha}+O(\alpha^{-2}).
```

Independent eigenvalue controls illustrate the asymptote:

| $\alpha$ | Outer-pair growth in original time | $\alpha$ times that growth |
|---:|---:|---:|
| 128 | 0.00979872 | 1.254237 |
| 2048 | 0.000610499 | 1.250303 |
| 32768 | 0.0000381476 | 1.250019 |

These are equilibrium eigenvalues, not trajectory Lyapunov exponents or evidence of an attractor. Small positive growth near the parabola is expected from the linear equations alone.

## The limits and initial scales do not commute

Several distinctions matter when interpreting a growing parameter atlas:

- Holding fast time $s$ fixed as $\alpha\to\infty$ sends the original time $t$ to zero. On fixed original-time intervals, the oscillator performs increasingly many cycles and its amplitude can change by order one. The frozen-action limit is not a uniform long-time approximation.
- The limiting phase oscillates ever faster in original time. An envelope, phase-aligned return map, or invariant measure is generally more informative than a pointwise time-trace limit.
- Lyapunov rates satisfy $\lambda_s=\lambda_t/\alpha$. A rate tending to zero in fast time need not tend to zero in original time. The constant linear state scaling affects finite-window conditioning but not asymptotic exponents at each fixed $\alpha$.
- Bounded original $z$ becomes $Z=z/\alpha\to0$. A usual start with $y=z=0$ collapses to zero limiting oscillator action. Near an outer equilibrium, the finite-$\varepsilon$ equilibrium has $Z_e=-\varepsilon x_e$, so an initial excitation of order $\varepsilon$ can take order $\log\alpha$ in original time to grow to order one when the growth rate stays positive. On $\kappa=10/49$, the weak rate above can instead lead to order $\alpha\log\alpha$ delays. These estimates describe appropriate local unstable departures, not every initial condition.
- Compactness must be checked in the rescaled variables. The singular-limit argument does not cover families whose rescaled amplitudes diverge or whose slopes also vary with $\alpha$.

Thus the present algebra leaves both finite termination and persistence of chaotic regimes at arbitrarily large finite $\alpha$ unresolved. A regular limiting fiber can coexist with singularly shrinking parameter windows or very long finite-parameter transients. Uniform attraction on a verified regular fiber branch could exclude chaos locally; approaching a nonuniform or singular region would invalidate that shortcut.

A useful next calculation is to locate the fixed points and multipliers of the scalar fiber maps, evaluate their averaged action drift, and compare finite-$\alpha$ return maps on matched nonzero amplitudes. The small frozen-fiber check below begins this calculation. The original fixed-$\kappa$ question should remain distinct from a boundary-layer study with $\kappa-10/49$ of order $1/\alpha$.

## Frozen fibers sampled from one finite-parameter scout

The [independent scalar probe](../research/probe_chua_fiber.py) takes all five saved observation-boundary amplitudes from scout row 2198, with $(\alpha,\kappa,m_0,m_1)=(2048,0.1425,-2,-5/7)$. It then sets $\varepsilon=0$ and solves the scalar equation with $y=R\sin(\sqrt\kappa s)$. The section is $y=0$ with increasing $y$, and the fast forcing period is approximately 16.64456. These are frozen problems; they are not another simulation of the full finite-$\alpha$ trajectory.

Here $x_e=4.5$. The exact outer affine response first touches a breakpoint at $R=(x_e-1)\sqrt{a^2+\kappa}=1.65699276$. At each saved amplitude, two independent starts $x(0)=\pm4.5$ were iterated to a periodic response. Additional controls solve the unique half-wave condition $x(T/2)=-x(0)$ and explicitly integrate the symmetry partner of the positive-well response. Root shooting also finds unstable periodic responses, which must not be called attractors.

The table uses the positive-well response; its negative half-wave partner has the same multiplier and averaged drift. The central response has half-wave symmetry and crosses both breakpoints. The scalar multiplier is $q=\exp[-\int_0^T g'(x(s))\,ds]$. The drift column is $\kappa H$, where $H=\langle y(x-y)\rangle$ is the period average; this is the candidate leading averaged action rate in original time on a selected branch.

| Frozen amplitude $R$ | Outer $q$ | Outer $\kappa H$ | Central $q$ | Central $\kappa H$ |
|---:|---:|---:|---:|---:|
| 1.430266 | 0.00860347 | +0.0400465 | 3.50687 | -0.299199 |
| 1.476419 | 0.00860347 | +0.0426727 | 1.93963 | -0.311715 |
| 1.527137 | 0.00860347 | +0.0456548 | 1.16643 | -0.323143 |
| 1.675276 | 0.0297122 | +0.0544876 | 0.441353 | -0.346974 |
| 1.714508 | 0.154060 | +0.0490613 | 0.368975 | -0.351642 |

At the two larger sampled amplitudes, three distinct attracting scalar periodic responses coexist: the two symmetry-related well responses and the central response. In particular, the negative initial value $-4.5$ reaches the central response, while $+4.5$ reaches the positive-well response. Thus initial sign alone does not identify a frozen basin. At the three smaller amplitudes, the computed central response is unstable while the well responses remain attracting. Its multiplier crosses one somewhere between the sampled amplitudes 1.52714 and 1.67528; this probe neither locates nor classifies that bifurcation and is not a complete fixed-point census.

The two attracting branch types have **opposite averaged action drift** where they coexist. This is a concrete reason not to assume one globally selected, uniformly attracting scalar response. It suggests checking branch selection and departure in a matched finite-$\alpha$ return map. It does not yet establish delayed switching, a relaxation cycle, or a chaotic mechanism in that full system. None of the frozen scalar responses is chaotic.

All 20 periodic integrations/shooting controls have independent Radau closure residuals below $1.5\times10^{-9}$, and DOP853 residuals below $3.6\times10^{-10}$. Breakpoint events determine residence times and multipliers. Event-split quadrature of the dense periodic trajectory agrees with a 4001-point Simpson check within $5.7\times10^{-10}$ in $H$. Reproduction takes approximately 21 seconds; [JSON settings and diagnostics](../results/chua-scaling/fiber-probe.json) and [sampled cycles and iterates](../results/chua-scaling/fiber-probe.npz) retain source identity, both starts, all residuals, and qualifications. These numerical checks are not a proof of completeness or scientific novelty.

## Numerical checks and literature scope

`python research/audit_chua_scaling.py` computes 256 seeded identity checks, six small eigenvalue problems, two one-period affine controls, one center-ellipse control, and seven scalar return-map values. It runs no finite-$\alpha$ chaos survey. Errors are below $2.9\times10^{-16}$ for the scaled field, $7.2\times10^{-15}$ for the action identity, and $4.0\times10^{-15}$ for the characteristic polynomial. Settings, values, access limits, and qualifications are in [structure-audit.json](../results/chua-scaling/structure-audit.json).

The primary literature establishes that these methods and mechanisms are longstanding. Chua, Komuro, and Matsumoto develop the piecewise-linear family and return maps [1]. Rossetto explicitly applies singular perturbation to Chua [2]; only the publisher abstract was recovered in this pass, so coverage of this exact quadratic parameter path was not verified. Euzébio and Llibre use averaging for a different six-parameter **cubic** Chua model and zero–Hopf scaling [3]. Their stated smoothness hypotheses cannot automatically be transferred to our piecewise-linear limit. Nonsmooth averaging has its own theory [4]. None of these bounded source checks establishes novelty or an absence of earlier results on our exact path.

## References

1. L. O. Chua, M. Komuro, and T. Matsumoto, *The Double Scroll Family*, **IEEE Transactions on Circuits and Systems 33**, 1073–1118 (1986). [Primary author-hosted paper](https://people.eecs.berkeley.edu/~chua/papers/Chua86.pdf), [DOI:10.1109/TCS.1986.1085869](https://doi.org/10.1109/TCS.1986.1085869). Introductory equations and the return-map/bifurcation scope were inspected; this note does not claim to reprove their chaos theorem.
2. B. Rossetto, *Chua's circuit as a slow-fast autonomous dynamical system*, **Journal of Circuits, Systems and Computers 3**, 483–496 (1993); reprinted in *Chua's Circuit: A Paradigm for Chaos*, pp. 711–724. [Publisher volume with chapter abstract](https://doi.org/10.1142/1997), [chapter DOI](https://doi.org/10.1142/9789812798855_0037). Full chapter retrieval failed; no specific theorem from it is asserted here.
3. R. D. Euzébio and J. Llibre, *Zero-Hopf bifurcation in a Chua system*, **Nonlinear Analysis: Real World Applications 37**, 31–40 (2017). [Primary author manuscript](https://files.cercomp.ufg.br/weby/up/1076/o/Euzebio_Llibre_chua.pdf), [DOI:10.1016/j.nonrwa.2017.02.002](https://doi.org/10.1016/j.nonrwa.2017.02.002). Model equation (1), rescaling, and averaging hypotheses in section 2 were read.
4. J. Llibre, A. C. Mereu, and D. D. Novaes, *Averaging theory for discontinuous piecewise differential systems*, **Journal of Differential Equations 258**, 4007–4032 (2015). [Author-institution abstract and record](https://portalrecerca.uab.cat/en/publications/averaging-theory-for-discontinuous-piecewise-differential-systems/), [DOI:10.1016/j.jde.2015.01.022](https://doi.org/10.1016/j.jde.2015.01.022). Abstract/metadata inspected; full hypotheses were not recovered or applied in this pass.
