# An order constraint on the smooth two-operator model

**Verified conclusion:** when $k,b>0$ and $b^2\ge4k$, the explicitly logistic-smoothed two-operator model is a strongly 2-cooperative cyclic negative-feedback system after an invertible linear change of coordinates. Every compact omega-limit set that **contains no equilibrium** is a periodic orbit. This is an application of an established theorem, not a new chaos result.

The equilibrium-free qualification is retained. The algebra needed for a stronger classical cyclic-feedback conclusion is favorable, but the original 1990 full theorem was not recovered in this audit. We do not infer a global no-chaos theorem from equilibrium uniqueness alone.

## Model and scope

The continuous-heat interpretation of Guo et al.'s two material operators [1] is recorded in [the research report](coupled-operators-research.md). Here the optical law is the explicitly assumed profile

```math
S_w(u)=\frac{1}{1+e^{-u/w}},\qquad w>0.
```

It is not a measured source beam profile. The six equations are

```math
\begin{aligned}
\dot x_-&=v_-,&\dot v_-&=k(y_--x_-)-bv_-,&
\dot y_-&=S_w(x_+-\delta_+)-y_-,\\
\dot x_+&=v_+,&\dot v_+&=k(y_+-x_+)-bv_+,&
\dot y_+&=S_w(\delta_--x_-)-y_+.
\end{aligned}
```

This audit applies to positive widths. It does **not** apply the smooth theorem to the ideal Heaviside relay, a numerical sigmoid clipped to constant plateaus, or an event-driven limit without verifying the relevant nonsmooth theory. Nor does it apply to an operator with underdamped mechanics through the coordinate change below.

## Positive mechanical stages when the poles are real

For $b^2\ge4k$, define positive real rates

```math
r_1+r_2=b,\qquad r_1r_2=k,\qquad
r_2=\frac{b+\sqrt{b^2-4k}}2,\qquad r_1=\frac{k}{r_2}.
```

Set $z_\pm=v_\pm+r_1x_\pm$. Then each mechanical subsystem becomes

```math
\dot z_\pm=ky_\pm-r_2z_\pm,
\qquad
\dot x_\pm=z_\pm-r_1x_\pm.
```

This is an exact change of coordinates, not removal of inertia or a slow-manifold approximation. It remains invertible at critical damping $r_1=r_2=\sqrt{k}$: the inverse is simply $v_\pm=z_\pm-r_1x_\pm$.

In cyclic order

```math
u=(y_+,z_+,x_+,y_-,z_-,x_-),
```

the equations read

```math
\begin{aligned}
\dot u_1&=S_w(\delta_--u_6)-u_1,\\
\dot u_2&=ku_1-r_2u_2,\\
\dot u_3&=u_2-r_1u_3,\\
\dot u_4&=S_w(u_3-\delta_+)-u_4,\\
\dot u_5&=ku_4-r_2u_5,\\
\dot u_6&=u_5-r_1u_6.
\end{aligned}
```

Each equation depends only on its own coordinate and its predecessor. The five forward links have positive derivatives; the closing link $u_6\to u_1$ has a negative derivative. All are strictly nonzero for finite states and $w>0$. The sign product is negative. The whole loop is not an ordinary cooperative system, despite the positive mechanical stages.

## Boundedness and the theorem's hypotheses

An explicit invariant box is

```math
0\le y_\pm,x_\pm\le1,
\qquad 0\le z_\pm\le r_1.
```

At a thermal boundary the sigmoid points inward; at $z=r_1$, the derivative is $k(y-1)\le0$; at $x=1$, it is $z-r_1\le0$. Lower boundaries are also inward. The previously used physical starts with thermal fractions and displacements in $[0,1]$ and zero velocities belong to this box because $z=r_1x$.

More generally, any finite initial state has a bounded forward orbit: each thermal equation is a stable first-order filter of a bounded sigmoid, followed by two stable first-order mechanical filters. The smooth field is globally Lipschitz for each fixed positive width. Thus omega-limit sets are nonempty and compact. These bounds concern the normalized model, not the validity of arbitrary actual cantilever deflections.

The transformed Jacobian has positive subdiagonal entries, one negative upper-right corner, and no other off-diagonal entries. Its second additive compound is Metzler (all off-diagonal entries nonnegative) and irreducible. One can view its graph as the allowed forward movements of two distinct indices around a six-cycle; all 15 index pairs communicate. The sign at the corner is precisely what makes its wedge-space off-diagonal contributions positive. This structure persists in the averaged Jacobian between two solutions on the convex state space.

Weiss and Margaliot [3], Corollary 4 and Theorems 5–6, give the corresponding strong 2-cooperativity criterion. Their Theorem 12 states that a compact omega-limit set without an equilibrium is a periodic orbit. The primary author manuscript and these specific statements were read. Consequently an equilibrium-free quasiperiodic or chaotic omega-limit set is excluded in this parameter region. This does not exclude finite-time expansion, transient irregularity, or a positive exponent on an unstable invariant solution.

## Equilibrium uniqueness and the remaining stronger-theorem check

At equilibrium, $y_\pm=x_\pm$ and $v_\pm=0$. Eliminating $x_+$ gives

```math
x_-=
S_w\big(S_w(\delta_--x_-)-\delta_+\big).
```

The right side is strictly decreasing, while the left side increases. The difference is positive at $x_-=0$ and negative at $x_-=1$, so there is exactly one equilibrium, with both displacements strictly between zero and one.

Let $a_-,a_+>0$ be the two illumination-slope magnitudes there. The exact characteristic polynomial and determinant are

```math
\det(\lambda I-J)
=[(\lambda+1)(\lambda^2+b\lambda+k)]^2+k^2a_-a_+,
\qquad
\det(-J)=k^2(1+a_-a_+)>0.
```

This determinant is nonzero even at a Hopf candidate; it is not a Hurwitz-stability test.

Mallet-Paret and Smith's original cyclic-feedback paper [2] was verified through its publisher record, but the PDF request returned subscription content. Hori, Hara, and Kim's accessible primary research report [4], Proposition 1, explicitly attributes a stronger alternative to that paper: for its monotone gene-chain setting with a unique equilibrium, bounded positive orbits have equilibrium, periodic, or homoclinic omega-limits, and the additional positive-determinant condition eliminates the homoclinic alternative. Its model, proposition, and determinant calculation were inspected.

Our six-cycle has the same relevant sign structure and a positive determinant, making that extension a well-supported mathematical lead. However, [4] states the corollary for its own gene-chain equations, and we have not directly checked the full original hypotheses needed to transfer every equilibrium-containing conclusion. Therefore the published conclusion of this audit stays at the fully verified equilibrium-free theorem above. In particular, uniqueness of our equilibrium is not presented as a substitute for the missing theorem check. No claim of uniqueness or attraction of a periodic orbit follows from the present audit.

## What this says about the existing numerical candidates

For the two existing references:

| Reference | $b^2/(4k)$ | Actual $b$ | Critical-damping $b=2\sqrt{k}$ |
|---|---:|---:|---:|
| A, scout 1026 | 0.0107317 | 1.78418 | 17.2229 |
| B, scout 3054 | 0.000646340 | 4.71772 | 185.568 |

Both are strongly underdamped and lie outside this order argument. Their numerical evidence is neither contradicted nor upgraded by the theorem.

The mechanical distinction has a direct interpretation. A cascade of the two real positive stages has a nonnegative impulse response. An underdamped second-order mechanical filter instead has an exponentially decaying sine impulse response and can overshoot or reverse its response to past heating. Thus real-pole mechanics allow the cyclic order argument; mechanical ringing removes this particular route to excluding complex recurrent dynamics.

Within the stated smooth family, underdamping is therefore necessary to evade this **equilibrium-free** exclusion. It has not been proved here to be necessary for every possible notion of chaos, and it is certainly not sufficient to produce chaos. Changing damping also changes equilibrium stability, transient times, and periodic orbits; the threshold $b^2=4k$ is not asserted to be an attractor bifurcation boundary.

## Reproduction

Run `python research/audit_coupled_order.py`. No GPU or trajectory integration is used. The [JSON results](../results/coupled-operators/order-audit.json) contain 128 seeded algebra checks, including 32 critical-damping cases, direct equilibrium controls, and exact reference parameters. They verify field conjugacy, Jacobian similarity, the positive invariant box, and the Metzler/irreducible second-compound pattern. Four matrix-exponential controls compare the additive-compound exponential with all corresponding minors of the original exponential.

Maximum relative field and Jacobian discrepancies are $2.5\times10^{-15}$ and $2.2\times10^{-14}$; the compound identity agrees within $1.2\times10^{-16}$. These are implementation checks of the derivations. They are not numerical evidence for the theorem itself and do not constitute a new trajectory survey.

## References

1. H. Guo et al., *Light-mediated communication in responsive materials ranging from individual self-oscillators to feedback-driven network*, **Nature Communications 16**, 11315 (2025). [DOI:10.1038/s41467-025-66395-3](https://doi.org/10.1038/s41467-025-66395-3). Model provenance and the continuous-heat interpretation are audited in [light-operators-audit.md](light-operators-audit.md); the logistic profile is our explicit modification.
2. J. Mallet-Paret and H. L. Smith, *The Poincare–Bendixson theorem for monotone cyclic feedback systems*, **Journal of Dynamics and Differential Equations 2**, 367–421 (1990). [Publisher record](https://link.springer.com/article/10.1007/BF01054041). Abstract/metadata verified; full original theorem not recovered.
3. E. Weiss and M. Margaliot, *A Generalization of Linear Positive Systems with Applications to Nonlinear Systems: Invariant Sets and the Poincare–Bendixson Property*, **Automatica 123**, 109358 (2021). [DOI:10.1016/j.automatica.2020.109358](https://doi.org/10.1016/j.automatica.2020.109358), [full primary manuscript, arXiv:1902.01630v4](https://arxiv.org/abs/1902.01630v4). Corollary 4, Theorems 5–6, and Theorem 12 inspected; the last explicitly retains the absence-of-equilibria condition.
4. Y. Hori, S. Hara, and T.-H. Kim, *Existence Criteria of Periodic Oscillations in Cyclic Gene Regulatory Networks*, University of Tokyo technical report METR 2010–11 (2010). [Full author/institutional report](https://www.keisu.t.u-tokyo.ac.jp/data/2010/METR10-11.pdf). Equation (1), Proposition 1, Lemmas 2–3, and Appendix B inspected. Its stronger determinant-based conclusion is recorded as a theorem-transfer lead, not silently generalized here.
