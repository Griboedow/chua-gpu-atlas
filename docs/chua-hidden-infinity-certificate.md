# A computer-assisted existence argument for infinite hidden periodic coexistence

## Statement and scope

Consider the dimensionless continuous three-segment Chua equations

```math
\begin{gathered}
\dot x=\alpha[y-x-h(x)],\quad\dot y=x-y+z,\quad\dot z=-\beta y,\\
h(x)=m_1x+\frac{m_0-m_1}{2}(|x+1|-|x-1|).
\end{gathered}
```

with γ = 0 and the exact decimal parameters

```math
\alpha=13.803545907139778,\quad m_0=-1.5179789642291142,\quad m_1=-0.10174607415658823.
```

**Theorem.** There is a parameter β* in (29.51661848407718, 29.51661848408720) at which a dissipative periodic saddle has a quadratic homoclinic tangency, unfolding nondegenerately with β. Arbitrarily close to β* there are individual β values for which the flow has infinitely many geometrically distinct, simultaneous, hidden attracting periodic orbits.

The argument combines interval computations, the analytic steps below and the classical parameter theorem. It is not a complete Lean formalization, an externally reviewed result, a decimal witness β∞, or a certificate of chaotic attraction. The numerical examples used to find the connection are not substituted for the exact-manifold argument.

**Meaning of the conclusion.** One fixed choice of α, m₀ and m₁ is used throughout. For every ε > 0 there is an individual β∞ with 0 < |β∞ − β*| < ε for which infinitely many attracting periodic flow orbits coexist. Each orbit has a basin avoiding an open neighborhood of every equilibrium. The claim concerns the number of simultaneous periodic attractors, not the number of connected surfaces in a global atlas.

The proof has three dependencies: **T**, a quadratic tangency of exact invariant manifolds (§§1–3); **H**, uniform equilibrium-neighborhood exclusion (§4); and **N**, the localized parameter theorem (§5). T alone does not prove hiddenness; H alone does not prove existence of an attracting cycle. The [review guide](chua-review-guide.md#3-infinite-hidden-coexistence) gives the corresponding data fields and executable checks.

| Notation | Meaning |
|---|---|
| q = (y,z) | Physical coordinates on the section x = 1 |
| ξ = (u,s) | Coordinates centered at the saddle, in a fixed approximate eigenbasis |
| Fβ | Ordinary first return, two affine flights |
| Pβ = Fβ² | Two first returns; the map to which the parameter theorem is applied |
| I | The tested tangency parameter box, β₀ ± 5 × 10⁻¹² |
| J | [29.516618475, 29.516618493], the larger interval on which hiddenness is uniform |
| R | The rational matrix of the Lyapunov quadratic form in §4 |

![Connection, certified signs and the implication for hidden coexistence](assets/chua-followup-hidden.png)

## 1. An analytic return branch

Let Fβ be the first return to x = 1, with positive x velocity, along an outer-right flight and then a central flight. All open outer flight intervals satisfy x > 1; all open central flight intervals satisfy 0 < x < 1. Interval Newton enclosures give positive crossing times and exclude zero crossing speed. Adaptive time covers certify the complete flights, including collars at their endpoints. The branch is therefore the intended first return and is analytic on an open neighborhood of the tested compact sets. The derivative jump of the vector field at the diode surface does not invalidate this branchwise analyticity.

The affine solutions are c + exp(At)(q − c). Their parameter and state derivatives, including derivatives of the implicit event time, are enclosed to second order. The [event-map implementation](../research/chua_interval_lobe.py) uses the same validated flight machinery as the earlier classical certificate, with section x = 1 and no symmetry inversion.

A parameter-uniform Krawczyk inclusion encloses a fixed saddle p(β) of Fβ. At the numerical center its section coordinates are approximately (y,z) = (0.282219168263034, 1.998952350199077). The multipliers of Fβ are near −1.76031507122 and −0.000210541402077. Write Pβ = Fβ². Conservative enclosures over the entire tangency box are

| Quantity for Pβ | Enclosing interval |
|---|---:|
| Unstable multiplier | [3.09870905, 3.09870925] |
| Stable multiplier | [4.4316 × 10⁻⁸, 4.4339 × 10⁻⁸] |
| Determinant | [1.3732 × 10⁻⁷, 1.3739 × 10⁻⁷] |

Pβ is a strongly dissipative, orientation-preserving local diffeomorphism. It represents two ordinary first returns; it is not a new physical model.

## 2. Exact local invariant graphs

Use ξ = V⁻¹(q − p(β)), with the fixed decimal matrix V in [connection-center-14-2.json](../results/chua-followup/hidden-search/connection-center-14-2.json). V need not be an exact eigenbasis. In the square |u|, |s| ≤ δ = 10⁻⁷, write Fβ = (U,S) in these centered coordinates. The certificate checks

```math
U_u\lt-1.7,\quad |U_s|\lt0.000012,\quad
|S_u|\lt0.000012,\quad |S_s|\lt0.000223.
```

All state Hessian entries have magnitude below M = 60. The centered parameter derivatives have magnitude below Q = 10⁻⁶. The derivative enclosure is intersected with its mean-value enclosure about ξ = 0, using the Hessian over the entire square. This tightens the bound without sampling the square or assuming exact diagonalization.

The graph-transform inequalities derived in the [earlier tangency proof](chua-infinity-certificate.md#3-bounds-on-the-exact-invariant-manifolds) apply also when both saddle multipliers of Fβ are negative. With

```math
L=10^{-5},\qquad K=100,\qquad \Gamma=2\times10^{-6},
```

they preserve the graphs s = ψ(u,β) and u = χ(s,β), through the origin, with slopes bounded by L, second derivatives by K and parameter derivatives by Γ. The checked upper bounds are:

| Graph-transform quantity | Unstable graph | Stable graph |
|---|---:|---:|
| Slope | 7.061 × 10⁻⁶ | 7.061 × 10⁻⁶ |
| Contraction in graph values | 0.000223001 | 0.588236 |
| Second derivative | 20.770 | 35.296 |
| Parameter derivative | 1.000457 × 10⁻⁶ | 1.764712 × 10⁻⁶ |

Every bound is strictly within its prescribed class. The expanding projection covers the whole local interval, and the stable image stays inside the local square. The uniform graph transforms and smooth stable-manifold theorem identify their limits with the exact parameter-dependent invariant manifolds.

## 3. The global connection and four strict signs

The numerical line-based solve supplies only the center, approximately

```math
a_0=6.973622953248532\times10^{-8},\qquad
\beta_0=29.516618484082188876.
```

The complete decimal centers are retained in the data. Test the box

```math
|a-a_0|\le10^{-15},\quad |\beta-\beta_0|\le5\times10^{-12},\quad |v|\le8\times10^{-13}.
```

Because L|a| < 8 × 10⁻¹³, this box contains the entire exact departure graph v = ψ(a,β). Compose 36 returns, giving 72 affine flights. There are 14,555 accepted global time cells in each working-precision run. At least one iterate leaves the local square, while the final stable coordinate lies strictly inside it. Thus the excursion is nontrivial and both infinite tails belong to the local manifold domains.

Let (U₃₆,S₃₆) denote the centered composite map. Define its exact normal gap by

```math
g(a,\beta)=U_{36}(a,\psi(a,\beta),\beta)
-\chi\bigl(S_{36}(a,\psi(a,\beta),\beta),\beta\bigr).
```

The chain rule retains all manifold slopes, curvatures and parameter derivatives within their proved bounds. In particular, neither ψ nor χ is set to zero in the certificate. Conservative outward-rounded results from the complete box are

| Quantity | Enclosing interval |
|---|---:|
| g on lower β face, every a | [−9.334 × 10⁻¹², −9.295 × 10⁻¹²] |
| g on upper β face, every a | [9.295 × 10⁻¹², 9.334 × 10⁻¹²] |
| gₐ on lower a face, every β | [−36.754, −22.344] |
| gₐ on upper a face, every β | [22.344, 36.754] |
| gₐₐ throughout the box | [2.2495 × 10¹⁶, 3.6604 × 10¹⁶] |
| gβ throughout the box | [1.862685, 1.862993] |

For example, the parameter derivative includes both moving manifolds:

```math
g_\beta=U_\beta+U_v\psi_\beta
-\chi_s(S_\beta+S_v\psi_\beta)-\chi_\beta.
```

All derivatives here belong to the centered composite map; its saddle center p(β) is also differentiated. Omitting that dependence or replacing the exact graphs by their tangent lines would test a different gap.

For each β, strict convexity and the two a-face signs give a unique interior stationary point a(β). The minimum g(a(β),β) changes sign between the two parameter faces. Its derivative is gβ > 0, so this minimum vanishes at exactly one β in the tested box. This uniqueness concerns the specified connection and departure interval. At that zero, g = gₐ = 0, gₐₐ > 0 and gβ > 0. This is a quadratic homoclinic tangency with a nondegenerate parameter unfolding. Since 36 is even, the connection also belongs to Pβ = Fβ².

### 3.1 Controlling interval overestimation

Direct composition in fixed physical coordinates could not enclose the complete excursion: the interval widths grew enough to make a central flight unresolved. The successful implementation transports first and second derivative enclosures in a new rational basis at each return. A floating QR factorization merely chooses each basis; Arb encloses its inverse and every subsequent change of coordinates. No interval is replaced by its midpoint in the certificate. Centered Taylor intersections retain the parameter dependence and reduce wrapping. The rejected fixed-coordinate attempt remains in the evidence.

This is distinct from numerical seed refinement, where midpoint iterations are explicitly permitted because that stage makes no existence claim. The final fixed-point and manifold enclosures restore all uncertainties.

## 4. Uniform exclusion of every equilibrium neighborhood

The hiddenness certificate uses a larger β ball centered at the complete β₀ above with radius 10⁻⁸. In particular it covers the closed interval J = [29.516618475, 29.516618493], with the full tangency box strictly inside J. The slopes satisfy m₀ < −1 and m₁ > −1. The equilibria are the origin and

```math
E_\pm=\pm(x_e,0,-x_e),\qquad x_e=\frac{m_1-m_0}{1+m_1}\gt1.
```

An interval-Newton enclosure for the central characteristic polynomial gives one positive real eigenvalue λ and a stable pair σ ± iω, uniformly in β. A parameter-dependent real eigenbasis has first row (1,1,0). In these coordinates the unstable scalar grows as exp(λt), while the stable plane rotates and contracts as exp(σt).

Take a stable-coordinate disk of radius 10⁻⁶. While |u| ≤ 1/2, the first state coordinate satisfies |x| ≤ 1/2 + 10⁻⁶, so the trajectory remains in the exactly linear central region. Any point in a sufficiently small origin ball either has u = 0 and tends to the origin, or exits at u = ±1/2 with stable coordinates still in that disk. The entire positive exit rectangle is enclosed through its first x = 1 crossing and then for 160 time units in the outer affine region; 2,702 time cells verify that this outer flight never returns through x = 1. The negative exit follows by odd symmetry.

A fixed rational symmetric matrix R is positive definite, and so is −(A₁ᵀR + RA₁), uniformly over J. Around E₊ choose

```math
\mathcal E_+=\{X:(X-E_+)^TR(X-E_+)\le\ell\},\qquad
\ell=\frac{(x_e-1)^2}{2(R^{-1})_{11}}.
```

This ellipsoid lies strictly inside x > 1. Its quadratic energy decreases, and the entire transported exit rectangle lies inside it. The mirrored ellipsoid has the same property around E₋. The calculation provides the following conservative Euclidean neighborhoods, at every β ∈ J:

| Center | Radius | Destination of every point |
|---|---:|---|
| Origin | 2.53 × 10⁻⁷ | Origin or one of the two outer equilibria |
| Either outer equilibrium | 0.00322 | That equilibrium |

Thus these three open balls are disjoint from the basin of every attracting periodic orbit. This is a neighborhood statement about all points, not an unstable-eigenvector sampling test. It is the hidden-attractor definition needed for the theorem. The certificates repeat at 80 and 100 digits and retain the rational R, its positive principal minors, the spectral enclosures and every path check. The JSON field storing R is `outer_lyapunov.P`; it does not denote the return map Pβ.

**Boundedness of the model.** The same matrix also gives a global forward-boundedness check, separately from dissipativity of the return saddle. Write the complete field as

```math
\begin{gathered}
\dot X=A_1X+b\mathrm{sat}(x),\qquad
b=-\alpha(m_0-m_1)e_1,\\
\mathrm{sat}(x)=\tfrac12(|x+1|-|x-1|).
\end{gathered}
```

For V(X) = XᵀRX and Q = −(A₁ᵀR + RA₁), positive definiteness gives

```math
\dot V\le-\lambda_{\min}(Q)\|X\|^2+2\|Rb\|\|X\|.
```

Since |sat(x)| ≤ 1 and Q is uniformly positive on the compact interval J, V decreases outside a bounded ellipsoid. Every forward trajectory is therefore bounded for these parameters. This uses the full three-segment field, rather than continuing the outer affine equation through the central region. It is a qualitative consequence of the already verified matrix inequalities, not an additional numerical attractor search.

## 5. From the certified tangency to infinite hidden coexistence

We use the parameter Newhouse–Robinson theorem: a sufficiently smooth fixed one-parameter family unfolding a quadratic homoclinic tangency of a dissipative surface saddle has arbitrarily nearby parameter intervals containing parameters with infinitely many simultaneous sinks. See [Robinson (1983)](https://doi.org/10.1007/BF01206892) and the formulation in Theorem 2.6 of [Catsigeras, Cerminara and Enrich, author manuscript](https://www.fing.edu.uy/~eleonora/Recopilacion/Archivos/InfinitosPozos_Premat-2008-108.pdf), [published version (2011)](https://doi.org/10.3934/dcds.2011.29.693).

The hypothesis correspondence is explicit:

| Parameter-theorem hypothesis | Verification for this family |
|---|---|
| A fixed C³ one-parameter family of surface diffeomorphisms | β varies with α and both slopes fixed; transverse affine flights give analytic local maps; localization is addressed below. |
| A saddle with positive stable and unstable multipliers and product below one | The enclosures for Pβ = Fβ² in §1. |
| Quadratic contact of the exact local manifolds after a nontrivial excursion | The graph bounds, domain containments and g = gₐ = 0, gₐₐ > 0 in §§2–3. |
| Generic unfolding in the given parameter | The minimum gap has derivative gβ > 0; parameter reversal is permitted to match the theorem's sign convention. |
| The resulting sinks belong to the original flow | Every iterate remains in the agreement region of the return map and its extension, as explained below. |

Analyticity on the validated branches supplies C³ regularity; a bound only on second derivatives would not by itself supply this hypothesis. The theorem is applied to the specified β family, rather than to an unspecified perturbation of the Chua equations.

The compact-surface formulation requires a localization argument for this partial return map. The saddle and complete homoclinic orbit form a compact set K: the two tails converge to the saddle, and only finitely many excursion points lie outside its neighborhood. The certified transverse flights and local saddle domain give analytic return branches and local inverses on a neighborhood of K. Pβ is injective there by uniqueness of the flow and of preceding section crossings, and preserves orientation.

Choose one closed disk containing the saddle and both sufficiently short tails, and finitely many pairwise disjoint disks around the remaining excursion points, with closures in this valid region. Injectivity makes their images disjoint disks. Their orientation-preserving restrictions extend across collars and the complementary punctured sphere to a smooth surface diffeomorphism. Isotopy extension supplies the small parameter variations while preserving agreement on smaller disks covering K. This is an extension of the map germ; it is not an assertion that the original global return is defined everywhere.

Fix a smaller neighborhood Ω of K contained in the agreement region. The local saddle returns and the finite transition can be chosen with all intermediate images inside Ω. The localized thick-horseshoe and sink construction uses only finite words in these pieces, so all its periodic orbits remain in Ω. This is the domain argument detailed in [the preceding coexistence proof](chua-basin-and-coexistence.md#32-localization-and-extension-argument); the extension outside Ω is never used by the resulting orbits. Agreement only at the saddle and tangency point would be insufficient: it is required along every intermediate iterate.

Shrink the parameter interval to remain inside J. The parameter theorem then gives infinitely many simultaneous sinks of the actual Chua return map at some individual β values arbitrarily near β*. A sink of Pβ gives an attracting periodic flow orbit. A flow orbit has finitely many crossings in the compact transverse return region, so at most finitely many Pβ cycles represent the same flow orbit. Infinitely many map cycles therefore yield infinitely many geometrically distinct attracting flow cycles. Section 4 excludes every equilibrium neighborhood from each basin; consequently all these cycles are hidden.

The infinite existence conclusion comes from the theorem, not from enumerating cycles. The exact β∞ is not computed. Neither infinite coexistence at β* itself nor an interval on which an infinite collection persists is claimed. Component tolerances and smooth modifications of the ideal diode require separate analysis.

## 6. Controls and remaining verification boundary

The 140- and 160-digit tangency certificates pass all local, graph, path, face-sign, curvature, unfolding and dissipativity conditions. Independent 90-digit mpmath matrix exponentials and event solves recover the numerical center. Finite differences with two step sizes give curvature approximately 2.9549291146 × 10¹⁶ and parameter derivative 1.86283893632, inside the interval bounds. Direct DOP853 integration checks selected full returns. Moving the whole β box by ±10⁻⁹ makes the required parameter-face sign fail. The original classical tangency certificate also passes after the shared map changes.

These are separate numerical implementations and negative controls, not independent mathematical review. Python/FLINT correctness, the analytic graph-transform argument, smooth localization and the cited parameter theorem remain part of the trust basis. The existing Lean development does not verify this new result. The [evidence index](../results/chua-followup/README.md) lists the exact commands and files.

The [record-level review](../results/chua-followup/review/record-check.json) recomputes the graph inequalities, displayed signs, Lyapunov positivity and published radii from the retained bounds. It also checks prediction scores and amplitude coverage. Fresh certificate replays and a second matrix-exponential calculation are recorded in the [revision audit](../results/chua-followup/review/README.md). A successful replay is evidence about the implementation; it does not replace review of the analytic implications in §§2–5.
