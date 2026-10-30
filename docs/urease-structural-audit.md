# Structural audit of the original four-species urease model

Audit date: 30 October 2026 UTC. This note concerns the deterministic four-species equations, before the additional quasi-steady reductions in Straube et al. (2023). It adds source verification and elementary structural results to [the candidate screening](enzyme-feedback-screening.md). The derivations below are ours; they are not claims of previously unpublished mathematics.

**Decision:** the published parameter choice has an exact attracting three-dimensional hyperplane. We have not established that all recurrent dynamics on that hyperplane are planar, or proved a global exclusion of chaos there. A stronger exclusion is available when all three transport rates are equal. A stiff numerical control of the full four-species equations remains appropriate at the published unequal transport rates.

## Exact model and source control

Let `S`, `H`, `P`, and `Q` denote urea, hydrogen ions, ammonia, and ammonium concentrations in mol/L. Time is in seconds. Equations S1a–d of the 2021 supporting information [1] give

```math
v(S,H)=\frac{v_{\max}S}{K_M+S}
\frac{H}{H+H^2/K_{E1}+K_{E2}},\qquad
j=k_2PH-k_{2r}Q,
```

```math
\dot S=-v+k_S(S_e-S),\qquad
\dot H=-j+k_H(H_e-H),
```

```math
\dot P=2v-j-kP,\qquad \dot Q=j-kQ.
```

The activity expression above is algebraically the source expression `1/(1+H/KE1+KE2/H)` for positive `H`, and has a smooth extension to `H=0`. It is useful for evaluating the boundary without division by zero. All parameters below are strictly positive.

| Parameter | Published value | Unit |
|---|---:|---|
| `vmax` | 1.85e-4 | mol/L/s |
| `KM` | 3e-3 | mol/L |
| `KE1` | 5e-6 | mol/L |
| `KE2` | 2e-9 | mol/L |
| `k2` | 4.3e10 | L/mol/s |
| `k2r` | 24 | 1/s |
| `kS` | 1.4e-3 | 1/s |
| `kH` | 9e-3 | 1/s |
| `k` | 1.4e-3, set equal to `kS` | 1/s |
| `Se` | 3.8e-4 | mol/L |
| `He` | 1.3e-4 | mol/L |
| `S(0)` | 5e-5 | mol/L |
| `H(0)` | 1e-5 | mol/L |

The 2021 main text and its two-page supporting information do **not specify the initial product concentrations `P(0)` and `Q(0)`**. Do not silently label a chosen product initialization as the published one. A reasonable reproduction can state its choice explicitly and compare the final periodic state after a long transient, including different product starts.

The source reports a deterministic period of **17.79 ± 0.01 minutes** and approximately **pH 4.2–8.3** in the reference example. These are useful controls for the full model. The stochastic period scatter and loss of regularity in small vesicles are not evidence of deterministic chaos. The 2023 article [2] analyzes reductions and their accuracy; an exact mapping inside one reduced model must not be mistaken for an exact reduction of these original four equations.

Source access: the 2021 arXiv v1 PDF contains the supporting information on its last two pages; the published article PDF was also inspected through the authors' institutional repository. The 2023 institutional record confirms the publication but provides no final PDF, so the equation audit used arXiv v1. Cache hashes and source URLs are recorded in [the machine-readable audit](../results/urease/structural-audit.json).

## Physical domain, boundedness, and interior dynamics

The nonnegative orthant is forward invariant. At a zero coordinate the corresponding derivative is nonnegative; in particular, at `S=0` and `H=0` the positive reservoir supply points inward. Negative concentrations or a manually imposed pH floor do not belong to this ODE model.

Two reaction-invariant linear combinations are useful:

```math
N=2S+P+Q,\qquad Z=H+Q.
```

Direct addition, with no quasi-steady assumption, gives

```math
\dot N=2k_SS_e-2k_SS-k(P+Q),
```

```math
\dot Z=k_HH_e-k_HH-kQ
=k_HH_e-kZ-(k_H-k)H.
```

Consequently,

```math
\dot N\leq 2k_SS_e-\min(k_S,k)N,\qquad
\dot Z\leq k_HH_e-\min(k_H,k)Z.
```

These scalar inequalities bound every concentration: `S <= N/2`, `P,Q <= N`, and `H <= Z`. They also rule out finite-time blow-up for finite nonnegative initial conditions. For each fixed positive parameter set the model is dissipative in this elementary boundedness sense. Boundedness alone does not rule out chaos.

There are also eventual positive lower bounds, although they can be extremely loose. Since `v <= (vmax/KM) S`,

```math
\dot S\geq k_SS_e-(k_S+v_{\max}/K_M)S.
```

Choose eventual strict upper bounds `Nbar,Zbar` from the inequalities above. Then

```math
\dot H\geq k_HH_e-(k_H+k_2\overline N)H.
```

Thus `S,H` eventually exceed positive constants. The continuous rate `v` has a positive minimum on the resulting compact rectangle. Then

```math
\dot P\geq 2v_{\min}-(k+k_2\overline Z)P,
\qquad
\dot Q\geq k_2P_{\min}H_{\min}-(k_{2r}+k)Q
```

give eventual positive bounds for `P,Q` as well. One can take slightly smaller constants than the scalar comparison limits. Compact recurrent dynamics therefore stay in the interior for each fixed strictly positive parameter set. This is a mathematical concentration statement: it does not ensure that the corresponding molecular copy numbers are large enough for a deterministic approximation.

The four-dimensional divergence is

```math
\mathrm{div}\,F=-v_S-k_S-k_H-k_{2r}-2k-k_2(P+H)<0.
```

Volume contraction is another useful consistency check, but negative divergence in four or three dimensions does not exclude a chaotic attractor.

## Exact reduction at equal product and substrate transport

For the source convention `k=kS`,

```math
\dot N=k(2S_e-N),\qquad
N(t)=2S_e+[N(0)-2S_e]e^{-kt}.
```

Every compact recurrent invariant set lies on `N=2Se`. This is an exact attracting three-dimensional hyperplane, with a decoupled normal perturbation equation `delta N'=-k delta N`. It is not yet a planar reduction.

On that hyperplane use coordinates `(S,H,Z)` and reconstruct

```math
Q=Z-H,\qquad P=2S_e-2S-Z+H.
```

The exact system is

```math
\dot S=-v(S,H)+k(S_e-S),
```

```math
\dot H=k_{2r}(Z-H)-k_2H(2S_e-2S-Z+H)+k_H(H_e-H),
```

```math
\dot Z=k_HH_e-kZ-(k_H-k)H.
```

At the published point `kH>k`, the last variable retains a finite memory of `H`. On a bounded complete trajectory it satisfies

```math
Z(t)=\frac{k_HH_e}{k}-(k_H-k)
\int_0^\infty e^{-ks}H(t-s)\,ds.
```

Eliminating this history integral does not produce a two-variable ordinary differential equation.

## A stronger exclusion when all transport rates agree

If `kH=kS=k`, then both

```math
\dot N=k(2S_e-N),\qquad \dot Z=k(H_e-Z)
```

are closed. Recurrent dynamics lie on the plane `N=2Se, Z=He`. In `(S,H)` the physical domain is the convex polygon

```math
S\geq0,\quad 0\leq H\leq H_e,\quad
P=2S_e-2S-H_e+H\geq0.
```

The planar divergence equals

```math
-v_S-2k-k_{2r}-k_2(P+H)<0.
```

Bendixson's criterion excludes a nonconstant periodic orbit in the interior of this simply connected domain. The reservoir-driven inward conditions preclude a recurrent boundary cycle. Planarity also excludes a chaotic limit attractor. This parameter slice is a useful negative control; it does **not** contain the published oscillatory point.

This result can be placed in a complete classification of closed affine reaction-invariant coordinates. The two internal stoichiometric directions, in state order `(S,H,P,Q)`, are

```math
\nu_v=(-1,0,2,0),\qquad \nu_j=(0,-1,-1,1).
```

Every linear combination that cancels both reactions has the form

```math
L=aN+bZ,\qquad w=(2a,b,a,a+b).
```

To satisfy a scalar relaxation law `L'=c-lambda L` for all states, this vector must also be an eigenvector of the diagonal transport matrix `diag(kS,kH,k,k)`. With positive, nondegenerate kinetics, the possibilities are:

| Transport relation | Available reaction-invariant scalar relaxation |
|---|---|
| `kS=k` | `N` |
| `kH=k` | `Z` |
| `kS=kH=k` | Every `aN+bZ`; two independent coordinates |
| Neither equality | None in this two-dimensional linear space |

Hence there is no second independent closure of this affine type at the source point. This does not classify nonlinear first integrals or nonlinear invariant surfaces.

## Fast protonation: a useful structure, not an exact finite-rate theorem

Freeze `S` and the reaction totals `N,Z`. The fast protonation subsystem reduces to

```math
\dot H_{\mathrm{fast}}=k_{2r}(Z-H)
-k_2H(N-2S-Z+H).
```

The physical interval is `max(0,Z-(N-2S)) <= H <= Z`, with `N-2S >= 0` and `Z >= 0`. At its endpoints the fast vector field points inward. Its derivative is

```math
\frac{\partial\dot H_{\mathrm{fast}}}{\partial H}
=-k_{2r}-k_2(P+H)\leq-k_{2r}<0.
```

There is exactly one physical fast equilibrium and no fold of this frozen fast equilibrium branch. At the reference point the frozen relaxation rate is at least `24 /s`, whereas the transport rates are much smaller. This supports investigating a normally attracting slow manifold in a suitable singular limit.

It does **not** establish a globally attracting two-dimensional invariant manifold at the finite published rates. Such a claim would require an explicit small-parameter formulation, compactness and regularity conditions, and adequate control of coupling derivatives and the relevant invariant set. The number `-24 /s` is a frozen subsystem bound, not a computed Lyapunov exponent of the full model. We did not prove that all recurrent dynamics are planar at `kH>k`.

Nor is setting the time derivative of a product, or of the product `P H`, to zero an identity of the full ODE. Those additional closures remove a genuine finite relaxation time. They may approximate the regular orbit very accurately while being insufficient for a universal no-chaos conclusion.

## Why the simple competitive-system argument does not transfer

For the exact three-dimensional coordinates `(S,H,Z)`,

```math
\frac{\partial\dot H}{\partial S}=2k_2H>0,\quad
\frac{\partial\dot H}{\partial Z}=k_{2r}+k_2H>0,\quad
\frac{\partial\dot Z}{\partial H}=k-k_H.
```

At the source point the last two derivatives have opposite signs. Their product is unchanged by diagonal sign flips or permutations, so this representation cannot be made wholly cooperative or wholly competitive by such transformations. Moreover `partial S'/partial H=-v_H` changes sign at `H=sqrt(KE1 KE2)=1e-7 M`; the reference cycle passes across this value.

Thus the straightforward orthant monotonicity argument used for the three-state Milster reactor is unavailable here. This is a limitation of that argument, not a proof that no other cone, change of variables, or specialized theorem could apply.

## Chemical validity and the next numerical check

The four-variable equations are already a chemical simplification. The 2021 supporting information starts from an earlier eight-species description, omits carbon dioxide, and eliminates hydroxide dynamics using a quasi-steady assumption with hydroxide exchange neglected. Its numerical comparison concerns the reference regime. It does not validate arbitrary parameter extrapolation to extreme pH, membrane transport ratios, or enzyme concentrations. The 2023 work separately examines water ionization in reduced models; that is not a license to add an independent hydroxide state without specifying the altered network.

The membrane supplies substrate and acid from a fixed reservoir. There is therefore no finite-fuel monotonicity exclusion analogous to the closed artificial-cell model. The model also assumes well-mixed concentrations and effective membrane transport; electrochemical and spatial effects are not additional hidden states in these equations. At nanoscale volumes, discrete stochastic reaction events matter. A deterministic attractor claim must be separated from experimentally observable molecular dynamics.

The next control should integrate all four concentrations with a stiff independent solver, state the product initialization, recover the roughly 17.79-minute regular cycle, verify the exact `N(t)` law and positivity, and compare multiple starts and tolerances. The equal-transport slice above offers a structural negative control. Only after that is a parameter search justified. This audit itself contains no numerical trajectory or positive Lyapunov result.

## References

1. A. V. Straube, S. Winkelmann, C. Schütte, and F. Höfling, “Stochastic pH oscillations in a model of the urea–urease reaction confined to lipid vesicles,” *J. Phys. Chem. Lett.* **12**, 9888–9893 (2021), [doi:10.1021/acs.jpclett.1c03016](https://doi.org/10.1021/acs.jpclett.1c03016). [Published PDF through the authors' institution](https://publications.imp.fu-berlin.de/2630/1/acs.jpclett.1c03016.pdf); [arXiv:2110.03465v1, including supporting information](https://arxiv.org/pdf/2110.03465v1).
2. A. V. Straube, S. Winkelmann, and F. Höfling, “Accurate Reduced Models for the pH Oscillations in the Urea–Urease Reaction Confined to Giant Lipid Vesicles,” *J. Phys. Chem. B* **127**, 2955–2967 (2023), [doi:10.1021/acs.jpcb.2c09092](https://doi.org/10.1021/acs.jpcb.2c09092). [Inspected arXiv version 2212.14503v1](https://arxiv.org/html/2212.14503v1); [institutional publication record](https://publications.imp.fu-berlin.de/2959/).
3. A. V. Straube, S. Winkelmann, F. Höfling, and M. Engel, “Unfolding the geometric structure and multiple timescales of the urea–urease pH oscillator,” [arXiv:2508.07275v1](https://arxiv.org/abs/2508.07275), used in the preceding screening for the explicitly planar model. An institutional [project publication list](https://www.zib.de/research/projects/synchronization-and-geometric-structures-stochastic-biochemical-oscillators) lists acceptance in *Journal of Nonlinear Science*; the final journal version was not inspected in this block.
