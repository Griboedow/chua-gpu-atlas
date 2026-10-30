# Late screening IV: a rigid microfluidic network

Audit: 30 October 2026 UTC. **The exact four-state model is ready for a reproducibility check, but the broad fluidic mechanism already has published chaos.** No trajectories were simulated in this block and no new-attractor claim is made.

**Subsequent check:** the [literal-model CPU study](microfluidic-research.md) did not recover the periodic source control. A [final-supplement audit](microfluidic-source-audit.md) confirms that the printed coefficients match the manuscript, and [structural checks](microfluidic-structure.md) identify physical-domain limits. The text below records the original extraction stage; its proposed simulation is now superseded by those results.

## Source and scope

Case, Angilella, and Motter [1] model a pressure-driven rigid network using five channel flow rates and one vortex size. Two flow-conservation constraints leave **four independent states**. The primary supplement contains the complete equations and fitted constants. Its bifurcation diagram, Fig. S5, contains periodic and multistable regimes. The separate, isolated chicane model in the main text is planar and must not be mistaken for the full network.

The later primary paper by Shao, Angilella, and Motter [2] already demonstrates chaotic flow-rate dynamics in a **different serialized network**, using two-dimensional Navier–Stokes simulations. At inlet pressures 105 and 70 Pa it reports an estimated largest exponent of 31.9 inverse seconds. This is Eulerian flow-rate chaos, distinct from chaotic tracer motion in an otherwise periodic flow.

Thus, the remaining question is narrow: does the **literal 2020 four-state closure** support physically admissible chaos? This audit found no verified answer for those exact equations, but does not establish that an answer is unpublished. The 2025 result sharply reduces priority for a first-chaos contribution.

## Exact source equations

All quantities below are dimensionless. Dots denote the source's scaled time. Source locations: supplement S1–S2 and S5–S13, PDF pages 18 and 22 of the combined author manuscript [1]. Define

```math
\eta(Q_1)=p+qQ_1,\qquad
\xi(Q_1)=\frac{(r^*)^3\eta(Q_1)-Q_3^*}{r^*},\qquad
R(r)=L_b+\gamma(r-r_b)^2.
```

The dependence of **both** coefficients on $Q_1$ follows from the two source equations. Freezing $\xi$ while varying $\eta$ would violate their stated zero-flow condition. The published equations are

```math
\dot Q_1=12\left(\frac{P_1^{\mathrm{in}}-P_3}{L_1}-Q_1\right),
\qquad
\dot Q_2=12\left(\frac{P_2^{\mathrm{in}}-P_4}{L_2}-Q_2\right),
```

```math
\dot Q_3=12\left(\frac{\zeta P_3-P_4}{L_3}
-\frac{R(r)}{L_3}Q_3+k_3\frac{Q_3Q_5}{Q_1}\right),
```

```math
\dot Q_4=12\left(\frac{P_4}{L_4}-\frac{\alpha+\beta Q_4}{12}Q_4\right),
\qquad
\dot Q_5=12\left(\frac{P_3}{L_5}-Q_5+k_5\frac{Q_3Q_5}{Q_1}\right),
```

```math
\dot r=\epsilon\left[Q_3-Q_3^*-\eta(Q_1)(r-r^*)^3+\xi(Q_1)(r-r^*)\right],
\qquad Q_1=Q_3+Q_5,\qquad Q_2=Q_4-Q_3.
```

The plus signs on the $k_3$ and $k_5$ terms are literal source signs, visually checked in S7 and S9. They should not be silently changed on the assumption that every term labeled a loss must have a negative sign.

## Our explicit elimination of the algebraic pressures

This is an algebraic reduction of the source equations, not an added physical assumption. Choose $x=Q_1$, $y=Q_2$, $z=Q_3$ and $r$ as independent states. Then $Q_4=y+z$ and $Q_5=x-z$. On the regular domain $x\ne0$, set

```math
m=\frac{z(x-z)}{x},\qquad
a_3=-\frac{R(r)z}{L_3}+k_3m,\qquad
a_4=-\frac{[\alpha+\beta(y+z)](y+z)}{12},\qquad
a_5=-(x-z)+k_5m.
```

Differentiating the two flow constraints gives a constant two-by-two linear solve:

```math
\begin{pmatrix}
L_1^{-1}+\zeta L_3^{-1}+L_5^{-1}&-L_3^{-1}\\
-\zeta L_3^{-1}&L_2^{-1}+L_4^{-1}+L_3^{-1}
\end{pmatrix}
\begin{pmatrix}P_3\\P_4\end{pmatrix}
=
\begin{pmatrix}
P_1^{\mathrm{in}}/L_1-x-a_3-a_5\\
P_2^{\mathrm{in}}/L_2-y-a_4+a_3
\end{pmatrix}.
```

After this solve, use the first three flow equations and the vortex equation. No pressure dynamics or extra storage variable is introduced. With $A_0=L_1^{-1}+L_5^{-1}$ and $B_0=L_2^{-1}+L_4^{-1}$, the matrix determinant is

```math
A_0B_0+\frac{A_0+\zeta B_0}{L_3}>0
```

for positive channel lengths and $\zeta>0$. Thus, this pressure elimination is nonsingular throughout that parameter domain. At the source constants its determinant is 0.003138888888888888 and its two-norm condition number is 4.75035469359741.

A 1,000-point **algebraic** check, using seeded random positive flow states and pressures, reconstructed all five derivatives and verified both differentiated conservation identities. The maximum residual divided by the largest derivative component was $2.103\times10^{-15}$. This check did not integrate the system or demonstrate boundedness.

Substituting the coefficient identity also gives a useful independent vortex formula:

```math
\dot r=\epsilon\left[z-r\left(\frac{Q_3^*}{r^*}
+\eta(x)(r-r^*)(r-2r^*)\right)\right].
```

In particular, at $r=0$ the derivative is exactly $\epsilon z$. This is a useful implementation check and a conditional positivity result when the transverse flow remains nonnegative.

## Parameters and dimensional interpretation

The source width is $w=5\times10^{-4}$ m, density $\rho=1000$ kg/m³, and kinematic viscosity $\nu=10^{-6}$ m²/s. Physical channel lengths in metres are $(0.001,0.006,0.001,0.010,0.005)$. Under the source scaling, the constants are

| Quantity | Dimensionless value |
| --- | --- |
| $L_1,L_2,L_3,L_4,L_5$ | $24,144,24,240,120$ |
| $\alpha,\beta$ | $31,0.2565$ |
| $\zeta,L_b,\gamma$ | $1.08,42,7.35$ |
| $r_b,r^{\ast},Q_3^{\ast}$ | $4.1,5.5,25$ |
| $\epsilon,k_3,k_5$ | $10.9,15.75,12$ |
| $p,q$ | $0.12,8\times10^{-6}$ |

The dimensional Forchheimer constants in Fig. S4 are $1.24\times10^8$ m⁻² and 513 m⁻¹; these yield the table's $\alpha$ and $\beta$. Do not combine the isolated-chicane Fig. 6 constants with the full-network S13 fit.

The conversions are

```math
t_{\mathrm{physical}}=0.25\tau\ \mathrm{s},\qquad
P_{\mathrm{physical}}=0.004P\ \mathrm{Pa},\qquad
r_{\mathrm{physical}}=\frac{w}{12}r,
\qquad Q_{\mathrm{physical}}=\nu Q.
```

Here physical $Q$ is area flow per unit out-of-plane depth. Numerically, one dimensionless flow unit is one microlitre per second per millimetre of depth. The source Reynolds convention is $\mathrm{Re}=2Q$ for a positive channel flow. Equal inlets are the source setting; 250 Pa corresponds to $P_1^{\mathrm{in}}=P_2^{\mathrm{in}}=62500$ and lies in the periodic window pictured in Fig. S5. This is a proposed first reproduction setting, **not** a cycle reproduced in this audit. No precise initial condition or tabulated source trace is supplied by the recovered supplement.

## Physical guards and priority limits

The source's branch-flow assumptions suggest monitoring $x>0$, $y\ge0$, $0\le z\le x$, and $r\ge0$. These conditions have not been proved forward invariant. In particular, the closure has a division by $x$ and its Forchheimer term is fitted for forward flow. Do not repair a flow reversal by clipping or replacing $Q_4^2$ with $Q_4|Q_4|$ and then call it the same model. Flag exits from the fit's operating regime separately from numerical blow-up.

The source measured the chicane resistance fit over approximately 60–400 micrometres of vortex size. Extrapolation outside that interval needs an explicit qualification. Varying fitted coefficients independently is a mathematical sensitivity study; it does not automatically correspond to a realizable new channel geometry. A physically meaningful initial exploration can vary constant inlet pressures while retaining the fitted geometry, and must first recover its published periodic behavior.

The model is autonomous and continuously powered by maintained inlet pressure. It has no prescribed oscillatory input or finite fuel budget. Nevertheless, neither positivity nor a dissipative-looking term proves an absorbing set or a strange attractor. A computed four-state chaotic orbit would also need validation against the hydrodynamic approximation; a three-dimensional plot would be a projection.

Forward-citation metadata returned eight works. The most relevant, [2], was read in full. It changes the geometry, removes the strongly nonlinear obstacle-laden branch, uses two-mode reduction for an elementary network, and demonstrates chaos after serialization using the full fluid solver. It therefore excludes novelty for the **mechanism**, but is not an established equivalence to the four-state S5–S12 closure. The other seven records were screened at metadata level; this is not a complete full-text citation audit. Exact-title, DOI, author/mechanism, and Lyapunov queries yielded no verified exact-S5–S12 chaos example in this bounded pass. That retrieval outcome is not proof of absence.

## References and provenance

1. D. J. Case, J.-R. Angilella, and A. E. Motter, *Spontaneous oscillations and negative-conductance transitions in microfluidic networks*, **Science Advances 6**, eaay6761 (2020). [DOI](https://doi.org/10.1126/sciadv.aay6761), [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7220308/), [author manuscript including complete supplement](https://arxiv.org/pdf/2006.09400). Equations were checked against rendered PDF pages 18, 21, and 22; source bifurcation diagram checked on page 23.
2. Y. Shao, J.-R. Angilella, and A. E. Motter, *Emergent oscillations and chaos in noncompliant microfluidic networks*, **Physical Review Fluids 10**, 054401 (2025). [DOI](https://doi.org/10.1103/PhysRevFluids.10.054401), [full primary author manuscript](https://arxiv.org/html/2505.00068v1), especially section IV and Figs. 6–7.
3. [OpenAlex forward-citation query](https://api.openalex.org/works?filter=cites:W3024896709&per-page=200), retrieved 30 October 2026. Index records may have incomplete coverage or inconsistent dates.

Source files remain under ignored `outputs/candidate_research/`; source prose, figures, and author code are not redistributed in this MIT repository. The companion [JSON audit](../results/candidate-audits/late-candidate-screening-4.json) records equations, parameters, citations, provenance, and the algebraic check. No integration, GPU job, author contact, or commit was made in this screening block.
