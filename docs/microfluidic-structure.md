# Structural audit of the literal four-state microfluidic closure

This note concerns the exact S5–S13 closure of Case, Angilella, and Motter [1], with two constant positive inlet pressures and the source's fitted coefficients. It retains the **positive** junction terms and the stated dependence of both vortex coefficients on $Q_1$. The results are algebraic: a counterexample to global forward-flow invariance, a conditional energy bound, and equilibrium controls. No trajectories were computed here and no new-attractor claim is made.

The broader fluidic mechanism already has published chaotic behavior in a different network [2]. Neither that result nor the planar isolated-chicane model is a theorem about this four-state closure. The equations and source provenance are recorded separately in [the source audit](late-candidate-screening-4.md).

Reproduce the calculations with `python research/audit_microfluidic_structure.py`. The [machine-readable record](../results/microfluidic/structural-audit.json) includes constants, equilibrium spectra, the boundary counterexample, root-search limitations, and the energy identity check.

## The forward-flow cone is not globally invariant

Write $x=Q_1$, $y=Q_2$, $z=Q_3$, $w=Q_5=x-z$, and $u=Q_4=y+z$. The proposed forward-flow domain is

```math
x>0,\qquad y\geq0,\qquad 0\leq z\leq x,\qquad r\geq0.
```

Let $A_0=L_1^{-1}+L_5^{-1}$, $B_0=L_2^{-1}+L_4^{-1}$, and let $\Delta$ be the positive determinant of the two-by-two pressure-elimination matrix. On the boundary $z=0$, define

```math
H_2(y)=\frac{P_2}{L_2}-y+\frac{(\alpha+\beta y)y}{12}.
```

The exact eliminated equations give

```math
\left.\frac{\dot z}{12}\right|_{z=0}
=\frac{\zeta B_0P_1/L_1-A_0H_2(y)}{L_3\Delta}.
```

Since $\beta>0$, this expression is negative for sufficiently large nonnegative $y$, despite positive inlet pressures. At $P_1=P_2=62500$, the sign changes at $y\simeq64.48200099$. A concrete admissible boundary state is

```math
(x,y,z,r)=(1000,100,0,4.1),\qquad
\dot z=-1442.4778761.
```

By continuity at this regular point, strictly positive $z$ sufficiently near this boundary also exits the forward-flow domain. Thus positive inlet pressures do not by themselves justify forward-flow invariance. The boundary $r=0$ has $\dot r=\epsilon z\geq0$ only as long as $z$ remains nonnegative. The term $z(x-z)/x$ is singular at $x=0$; no extension, clipping, or reversal repair has been introduced.

## An energy bound conditional on staying in that domain

An exact weighted kinetic expression cancels both internal node pressures:

```math
V=\frac{1}{24}\left[
\zeta L_1x^2+L_2y^2+L_3z^2+L_4(y+z)^2+\zeta L_5(x-z)^2
\right].
```

Using all five source momentum equations and the differentiated flow constraints gives

```math
\begin{aligned}
\dot V={}&\zeta P_1x+P_2y-\zeta L_1x^2-L_2y^2-R(r)z^2\\
&-\frac{L_4}{12}(\alpha+\beta u)u^2-\zeta L_5w^2
+L_3k_3zm+\zeta L_5k_5wm,
\end{aligned}
\qquad m=\frac{zw}{z+w}.
```

The channel weights are $(\zeta,1,1,1,\zeta)$; omitting these weights would leave an uncancelled contribution from the source's $\zeta P_3$ term. The two positive junction terms are retained explicitly.

While $w,y,z\geq0$ and $x=z+w>0$, one has $0\leq m\leq z$. Put

```math
a_w=\zeta L_5,\qquad B_w=\zeta L_5k_5,\qquad
C_z=L_3k_3+\frac{B_w^2}{2a_w},\qquad d=\frac{L_4\beta}{12}>0.
```

Dropping nonpositive terms and using $B_wwz\leq a_ww^2/2+B_w^2z^2/(2a_w)$ yields

```math
\dot V\leq \zeta P_1(w+z)+P_2y
-\frac{a_w}{2}w^2-L_2y^2+C_zz^2-dz^3.
```

Young's inequality absorbs the linear $w,y$ forcing. The remaining cubic polynomial in $z$ has a finite maximum, so explicitly

```math
\dot V\leq K-\frac{a_w}{4}w^2-\frac{L_2}{2}y^2-z^2,
```

```math
K=\frac{(\zeta P_1)^2}{a_w}+\frac{P_2^2}{2L_2}
+\max_{z\geq0}\{\zeta P_1z+(C_z+1)z^2-dz^3\}.
```

The maximum is attained at the positive root of $3dz^2-2(C_z+1)z-\zeta P_1=0$. Write $V=v^{\mathsf T}Gv/24$ for $v=(w,y,z)$. The positive-definite Gram matrix is

```math
G=\begin{pmatrix}
\zeta(L_1+L_5)&0&\zeta L_1\\
0&L_2+L_4&L_4\\
\zeta L_1&L_4&\zeta L_1+L_3+L_4
\end{pmatrix}.
```

It follows that $\dot V\leq K-kV$, where

```math
k=\frac{24\min\{a_w/4,L_2/2,1\}}{\lambda_{\max}(G)}>0.
```

Thus $V(t)\leq V(0)e^{-kt}+(K/k)(1-e^{-kt})$ on every interval where the trajectory stays in the regular forward-flow domain. For a trajectory that remains there for all future time, this also supplies a parameter-dependent asymptotic flow bound.

The exact vortex equation can be written

```math
\dot r=\epsilon\left[z-r\left(\frac{Q_3^*}{r^*}
+(p+qx)(r-r^*)(r-2r^*)\right)\right].
```

If $z\leq Z$, then for $r\geq4r^*$ its derivative is at most $\epsilon[Z-3pr^3/8]$. A corresponding bound on $r$ follows. This uses the fixed source values $p,q>0$ and $x>0$.

These estimates are conservative: at equal dimensionless pressures 62500, the calculated asymptotic bound on each of $w,y,z$ is about $1.89\times10^5$, far above the typical source flows. They establish conditional finiteness, not realistic amplitude predictions. Since the forward-flow cone is not invariant and $x=0$ is singular, this is **not** a global absorbing-set theorem for the unrestricted four-state vector field. It also does not establish that the fitted resistance remains valid over the bounded region.

## Independent physical-equilibrium reduction

At an equilibrium, $P_3=P_1-L_1x$ and $P_4=P_2-L_2y$. The fifth-channel equation implies

```math
(L_1+L_5)x^2-[P_1+L_5(1+k_5)z]x+L_5k_5z^2=0.
```

The fourth-channel equation, with $u=y+z$, becomes

```math
\frac{L_4\beta}{12}u^2+
\left(L_2+\frac{L_4\alpha}{12}\right)u=P_2+L_2z.
```

It has a unique positive root for each $z\geq0$. The condition $y\geq0$ further requires

```math
0\leq z\leq z_{\max}(P_2),\qquad
\frac{L_4}{12}(\alpha+\beta z_{\max})z_{\max}=P_2.
```

For the source pressure pair, $z_{\max}=65.40787672$. Since this is much less than $P_1/L_1$, the smaller quadratic root for $x$ lies below $z$ and is inadmissible. The larger root is used; the script explicitly refuses this simplification when that domain condition fails.

The transverse-channel equation now specifies the resistance required at each admissible $z$:

```math
R_{\mathrm{req}}(z)=
\frac{\zeta(P_1-L_1x(z))-[P_2-L_2y(z)]+L_3k_3m(z)}{z}.
```

Where $R_{\mathrm{req}}\geq L_b$, test both branches $r=r_b\pm\sqrt{(R_{\mathrm{req}}-L_b)/\gamma}$, retaining $r\geq0$, and impose the exact vortex nullcline. This reduces the equilibrium search to scalar brackets without integrating trajectories or changing any coefficient.

At $P_1=P_2=62500$ it finds the physical equilibrium

```math
(x,y,z,r)=(694.316102856,\ 36.4295124299,\ 24.8474262027,\ 5.29726846417),
```

with eigenvalues in source inverse-time units

```math
-11.85304726,\qquad -52.39025782,\qquad
-8.37345678\pm18.32414871i.
```

It is strongly linearly stable. At 25 equally spaced equal-pressure settings from 40000 to 100000, corresponding to 160–400 Pa, the same bracketed search finds one physical equilibrium at each setting and all computed eigenvalues have negative real parts. Their largest real parts range approximately from $-9.44$ to $-7.96$. No nearby Hopf crossing was found in these sampled controls. This finite bracketing scan does not certify that all roots were found, exclude a narrow intervening bifurcation, or exclude a coexisting cycle or attractor.

## Vortex contraction and the unresolved source reproduction

The literal coefficient identity gives

```math
\xi(x)=(r^*)^2(p+qx)-\frac{Q_3^*}{r^*}.
```

With the source constants, $\xi>0$ requires $x>3782.87002254$. At the source-pressure equilibrium above, $\eta=0.1255545288$ and $\xi=-0.7474300486$. Thus the vortex equation, with the other variables held fixed, is strictly contracting in $r$ there:

```math
\frac{\partial\dot r}{\partial r}
=\epsilon\{\xi(x)-3\eta(x)(r-r^*)^2\}<0.
```

This quantitative consequence should be checked against the source's coefficient convention and stated periodic example. A negative $\xi$ is not itself a four-state no-oscillation theorem. We do not increase $q$, freeze $\xi$, reverse the junction signs, or otherwise modify the closure to obtain a desired cycle.

No further dimensional reduction has been established. The two flow constraints have already been used to obtain four states. The derivative of the transverse-flow equation with respect to vortex size changes sign across $r=r_b$, while the vortex is driven by $z$. Therefore, a fixed-sign monotone-feedback argument cannot be applied globally in this domain without additional work. The conditional bound and local stability do not settle global chaos.

## Checks and references

The source-momentum energy identity was checked independently against the gradient of $V$ and the eliminated vector field at 128 seeded forward-domain states; the maximum scaled residual was $6.8\times10^{-12}$. Equilibrium Jacobians were evaluated with complex-step differentiation of this independent field implementation. The scalar-root residual at the source equilibrium is approximately $1.0\times10^{-11}$ in the original source rates.

1. D. J. Case, J.-R. Angilella, and A. E. Motter, *Spontaneous oscillations and negative-conductance transitions in microfluidic networks*, **Science Advances 6**, eaay6761 (2020). [Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7220308/), [author manuscript with complete SI](https://arxiv.org/pdf/2006.09400), [DOI](https://doi.org/10.1126/sciadv.aay6761).
2. Y. Shao, J.-R. Angilella, and A. E. Motter, *Emergent oscillations and chaos in noncompliant microfluidic networks*, **Physical Review Fluids 10**, 054401 (2025). [Primary author manuscript](https://arxiv.org/html/2505.00068v1), [DOI](https://doi.org/10.1103/PhysRevFluids.10.054401). Its serialized geometry and full fluid simulations are distinct from the closure audited here.
