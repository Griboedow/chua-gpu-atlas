# Collective molecular motors: an exact old-model revisit

Status: 30 October 2026. This is a source and equivalence audit, with a separate reproducible periodic-control calculation. It is **not a new model, mechanism, or attractor claim**.

The starting lead was the two-row flagellar activation model of Alouges et al. (2025). Its five-state first-harmonic equations contain an exact closed three-state core. That core is linearly equivalent, after time and amplitude scaling, to Oriola et al. (2014), Eq. (10). The older paper explicitly tested for chaos and instead reported complex periodic and subharmonic orbits. This makes a further calculation a narrowly defined revisit of an existing question.

## Exact equations suitable for independent implementation

For positive parameters $\rho,r$, and dimensionless time $\tau$, use the state $(h,Y,Z)$:

```math
\begin{aligned}
h'&=[-\rho+r(1-Z)]h-rY,\\
Y'&=h(1-Z)-Y,\\
Z'&=hY-Z.
\end{aligned}
```

The analytic Jacobian is

```math
J=\begin{pmatrix}
-\rho+r(1-Z)&-r&-rh\\
1-Z&-1&-h\\
Y&h&-1
\end{pmatrix}.
```

[motor_pair.py](../research/motor_pair.py) provides the field, Jacobian, an independent reconstruction from the original equations, both coordinate maps, and the periodic-control command. No clipping, extra forcing, delays, stochastic terms, or additional feedback are introduced.

### Connection to Oriola et al. (2014)

The source uses $(\bar a_1,\bar b_1,\bar X)$, time $\tau=\alpha t$, and

```math
\begin{aligned}
\bar a_1'&=-\bar a_1-2\pi\rho\bar X\bar b_1-\mu\bar a_1\bar b_1,\\
\bar b_1'&=-\Omega-\bar b_1+2\pi\rho\bar X\bar a_1+\mu\bar a_1^2,\\
\bar X'&=-\rho\bar X-\frac{\mu}{2\pi}\bar a_1.
\end{aligned}
```

The invertible transformation, for $\rho,\Omega>0$, is

```math
h=2\pi\rho\bar X+\mu\bar a_1,
\qquad Y=\frac{\bar a_1}{\Omega},
\qquad Z=\frac{\bar b_1+\Omega}{\Omega},
\qquad r=\mu\Omega.
```

Consequently only two parameters control this normalized vector field. Changing $\mu$ and $\Omega$ at fixed product rescales motor amplitudes; initial conditions must also be transformed when comparing basins. It does not create a third independent dynamical control parameter.

The source transition rates obey $\omega_1=(\beta/2)[1+\cos(2\pi\xi/\ell)]$ and $\omega_2=\alpha-\omega_1$. Their nonnegativity requires $0\leq\Omega=\beta/(2\alpha)\leq1/2$; the nonsingular normalization above uses $0<\Omega<1/2$. The source definitions are $\rho=K/(\lambda_f\alpha)$ and $\mu=\pi^2U/(\lambda_f\ell^2\alpha)$, where $\lambda_f$ denotes mechanical friction. A constant chemical activity is a maintained input of this model, not a simulated finite ATP stock.

The closed first-harmonic equations are exact for the specified sinusoidal potential and constant total transition rate. They are not an exact reduction of every molecular-motor model, or of a whole bending flagellum in fluid.

### Exact connection to Alouges et al. (2025)

Use its Eq. (33), fixing chemical activity. Write $a=a_0>0$, $\zeta=2\pi k/(\eta\ell)$, $\lambda_m=2\pi^2U/(\eta\ell)$, and

```math
c=\frac{p_1^c+q_1^c}{2},\qquad
s=\frac{p_1^s-q_1^s}{2},\qquad
H=\zeta x+\lambda_m s.
```

Direct addition and subtraction give the closed subsystem

```math
\dot c=-ac+Hs+a_1/\ell,\qquad
\dot s=-as-Hc,\qquad
\dot x=-\ell H/(2\pi).
```

For the source's active sign $a_1<0$, set

```math
\tau=at,\quad h=H/a,\quad
\Omega=-a_1/a,\quad
Y=\ell s/\Omega,\quad
Z=(\ell c+\Omega)/\Omega,
\quad \rho=\frac{k}{\eta a},\quad
r=-\frac{\lambda_m a_1}{a^2\ell}.
```

These produce exactly the normalized field above. The remaining combinations $d=(p_1^c-q_1^c)/2$ and $e=(p_1^s+q_1^s)/2$ obey

```math
\dot d=-ad+He,\qquad
\dot e=-ae-Hd+b_1/\ell.
```

They do not feed back on the core. Two solutions of this passive subsystem over the same core trajectory have squared separation decaying exactly as $\exp(-2at)$. Thus the five-state notation does not provide two extra independent sources of instability.

The source rates are nonnegative when $\sqrt{a_1^2+b_1^2}\leq a/2$. For the illustrated choice $b_1=a_1$, this imposes $\Omega\leq1/(2\sqrt2)$. The core itself is independent of $b_1$; the admissible physical interpretation is not.

## Mathematical checks before a survey

Let $R^2=Y^2+(1-Z)^2$. The equations give

```math
\frac12(R^2)'=-R^2+(1-Z)\leq-R^2+R.
```

The disk $R\leq1$ is forward invariant. It is a useful smaller attracting region, not a test declaring every transient state outside it unphysical. Put $q=h-rY$. Then

```math
q'=-\rho(q+rY).
```

Inside the disk, $|Y|\leq1$, so $|q|\leq r$ is also invariant. Their product defines a compact invariant solid cylinder in $(q,Y,Z)$, with $|h|\leq2r$. More generally, the radial inequality bounds $R$ and the stable scalar filter bounds $q$: the smooth field has no finite-time escape from finite initial states. These bounds do not prove chaotic recurrence.

There is one equilibrium, $(0,0,0)$. Its characteristic polynomial is

```math
(s+1)\,[s^2+(\rho+1-r)s+\rho].
```

It is asymptotically stable for $r<\rho+1$. The Hopf boundary is $r=\rho+1$, with angular frequency $\sqrt\rho$. The 2025 source supplies a local supercritical-Hopf result. This says nothing about the global dynamics far above threshold. The symmetry $(h,Y,Z)\mapsto(-h,-Y,Z)$ is inherited from the source; two visually prominent lobes alone are therefore weak evidence of a new geometric type.

Two additional elementary deductions help constrain a search. They are our calculations from the verified equations; no mathematical-priority claim is made.

First, the coercive function

```math
V=\frac{q^2}{2\rho r}+\frac{Y^2+Z^2}{2}
```

has the exact derivative

```math
V'=-\frac{q^2}{r}+(r-1)Y^2-Z^2.
```

For $0<r<1$ it is negative definite. At $r=1$ the zero-derivative set is $q=Z=0$, and its only invariant subset is the origin, because $q'=-\rho Y$ there. LaSalle's argument on compact sublevel sets therefore gives global attraction to the origin for **all finite initial states when $0<r\leq1$**. This is a useful guaranteed negative control. The interval $1<r<\rho+1$ has a locally stable equilibrium, but the preceding global argument does not cover it.

Second, when $r>\rho+1$, the two roots of the quadratic factor have positive real parts, whereas the remaining eigenvalue is $-1$. The invariant axis $h=Y=0$, on which $Z'=-Z$, is tangent to the stable eigenspace. Uniqueness of the local stable manifold identifies it with this axis; backward uniqueness then gives the complete stable manifold. Every nonzero point on the axis diverges as time tends to minus infinity. Thus **there is no homoclinic orbit to the origin above the Hopf threshold**. The unique equilibrium also excludes a heteroclinic network between distinct equilibria. This rules out those particular equilibrium-based mechanisms, not chaotic sets or crises involving periodic orbits.

## Published controls and parameter interpretation

| Source case | $\rho$ | $\mu$ | $\Omega$ | $r=\mu\Omega$ | Published behavior |
| --- | ---: | ---: | ---: | ---: | --- |
| 2014 Fig. 2 | 0.05 | 200 | 0.2 | 40 | Complex periodic orbit; explicit Lyapunov check |
| 2014 Fig. 3a | 0.7 | 1,000 | 0.4 | 400 | Periodic cusped oscillation |
| 2014 Fig. 3b | 4.5 | 1,000 | 0.4 | 400 | Periodic oscillation |
| 2014 Fig. 4a, left | 6 | 63,096 | 0.425 | 26,815.8 | Subharmonic periodic oscillation |
| 2014 Fig. 4a, right | 3 | 6,500 | 0.4 | 2,600 | Subharmonic periodic oscillation |
| 2014 Fig. 4d, left | 6 | 19,000 | 0.4 | 7,600 | Subharmonic periodic oscillation |

The source's physical estimates are roughly $\rho=10^{-2}$ to $10^2$ and $\mu=10^3$ to $10^5$. These arise from ranges of stiffness and friction, with shared dependencies. A rectangular mathematical survey is not automatically a calibrated experimental parameter range. In particular the complete-PDE comparison labelled $\Omega=4$ uses a different rate law; **do not transfer 4 into the reduced model whose admissible $\Omega$ is below 1/2**.

The Fig. 5 phase-diagram initial state is $(\bar a_1,\bar b_1,\bar X)=(0,0,0.074)$, hence $(h,Y,Z)=(2\pi\rho\,0.074,0,1)$. The Fig. 2 initial state was not recovered. Our independent control instead uses equilibrium motor harmonics with a stated displacement perturbation, $(h,Y,Z)=(2\pi\rho\,0.01,0,0)$. Coexisting cycles mean a different start need not reproduce the exact plotted orbit.

For the 2025 source, Table 1 gives $\ell=10$ nm, $k_BT=0.0042668$ nN nm, $\eta=10^{-7}$ kg/s, $k=9.5\times10^{-5}$ kg/s², $U=10k_BT$, $c=1000$ in $a_0=c\sqrt{\Omega_{\rm src}}$, and $\Omega_0=15k_BT$. Inserting the numerical energy $\Omega_0=0.064002$ in nN nm recovers its printed $d=-56.4588$ in $a_1=b_1=d\Omega_{\rm src}$. Treating $\Omega_0$ as the bare number 15 does not. This audit keeps the source activity and the 2014 normalized $\Omega$ distinct. The 2014 dimensionless controls avoid this unit convention.

## Verification and priority status

The implementation has an independent Eq. (7) reconstruction and analytic-Jacobian checks. The reproducible control command is

```bash
python research/motor_pair.py --audit-output outputs/motor-algebra-audit.json
python research/motor_pair.py --control-output outputs/motor-periodic-control.json
python research/motor_pair.py --refine-control outputs/motor-periodic-control.json
```

Its output records two-solver agreement, long section returns, full-state shooting closure, an independent Radau closure/monodromy, and transverse Floquet multipliers. The actual run is retained in [source-audit.json](../results/motor-pair/source-audit.json). A stable reached cycle is a positive implementation control, not a claim that every basin is periodic.

The first control reaches a cycle of period approximately 141.492 in source time. DOP853 and Radau agree over the first 20 time units within $2.7\times10^{-9}$ in $h$ and $1.3\times10^{-10}$ in the other coordinates. After 2,000 time units, return intervals have stabilized to a few parts in $10^8$. The initial shooting has closure residuals about $1.8\times10^{-7}$ for DOP853 and $6.6\times10^{-6}$ for independent Radau; its transverse multiplier modulus is about 0.066, corresponding to exponent about $-0.0192$. These are numerical accuracy limits, not exact closure. The section is transverse, with $h'\simeq9.79$; long-orbit conditioning remains relevant.

A tighter unconstrained root-refinement attempt was manually stopped after roughly six minutes without replacing the saved control. The implementation now makes exactly two Newton corrections in the algebraically equivalent memory coordinate and records the independent residuals. Neither a solver-success message nor a very small requested tolerance is treated as proof of achieved accuracy. The initial diagnostics remain preserved alongside any refinement.

The bounded check finished: it **did not improve closure**, with residuals approximately $6.0\times10^{-6}$ (DOP853) and $2.1\times10^{-5}$ (Radau). Both solvers still give transverse multipliers close to $-0.00422\pm0.06581i$, agreeing within about $1.6\times10^{-6}$, and corresponding exponents approximately $-0.01922$. The robust conclusion is a numerically supported stable periodic control with the stated accuracy bounds. No further precision-chasing integrations were performed.

The source audit recovered the 2014 author-uploaded full paper, the matching detailed thesis chapter, the complete 2025 arXiv v2, and the 2025 multiscale follow-up. Bounded forward discovery returned 11 records for the 2014 DOI and four for the 2025 DOI, including duplicate preprint/journal records. These counts are retrieval metadata, not an exhaustive literature census. Closest later branches include spatial dynein/flagellar models, a three-state power-stroke model, and the already audited DCR motor-contact network; none was shown equivalent to the core by those titles alone.

The 2014 paper's negative chaos check is direct prior work. No exact-core later positive-chaos demonstration was verified in this bounded search. That is an unresolved priority check, not evidence that a positive result would be unpublished. The 2026 cilia/photochemical-rod follow-up was inspected through primary publisher text and a partially retrievable full-paper text; its complete coefficient audit remains unfinished.

## References

1. D. Oriola, H. Gadêlha, C. Blanch-Mercader, J. Casademunt, “Subharmonic oscillations of collective molecular motors,” *EPL* **107**, 18002 (2014). [DOI](https://doi.org/10.1209/0295-5075/107/18002); [author-uploaded full text](https://www.researchgate.net/publication/263668359_Subharmonic_oscillations_of_collective_molecular_motors).
2. D. Oriola, *Self-organization and cooperativity*, doctoral thesis, Chapter 6, especially pp. 115–124. [Institutional full text](https://diposit.ub.edu/bitstreams/9d37fdc6-1352-492b-940f-bee5ddcbe162/download).
3. F. Alouges, I. Anello, A. DeSimone, A. Lefebvre-Lepot, J. Levillain, “Some mathematical models for flagellar activation mechanisms,” *Mathematical Models and Methods in Applied Sciences* **35**, 2395–2424 (2025). [DOI](https://doi.org/10.1142/S0218202525500423); [full arXiv v2](https://arxiv.org/html/2409.03506v2).
4. I. Anello, F. Alouges, A. DeSimone, “Beating of eukaryotic flagella via Hopf bifurcation of a system of stalled molecular motors,” *European Journal of Mechanics A/Solids* **114**, 105729 (2025). [DOI](https://doi.org/10.1016/j.euromechsol.2025.105729); [full preprint](https://arxiv.org/html/2412.06067v1).
5. F. Jülicher, J. Prost, “Spontaneous oscillations of collective molecular motors,” *Physical Review Letters* **78**, 4510–4513 (1997). [DOI](https://doi.org/10.1103/PhysRevLett.78.4510).
6. D. Oriola, H. Gadêlha, J. Casademunt, “Nonlinear amplitude dynamics in flagellar beating,” *Royal Society Open Science* **4**, 160698 (2017). [Full primary text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5383814/).
7. D. Agostinelli et al., “Spontaneous oscillations in eukaryotic cilia and photo-responsive rods,” *Journal of the Mechanics and Physics of Solids* **210**, 106529 (2026). [DOI and primary preview](https://doi.org/10.1016/j.jmps.2026.106529).

All source files remain in the ignored research cache. Original articles and theses retain their own copyright/licensing; the independent implementation and our derivations do not relicense those sources.
