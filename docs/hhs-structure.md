# Exact single-cell HHS: boundedness and equilibrium controls

This note audits the unsaturated, autonomous five-variable single-cell model in equation S2 of [1], with the source's explicit choice of external force equal to zero. It proves an initial-condition-dependent global bound and supplies independent equilibrium calculations. It does not prove a chaotic attractor, exclude chaos at nonzero coupling, or establish scientific novelty. The source already reports spatiotemporal chaos in its spatial model; that model contains additional terms and is not this five-state ODE.

Reproduce the algebraic checks with `python research/audit_hhs_structure.py`. The script computes no trajectories and needs only NumPy and SciPy. Its complete numerical record, source-cache hash, and qualifications are in [structural-audit.json](../results/hhs/structural-audit.json).

## Exact normalization and assumptions

For positive mechanical coupling scale $\alpha$, set $x=(1-L)/\alpha$, $y=(1-L_0)/\alpha$, and $\mu=\alpha\beta$. In state order $(M,E,D,y,x)$, equation S2 becomes

```math
\dot M=a-(b+1)M+cM^2E,\qquad
\dot E=bM-cM^2E-DE,
```

```math
\tau_D\dot D=D_0-D+\mu Dx,\qquad
\tau_L\dot y=E-E_r-y,\qquad
\tau_R\dot x=y-x.
```

Here $E_r$ is a fixed reference concentration, distinguished from the initial value of $E$. All $a,b,c,D_0,E_r,\tau_D,\tau_L,\tau_R$ are positive, and $\mu\geq0$. The proof below permits arbitrary finite real initial $x,y$, nonnegative initial $M,E$, and strictly positive initial $D$.

This affine change of coordinates is exact. For fixed chemical parameters and time constants, the normalized dynamics depend on $\alpha,\beta$ only through their product. The physical lengths still require a choice of $\alpha$: $L=1-\alpha x$ and $L_0=1-\alpha y$. Positive concentrations do not imply positive lengths, and linear-elastic validity is a further physical restriction. This reduction does not automatically apply after adding the continuum model's cubic regularizers. The case $\alpha=0$ is singular for this coordinate change and should be treated in the original coordinates.

Time can be expressed in units of any one time constant, but doing so multiplies the chemical rates as well as changing the remaining time-constant ratios. One cannot change $\tau_R$ and hold all normalized chemical rates fixed while claiming an equivalent model.

## A global bound for the literal S2 model

The following is our derivation. It is a bound for each finite admissible initial state, not a claimed parameter-only absorbing set.

**1. Positivity and lower bounds.** At the chemical boundaries, $\dot M=a$, $\dot E=bM$, and $\tau_D\dot D=D_0$. Thus $M,E\geq0$ and $D>0$ persist. Define

```math
h=\min\{x(0),y(0),-E_r\},\qquad
C_h=1-\mu h,\qquad
\delta=\min\{D(0),D_0/C_h\}>0.
```

The two stable first-order filters imply $x,y\geq h$. Consequently, $\tau_D\dot D\geq D_0-C_hD$, and scalar comparison gives $D\geq\delta$.

**2. Bounds on the chemistry and filters.** For $N=M+E$,

```math
\dot N=a-M-DE\leq a-\min\{1,\delta\}N.
```

Therefore, with

```math
B=\max\left\{N(0),\frac{a}{\min\{1,\delta\}}\right\},
\qquad H=\max\{x(0),y(0),B-E_r\},
```

we have $0\leq M,E\leq B$ and $h\leq x,y\leq H$. These inequalities hold on any interval of local existence. The remaining equation satisfies $\tau_D\dot D\leq D_0+(-1+\mu H)D$, so $D$ cannot diverge in finite time. The locally Lipschitz vector field therefore has a global forward solution. The next step strengthens its finite-time bound to a time-uniform bound.

**3. A large-D barrier.** Put $C=1+\mu E_r$ and choose any $\rho>0$ such that $\rho bB<C/(4\tau_D)$. Define

```math
V=\log D+\frac{\mu}{\tau_D}(\tau_Ly+\tau_Rx)+\rho E.
```

Direct differentiation cancels both filter terms and gives the exact identity

```math
\dot V=\frac{D_0}{\tau_DD}-\frac{C}{\tau_D}
+\left(\frac{\mu}{\tau_D}-\rho D\right)E
+\rho bM-\rho cM^2E.
```

For

```math
D\geq D_*:=\max\left\{\frac{4D_0}{C},
\frac{\mu}{\rho\tau_D}\right\},
\qquad \dot V\leq-\frac{C}{2\tau_D}.
```

Meanwhile the already bounded correction $Z=V-\log D$ obeys

```math
Z_-:=\frac{\mu(\tau_L+\tau_R)}{\tau_D}h
\leq Z\leq
Z_+:=\frac{\mu(\tau_L+\tau_R)}{\tau_D}H+\rho B.
```

Whenever $V\geq\log D_{\ast}+Z_+$, we necessarily have $D\geq D_{\ast}$ and hence $\dot V<0$. A scalar barrier argument yields the explicit finite bound

```math
D(t)\leq
\exp\!\left(\max\{V(0),\log D_*+Z_+\}-Z_-\right)
\qquad(t\geq0).
```

Together these inequalities bound all five states. The bound is intentionally conservative: for the example $\mu=1.5$ and initial state $(1.05,3.43,0.1,0,0)$, the calculated bound on $\log D$ is about $209.96$. Such a bound establishes finiteness but is not an estimate of the observed amplitude. The function $V$ need not decrease at small $D$; this is not a global stability proof for an equilibrium.

For the common initial choice $x=y=0$, one has $h=-E_r$ and $\delta=\min\{D(0),D_0/(1+\mu E_r)\}$. No spatial saturation, numerical clipping, or replacement chemical law was used.

## Equilibrium polynomial and nearby distinct events

At a positive equilibrium,

```math
\frac{a}{b+1}<M<a,\qquad
E=\frac{(b+1)M-a}{cM^2},\qquad
D=\frac{a-M}{E},\qquad x=y=E-E_r.
```

The scalar equation $D[1-\mu(E-E_r)]=D_0$ reduces to

```math
cCM^3-[cCa+\mu(b+1)]M^2
+[(b+1)D_0+\mu a(b+2)]M-aD_0-\mu a^2=0,
\qquad C=1+\mu E_r.
```

The source defaults are $a=c=1$, $b=3$, $D_0=0.1$, $\tau_D=1.1$, $\tau_L=1.2$, and $\tau_R=1$. The uncoupled reference root and concentration are

```math
M_r=0.6137372763943942,\qquad E_r=3.8626272360560576.
```

The equilibrium $(M_r,E_r,D_0,0,0)$ persists for **every** $\mu$. Factoring its known root out of the cubic leaves a quadratic, which avoids inaccurate extraction of almost coincident roots. Independent polynomial and full-Jacobian controls distinguish these events:

| Event | Coupling $\mu$ | Linear result at source time constants |
|---|---:|---|
| Birth of two additional equilibria | 1.05624167401125 | Double root near $M=0.60898536$; one zero eigenvalue, but an unstable complex pair remains |
| Additional branch meets persistent reference equilibrium | 1.05664228907977 | Zero eigenvalue at the reference state; another pair has real part about $0.05967$ |
| Compressed branch gains linear stability | 1.08611380421898 | Pair crosses the imaginary axis near $\pm0.69146390i$ |
| Compressed branch loses linear stability | 2.42673506600070 | Pair crosses near $\pm0.73880131i$ |

The last two are **Hopf candidates from linear analysis**. Transversality was checked numerically, but no first Lyapunov coefficient or nonlinear periodic-branch continuation was computed here. The exact source-default equilibrium is stable at sampled $\mu=1.1,1.5,2$, and unstable at $\mu=3$. No classification of the resulting nonstationary motion follows from this alone.

The source's rounded value $\mu_c\simeq1.0566$ corresponds to its zero-determinant calculation at the persistent reference equilibrium. Our factorization resolves a slightly earlier saddle-node; these events are separated by approximately $4.00615\times10^{-4}$. The equilibrium locations and these zero-eigenvalue events do not depend on the positive time constants, although Hopf crossings and stability generally do.

## Primary-source reproduction questions

The inspected source is the full public **v3** manuscript and SI [1], not just the earlier preprint title. The source reports the spatial chaotic regime itself. The finite-dimensional model's exact implementation still deserves its own controls.

There is a checkable discrepancy between the **displayed** determinant expansion in SI S7–S8 and the matrix in S6. At the reference state, define $A=-(b+1)+2cM_rE_r$, $B=cM_r^2$, and $Q=b-2cM_rE_r$. Direct expansion of the displayed S6 matrix gives

```math
\tau_D\tau_L\tau_R\det J
=A(B+D_0+\mu E_rD_0)+BQ.
```

The displayed S7–S8 expression instead has $-BQ$. Taking that written expression literally gives a zero at $\mu\simeq-3.524773$, whereas the matrix gives the reported positive value $1.056642289$. The independent script checks the corrected expansion against numerical determinants. This locates a discrepancy between displayed formulas; it does **not** establish how the authors' numerical implementation was coded.

For fixed $E_r$ and source chemical parameters, at $\mu=0.8$ there is exactly one equilibrium, with $E=E_r$. Distinct exact constant asymptotes at that gain cannot all belong to this same vector field. Source Fig. S1 and its prose use $E_0$ also for sampled initial concentrations, while the model defines a reference value. This is a concrete reproduction question about reference versus initial concentration, the coupling convention, and whether plotted traces represent asymptotic constants. It is not resolved by reading several apparent plateau heights from a figure. No alternate reference rule has been silently substituted in our equations.

The source SI X discusses QR-based exponents in the spatial context; its written description alone does not establish a full asymptotic five-state spectrum for S2. An independent calculation should evolve the full trajectory-dependent Jacobian, check the time-averaged divergence against the sum of exponents, and check the neutral flow direction on bounded nonstationary recurrent trajectories. At zero feedback, the block-triangular system supplies an especially strong control.

## What is and is not excluded

With $\mu=0$ and $D=D_0$, the autonomous chemistry is planar; the other variables are stable driven filters. A stable chemical cycle has the five-state spectrum consisting of its neutral and radial exponents plus $-1/\tau_D,-1/\tau_L,-1/\tau_R$. Arbitrary positive initial $D$ relaxes exponentially to $D_0$.

For $\mu>0$, we have not proved a planar reduction or a monotone-system exclusion. The derivative $\partial\dot E/\partial M=b-2cME$ changes sign in the positive domain, and the $M,E$ chemistry adds a side loop to the negative feedback chain $E\to y\to x\to D\to E$. This blocks a direct invocation of the elementary monotone cyclic-feedback argument; it does not prove that every more sophisticated exclusion is impossible. Boundedness leaves periodic, quasiperiodic, transient, and chaotic behavior to be distinguished by further analysis.

The script verifies the Jacobian by complex-step differentiation and checks the barrier derivative, equilibrium polynomial, and determinant identities at 64 reproducibly selected finite states. Maximum absolute errors are approximately $1.4\times10^{-14}$, $1.4\times10^{-14}$, $9.1\times10^{-13}$, and $4.4\times10^{-16}$, respectively. These numerical identity checks support transcription accuracy; the boundedness argument itself is analytic.

## Reference

1. S. Mondal, P. Dewan, L. S. Kumar, and S. Sarkar, *Symmetry-protected phases in a 1D active solid with mechanochemical feedback*, [arXiv:2504.19655v3, full manuscript and SI](https://arxiv.org/html/2504.19655v3); *Physical Review Letters* **137**, 078401 (2026), [DOI:10.1103/p1bw-v87g](https://doi.org/10.1103/p1bw-v87g). Equations S2 and S5–S8, SI I.1, II, III, X, and Fig. S1 were examined. The public author manuscript was used; no access restriction was bypassed.
