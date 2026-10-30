# Motor rows as finite-memory motion

The existing motor-row equations admit an exact memory representation and an energy balance. These clarify what drives the motion and give a sufficient region where all trajectories in the stated physical invariant domain converge to rest. They do not introduce another model or establish scientific priority for the derivations.

The [source audit](motor-rows-source-audit.md) defines the normalized variables and source assumptions. Write the complex population coefficient as $z_j=c_j+i s_j$, and let $V_j=X_j'$. The population equations become

```math
z_j'=1+i b+(-1+i V_j)z_j.
```

Multiplication by the integrating factor gives an exact finite-past solution:

```math
z_j(t)=e^{-(t-t_0)+i[X_j(t)-X_j(t_0)]}z_j(t_0)
+(1+i b)\int_{t_0}^t e^{-(t-u)+i[X_j(t)-X_j(u)]}\,du.
```

Its imaginary part contains the exponentially fading history of the row's own displacement:

```math
s_j(t)=\mathrm{Im}\,H_j(t,t_0)
+\int_{t_0}^t e^{-(t-u)}
\left[\sin\Delta X_j(t,u)+b\cos\Delta X_j(t,u)\right]du,
```

where $H_j$ is the first term in the previous expression and $\Delta X_j(t,u)=X_j(t)-X_j(u)$. The force equation is still exactly

```math
X_j'=-\rho X_j+r(s_j-\bar s),\qquad \sum_jX_j=0.
```

Thus each row responds to its own fading displacement history, coupled through the subtraction of the common response. On a bounded complete trajectory one may take $t_0\to-\infty$ and the homogeneous term vanishes. For an arbitrary finite initial condition, dropping that term would change the problem. This is the collective counterpart of the [two-row memory equivalence](motor-memory-priority-audit.md); it does not identify the full ensemble with an already studied scalar trapped-particle model.

## An exact energy balance

Differentiate the force equation and use the population equations. With $P=I-\mathbf{1}\mathbf{1}^{T}/N$ and componentwise multiplication denoted by $\odot$,

```math
V'=-(\rho+1)V-\rho X+rP(c\odot V).
```

The mean term does no work because $\sum_jV_j=0$. For the nonnegative auxiliary function

```math
E=\frac12\sum_j\left(V_j^2+\rho X_j^2\right),
```

the derivative is

```math
E'=\sum_j\left[r c_j-(\rho+1)\right]V_j^2.
```

This is an auxiliary mathematical energy, not a claim about the dimensional thermodynamic energy of the source apparatus. The population coefficient changes the instantaneous damping of row motion. Sustained motion requires the weighted feedback input to balance damping on average; its presence alone does not imply chaos.

## A sufficient convergence region

The source audit establishes the positively invariant row disks

```math
c_j^2+s_j^2\leq R^2,\qquad R=\sqrt{1+b^2}.
```

Within these disks, if $rR\lt\rho+1$, then

```math
E'\leq-\left(\rho+1-rR\right)\sum_jV_j^2.
```

All trajectories in this domain are bounded: the disks bound the population coefficients, and the stable linear displacement filter bounds each $X_j$. A trajectory therefore remains in a compact positively invariant set. LaSalle's invariance principle restricts its limit set to the largest invariant subset of $V=0$, implying $V\to0$. Since $z_j$ stays bounded, the forcing $iV_jz_j$ in the stable linear equation $z_j'+z_j=1+ib+iV_jz_j$ tends to zero. Consequently $z_j\to1+ib$, and the force equation then gives $X\to0$. All such trajectories converge to the unique resting equilibrium. The invariant subset of $V=0$ itself includes transient population relaxation; it must not be identified with the equilibrium alone.

This is a sufficient condition in the specified invariant domain, with a strict inequality. It is not a sharp chaos boundary, and no conclusion at equality or outside the region follows from this argument. For $b=1$, the sufficient condition is $r\sqrt{2}\lt\rho+1$, stricter than the local Hopf threshold $r=\rho+1$. The independently checked three-row candidate lies outside the exclusion region. Local equilibrium stability and global exclusion must not be confused.

## Reproduction and references

[The algebra checker](../research/audit_motor_rows_memory.py) compares the complex, velocity and energy identities against the independently implemented source force-constraint equations at 256 seeded states for two, three, four and eight rows. It performs no trajectory integration. [Recorded errors](../results/motor-rows/memory-structure.json) are floating-point checks of the derivations, not a substitute for their assumptions.

```bash
python research/audit_motor_rows_memory.py
```

1. F. Alouges, I. Anello, A. DeSimone, A. Lefebvre-Lepot and J. Levillain, “Some mathematical models for flagellar activation mechanisms,” *Mathematical Models and Methods in Applied Sciences* **35**, 2395–2424 (2025). [DOI](https://doi.org/10.1142/S0218202525500423); [primary full text](https://arxiv.org/html/2409.03506v2). Source equations and their physical scope.
2. The [exact-model priority audit](motor-memory-priority-audit.md) compares the two-row reduction with published finite-memory particle models and records access limits. Similar memory language alone does not establish equivalence of different systems.
