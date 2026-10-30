# Molecular motors and a trapped particle with sinusoidal memory

Status: 30 October 2026. Equation-level literature audit, without new integrations. **The motor core has an exact interpretation as the zero-inertia limit of an existing trapped pilot-wave model.** This audit did not identify a published chaotic attractor specifically in that three-state limit; it does not establish its absence or scientific novelty.

## Exact coordinate change

Starting from the [motor core](motor-pair-source-audit.md), with positive parameters, define

```math
x=\frac{h-rY}{\rho},\qquad c=1-Z,\qquad s=-Y.
```

Direct substitution gives

```math
v=x'=-h=-\rho x+rs,\qquad
c'=1-c-vs,\qquad s'=-s+vc.
```

The inverse is explicit:

```math
h=\rho x-rs,\qquad Y=-s,\qquad Z=1-c.
```

Thus this is an invertible affine change of the same motor vector field, not an extra feedback law. Setting the complex memory variable to $z=c+is$ gives

```math
z'=(-1+iv)z+1.
```

The independent [algebra checker](../research/audit_motor_memory.py) verifies the inverse, transformed motor field, complex equation and zero-inertia parent field on 128 declared random cases. Maximum scaled discrepancies are below `8.6e-13`; the [saved check](../results/motor-pair/memory-equivalence-check.json) records the seed and normalization. This checks implementation consistency alongside the explicit derivation and performs no trajectory integration.

For any finite starting time, its exact integrating-factor formula is

```math
z(t)=e^{-(t-t_0)+i[x(t)-x(t_0)]}z(t_0)
+\int_{t_0}^{t}e^{-(t-u)+i[x(t)-x(u)]}\,du.
```

For a compatible complete past history, the homogeneous term disappears and therefore

```math
x'=-\rho x+r\int_{-\infty}^{t}
e^{-(t-u)}\sin[x(t)-x(u)]\,du.
```

This is overdamped motion in a harmonic trap, driven by exponentially fading sinusoidal path memory. A complete-history memory satisfies $|z|\leq1$. Arbitrary ODE initial memories retain the explicitly displayed decaying homogeneous term; an arbitrary initial condition must not silently be identified with a prescribed past trajectory.

## Primary comparison

Valani's 2022 paper, Sec. V.A, Eq. (12), explicitly gives the **finite-inertia** sinusoidal-wave model in a harmonic potential [1]. Renaming its variables to avoid collisions with the motor coordinates, the source equations are

```math
q'=V,\qquad
\kappa V'=P-V-kq,\qquad
P'=-P+\beta V-VQ,\qquad
Q'=-Q+VP.
```

Here the source has $\sigma=1/\kappa$, $r=\beta$, and $b=1$. Its parameters represent inertia relative to drag, wave forcing relative to drag, and trap strength. The exact algebraic zero-inertia reduction uses

```math
\kappa=0,\quad k=\rho,\quad \beta=r,\quad
q=x,\quad P=rs,\quad Q=r(1-c),\quad V=P-kq.
```

Substitution reproduces all three equations above. **This limiting identification is our derivation from the published equations.** The inspected source presents a four-state system at finite inertia; it does not demonstrate a chaotic attractor in this zero-inertia reduction. A nearby discussion of hyperchaotic Lorenz systems is not such a demonstration. Likewise, its illustrative unrestricted Lorenz plots with $b=8/3$ do not validate the physical sinusoidal-memory slice $b=1$ [1].

Durey (2020), Eqs. (1)-(6), derives the same cosine-wave memory construction without a trap [2]. Its three dynamical states are velocity and two wave variables. Our three states instead comprise position and two memories after eliminating velocity. Equal dimension does not make those two systems equivalent.

Molacek's 2013 thesis, Appendix E, Eqs. (E.4)-(E.6), already studies a harmonically trapped sinusoidal-memory walker and reports irregular regimes [3]. The inspected thesis transcription explicitly keeps position, velocity and two memory integrals. Its survey uses very large memory, not elimination of inertia. The MIT record was verified; its downloadable PDF returned HTTP 429, so the appendix was checked in a public transcription rather than a newly verified page image. This limits the strength of exact typographic checks on that source.

Valani and Dandogbessi (2024), Eqs. (1)-(2), retain acceleration and apply constant bias, rather than harmonic confinement [4]. Valani's infinite-memory work also permits an interpretation as an overdamped particle driven by a three-state internal engine [5]. That wording does not eliminate the engine's velocity state and is not the same operation as imposing the algebraic relation above.

## Consequences for priority and further work

- The exponential complex-memory closure and harmonic-trap parent model are established literature. They must be acknowledged in any account of the motor calculation.
- The equivalence is to a singular parameter boundary. Chaos at finite inertia, in an infinite-memory limit, in a spatially extended wave field, or in an externally imposed internal Lorenz engine does not by itself prove chaos on that boundary.
- This bounded search found no equation-level demonstration of chaos specifically for the finite-memory, zero-inertia, harmonically trapped three-state model. This is a search result, not proof of novelty. It did not exhaust all citations or equivalent circuit models.
- A new plot of this vector field, even if chaotic, would require comparison with both molecular-motor and trapped-memory literature. No new attractor or mechanism is claimed here.

The bounded search used Durey's cited precursor, Valani's explicit harmonic-trap equations, later memory-driven-particle papers, and focused searches combining overdamped, sinusoidal memory, harmonic potential and chaos. Source PDF equation pages for [1] and [2] were rendered and visually inspected locally. Third-party PDFs and screenshots remain outside the published MIT source tree.

## References

1. R. N. Valani, *Lorenz-like systems emerging from an integro-differential trajectory equation of a one-dimensional wave-particle entity*, Chaos **32**, 023129 (2022). [DOI](https://doi.org/10.1063/5.0076162); [inspected author preprint, Sec. V.A, Eq. (12)](https://arxiv.org/pdf/2110.09754).
2. M. Durey, *Bifurcations and chaos in a Lorenz-like pilot-wave system*, Chaos **30**, 103115 (2020). [DOI](https://doi.org/10.1063/5.0020775); [author PDF](https://eprints.gla.ac.uk/276681/1/276681.pdf).
3. J. Molacek, *Bouncing and walking droplets: towards a hydrodynamic pilot-wave theory*, PhD thesis, MIT (2013), Appendix E, pp. 175-196. [Institutional record](https://hdl.handle.net/1721.1/83698); [inspected public transcription](https://studylib.net/doc/10674740/bouncing--and--walking--droplets--jan--molacek-archnie).
4. R. N. Valani and B. S. Dandogbessi, *Asymmetric limit cycles within Lorenz chaos induce anomalous mobility for a memory-driven active particle*, Physical Review E **110**, L052203 (2024). [DOI](https://doi.org/10.1103/PhysRevE.110.L052203); [author text](https://arxiv.org/html/2408.00168v2).
5. R. N. Valani, *Infinite-memory classical wave-particle entities, attractor-driven active particles, and the diffusionless Lorenz equations*, Chaos **34**, 013133 (2024). [DOI](https://doi.org/10.1063/5.0171007); [author text](https://arxiv.org/html/2308.03062v1).
