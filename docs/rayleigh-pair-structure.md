# Structural controls for the quadratic-drag Rayleigh pair

Audit date: 30 October 2026 UTC. This note concerns exactly the four-state canonical reduction of Li et al. [1] recovered in [the source screening](late-physical-candidates-2.md):

```math
x_i''=-x_i-\kappa(x_1+x_2)+\mu(1-|x_i'|)x_i',
\qquad \mu>0,\quad\kappa\geq0.
```

Time is the source-derived canonical variable `sigma`; it is not physical seconds. The drag is quadratic with an absolute value. Replacing it by cubic Rayleigh drag changes the model.

**Outcome:** an absorbing-set argument proves global mathematical boundedness. The two modal planes admit exact transverse Floquet controls, and both source-derived modal cycles are transversely stable. We did not establish a global theorem excluding chaos in the active four-dimensional system. No off-modal attractor trajectory or parameter survey was computed in this block.

Reproduce the algebra and seven small modal controls with `python research/audit_rayleigh_pair.py`. The [results record](../results/candidate-audits/rayleigh-pair-structure.json) contains both DOP853 and Radau calculations. These are numerical controls, not rigorously validated orbit-existence proofs.

## Reflection, symmetry, and modal equations

Reflect the second coordinate and velocity, `y1=x1, y2=-x2`. The elastic terms become the ordinary difference coupling `-kappa*(yi-yj)`. This is an exact linear conjugacy, so the source's sum coupling does not create a distinct mechanism through its sign alone.

The original system is invariant under exchanging the two oscillators and under simultaneous sign reversal. Set

```math
s=(x_1+x_2)/2,\qquad a=(x_1-x_2)/2,
\qquad g(v)=(1-|v|)v.
```

Then

```math
s''=-(1+2\kappa)s+
\frac{\mu}{2}[g(s'+a')+g(s'-a')],
```

```math
a''=-a+\frac{\mu}{2}[g(s'+a')-g(s'-a')].
```

The symmetric plane `a=a'=0` carries a scalar oscillator of linear frequency `sqrt(1+2*kappa)`. The antisymmetric plane `s=s'=0` carries one of frequency 1. Either is planar and cannot itself contain a strange attractor. At zero coupling the two identical oscillators decouple, providing another non-chaotic control. Initializing exactly on a modal plane cannot test whether a separate off-modal attractor exists.

The scalar modal velocity also satisfies the Liénard equation

```math
v''+\mu(2|v|-1)v'+\omega_b^2v=0.
```

The first derivative of `g` is continuous, including at zero: `g'(v)=1-2|v|`. Ordinary variational equations apply; there is no impact, reset, or saltation matrix. This planar observation does not reduce the complete pair to two dimensions.

## Exact transverse damping equation

Let a modal orbit have velocity `v_b(t)`, base frequency `omega_b`, and period `T`. Its transverse displacement obeys

```math
\eta''+c(t)\eta'+\omega_t^2\eta=0,
\qquad c(t)=\mu(2|v_b(t)|-1).
```

For the symmetric orbit, `omega_t=1`; for the antisymmetric orbit, `omega_t=sqrt(1+2*kappa)`. The longitudinal variational equation has the identical damping function and stiffness `omega_b^2`. Therefore its flow multiplier is one, and the other longitudinal multiplier is `exp(-integral(c,t=0..T))`. The transverse determinant has exactly that same exponential value. This is a useful independent implementation check.

Positive mean damping alone does not establish transverse stability. The substitution

```math
\eta(t)=\exp\left[-\frac12\int_0^t c(u)\,du\right]z(t)
```

gives the Hill equation

```math
z''+\left[\omega_t^2-\frac{c'}2-\frac{c^2}4\right]z=0.
```

The transverse growth rate equals the corresponding Hill growth rate minus half the mean damping. The derivative `c'` is understood piecewise at velocity zero; `c` itself is continuous, so this transformation introduces no impulse. A positive transverse exponent would establish instability of that modal cycle, not a chaotic attractor elsewhere.

## Weak-damping interpretation

At leading order for small positive `mu`, the modal velocity is sinusoidal with amplitude `V=3*pi/8`. This follows by balancing average injection and drag:

```math
\langle v^2\rangle=\langle|v|^3\rangle,
\qquad V^2/2=4V^3/(3\pi).
```

Thus mean transverse damping is approximately `mu/2`, and away from resonance both transverse real exponents are approximately `-mu/4`. More explicitly, choosing the leading velocity as `V*cos(omega_b*t)` gives

```math
c(t)\simeq\mu\left[\frac12+
3\sum_{n=1}^{\infty}\frac{(-1)^{n-1}}{4n^2-1}
\cos(2n\omega_b t)\right].
```

At a primary first-order parametric resonance `omega_t≈n*omega_b`, the leading maximal real exponent is approximately `mu*(|3/(4n^2-1)|-1)/4` at exact tuning. It is negative for `n>=2`; the `n=1` result is marginal at this order. In particular, generic weak damping does not make every rational internal frequency ratio an instability. Higher-order effects, small-coupling degeneracy near `kappa=0`, and finite or large `mu` require their own analysis. This expansion is not a uniform global theorem.

## A global absorbing set

The elementary energy derivative from the source screening proves finite-time non-blow-up. A saturating cross term strengthens it to ultimate boundedness, including large displacement with small instantaneous velocity.

Let `M=I+kappa*11^T`, `y=x'`, `z=-Mx`, and `h_i(y)=mu*(|yi|-1)*yi`. Since `M` is positive definite,

```math
y'=z-h(y),\qquad z'=-My.
```

Put `r=||y||`, `phi(y)=y/sqrt(1+r^2)`, and choose `0<epsilon<=1/(2*mu)`. Define

```math
W(y,z)=\frac12\|y\|^2+\frac12z^TM^{-1}z
-\varepsilon\phi(y)^Tz.
```

This function is proper: its quadratic terms dominate the cross term, whose magnitude is at most `epsilon*||z||`. A constant can be added to make it positive. Direct differentiation gives

```math
W'=-y^Th-\varepsilon z^TD\phi\,z
+\varepsilon h^TD\phi\,z+\varepsilon\phi^TMy.
```

The symmetric matrix `Dphi` has smallest eigenvalue `(1+r^2)^(-3/2)` and largest eigenvalue `(1+r^2)^(-1/2)`. Use Young's inequality in its quadratic form, `sum|yi|^3>=r^3/sqrt(2)`, `||h||<=mu*r*(1+r)`, and `lambda_max(M)=1+2*kappa`. They yield

```math
W'\leq
\left(\mu+\frac{\varepsilon\mu^2}{\sqrt2}\right)r^2
-\frac{\mu}{2\sqrt2}r^3
+\varepsilon(1+2\kappa)r
-\frac{\varepsilon\|z\|^2}{2(1+r^2)^{3/2}}.
```

For sufficiently large `r` the cubic negative term dominates. On the remaining bounded interval of `r`, sufficiently large `||z||` makes the last term dominate. Hence `W'<-1` outside a compact set. A sufficiently large sublevel set of the proper function `W` is absorbing and forward invariant. Together with global existence, this proves an ultimate bound for every finite initial state and fixed `mu>0,kappa>=0`.

This argument establishes a bounded mathematical dynamical system. It neither makes the original mechanical energy monotone everywhere nor excludes recurrent or chaotic dynamics inside the absorbing set. The bound is also not a guarantee of physical fiber lengths, tension, or validity of the thermal approximation.

## Small CPU controls

Each modal cycle was located by a one-dimensional section return solve. Its longitudinal and transverse variational equations were then integrated independently by DOP853 and Radau. No generic four-dimensional orbit was searched.

| `mu` | `kappa` | Mode | Period in `sigma` | Largest transverse multiplier modulus |
|---:|---:|---|---:|---:|
| 0.00880947943 | 0.557378006 | Symmetric | 4.32066008 | 0.99052943 |
| 0.00880947943 | 0.557378006 | Antisymmetric | 6.28319652 | 0.98625736 |
| 0.2 | 0 | Uncoupled control | 6.28896156 | 1, to numerical error |
| 1 | 1.5 | Symmetric | 3.15963204 | 0.45234836 |
| 1 | 1.5 | Antisymmetric | 6.42708568 | 0.28020196 |
| 3 | 1.5 | Symmetric | 3.30238628 | 0.07641236 |
| 3 | 1.5 | Antisymmetric | 7.49103923 | 0.00980655 |

The first two controls use the source-derived canonical parameters, rounded as recorded in the source screening. Their transverse exponents are approximately `-0.00220237`, in agreement with `-mu/4`. Their computed Floquet spectra indicate local attraction in the full pair. That does not establish that these basins exhaust all initial conditions. The four larger-`mu` examples are mathematical controls, not validated physical LCE operating points. No transverse instability was found in this small set.

The 256 random-state algebra checks verify transcription of the modified-energy derivative and its upper bound; the analytical proof does not rely on sampling. Numerical residuals, complete multipliers, and integrator comparisons are saved in the results record.

## Stronger damping: two transverse instabilities

A subsequent set of twelve hand-chosen modal controls used both planes at `mu=10, kappa=0.01,1.5,4` and `mu=30, kappa=0.1,1.5,4`. Reproduce them with `python research/rayleigh_modal_followup.py`; the separate [follow-up record](../results/rayleigh-pair/modal-followup.json) preserves the original seven controls unchanged. All twelve completed in about 75 seconds, using DOP853 and Radau independently.

Two antisymmetric cycles were transversely unstable:

| `mu` | `kappa` | Period in `sigma` | Largest transverse multiplier | Positive transverse exponent |
|---:|---:|---:|---:|---:|
| 10 | 0.01 | 13.85173814 | 1.01904760 | 0.00136217 |
| 30 | 0.1 | 34.34908472 | 1.15097732 | 0.00409360 |

The symmetric cycles at those same points remained stable, with transverse multiplier moduli approximately `0.98013157` and `0.78658473`. The other eight tested modal cycles were stable. Both unstable multipliers are positive real. This establishes numerical growth of a transverse perturbation to those modal cycles, not a period-doubling finding or a chaotic attractor. The corresponding growth times are long: one e-fold takes approximately 734 and 244 canonical time units. A trajectory seeded extremely close to an unstable modal plane can consequently look persistent for a long time before leaving it.

The dominant multipliers agree between the two solvers. At large damping, the tiny contracting multiplier cannot be recovered reliably from a direct double-precision monodromy determinant: for the second unstable case, its true logarithm is hundreds of units below zero. Unresolved direct exponents are marked missing. When the leading multiplier is resolved, the record separately infers the other logarithmic multiplier from the integrated trace; that is explicitly labeled an inference, not a direct eigenvalue measurement. The largest full augmented-state difference between solvers among these twelve controls was approximately `2.3e-9`.

These stronger-damping examples are mathematical bifurcation leads. They do not establish physical validity of the reduced LCE constitutive approximation at those parameters, nor suggest that all trajectories avoid the competing stable symmetric cycle.

## Scope of the exclusion

The exact modal planes and the uncoupled limit exclude chaos in those restrictions. The passive parameter regime discussed in the source screening also converges to equilibrium. We have not found an applicable theorem that excludes autonomous chaos throughout the active coupled four-state model. Opposite signs in the displacement/velocity reciprocal derivatives obstruct the simplest orthant cooperative/competitive sign-flip argument, but this is not a classification of all possible coordinate transformations or monotone cones.

A subsequent [full-system survey](rayleigh-pair-research.md) uses starts away from both modal planes: 3,081 scouts and 468 completed tighter follow-ups retain no chaos candidate. Independent calculations also close an off-modal saddle cycle at weak coupling. These unstable periodic controls are useful bifurcation leads, but are not chaotic attractors. No literature absence or priority conclusion follows from this structural audit or the finite survey.

## Reference

1. K. Li, H. Wu, B. Zhang, Y. Dai, and Y. Yu, “Heat-Driven Synchronization in Coupled Liquid Crystal Elastomer Spring Self-Oscillators,” *Polymers* **15**, 3349 (2023). [Primary article](https://doi.org/10.3390/polym15163349); [full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10458843/fullTextXML). The exact Equation 14 reduction, parameter conversion, physical limits, and bounded prior-work audit are documented in [the source-screening note](late-physical-candidates-2.md).
