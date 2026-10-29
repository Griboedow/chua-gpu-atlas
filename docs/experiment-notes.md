# Notes from returning to the Chua system

These are interpretations of the supplied calculations. They are not claims of new general theorems about Chua's circuit.

## The surface is a question, not an answer

The first three-dimensional picture looked like a curved solid. A grid renderer can easily make that impression stronger than the evidence warrants. In fact it marks sampled parameter tuples where at least one tested initial condition produced bounded, moving, positive-Lyapunov behavior during a finite observation.

The original surfaces hit their scan boundaries. Their apparent edges were partly edges of the experiment. Going from 51³ to 81³ improved the local shape, but could not answer where the whole family ended.

The wide survey consequently uses points rather than a filled surface. Across 32,768 logarithmically distributed Sobol samples, 2,249 candidates survived a second integration with smaller steps and longer observation. They reach almost to α=256 and m₀=−8. We need a larger search or an analytic argument to decide global boundedness. We also have not established that every distant cluster connects continuously to the original region.

## Bounded motion and an unbounded parameter set are different questions

An attractor can fit inside a finite region of state space (x,y,z) while a family of such attractors exists for arbitrarily large parameter values. Conversely, a curved parameter cloud can eventually terminate. The boundedness of one trajectory does not decide the extent of the parameter family.

Even for one parameter tuple, staying below the escape cutoff for a finite interval only demonstrates finite-time boundedness. It does not prove that the trajectory stays there forever.

## A narrow window can hide between neighboring grid points

The first classic grid used Δβ=0.16. A finer transect found a regular window around β=14.31–14.33. At β=14.32 the retained maxima form a repeating five-peak pattern; adaptive DOP853 independently reproduces it.

The local zoom uses Δα=0.003, Δβ=0.004 and Δm₀=0.0003. What looks like a small gap in one slice becomes an inclined sheet in three dimensions. This is a practical reason to scan three parameters: a window can move between slices and appear or disappear in a two-dimensional view.

The pink surface uses near-zero/nonpositive exponents at all three starts. Only selected examples were checked for peak recurrence. It should not be read as a rigorously determined periodic-orbit bifurcation surface.

## One parameter tuple does not necessarily select one attractor

At (α,β,m₀,m₁)=(15.6,32,−1.143,−0.714), starts at (±0.7,0,0) approach opposite single-scroll chaotic regimes. A basin slice at z₀=0 assigns 11,400 sampled starts to each side; the exact origin remains an equilibrium.

The equal counts are a consequence of central symmetry and the chosen symmetric sampling. They are not a universal 50/50 probability. Changing the distribution of initial conditions changes such fractions.

This also explains why “chaos exists at this parameter” and “every tested start is chaotic” are different quantities. Both union and intersection masks are retained in the datasets.

## Why the cloud curves

For the studied outer slope m₁=−5/7 and positive α,β, the outer-branch Jacobian has characteristic polynomial

```text
λ³ + A λ² + B λ + C
A = 1 + α(1+m₁)
B = β + αm₁
C = αβ(1+m₁).
```

The cubic stability condition A B > C becomes

```text
β > −αm₁[1 + α(1+m₁)] = 5α/7 + 10α²/49.
```

The other positivity conditions also hold above this threshold for the stated slope and positive α. The threshold depends on m₁, not on m₀, because the outer branch has a constant derivative. In the explored m₀<−1 family, the outer equilibria are indeed on this branch.

All repeated wide-survey candidates lie below the surface. This gives a useful reference scale for the bend. Local equilibrium stability does not by itself determine all global attractors or exclude coexistence elsewhere.

## A positive Lyapunov exponent needs context

The largest exponent estimates the average exponential growth rate of an infinitesimal perturbation. Positive means nearby trajectories initially separate; negative means contraction along that perturbation after alignment. The estimate here is finite-time, obtained by jointly evolving and repeatedly normalizing a tangent vector.

A positive exponent at an unstable equilibrium is not evidence that a stationary trajectory is chaotic. A positive finite-time value can also come from a transient. Conversely, the largest exponent of a nonstationary regular autonomous flow is often near zero rather than negative. This is why the screen also checks motion, escape, several observation windows and selected peak patterns.

## Faster calculations expose numerical problems sooner

At α≈255.69, one fixed-step GPU estimate is about 0.93, while the independent event-resolved estimate is about 0.18. Calling this agreement would be misleading. Both remain positive in their four windows, but quantitative convergence is unresolved.

The vector field is continuous at x=±1, while its Jacobian changes abruptly. A fixed RK step may cross a kink without locating it. The independent solver stops at crossings and restarts on the appropriate affine branch. Observation lengths also differ. Both effects require further study before using distant-point exponents as precise measurements.

A useful next experiment would compare matched post-transient windows with systematic step refinement against event-resolved integration, then repeat the wide search where the classification itself changes. More points alone will not resolve an integration bias.

## Questions still open in this project

- Does the explored chaotic family eventually terminate as α increases or m₀ decreases?
- Which distant points belong to the same connected parameter component?
- How do the thin regular sheets bifurcate and intersect away from the local zoom?
- Which apparent gaps are regular windows, missed basins, long transients or numerical artifacts?
- How does changing the fixed outer slope m₁ reshape the three-dimensional sections?

For bibliographic context, see the numbered references in the [README](../README.md#references). Full settings and raw evidence accompany the [results release](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.0.0).
