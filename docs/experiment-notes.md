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

## Going four times farther

The expanded survey increases the α limit to 1024 and the depth −1−m₀ limit to 28. Sampling β relative to its quadratic stability scale raises the upper β limit to about 279,144. Among 16,384 new tuples, 1,832 retain a chaotic start in two adaptive GPU runs. The cloud reaches α≈1023.64, β≈173,658 and m₀≈−28.91; these are separate extrema.

The tail still reaches the imposed limits in α and slope depth. This strengthens the evidence that the old search box was a small fragment, but does not establish an unbounded family. The sample has 178 candidates at α>768 and 95 at m₀<−22. All repeated candidates remain below the outer-equilibrium stability surface.

The overlay uses a different integrator, escape cutoff and sampling measure. It answers a geometric question about extent; comparing raw node counts would not measure a change of chaotic volume.

## Faster calculations expose numerical problems sooner

The initial v1.0 investigation reported about 0.93 from fixed-step GPU integration and about 0.18 from an affine-branch CPU event routine at α≈255.69. The expanded study uncovered a failure in that reference strategy: if a crossing is missed, the solver can keep extending the wrong affine field. At the maximum-α example, a diagnostic rerun found an accepted state at x≈−1.3×10⁸ while the solver still assumed the inner branch |x|≤1. The corrected full-field reference stays bounded below 65 in the same 100-time-unit check. The exact parameters, starting state and historical commit are retained in `results/expanded/reference_failure_diagnostic.json`.

The old CPU records remain public for provenance, but they are not definitive reference measurements. The corrected DOP853 implementation evaluates the actual piecewise field at every stage. Events record crossings without selecting the active branch, and tangent normalization uses a scaled norm. A missed diagnostic crossing therefore cannot force the wrong dynamics.

The new adaptive GPU backend controls state and tangent error together. Short states agree with independent DOP853 at α=9, 255.69 and 1024; a linear exact-solution control and stationary-equilibrium rejection also pass. Selected distant examples receive additional full-field CPU checks at matched intervals. These checks support the selected finite-time behavior, not precise exponent convergence over the entire cloud.

A useful next experiment would systematically refine tolerances and observation lengths at fixed parameter tuples, then reassess any classification changes. More sampling alone cannot resolve integration bias or long transients.

## Questions still open in this project

- Can we find a robust chaotic attractor with geometry substantially different from the familiar scroll, butterfly, funnel and toroidal examples? Treat visual distinctness, a different dynamical mechanism, and literature-supported novelty as separate questions. A rotated, stretched, nonlinearly remapped or unusually projected known attractor does not by itself meet the research objective. Compare multiple projections and Poincare sections, Lyapunov spectra, symmetries, and the equations that generate the dynamics. For higher-dimensional candidates, label three-dimensional views as projections. Preserve reproducible trajectories and prepare a rotation GIF for any substantive candidate.
- Does a chaotic attractor family persist as α grows without bound along β = κ α², with κ, m₀ and m₁ held fixed? After suitable time and state rescaling, does a limiting dynamical model or a limiting structure of the chaotic parameter region emerge? The apparent bend in a finite cloud motivates this question but does not establish such a limit; its novelty also needs a literature check.
- Does the explored chaotic family eventually terminate as α increases or m₀ decreases?
- Which distant points belong to the same connected parameter component?
- How do the thin regular sheets bifurcate and intersect away from the local zoom?
- Which apparent gaps are regular windows, missed basins, long transients or numerical artifacts?
- How does changing the fixed outer slope m₁ reshape the three-dimensional sections?

For bibliographic context, see the numbered references in the [README](../README.md#references). Full settings and raw evidence accompany the [results release](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.1.0).
