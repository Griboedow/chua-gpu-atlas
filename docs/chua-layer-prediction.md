# Predicting Chua layers from an organizing periodic orbit

Research in progress. The theory and forecast protocol below are written before searching for the 75- and 85-peak members. Completed checks will be reported separately from these forecasts.

## The question

Can a measured family of thin layers be explained by a mechanism that predicts its motions, parameter locations and widths? The target is the central sequence with 15, 25, 35, 45, 55 and 65 maxima, at alpha = 9, m0 = -8/7 and m1 = -5/7. This is one sequence inside the atlas, not a proposed classification of all periodic windows.

## A layer needs more than a left/right word

Let L, C and R denote the affine regions x < -1, -1 < x < 1 and x > 1. The successive region visits form a useful itinerary, but this partition is not a demonstrated generating partition. In particular, S08 and S11 have the same half-cycle itinerary C R C R C, although their periods are about 10.89 and 14.50. Different return times solve the same region-word closing problem. A word alone neither establishes orbit identity nor guarantees that an attracting solution exists.

The appropriate description retains:

1. The primitive full orbit and its inversion symmetry, with phase shifts and inversion partners treated explicitly.
2. The region word and residence times, together with a starting section state. The affine equations reconstruct the orbit from these data.
3. The period, ordered section returns and extrema, to distinguish repeated-looking passes.
4. The transverse multipliers and the branch's known continuation ancestry.

Across the adding sequence, the half-cycle outer-region word is `(RRLL)^n RR`, for n = 1,...,6. Inversion supplies the other half. Its full outer-region word is `(RRLL)^(2n+1)`. Thus the coarse word has a shorter symbolic period than the actual orbit: unequal residence times and return coordinates carry the longer recurrence. The forecast for n = 7,8 is 75,85 maxima and 15,17 upward section returns. Counts must be checked independently; a formal word is only a candidate itinerary.

## A local-to-global mechanism to test

The short five-peak symmetric cycle continues as a saddle cycle in this parameter range. Denote its full period by Tp, its unstable transverse multiplier by Lambda > 1 and its stable multiplier by lambda, with abs(lambda) < 1. A long candidate follows the saddle for additional returns before completing a global excursion. Adding one full saddle return to each half of an inversion-symmetric orbit predicts

```math
T_{n+1}-T_n\longrightarrow 2T_p.
```

The proposed organizing event is a nondegenerate quadratic tangency of the relevant invariant manifolds. This event has **not** been proved for this branch. Under that hypothesis, a local expanding passage multiplies displacement by Lambda^n. A quadratic global return therefore has curvature of order Lambda^(2n). After centering and rescaling, its effective one-dimensional fixed-point problem has the form

```math
X_{k+1}=M-X_k^2+\text{small corrections}.
```

For this quadratic map, a fixed point satisfies M = X + X^2 and its multiplier is -2X. The fold occurs at X = -1/2, M = -1/4; the flip at X = 1/2, M = 3/4. The interval in M has width one. For the inverted half return used here, the flip corresponds to symmetry breaking of the flow cycle, rather than a doubling of its full period.

If the unfolding derivative with respect to beta is nonzero and the strongly contracting direction remains negligible in the rescaling, this gives two different scales:

```math
\beta_n-\beta_\infty\sim A\Lambda^{-n},\qquad W_n\sim B\Lambda^{-2n}.
```

Here Wn is the symmetric-core width from fold to symmetry breaking. The signs and constants of the position law depend on the global return. Neither the full window including descendants nor all families share a universal numerical constant. In particular, width should shrink approximately as the **square** of the factor governing the distance to accumulation.

For an initial finite-order correction to locations, we freeze the ansatz

```math
\beta_n=\beta_\infty+A r^n+(B+Cn)r^{2n},\qquad r=\Lambda^{-1}.
```

The n-dependent term allows the first correction from parameter-dependent local expansion. This is an asymptotic model, not an exact identity. Its constants are calibrated using existing layers, while Lambda comes from the short saddle cycle. Forecasts retain a simpler leading-order alternative so the effect of the correction is visible.

## What this would say about the number of layers

An infinite sequence is a theorem-level conclusion only if the required global connection, quadratic nondegeneracy, dissipativity, admissible itineraries and transverse parameter unfolding are established. A finite fit, however accurate, does not verify those hypotheses for all n. Counting unstable periodic orbits also does not count attracting windows.

Conditional on Wn approximately B Lambda^(-2n), the number of cores wider than a resolution epsilon grows only logarithmically with 1/epsilon:

```math
n_{\max}\approx\frac{\log(B/\epsilon)}{2\log\Lambda}.
```

An additive index convention and finite-order corrections remain. This is a count of resolvable members of this sequence, not the total atlas count and not the number of attractors coexisting at a single parameter tuple. A finite experimental resolution can therefore see only a modest number even if the mathematical hierarchy is infinite.

## Forecast and falsification protocol

- Use only the previously published 15–65-peak sequence for calibration. Existing values are not called blind tests.
- Freeze a machine-readable prediction file, its input hashes, target itineraries and tolerances in Git before solving for 75 or 85 peaks.
- Predict peak count, number of section returns, symmetry, fold location, symmetric-core width and period. Publish both the simple and corrected location forecasts.
- Construct the candidate by inserting the predicted saddle return; this is a targeted search, not an exhaustive scan. Compare the correction distance with the window's width, not only with beta itself.
- Use higher precision if floating-point error approaches the predicted width. Preserve unsuccessful trials. Independently integrate the full piecewise vector field and check region validity and attraction.
- Keep failed forecasts, including a location forecast that misses a very thin stable core. Re-fitting after discovery does not repair a frozen forecast.

## Relation to prior work

Symbolic descriptions and one-dimensional approximations to Chua return maps have a substantial history. Chua and Tichonicky's [1989 Berkeley report](https://www2.eecs.berkeley.edu/Pubs/TechRpts/1989/ERL-89-108.pdf) uses kneading information to organize periodic windows and relate map cycles to flow orbits. Its coding is richer than merely labeling our three affine regions. The present coarse words must not be mistaken for a proved symbolic conjugacy.

The local/global return construction near a homoclinic tangency and the creation of periodic sinks are also established theory; see Gonchenko, Shilnikov and Turaev, [Nonlinearity 21 (2008), 923–972](https://doi.org/10.1088/0951-7715/21/5/003). Applying that theory requires checking its hypotheses for the particular return map. The continuous piecewise-linear vector field is not globally smooth, but fixed transverse crossing itineraries define smooth local flight maps; grazing would require a separate treatment.

The intended contribution is a documented, falsifiable prediction for this exact numerical hierarchy and a clearer description of what the atlas labels mean. General symbolic dynamics, period addition and exponential window shrinkage are not new discoveries.
