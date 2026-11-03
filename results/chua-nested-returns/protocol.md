# Prospective test of a second level of periodic returns

This protocol precedes searches for the proposed second-level cycles. The first-level studies found long attracting cycles that spend extra circuits near a short saddle. The present hypothesis is that a continued branch of one such long cycle can itself organize a further sequence of attracting cycles.

## Selected parents and section

Use the existing 35-peak S06 cycle (organized at the first level by S08) as the primary case, and the existing 35-peak S11-n3 cycle as a second case. Fix alpha = 9, m0 = -8/7, m1 = -5/7 and gamma = 0; vary beta. Both long parents have seven inversion-normalized x = 0 transitions during a symmetric half-period.

The parent block return follows the parent's half-itinerary and applies central inversion. It is the seventh iterate of the elementary signed section return. Continue both the originally attracting branch and its fold-connected saddle companion where possible. A continuation through the saddle-node is an allowed change from attracting to saddle motion within the same periodic-orbit family. Which saddle actually organizes the proposed descendants must be determined from the data, not assumed from its branch name.

## Forecasts before target searches

Use only previously known parent roots, their fold and symmetry endpoints, and their companions to define the initial quadratic approximation. In the parent-block coordinate, Q_M(X) = M - X^2 has its fixed-point fold at M = -1/4 and its flip at M = 3/4. Define M = -1/4 + (beta_fold - beta)/W, with W = beta_fold - beta_symmetry for these two parents. The section secant between the existing attracting root (half multiplier 1/2) and its companion fixes the initial state scale: X = -1/4 - s/2.

For odd block counts k = 3, 5 and 7, select the quadratic-map periodic window whose superstable critical orbit has one positive iterate followed by negative iterates before returning to zero (the window nearest M = 2 with this itinerary). Compute its fold, multiplier-1/2 representative and flip before searching the Chua targets. Publish numerical beta seeds, width forecasts and period forecasts. The period forecast uses k times the continued 35-peak saddle-companion period at the predicted representative beta. Coefficients must not be fitted to any new child cycle.

Predicted full-period peak counts are 105, 175 and 245. Record beta-location error in measured window widths, width relative error and period absolute error. Initial engineering gates are one measured width for location, 30% relative width error and 0.5 time units for period; these are not confidence bounds. Search corrections may exceed a gate, but the forecast must remain a miss. If the first quadratic extrapolation fails to seed a target, keep the failure and record any subsequent continuation separately.

## What counts as evidence for nesting

1. The child must be a valid primitive attracting Chua orbit, not a multiple traversal of the parent. Verify phase-space closure, region itinerary, independent full-field peak counts and Floquet attraction. Repeat representative and endpoint solves at higher precision.
2. Group the child's elementary section crossings into seven-crossing blocks. Show distinct block-start points and the repeated internal 35-peak parent motif. Merely concatenating an itinerary or obtaining a peak count divisible by 35 is insufficient.
3. Continue the proposed 35-peak organizing saddle to each child's identical beta. Measure shadowing of its seven section points in order and residence over complete parent blocks. Compare with the continuation of the other 35-peak branch and with the original five-peak saddle.
4. Across increasing k, test whether extra residence is near the long saddle and whether period increments approach twice its full period. This distinguishes a second-level mechanism from one additional circuit near the original short saddle.
5. Use the previously computed ordinary 105-peak S08-sequence orbit as a negative control. Its identical peak count must not make it pass the 35-peak-block ancestry check automatically.

Sample actual parent-block returns across the relevant state and parameter range. The local quadratic curve already measured near multiplier 1/2 does not establish the much wider range needed here. Inspect the return's critical region, the two long-parent branches, and finite-period correction errors.

## Scope

Successful computations would support one further level of organization for specified branches. They would not prove arbitrary recursive depth, a complete hierarchy for the atlas, a unique parent for every orbit, infinitely many descendants of either new saddle, or new knot types. Failure to find targets in this experiment is not a global nonexistence result. Preserve all search failures and forecast misses, and distinguish numerical controls from interval or formal proofs.
