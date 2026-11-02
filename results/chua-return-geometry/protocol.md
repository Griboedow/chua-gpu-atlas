# Return geometry and transient amplification: prospective tests

This protocol is recorded before computing the new transverse amplification and companion-saddle data. It uses already published periodic roots and window widths. It is a new test of their dynamical relationship, not a new blind test of those roots or widths. Work is local; original evidence files remain unchanged.

## Questions

1. How much can a perturbation transverse to a stable cycle grow during one period, even when the full-period Floquet spectral radius is fixed?
2. Does this amplification grow with the extra residence time near the organizing saddle in the eight-member period-adding sequence?
3. Does the same quadratic return geometry predict the separation between the attracting cycle and its saddle companion at identical parameters?
4. Which conclusions survive comparison with the other reference families?

## Definitions fixed before measurement

Use dimensionless coordinates (x,y,z) with their Euclidean norm. For the variational propagator Phi(t+s,t), let B(t) be an orthonormal basis perpendicular to the vector field along the orbit. The transverse propagator is B(t+s)^T Phi(t+s,t) B(t). It removes a displacement along the orbit, which represents a phase shift. The field is continuous at diode breakpoints, so the variational solution has no saltation jump there.

For each sampled initial phase, G(t) is the largest singular value encountered over sampled durations 0 <= s <= T. G_max is the maximum over sampled phases. These are finite-time linear quantities, in a stated norm, not basin volumes, finite escape thresholds or energy gains. Sampled maxima are lower estimates of the continuous supremum. All 34 reference cycles are measured; they use their own reference beta values.

## Predictions and possible failures

The established local/global return model gives T_(n+1)-T_n -> 2 T_p and W_(n+1)/W_n -> Lambda^(-2), with Lambda approximately 2.71343178. Here W is the symmetric interval from fold to symmetry breaking. Those are existing predictions, not discoveries of this study.

The new amplification hypothesis is G_(n+1)/G_n -> Lambda, so W_n G_n^2 tends to a nonzero constant. This requires the dominant sampled amplification to follow the added local saddle passage; different global excursions or the maximizing phases can invalidate a simple finite-n law. Retain the result even if the last two gain ratios differ from Lambda by more than 20%, or the last four W G^2 values have max/min greater than 1.5. These gates are exploratory engineering tolerances, not confidence intervals or rigorous asymptotic bounds.

For Q_M(X)=M-X^2 at the normalized attracting multiplier +1/2, M=-3/16. The two fixed points are -1/4 and -3/4; the companion multiplier is +3/2. A Chua companion on the same half-itinerary is predicted to approach this multiplier, giving full-period multiplier 9/4. The section separation is expected to have scale Lambda^(-n) in a consistently aligned section. The test must match phase: arbitrary first crossings of a long cycle need not supply a consistent coordinate for comparing families. A saddle periodic point is not by itself a measurement of distance to the basin boundary.

Potential geometric evidence includes increasing numbers of section returns close to the short saddle, residence times approaching its half-period, and convergence after rescaling. Conclusions about all other families require separate measurements. Neither the topology of all global layer connections nor completeness of the atlas is inferred.

## Controls

- Reconstruct affine arcs from the published 70-digit roots, and retain source hashes.
- Verify full-period transverse multipliers and tangent transport.
- Repeat selected gains with a finer phase mesh and higher precision where needed.
- Independently integrate the full field and variation on selected maximizing intervals, projecting at intermediate points to avoid neutral-direction contamination.
- Check finite perturbations against linear predictions at multiple small amplitudes on selected intervals.
- Solve companions at fixed beta with distinct-root and itinerary checks; repeat the longest tested roots at higher precision.
- Keep unsuccessful predictions and solver failures. Separate algebraic derivations, numerical evidence and existing interval theorems in the report.

## Relation to earlier work

Homoclinic rescaling and periodic-sink cascades are established theory. Phase-dependent transient amplification despite Floquet stability is also established; see Kiorpelidis et al., *Transient amplification in stable Floquet media*, Physical Review B 110, 134315 (2024), https://doi.org/10.1103/PhysRevB.110.134315. This experiment concerns the specific Chua reference cycles, their shrinking parameter windows and the connection to their organizing saddle. No priority claim is made.
