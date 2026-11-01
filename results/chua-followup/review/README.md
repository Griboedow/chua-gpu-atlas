# Revision audit: claims, inequalities and reproduction

This revision makes the three conclusions separately reviewable: local prediction of periodic windows, excitation across a certified basin boundary, and infinite hidden periodic coexistence. The [review guide](../../../docs/chua-review-guide.md) maps each conclusion to the evidence and explains the remaining mathematical dependencies.

## Repeated checks

| Check | Retained output | Interpretation |
|---|---|---|
| Scores recomputed from frozen forecasts and solved endpoints; graph inequalities, face signs, coverage and Lyapunov inequalities recomputed from saved bounds | [record-check.json](record-check.json) | Checks the numerical claims without accepting stored success flags as their justification. |
| Complete hidden-region tangency calculation at 160 decimal working digits | [tangency-160.json](tangency-160.json) | Repeats the saddle inclusion, local graph bounds, 36 returns, quadratic gap and parameter unfolding. |
| Uniform equilibrium exclusion at 100 digits | [equilibria-100.json](equilibria-100.json) | Repeats the spectral enclosure, complete exit transport and invariant capture. |
| Seven local impulse-ray crossings at 100 digits | [thresholds-100.json](thresholds-100.json) | Repeats the all-time flight enclosures, local-domain containment and transverse signs. |
| Independent matrix exponentials, finite differences and direct-field ODE checks | [independent-tangency-checks.json](independent-tangency-checks.json) | Repeats the separate numerical control, including two deliberately shifted parameter boxes. |
| Complete lower-y-ray amplitude coverage at 80 digits | [below-certificate-80.json](below-certificate-80.json) | Repeats the central segment and every captured interval; completed parts are joined only after coverage checks. |

The [execution record](execution.json) records input and output hashes, software versions, replay results and the merger controls. Replays use the published seeds and the same interval implementation. They are repeat executions, not external review or a second proof implementation.

## Findings and corrections

The quantitative conclusions remain unchanged in these checks: 32 raw and 48 corrected midpoint hits among 48 held-out targets; seven local crossing brackets; the y-ray excitation infimum; and the interval hypotheses used in the hidden-coexistence argument. No violated mathematical inequality was found in the checks listed above.

The revision corrects or makes explicit the following points:

- **Different counts.** Sampled family sheets, windows along a transect and simultaneous attractors at one parameter are different quantities. The infinite-coexistence conclusion does not determine the global connected-sheet count.
- **Different parameter sets.** The slope-patch predictors use serialized floating-decimal parameters. The earlier exact-manifold certificate is for the rational classical slice. It is not a certificate for every node of either patch.
- **Different excitation statements.** All seven directions have local crossings; only the y ray has complete exclusion of smaller amplitudes. Its critical amplitude is an infimum, not an attained minimum or a minimum over all directions.
- **Interval display precision.** The short `ball` field can be too coarse to recover a strict inequality. The record checker encloses the separately stored endpoints instead of treating a rounded midpoint as exact. The underlying central-flight bound remains strictly below one.
- **Return-map notation.** The proof now distinguishes physical z, local ξ, the return map P and the Lyapunov matrix R. Legacy JSON names are mapped explicitly to these objects.
- **Reproduction behavior.** Review commands write to separate files. A failed tangency or equilibrium-exclusion CLI run returns a nonzero exit status. The independent negative controls no longer overwrite the positive tangency output. The coverage merger verifies the prefix, ordered band metadata, central calculation and local certificate before joining chunks; the legacy prefix is accepted only with its published digest.
- **References.** The coexistence question is located in both the 2017 and 2023 papers. The 2008 author-manuscript title and theorem numbering are distinguished from the 2011 publication title.

The merger reproduces the original 662-interval certificate exactly from its retained chunks and rejects an omitted chunk. Separate record-level controls reject an omitted capture interval and a curvature enclosure containing zero even when the original success flags remain unchanged. Both shifted tangency boxes are also rejected.

## Verification boundary

The analytic review checked the gap chain rule with moving saddle and invariant-graph derivatives, the distinction between F and F², the interval-uniform hiddenness argument, and the localization needed for the parameter theorem. These arguments remain written mathematics. In particular, no program here proves the localized Newhouse–Robinson theorem, certifies historical priority or supplies a decimal β∞ with an enumerated infinite collection. The existing Lean project does not formalize this result.

The revised README, report, proof and review guide are also inspected as rendered GitHub pages. The browser record identifies the inspected revision and retains screenshot hashes separately from mathematical evidence.
