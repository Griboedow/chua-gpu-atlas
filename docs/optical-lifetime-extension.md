# Extending the optical cohort beyond 20,000 time units

**Completed observation horizon: 100,000 thermal-time units.** Both narrow-width cohorts eventually reach checked periodic cycles, while the 26 broad-width survivors remain irregular with positive finite-time estimates throughout the saved windows. This is stronger persistence evidence for the broad-width candidate, not a proof of a chaotic attractor.

This follow-up continues the **same initial conditions and tangent vectors** used in the [original lifetime experiment](coupled-operators-research.md#long-runs-reveal-periodic-capture-after-chaotic-transients). It addresses unresolved long chaotic transients; it is not a new parameter survey.

![Extended cohort and finite-time expansion](assets/optical_lifetime_extension.png)

## What changed when we waited longer?

| Optical width | Periodic captures at 20,000 | Additional verified captures | Final periodic captures | Still observed at 100,000 | Latest first qualifying checkpoint |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.003 | 19 / 32 | 13 | 32 / 32 | 0 / 32 | 73,500 |
| 0.005 | 29 / 32 | 3 | 32 / 32 | 0 / 32 | 25,500 |
| 0.04 | 6 / 32 | 0 | 6 / 32 | 26 / 32 | 1,000, in the original run |

All **16 additional periodic captures** passed independent CPU endpoint checks. Maximum scaled section mismatch is `1.92e-11`, maximum period mismatch `9.68e-12`, and maximum scaled one-period closure error `8.86e-11`. The original 54 captures retain their previous independent audit. Integration stops at confirmed capture, so captured trajectories are not each advanced for the entire 100,000-unit horizon.

The latest narrow-width example is especially instructive. At `w=0.003`, cohort ID 16 has **583 consecutive 125-unit windows above 0.01**, covering 72,875 time units. Its smallest estimate in that prefix is still `0.11128`, and the prefix mean is `0.26483`. Nevertheless, it subsequently reaches the checked periodic cycle: its first qualifying capture checkpoint is 73,500 and confirmation is at 74,000. This example is explicitly selected **retrospectively as the latest confirmed capture**, not as an unbiased typical lifetime. At `w=0.005`, the latest capture, ID 36, similarly has 198 positive windows covering 24,750 units before settling.

At `w=0.04`, **all 26 survivors have positive estimates in every one of their 800 saved windows**, including the original initial observation. Across all 20,800 estimates the minimum is `0.16932`. Their whole-history mean estimates range from `0.32908` to `0.33368`; their last 500-unit block means range from `0.29984` to `0.36905`. These are finite-time estimates for trajectories selected by survival, not confidence intervals or an asymptotic spectrum. The six periodic captures at this width were already present in the original run; no further capture occurred in the extension.

The figure's upper-right examples are selected by latest capture at each narrow width. Shading shows the four window estimates within each block. The lower-right panel summarizes the same 26 broad-width survivors over the full history; its band is the observed range, not a statistical confidence band. Neither the figure nor the finite cohort establishes convergence for every initial condition or the existence of a separate attracting chaotic invariant set.

## Model and scope

The six-state model couples two damped thermo-mechanical operators through opposing optical feedback. Its parameters are unchanged from reference A:

| Parameter | Value |
| --- | ---: |
| `k` | 74.1567393414624 |
| `b` | 1.784180421282182 |
| `delta_minus` | 0.4139517339970917 |
| `delta_plus` | 0.25683679452165964 |
| Assumed logistic width `w` | 0.003, 0.005 or 0.04 |

The logistic illumination law is an **explicit assumed modification of the ideal switching model**, not a measured beam profile of the source apparatus. Time is measured in thermal relaxation units. No forcing, noise, delay, new state, or changed parameter is introduced in this continuation. The [source audit](light-operators-audit.md) and [model report](coupled-operators-research.md) explain the physical interpretation and its limitations.

At the starting checkpoint, 54 of the original 96 trajectories had already been classified near independently checked periodic orbits. Their endpoints, tangents, IDs and capture records remain untouched. Only the other 42 trajectories are advanced: 13 at `w=0.003`, three at `w=0.005`, and 26 at `w=0.04`. The original random seed, initial conditions and all 40 original history rows are retained exactly.

## Continuation and independent endpoint checks

The GTX 1080 executes the existing FP64 adaptive Dormand–Prince implementation, with `rtol=1e-9`, `atol=1e-12`, and the original halved mechanical/optical step caps. A normalized tangent in velocity-balanced coordinates is carried between 500-unit blocks. Each block saves four separate 125-unit finite-time leading-exponent estimates. Restarting a block does not discard another transient or reset the tangent direction. Local clocks reset because the ODE is autonomous.

The periodic-capture diagnostic is unchanged. At two consecutive checkpoints, the one-period return residual must be below `1e-7` in the velocity-balanced maximum norm, and the magnitudes of all four exponent estimates must be below `0.01`. The stored capture time is the **first qualifying checkpoint**, not the exact time of basin entry or escape from a chaotic transient. Confirmation occurs one block later. The known periods are 1.8368786709, 1.8356888452 and 1.8879301309 for the three respective widths.

Every new classified endpoint is compared independently with the corresponding cycle using CPU DOP853 (`rtol=1e-12`, `atol=1e-14`, maximum step `0.005`). The check integrates four periods and compares upward section positions, return periods and one-period flow closure. Previously checked endpoints are reused only when their parameter bytes, endpoint bytes and canonical periodic-control content agree; the 54 original captures are not recomputed.

The extension reuses the previously validated GPU integrator. Its full 100,000-unit irregular histories are not independently reintegrated on the CPU here; the new independent CPU work specifically verifies captured endpoints. The earlier report contains shorter independent full-spectrum checks for the chaotic numerical regimes.

The calculation checkpoints every complete block. The first segment stops cleanly at total time 39,000 and is archived with its exact driver version. The second segment resumes that saved state and reaches 100,000. **All 160 extension blocks complete, with no integration or probe failure and no incomplete block.** The two GPU segments take 500.285 and 1,101.932 seconds, totaling **1,602.217 seconds (26.70 minutes)**. They perform 587,423,710 accepted and 31,857,877 rejected adaptive steps across the surviving lanes. A partial block would retain its actual elapsed time and unresolved status, rather than being counted as a complete common observation interval.

Across accepted steps in the extension, thermal fractions stay within `0.12652..0.70419`. The displacement ranges are `0.18870..0.80249` and `0.08578..0.63493`. These are monitored numerical ranges, not a rigorous enclosure. The model also has the general boundedness argument given in the original report: thermal inputs remain bounded and both mechanical filters are stable for `k,b>0`.

## Reproducibility and provenance

The original `results/light-operators/lifetimes/` directory remains read-only. New results are in `results/light-operators/lifetimes-extension/`:

- `ensemble.npz` preserves the original cohort and history, adds later checkpoints, and carries endpoint times, tangents, pending capture confirmations, integration/probe statuses, accepted/rejected step counts and state extrema.
- `ensemble.json` records numerical settings, wall/step budgets, completed segments and source hashes.
- `segments/` preserves both completed segments and the exact scripts used for continuation. The final driver additionally supports alternate output directories and the final figure; those interface changes do not alter the archived runs.
- `capture-cpu-audit.json`, `integrity-audit.json` and `window-persistence-audit.json` separate endpoint verification, data integrity and finite-time growth evidence.
- `provenance-manifest.json` identifies the binary snapshots, archived numerical drivers and final analysis drivers with declared raw/canonical hash rules.

Binary NPZ files use exact-byte SHA-256 checksums. Historical raw text hashes are retained; additional canonical text hashes replace CRLF by LF before hashing, so Git line-ending normalization does not masquerade as a scientific change. The continuation audit checks exact preservation of the original and first-segment histories, unchanged old captures, matching resume input, tangent normalization and successful numerical statuses.

From the repository root, a fresh continuation from the provided original 20,000-unit checkpoint can be written to a separate output directory:

```bash
python research/extend_optical_lifetimes.py --run --target 100000 --budget 2400 --output outputs/lifetime-reproduction
python research/extend_optical_lifetimes.py --audit --plot --output outputs/lifetime-reproduction
python research/audit_optical_extension_checkpoint.py --root outputs/lifetime-reproduction --integrity
python research/audit_optical_extension_checkpoint.py --root outputs/lifetime-reproduction --windows
```

The driver refuses to overwrite an existing extension silently. To continue a completed checkpoint, add `--resume` to its run command and keep the same output directory. To regenerate the supplied figure without integration, use `python research/extend_optical_lifetimes.py --plot`. The archived run used separate budgets for two segments; stopping between completed blocks at the planned wall budget is distinct from a numerical failure. Individual chaotic escape checkpoints can change under altered tolerances, hardware or roundoff, so exact cross-platform reproduction of those times is not claimed.

## Interpretation limits

Long positive finite-time exponents demonstrate sustained numerical expansion on the observed intervals. Together with a bounded trajectory, they remain insufficient to establish an asymptotically attracting chaotic invariant set. A very long transient near a nonattracting chaotic set can produce the same finite observations. Conversely, observed captures into a stable cycle do not exclude another basin containing a chaotic attractor at the same parameters.

Counts refer to this fixed cohort of 32 initial conditions per width. They are not basin-volume estimates or an escape-rate fit. Surviving trajectories are right-censored at the observation horizon. No new topological attractor class, experimental realization or scientific priority is claimed.

## References

1. Guo, H., Li, K., Yang, J., Li, D., Liu, F., and Zeng, H. (2025). Light-mediated communication in responsive materials ranging from individual self-oscillators to feedback-driven network. *Nature Communications* **16**, 11315. [doi:10.1038/s41467-025-66395-3](https://doi.org/10.1038/s41467-025-66395-3). Published physical system and supplement; the logistic profile used here is our explicit modeling assumption.
2. [Original numerical model and lifetime report](coupled-operators-research.md). Initial cohort, periodic controls, independent spectra and the original 20,000-unit observation.
3. [Order structure of the smooth coupled operators](coupled-operators-order-audit.md). Conditional order-based restrictions in the overdamped parameter regime; the present reference A is strongly underdamped.
