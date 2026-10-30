# Recorded evidence

These JSON files are copies of the computation settings and validation summaries supplied in the complete results release. They allow checking the headline claims without downloading the large arrays.

- `original/`: two 51³ scans, reference trajectories and eight numerical validation checks.
- `detailed/`: two 81³ grids, 41³ local zoom, beta transects, example regimes, connectivity summaries and independent validation.
- `expanded/`: adaptive solver controls, expanded-scout settings, coordinate extrema and corrected full-field DOP853 checks.
- `extent/`: wide-scout settings, repeated-candidate interpretation, per-parameter seed mapping and historical v1.0 affine-branch event checks (superseded; see correction below).

The large basin metadata, complete NPZ arrays, all offline interactive reports and browser-check records are in [release v1.1.0](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.1.0). Extract the archive into the repository root. Its internal `RESULTS_MANIFEST.json` provides per-file SHA-256 hashes; `SHA256SUMS.txt` verifies the ZIP itself.

These are finite-time numerical experiments. The old affine-branch CPU reference in `extent/event_resolved_checks.json` has been superseded after a missed-crossing failure at a distant point. Those files are historical evidence, not definitive reference values; use the new full-field checks under `expanded/`. Quantitative convergence over the entire cloud is not established. See the main README and methods for interpretation.

## Exploratory searches beyond Chua

These later datasets are small enough to include directly, along with their settings and original analysis code:

- `candidate-research/`: quantum-neuron surveys, independent controls, negative candidate follow-ups, and literature records.
- `synchronematic/`: phase/orientation cluster scouts, long reruns, full CPU spectra, robustness ensembles and a source-coefficient audit. The main positive evidence uses the literal article coefficients; the article and author code disagree on one sign. See [the notebook](../docs/synchronematic-research.md).
- `photofeedback/`: an explicitly phenomenological response model, a 3,072-trajectory survey, 12 unsuccessful follow-ups, linear stability roots and independent trajectories. This is not a fitted reconstruction of the experimental device. See [the notebook](../docs/photofeedback-research.md).
- `lce-slider/`: two 8,192-trajectory event-map surveys, long follow-ups, independent analytic/DOP853 spectra and 2,560 local sensitivity trajectories. The positive cases extend beyond the source paper's numerical stiffness/damping ranges. See [the notebook](../docs/lce-slider-research.md).
- `lce-slider/smooth*`: five logistic switching widths, a tighter/longer chaos check, periodic-orbit shooting/Floquet controls and a local equilibrium boundary. This explicitly changes the ideal switch; see [switching sensitivity](../docs/lce-slider-switching.md).
- `thermomechanical-sheet/`: an independent coordinate/Jacobian/equilibrium audit and forward-citation leads for Shen's 2025 model. No chaos survey yet; see [the audit](../docs/thermomechanical-sheet-audit.md).

All angles, time scalings, selection thresholds and source variants are described in their notebooks. Finite-time candidate membership is not a proof of an asymptotic attractor or publication novelty.

- `nifte/`: exact five-/nine-state thermofluidic equations, recovered stable/unstable source cycles, 14,336 scouts, tighter follow-ups, a full-spectrum weak-positive cycle control, source/physical audits and independent mode-interaction algebra. See [the notebook](../docs/nifte-research.md).
- `rayleigh-pair/`: exact four-state quadratic-drag spring pair; 3,081 scouts, 468 completed tighter checks with auditable resumptions, three very long weak-signal checks, full CPU spectrum, stable and saddle periodic cycles. See [the notebook](../docs/rayleigh-pair-research.md). No chaos candidate retained.
- `hhs/`: literal fixed-reference five-state mechanochemical single cell; exact uncoupled spectrum control, 3,267 scouts, 187 selected follow-ups, independent cycle checks, boundedness/equilibrium algebra and a source-figure reproduction discrepancy. See [the notebook](../docs/hhs-research.md). No retained chaos candidate among 186 completed long GPU intervals; one interval remains incomplete.
- `microfluidic/`: literal four-state fixed-geometry flow/vortex model; analytic-Jacobian and conservation checks, three source-pressure CPU controls, an independent Radau reconstruction, equilibrium and conditional energy calculations. [The periodic source control was not recovered](../docs/microfluidic-research.md); no mass GPU search was run while the printed-parameter discrepancy remains unresolved.
- `core-shell-gel/stiff-followup-2882.*`: independent implicit continuation of one saved budget-limited endpoint, followed by Newton/Radau and BDF closure to a stable cycle. [The follow-up](../docs/core-shell-gel-stiff-followup.md) preserves the original incomplete GPU interval and leaves other stiff cases unclassified.
- `candidate-audits/mineral-protocell-screening.json`: primary-source equations and scoped dimensional/monotonicity exclusions for fixed and growing mineral membranes. [No trajectories were simulated](../docs/mineral-protocell-screening.md).
