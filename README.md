# Chua GPU Atlas

**Revisiting an old coursework problem: from a two-parameter scan that took months to a three-dimensional numerical atlas of chaos.**

Python + OpenCL · actual GPU calculations · reproducible data · MIT license

![The expanded parameter cloud rotating about three axes](docs/assets/expanded-rotation.gif)

*Each dot is a parameter tuple (α, β, m₀). These are 1,832 finite-time chaos candidates from the expanded adaptive GPU survey. Gold outlines the previous search domain. The camera rotates about all three coordinate axes; the data do not change. The axes have different scales, and the cloud is not a filled or proven existence region.*

## Why I came back to this

My university coursework involved studying the Chua system by scanning two parameters. I spent months running calculations. I did not understand the subject very well, and my computer was weak. Much of the experience was simply waiting for another calculation to finish.

I decided to return to the problem and see whether there had been anything interesting hiding in it. This time I could use a GPU, vary three parameters, and look inside the shapes instead of just waiting for them to appear.

There was plenty to look at: thin regular windows inside chaotic regions, different attractors at exactly the same parameters, and a parameter cloud that keeps extending beyond the original search boxes. This repository is a record of that return, including the numerical checks and the places where the evidence is still incomplete.

The old “months” is a personal recollection, **not a speedup benchmark**. The old and new calculations are not matched experiments.

## What is the Chua system?

Chua's circuit is a small nonlinear electrical circuit with two capacitors, an inductor, a resistor, and a nonlinear resistor often called the Chua diode. Its dimensionless model has three state variables and can produce the famous **double-scroll attractor**: motion that keeps folding around two lobes without settling into a repeating cycle [1, 2].

```text
dx/dt = α [y − x − h(x)]
dy/dt = x − y + z
dz/dt = −βy

h(x) = m₁x + ½(m₀ − m₁)(|x + 1| − |x − 1|)
```

Here `m₀` and `m₁` are the inner and outer diode slopes. Some references absorb the `x` term into the nonlinear function, shifting both slopes by one; parameter conventions matter.

This project studies **two different 3D spaces**:

| Space | Coordinates | What a point means |
|---|---|---|
| State space | x, y, z | One instantaneous state of one system |
| Parameter space | α, β, m₀ | A different system, with m₁ held fixed |

![The classic double-scroll trajectory rotating about x, y and z](docs/assets/trajectory-rotation.gif)

*One finite trajectory at α=9, β=100/7, m₀=−8/7, m₁=−5/7. Only the camera moves.*

## What turned out to be interesting?

**A tiny parameter change can remove the observed chaos.** At the classic slopes and α=9, changing β from 100/7 to 14.32 is only a 0.24% change. The new trajectory has a repeating five-peak pattern and an exponent near zero. The initial coarse grid missed the window; a local 41³ zoom reveals thin inclined sheets of regular behavior.

![Rotation of the local regular-window surface](docs/assets/window-rotation.gif)

*Pink: an unsmoothed surface extracted from nodes with near-zero/nonpositive exponents for all three starts. It is clipped at the search box. A five-peak pattern was checked at the highlighted example, not at every surface vertex.*

**The same parameters can support two different chaotic outcomes.** At α=15.6, β=32, m₀=−1.143 and m₁=−0.714, symmetric initial states approach two symmetric single-scroll attractors. A 151² initial-condition slice shows their different basins. At β=29 in the same section, the calculation instead finds a six-peak regular window.

![Lyapunov transect, peak diagram, coexisting scrolls and basin map](docs/assets/findings.png)

**The original region was only a small piece.** A wider Sobol survey finds repeated candidates near α=256 and m₀=−8, with β reaching about 12,022. The search still reaches its imposed limits. We have not established whether the complete parameter set is bounded, or whether all distant points connect to the original component.

**There is an analytic clue to the bend.** For m₁=−5/7, the outer equilibria change linear stability at

```text
β = 5α/7 + 10α²/49.
```

All 2,249 repeated wide-survey candidates lie below this surface. It provides a useful growing β scale; it is not a proven global boundary of chaos.

These types of Chua phenomena are known in the literature [1, 2]. The contribution here is a reproducible numerical exploration of particular sections, not a claim to have discovered periodic windows or coexistence for the first time. More interpretation and open questions are in [the experiment notes](docs/experiment-notes.md).

## Four times farther: the expanded survey

The next experiment extends α from 256 to **1024**, and the slope depth −1−m₀ from 7 to **28** (m₀ down to −29). Sampling β relative to its quadratic stability scale extends the upper search limit to about **279,144**.

Of **16,384** new Sobol parameter samples, **1,832** remain candidates after a second calculation with tighter tolerances and longer observation. Their separate coordinate extrema are:

| Coordinate | Previous wide cloud | Expanded cloud |
|---|---:|---:|
| Largest α | 255.69 | **1023.64** |
| Largest β | 12,022.30 | **173,658.02** |
| Most negative m₀ | −7.9912 | **−28.9107** |

The expanded cloud contains 732 candidates beyond α=256 and 440 beyond m₀=−8; 981 lie outside the previous search domain in at least one coordinate constraint. There are still candidates near both newly extended limits, so **a natural outer cutoff has not been located in those directions**. Finite sampling does not establish unboundedness or connectedness.

The new FP64 GPU backend uses adaptive Dormand–Prince 5(4), controlling both the state and tangent vector. Scout and candidate verification took **281.5 seconds** in total on the GTX 1080, excluding independent CPU checks and plotting. The escape cutoff, sampling measure and integration settings differ from the previous survey, so candidate counts must not be compared as volumes.

![Old search limits inside the expanded cloud, in linear and logarithmic views](docs/assets/expanded_extent.png)

Four selected distant examples, including the largest α, largest β and most negative m₀, retain positive exponents in every window of independent full-field DOP853 checks. Their exponent magnitudes are not uniformly converged; the new report shows both solver estimates.

The new report adds a view of β/βcrit, with βcrit=5α/7+10α²/49, to separate the quadratic bend from the remaining shape. All repeated candidates in this sample remain below that stability surface. Small original regions and thin windows are still easy to miss in a coarse global survey.

The previous wide-map animation remains available for comparison:

<details><summary>Earlier cloud, before the fourfold expansion</summary>

![Original wide parameter cloud](docs/assets/extent-rotation.gif)

</details>

## Results and downloads

[**Download the complete results from release v1.1.0**](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.1.0).

The release contains `chua-results-v1.1.0.zip` and a SHA-256 checksum. Extract it into the repository root to obtain `outputs/`. Then open:

- `outputs/index.html` — the original atlas;
- `outputs/detailed/index.html` — refined surfaces, slices, bifurcation diagrams and basins;
- `outputs/extent/index.html` — the earlier wide point cloud;
- `outputs/expanded/index.html` — the fourfold expansion, old/new overlay and normalized β view.

Reports are **in English and work offline**, with a bundled Plotly runtime. They include links to NPZ arrays and JSON settings/validation records. Large raw arrays live in the release instead of Git history. The repository includes small illustration datasets, all GIFs and selected [evidence records](results/README.md).

| Experiment | Parameter samples | Initial states | Measured GPU scan time |
|---|---:|---:|---:|
| Original two local boxes | 2 × 51³ | 3 each | 34.4 + 36.7 s |
| Refined two local boxes | 2 × 81³ | 3 each | 137.2 + 147.1 s |
| Thin-window zoom | 41³ | 3 each | 83.4 s |
| Wide scout + candidate reruns + selected long reruns | 32,768 scout points | 3 each | about 15.8 min |
| Expanded adaptive scout + candidate reruns | 16,384 scout points | 3 each | 88.2 + 193.2 s |

Hardware: **NVIDIA GeForce GTX 1080, 8 GB**, FP64 OpenCL. Timings exclude plotting, archive creation and independent CPU checks; workloads use different integration settings. Exact metadata is retained with each experiment.

## Run it

Use **64-bit Python 3.11+** and a GPU driver exposing OpenCL with FP64 support. The recorded environment used Python 3.14.3. A separate CUDA toolkit is not required. The scanner requires an actual GPU; there is no silent CPU fallback.

```bash
python -m venv .venv
# Activate: source .venv/bin/activate (Linux/macOS)
# Activate: .\.venv\Scripts\Activate.ps1 (PowerShell)
python -m pip install -r requirements.txt
python run.py devices

# Small exploratory run
python run.py all --grid 21 --transient 100 --duration 200 --output outputs/quick

# Reproduce the original atlas and its numerical checks
python run.py all --grid 51
python validate.py
python run.py report

# Re-render all four GIFs from included illustration data; no GPU required
python render_animations.py
```

The [methods and reproduction guide](docs/methods.md) covers refined scans, the wide survey, validation, classification codes and array layouts. `requirements-lock.txt` records exact versions from the supplied run; `requirements.txt` gives portable version ranges. `requirements-render.txt` is sufficient if you only want to rebuild the GIFs.

Published starting points:

| Section | α | β | m₀ | Fixed m₁ | Source |
|---|---:|---:|---:|---:|---|
| Classic | 9 | 100/7 | −8/7 | −5/7 | Muthuswamy, Eqs. 11–14 [3] |
| Chua Circuits example | 15.6 | 28 | −1.143 | −0.714 | Published dimensionless example [4] |

These are separate sections of a four-parameter family. The slightly different outer slopes are deliberately preserved.

## How strong is the evidence?

The GPU integrates the state and a tangent vector with RK4. Repeated tangent normalization estimates the largest **finite-time Lyapunov exponent**: positive values indicate growth of small perturbations. To become a candidate, a trajectory must also stay inside the observation cutoff, keep moving, and pass a positive-exponent criterion in all four observation windows.

The original study passes eight numerical checks, including adaptive DOP853 comparisons, a linear exact-solution control, and rejection of a stationary unstable equilibrium. Separate refinement samples preserve candidate membership at 98.83% and 99.02% of sampled nodes. The thin-window check preserves 252/256 sampled categories. These are sample diagnostics, not confidence levels for the entire boundary.

**A reference-solver correction is part of the result.** During the expansion, the original v1.0 CPU routine failed at a distant example: after a missed crossing, it could keep extrapolating the wrong affine branch. Its historical values, including the earlier 0.18 estimate at α≈255.69, should not be treated as definitive reference measurements. The corrected DOP853 control evaluates the actual piecewise field at every stage; crossing events are diagnostics, not branch selectors. New independent checks are stored under `results/expanded/`. The original records remain available for provenance.

The adaptive backend also passes short-trajectory comparisons at α=9, 255.69 and 1024, an exact linear matrix-exponential control, and rejection of a stationary unstable equilibrium. Quantitative convergence of Lyapunov estimates across the full cloud is still not established. Cloud colors indicate m₀, not a claimed intensity of chaos.

A finite-time positive exponent does not prove an asymptotic strange attractor. Three initial conditions do not exhaust all basins. A ragged mesh does not prove fractality, and a finite search cannot settle global boundedness. These limitations matter when reading the pictures.

## Looking beyond Chua

A [separate research log](docs/candidate-research.md) screens unusual related systems against prior publications and investigates a quantum memristive neuron with Python and an adaptive FP64 GPU solver. Across 36,864 survey trajectories and 58 long follow-ups, **no new chaotic attractor was confirmed**. The log includes failed candidates, independent numerical controls, an analytic boundedness argument, and an energy-based exclusion for an uncoupled physical meminductor. Code and compact raw data are included; publication novelty remains unverified.

The next pass found a stronger candidate: **four synchronematic colloid oscillators** in the deterministic reduced equations printed by Leyva et al. (2026). [The new notebook](docs/synchronematic-research.md) includes a full Lyapunov spectrum with a positive leading exponent, smaller-step GPU reruns, geometry and coefficient sensitivity, and an analytic gradient-flow control. An audit also found a coefficient-sign disagreement between the article and the authors' code, which is investigated explicitly. The result is evidence for the literal printed Fourier approximation, not experimental confirmation or a certified first discovery.

![Independent numerical evidence for the synchronematic candidate](docs/assets/synchronematic_evidence.png)

An [optical-feedback notebook](docs/photofeedback-research.md) develops an explicitly assumed response model for the physical meminductor question. It identifies linear instability boundaries and sustained oscillation, but none of its 12 preliminary chaos candidates survive long follow-up. Both positive and negative results are reproducible from the included code and compact data.

A [magnetic energy audit](docs/meminductor-energy-audit.md) excludes recurrent chaos under fixed DC drive in a reciprocal moving-coil model even when electrical response is retained. The source's printed force/flux pair fails that reciprocity test, so this conditional result does not settle its optical-feedback experiment. The notebook distinguishes measured-device questions from illustrative model completions.

A [periodic-window study](docs/synchronematic-bifurcations.md) resolves fold candidates, a stable period-doubled branch and competing numerical outcomes in the four-oscillator control family. A separate [two-oscillator first integral](docs/synchronematic-pair-integral.md) explains an additional dimensional constraint in the pair model. These are mathematical model analyses; the control parameter changes the interaction law.

The [candidate queue](docs/candidate-queue.md) records additional physical systems, known prior results, dimensional exclusions and missing source material. An absence of search hits is never treated as proof of novelty.

A second positive lead is a **light-switching elastomer slider**. Its [research notebook](docs/lce-slider-research.md) reduces the published equations to four effective parameters and uses analytic flight segments, GPU return maps and independent DOP853 checks. An expanded search finds seven parameter tuples with persistent positive exponents, beyond the original paper's numerical stiffness/damping ranges. The representative trajectory remains small in strain and molecular concentration; realistic switching and publication novelty still need investigation.

![Independent event-resolved evidence for the slider candidate](docs/assets/lce_slider_evidence.png)

The [switching-width follow-up](docs/lce-slider-switching.md) finds that chaos survives narrow smooth transitions, while wider examples give stable periodic orbits and then rest. It supplies orbit-closure/Floquet controls, an equilibrium stability criterion, and a comparison with older relay equations. Physical switching resolution remains a substantive limitation.

The [weak-light limit](docs/lce-slider-relay-limit.md) reduces the ideal slider to a classical biased relay around a stable third-order plant. Independent CPU/GPU calculations retain positive exponents in the rescaled limit; a 102-offset survey finds chaotic samples away from a stable symmetric periodic orbit. The physical amplitudes shrink with intensity and the threshold must be retuned, so this does not predict motion in darkness. The connection to older relay theory narrows the possible novelty claim.

A separate [electrical rheostat oscillator](docs/rheostat-research.md) yields a negative chaos search across 4,096 scouts and 24 long checks. Its exact three-parameter reduction exposes a finite thermal-rate instability window. Independent periodic-orbit closure resolves the strongest weak positive exponent, and near-onset amplitudes agree with the quadratic-drag scaling.

A [16-cell epithelial-chain audit](docs/epithelial-chain-audit.md) separates mathematical expansion from physical validity: all 343 preliminary positive-exponent candidates have negative cell lengths. A source audit of [light-mediated operators](docs/light-operators-audit.md) identifies an exact equivalence between another published single photothermal oscillator and our existing relay limit.

The [six-state coupled optical operators](docs/coupled-operators-research.md) now provide a third independently checked numerical lead: 4,096 scouts, 436 long follow-ups, independent CPU spectra and a positive smooth optical-profile control. Their physical projections, Poincare sections and rotation GIF are included. The result concerns a continuous-heat interpretation of the published equations at changed parameters; novelty and experimental realization remain unconfirmed.

A [finite optical-width survey](docs/coupled-operators-research.md#smoothing-is-not-a-monotone-chaos-filter) finds that smoothing can remove and then restore a positive-exponent regime in the six-state model. Independent CPU spectra and a stable periodic control strengthen the evidence beyond a sharp-switch simulation; asymptotic chaos and experimental calibration remain open. The [one-transistor OECT model](docs/oect-planar-exclusion.md) is instead excluded from this autonomous chaos search by its two-state structure.

A [20,000-unit follow-up](docs/coupled-operators-research.md#long-runs-reveal-periodic-capture-after-chaotic-transients) exposes a concrete limit of short chaos scans: at optical width `0.005`, 29 of 32 new initial conditions eventually approach a stable cycle, sometimes after thousands of chaotic time units. At width `0.04`, 26 of 32 remain chaotic throughout the observation; six approach a coexisting cycle early. All 54 periodic endpoints across three widths were independently checked on the CPU. These are finite-time observations, not proofs of asymptotic attraction.

![Long optical transients and persistent candidates](docs/assets/coupled_operators_lifetimes.png)

An additional [geometry search](docs/attractor-geometry.md) asks whether a candidate has a substantially different shape and dynamical organization. A comparison figure, rotation GIF and finite-spectrum dimension estimates establish a baseline from existing data; they do not certify a new attractor type.

Further screening includes a [thermomechanical sheet study](docs/thermomechanical-sheet-audit.md): 8,192 scouts and 84 follow-ups reveal strong numerical chaos only outside the assumed physical domain, while a weak positive is independently identified as a stable periodic orbit. The study preserves unresolved cases and corrects a source equilibrium formula. A [single-macrospin exclusion](docs/macrospin-exclusion.md) supplies a different outcome: its deterministic axial dynamics reduces to a scalar tilt equation, so three plotted magnetization components do not provide an autonomous chaos mechanism.

## License

Project code, original text, generated figures and numerical datasets are provided under the [MIT License](LICENSE). Dependencies retain their own licenses. The results archive bundles an unmodified Plotly.js runtime; its MIT notice is retained and included in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Referenced publications are linked, not redistributed.

## References

1. **Chua, L. O.** (2007). Chua circuit. *Scholarpedia*, **2**(10), 1488. [doi:10.4249/scholarpedia.1488](https://doi.org/10.4249/scholarpedia.1488).
2. **Kennedy, M. P.** (1993). Three steps to chaos. II. A Chua's circuit primer. *IEEE Transactions on Circuits and Systems I: Fundamental Theory and Applications*, **40**(10), 657–674. [doi:10.1109/81.246141](https://doi.org/10.1109/81.246141). [Author-hosted PDF](https://people.eecs.berkeley.edu/~chua/papers/Kennedy93b.pdf).
3. **Muthuswamy, B.** (2012). Chua's Equation was Proved to be Chaotic in Two Years, Lorenz Equation in Thirty Six Years. Author-hosted chapter manuscript, dated March 25, 2012; Eqs. 11–14, printed p. 6. [PDF](https://www.bharathm.org/pubs/MuthuswamyWSPCSampleChapter-ChuaCircuit.pdf).
4. **Chua Circuits.** (n.d.). MATLAB: dimensionless Chua's circuit. [Online example](https://chuacircuits.com/matlab/). Accessed October 29, 2026.
5. **Matsumoto, T.** (1984). A chaotic attractor from Chua's circuit. *IEEE Transactions on Circuits and Systems*, **31**(12), 1055–1058. [doi:10.1109/TCS.1984.1085459](https://doi.org/10.1109/TCS.1984.1085459). [Berkeley technical-report version](https://www2.eecs.berkeley.edu/Pubs/TechRpts/1984/307.html).
6. **PyOpenCL contributors.** PyOpenCL documentation. [Runtime documentation](https://documen.tician.de/pyopencl/). Accessed October 29, 2026.
7. **SciPy contributors.** `scipy.integrate.solve_ivp`: initial-value ODE solvers, including DOP853. [API documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). Accessed October 29, 2026.
