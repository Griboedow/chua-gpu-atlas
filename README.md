# Chua GPU Atlas

**Revisiting an old coursework problem: from a two-parameter scan that took months to a three-dimensional numerical atlas of chaos.**

Python + OpenCL · actual GPU calculations · reproducible data · MIT license

![Rotation of the wide parameter cloud about three axes](docs/assets/extent-rotation.gif)

*Each dot is a parameter tuple (α, β, m₀). These are 2,249 finite-time chaos candidates surviving two GPU calculations. The camera rotates about all three coordinate axes; the data do not change. The axes have different scales, and the cloud is not a filled or proven existence region.*

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

## Results and downloads

[**Download the complete results from release v1.0.0**](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.0.0).

The release contains `chua-results-v1.0.0.zip` and a SHA-256 checksum. Extract it into the repository root to obtain `outputs/`. Then open:

- `outputs/index.html` — the original atlas;
- `outputs/detailed/index.html` — refined surfaces, slices, bifurcation diagrams and basins;
- `outputs/extent/index.html` — the wider point cloud in linear and logarithmic coordinates.

Reports are **in English and work offline**, with a bundled Plotly runtime. They include links to NPZ arrays and JSON settings/validation records. Large raw arrays live in the release instead of Git history. The repository includes small illustration datasets, all GIFs and selected [evidence records](results/README.md).

| Experiment | Parameter samples | Initial states | Measured GPU scan time |
|---|---:|---:|---:|
| Original two local boxes | 2 × 51³ | 3 each | 34.4 + 36.7 s |
| Refined two local boxes | 2 × 81³ | 3 each | 137.2 + 147.1 s |
| Thin-window zoom | 41³ | 3 each | 83.4 s |
| Wide scout + candidate reruns + selected long reruns | 32,768 scout points | 3 each | about 15.8 min |

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

# Re-render all three GIFs from included illustration data; no GPU required
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

**An unresolved numerical issue is part of the result.** At one distant point near α=255.69, a fixed-step GPU estimate is about 0.93 while an independent solver that locates every x=±1 crossing gives about 0.18. Both are positive, but their magnitudes have not converged. The wide cloud is therefore colored by m₀, not by a claimed intensity of chaos. All settings, observation lengths and checks are public.

A finite-time positive exponent does not prove an asymptotic strange attractor. Three initial conditions do not exhaust all basins. A ragged mesh does not prove fractality, and a finite search cannot settle global boundedness. These limitations matter when reading the pictures.

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
