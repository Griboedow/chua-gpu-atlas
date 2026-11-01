# Methods and reproduction guide

Run commands from the repository root. See the [main README](../README.md) for the study summary, figures and links to the reports and datasets.

Python code, actual GPU results, and an offline interactive report for two standard Chua double-scroll parameter sets. The parameter survey runs through **PyOpenCL on a GPU**. There is no silent CPU fallback. OpenCL is supported by the installed NVIDIA driver; a separate CUDA toolkit is not required.

## Open the results

Open `outputs/index.html` in a modern browser. Keep `outputs/plotly.min.js` next to it. Everything is local; no CDN or server is needed. Rotate the phase portrait and parameter surface, switch between the three-seed intersection and union, inspect Lyapunov-colored points, and view categorical slices.

- `outputs/overview.png` / `overview.pdf`: static overview.
- `outputs/*_scan.npz`: every sampled parameter, per-seed classification and Lyapunov diagnostics.
- `outputs/*_scan.json`: GPU identity, precise settings, timings, and class definitions.
- `outputs/*_trajectory.npz`: trajectory and relative time.
- `outputs/validation.json`: independent numerical checks and refinement sensitivity.
- `outputs/*_refinement.npz`: sampled points recomputed at half the time step and twice the transient and measurement duration.

## Installation

Use **64-bit Python 3.11 or newer**, plus a GPU driver exposing OpenCL. FP64 support is required by the default setting. The GTX 1080 used for the supplied results supports it.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run.py devices
```

## Run

```powershell
# Full calculation, plots, and report (two regions, three starts at each node)
.\.venv\Scripts\python.exe run.py all --grid 51

# Independent validation and sampled convergence check
.\.venv\Scripts\python.exe validate.py
.\.venv\Scripts\python.exe run.py report

# Quick experiment in a separate directory
.\.venv\Scripts\python.exe run.py all --grid 21 --transient 100 --duration 200 --output outputs/quick

# A denser survey, or just one center
.\.venv\Scripts\python.exe run.py all --grid 81 --center classic --output outputs/dense

# Rebuild plots from saved arrays without a GPU or reintegration
.\.venv\Scripts\python.exe run.py report
```

The main options are `--grid` (odd integer), `--dt`, `--transient`, `--duration`, `--threshold`, `--precision float64|float32`, `--device` (name substring), and `--config`. Defaults are a 41³ grid, dt=0.01, transient=300 and observation=800 in dimensionless time. **The delivered survey uses 51³**, with all other defaults. Increasing `--grid` from N to 2N roughly multiplies the work by eight. Smaller dt and longer observation improve the screen but do not prove asymptotic chaos.

## Equations and parameter conventions

```text
dx/dt = alpha * (y - x - h(x))
dy/dt = x - y + z
dz/dt = -beta * y
h(x)  = m1*x + (m0-m1)/2 * (abs(x+1)-abs(x-1))
```

`m0` is the **inner diode slope**, `m1` the outer diode slope. Some sources absorb the `x` term into their nonlinear function; their slopes are larger by one. Do not mix these conventions.

Three parameters vary independently: **alpha, beta, m0**. The fourth, **m1**, is fixed separately for each local parameter box. Edit `centers.json` to change the centers, box half-widths, or fixed slope. A center has four values in the order `[alpha, beta, m0, m1]`; half-widths contain three values for the varied axes. The selected literature examples are:

| Center | alpha | beta | m0 | fixed m1 | Box half-widths |
|---|---:|---:|---:|---:|---|
| Classic | 9 | 100/7 | -8/7 | -5/7 | 2, 4, 0.16 |
| Chua Circuits example | 15.6 | 28 | -1.143 | -0.714 | 3, 6, 0.16 |

These are two separate 3D sections of a four-parameter family. The slightly different fixed outer slopes are preserved from their respective sources; the two boxes must not be treated as one section with a single fixed m1.

## Numerical method and interpretation

Each GPU work item evolves one parameter tuple and one initial condition. RK4 jointly advances the ODE and `dv/dt = J(x) v`. The Jacobian uses the appropriate inner/outer slope at every RK stage. The vector field is continuous at x=±1; its Jacobian changes there. A fixed step does not locate the crossing exactly, so convergence checks are necessary. Tangent vectors are normalized every ten steps, including during the discarded transient. No fast-math flags are used.

Three starts are `(0.7,0,0)`, `(-0.2,0.1,0)`, `(1.5,0.02,-1.5)`. None is deliberately placed at an equilibrium. After the transient, the largest finite-time Lyapunov exponent is estimated from accumulated tangent growth in four equal observation windows.

Classifications:

| Code | Meaning |
|---:|---|
| 0 | A coordinate exceeded the escape cutoff 100, or the state became non-finite. This is an operational escape rule, not a theorem of unboundedness. |
| 1 | Late x-range < 1e-5 and negative LLE: equilibrium candidate. |
| 2 | Mean and final-window LLE ≤ 0.01: negative or near-zero candidate; periodicity is not proved. |
| 3 | Uncertain: inconsistent windows or other unresolved behavior. |
| 4 | Bounded in the observed interval, mean LLE > 0.01, all four windows > 0.005, late x-range > 1e-4: chaos candidate. |
| 5 | Tangent-vector numerical failure. |

A parameter node belongs to the displayed region if **at least one tested start** is class 4. The alternative view requires all three starts. This includes single-scroll and double-scroll chaos; the `visits_both_outer_regions` flag only records x<-1 and x>1 at some observed times and is not a topological classification.

The isosurface interpolates the **binary candidate mask** at 0.5. It does not draw a convex hull, smooth away holes, or extrapolate outside the scanned box. The node fraction is a grid statistic, not an exact parameter-space volume. Boundary-face flags identify truncation by the search box. Tiny periodic windows, long chaotic transients, coexisting attractors, and unsampled basins can all change a finer/longer survey. Absence of a candidate does not establish absence of an attractor. Finite-time positive exponents alone are not a mathematical existence proof.

## Validation

### Refined atlas and interpretation

`outputs/detailed/index.html` contains the follow-up study: two 81³ parameter surveys, a 41³ local zoom of a thin recurrent window, two 1,601-point beta transects, local-maxima bifurcation diagrams, four example regimes, and a 151² slice of the initial-condition basins at beta=32. The general grids use the original dt=0.01/transient=300/observation=800 settings. The local zoom and transects use dt=0.005/transient=600/observation=1600; example orbits use longer intervals. Their JSON metadata records exact settings.

```powershell
.\.venv\Scripts\python.exe run.py scan --grid 81 --output outputs/detailed
.\.venv\Scripts\python.exe detail_analysis.py
.\.venv\Scripts\python.exe run.py scan --config window_zoom.json --grid 41 --dt 0.005 --transient 600 --duration 1600 --output outputs/detailed
.\.venv\Scripts\python.exe validate_details.py
.\.venv\Scripts\python.exe detail_report.py
```

The local zoom finds sheet-like regular windows inside the chaotic region. At the classical slopes and alpha=9, beta=14.32 has a repeating five-peak pattern; it differs from the original beta=100/7 by only 0.24%. At the second reference's slopes and alpha=15.6, beta=29 has a repeating six-peak pattern, while beta=32 supports two numerically persistent symmetric single-scroll chaotic attractors. DOP853 independently confirms the five- and six-peak patterns. These are recurrence diagnostics, not rigorous mathematical period proofs.

`validate_details.py` also checks the GPU peak ring buffer against a retained trajectory, verifies the basin slice's central symmetry, and rechecks 256 stratified zoom nodes at dt=0.0025 with doubled durations: 252 retain their original category. `analysis_summary.json` records grid connectivity and seed robustness. Those statistics do not imply a fractal boundary, exhaustive basin coverage, or a rigorous existence region.

`*_mesh.npz` exports triangular surfaces with vertex coordinates `(alpha,beta,m0)`. Marching cubes uses every grid node, no smoothing, no convex hull, and no artificial closure at the scan boundary. The local regular-window mesh uses a near-zero-LLE criterion for all starts; peak periodicity is verified on the central transect examples, not every mesh vertex. `detail_analysis.py` adds a small OpenCL event kernel to retain local maxima without transferring complete trajectories.

### Wide exploration of the outer extent

`outputs/extent/index.html` shows a coarse point cloud well beyond the initial boxes, with both linear and logarithmic views. This section fixes m1=-5/7 throughout; it does not merge the slightly different m1 sections of the two original examples.

```powershell
.\.venv\Scripts\python.exe explore_extent.py
.\.venv\Scripts\python.exe validate_extent.py
.\.venv\Scripts\python.exe extent_report.py
```

The 32,768 scrambled Sobol points sample log(alpha), log(beta/alpha), and log(-1-m0), with alpha in [2,256], beta/alpha in [0.25,64], and -1-m0 in [0.003,7]. Three initial conditions include parameter-dependent neighborhoods of the outer equilibrium. The time step is grouped by the local rates; candidates are independently rerun at half the time step with a longer observation. The escape cutoff is 1000 because far-field attractors can be larger. The point cloud contains 2,249 candidates surviving both runs; it is not an interpolated solid or a guarantee that all sampled regions connect.

The new `GPU.scan(..., initial_states=array)` option accepts a finite `(N,S,3)` array for parameter-dependent starts; ordinary shared seeds remain the default. A DOP853 check of the point-to-seed mapping is recorded in `paired_initial_state_check.json`.

The original v1.0 `event_resolved_checks.json` files used a CPU routine that extrapolated one affine branch until a detected event. During the later expansion it failed at a distant point after a missed crossing. These records are retained for provenance, but their exponent values are not definitive references. `validate_extent.py` now uses the corrected full piecewise field and writes `full_field_checks.json`; the new expanded study has its own independent checks. No convergence certificate for the old full cloud is claimed.

The outer-equilibrium stability surface for m1=-5/7 is beta=5*alpha/7+10*alpha²/49. Every repeated candidate in this sample lies below it. This is a local linear stability result, not a proof that this surface exactly bounds global chaos. Candidates occur close to alpha=256 and m0=-8; no global closure or unboundedness follows from finite sampling. Only the m0<-1 continuation was explored.

### Fourfold expansion with adaptive GPU integration

The new section preserves the old data and uses `outputs/expanded/`. Its 16,384 scrambled Sobol points are uniform in log α, log(β/βcrit), and log(−1−m₀), where βcrit=5α/7+10α²/49. Ranges are [2,1024], [0.01,1.3], and [0.003,28], with m₁=−5/7. A rate-based sorting permutation groups similarly expensive trajectories on GPU lanes; it does not change the sampled points.

```bash
python validate_adaptive.py
python explore_expanded.py
python validate_expanded.py
python expanded_report.py
python render_animations.py --kind expanded
```

`explore_expanded.py --resume` reuses completed scout/verification files after verifying that their parameter indices match the requested Sobol sample. It does not resume an interrupted individual integration.

`adaptive_gpu.py` and `adaptive_kernels.cl` implement Dormand–Prince 5(4), with a separate clock for each trajectory and maximum scaled local error over all six state/tangent components. The tangent vector is normalized after every accepted step. Each kernel submission attempts at most 128 steps per trajectory. No fast-math flags are used. Steps are capped at min(0.1,0.5/√β); small accepted steps can cross diode kinks, so independent controls remain necessary.

The scout uses rtol=10⁻⁶, atol=10⁻⁹, transient=100 and observation=200. All 1,950 candidate parameter tuples are recomputed from the same three original starts at rtol=10⁻⁸, atol=10⁻¹¹, transient=200 and observation=600; 1,832 retain a class-4 start. The escape cutoff is 10⁶. The largest recorded coordinate among accepted verification runs is about 17,466, comfortably below that operational cutoff. No numerical-failure statuses were returned.

The class-4 criterion is unchanged: mean LLE>0.01, all four windows>0.005, late x-range>10⁻⁴ and no escape. Other codes match the original scanner. `stats` sums are over adaptive step endpoints; unlike a uniform grid, they cannot be interpreted as time averages. `step_counts` records accepted and rejected steps for each parameter/start pair.

`independent_reference.py` uses SciPy DOP853 on the full piecewise state and variational equations, rtol=2×10⁻⁹ and atol=2×10⁻¹¹. Tangent normalization every 0.25 time units uses a scaled norm to avoid squaring overflow. Nonterminal x=±1 events count observed crossings but never select the vector field. This corrects a failure mode in the older affine-branch reference; missed diagnostic events cannot trap the dynamics on the wrong branch.

For the selected new extremes, an independent check starts at the end of the verified GPU run, discards 20 time units and observes 80 more. An additional GPU calculation uses the identical starting state and intervals with tolerances 10⁻⁹/10⁻¹². Their values and per-window signs are retained together. This is validation of selected examples, not every point or the asymptotic parameter boundary.

### Original survey results

The supplied run evaluated 132,651 parameter nodes per box and three initial conditions per node: 795,906 trajectories in total. The two scans took 34.4 and 36.7 seconds on the GTX 1080 (FP64), excluding setup, plotting, and validation. The candidate union contains 44,273 classic-box nodes and 60,238 second-box nodes. Reference-point LLE estimates are 0.3380 and 0.4235. Eight independent checks passed. Candidate membership agreement in the separate 512-node-per-box refinement samples was 98.83% and 99.02%. These percentages describe that stratified sample, not a guaranteed accuracy of the entire boundary.

`validate.py` checks short GPU trajectories against adaptive SciPy DOP853, including time-step refinement; compares the LLE sign and approximate magnitude with independent adaptive variational integration; checks a smooth linear case against its exact matrix exponential and eigenvalues; rejects an exactly stationary unstable equilibrium as chaos; and confirms escaped runs do not retain misleading finite LLEs.

It also re-evaluates 256 candidate and 256 other nodes per box, using dt/2, twice the transient, and twice the observation duration. This stratified sample is a reproducible sensitivity diagnostic. It is not a statistical confidence interval or a convergence guarantee for the complete grid. The main surface uses the original grid; refined sample results are stored separately.

## Data layout

Arrays in `*_scan.npz` use flattened `ij` indexing, with m0 changing fastest. `parameters` has shape `(N³,4)`, `classes` and `lle` `(N³,3)`, and `lle_windows` `(N³,3,4)`. `candidate_any` and `candidate_all` have shape `(N³,)`. Reshape these to `(N,N,N)` for alpha/beta/m0 slices. `stats` per seed stores min x, max x, sum x, sum x², maximum absolute coordinate, reserved zero, last-window min x, last-window max x. Escaped/incomplete runs have NaN Lyapunov values and must be masked using `status` before interpreting their statistics.

```python
import numpy as np
d = np.load('outputs/classic_scan.npz')
region = d['candidate_any'].reshape(51, 51, 51)
center_lle = d['lle'][len(d['lle']) // 2]
```

`kernels.cl` contains the small OpenCL C kernels executed by the GPU. Python handles configuration, launching, validation, storage, and drawing. Short kernel submissions avoid running a single multi-minute kernel on a Windows display GPU.

## Sources

- [Leon O. Chua, “Chua circuit”, Scholarpedia (2007)](https://www.scholarpedia.org/article/Chua_circuit): circuit and mathematical context.
- [Bharathwaj Muthuswamy, “Chua’s Equation was Proved to be Chaotic in Two Years, Lorenz Equation in Thirty Six Years”](https://www.bharathm.org/pubs/MuthuswamyWSPCSampleChapter-ChuaCircuit.pdf): equations 11–14 and classical parameters on page 6; note the shifted slope convention in equations 16–19.
- [Chua Circuits, dimensionless simulation](https://chuacircuits.com/matlab/): second reference point, with the published rounded slopes.
- [PyOpenCL documentation](https://documen.tician.de/pyopencl/): GPU runtime API.
