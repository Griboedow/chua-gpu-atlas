# Chua GPU Atlas

## Why I came back to this

My university coursework involved studying the Chua system by scanning two parameters. I spent months running calculations. I did not understand the subject very well, and my computer was weak. Much of the experience was simply waiting for another calculation to finish.

I decided to return to the problem and see whether there had been anything interesting hiding in it. This time I could use a GPU, vary three parameters, and look inside the shapes instead of just waiting for them to appear.

The months are a personal recollection; those calculations and today's scans are not a matched speed benchmark.

## What this study found

A study of periodic layers, coexisting attractors and thin stability windows in the continuous piecewise-linear Chua system. GPU surveys locate candidate regions; periodic-orbit continuation resolves their structure; interval arithmetic checks a mechanism that generates infinitely many windows.

**Main result:** at α = 9, m₀ = −8/7 and m₁ = −5/7, infinitely many attracting periodic windows accumulate near **β = 14.30635909996146**. An executable interval certificate establishes a quadratic homoclinic tangency and the hypotheses of the classical sink-cascade theorem. The finite atlas contains **34 reference identities**, of which 32 have sampled three-dimensional sheets.

[Research guide](docs/research-guide.md) · [Proof and interval certificate](docs/chua-infinity-certificate.md) · [Motion gallery](docs/chua-family-gallery.md) · [Data and reproduction](results/README.md)

![Measured periodic windows and the certificate establishing an infinite cascade](docs/assets/chua-infinity-overview.png)

The infinite sequence concerns different periodic motions at different parameter values. It does not count globally disconnected surfaces, and it does not assert infinitely many attractors at the limiting parameter itself. The certificate encloses the tangency in **14.30635909996140 < β∗ < 14.30635909996152**. Its assumptions, derivation, positive checks and rejected controls are given in the [proof report](docs/chua-infinity-certificate.md).

## What a layer represents

A layer is a region in **parameter space** supporting a family of attracting periodic motions. Following that family changes its orbit continuously. Different families can differ in their complete repeating signal, period, visits to the two scrolls and symmetry. The corresponding orbit is a curve in **state space**.

![The 32 sampled families highlighted individually at a fixed camera angle](docs/assets/chua-readme-layers.gif)

The camera and coordinates stay fixed while individual families are highlighted. This view uses **5,447 checked periodic-cycle samples**, including 3,975 added parameter locations; 4,126 samples lie in the displayed box. Dots are closed attracting cycles. Facets show sample adjacency; their drawn thickness does not represent a window's physical width. The animation is the 32-family mesh snapshot. The later 75- and 85-peak entries have checked reference cycles but have not been meshed.

| Difference in motion | Example |
|---|---|
| Repetition time and waveform | [S08 and S11](docs/chua-family-gallery.md#s08) both have five peaks but periods of about 10.89 and 14.50. |
| Number of peaks before the motion repeats | The measured thin sequence has 15, 25, …, 85 peaks per full cycle. |
| Scroll visits | [E06](docs/chua-family-gallery.md#e06) stays near one outer scroll; the other mesh representatives visit both. |
| Symmetry | A symmetric cycle can lose symmetry and produce two inversion-related attracting cycles. |

Peak counts and orbit silhouettes alone do not identify a family. The catalogue follows normalized periodic roots and their continuation ancestry. The [gallery](docs/chua-family-gallery.md) shows an orbit, complete signal and peak sequence for each original family; [two further portraits](docs/chua-layer-prediction.md) show the 75- and 85-peak motions. The [refined atlas](docs/chua-refined-atlas.md) documents the identity convention and numerical checks.

## Prediction, search and proof

A short unstable periodic orbit organizes the measured 15 → 25 → … sequence. Longer cycles spend additional returns near this saddle before making a global excursion. Its expansion factor predicts the rate at which window spacing and width decrease, while its period predicts the additional repetition time.

Forecasts for the 75- and 85-peak members were [published before the search](https://github.com/Griboedow/chua-gpu-atlas/commit/0fc47ea). Both motions were subsequently found and checked independently. Their symmetric-core widths differed from the forecasts by **1.89% and 2.70%**.

![Published forecasts compared with the measured windows at their true width scale](docs/assets/chua-layer-prediction-proof.png)

The 75-peak forecast overlaps its core. The 85-peak location forecast misses its much narrower core and requires orbit correction. This miss remains part of the result. That core is about **6.65 million times narrower** than the old β grid spacing, which explains why a uniform scan can miss whole families.

The [prediction report](docs/chua-layer-prediction.md) preserves the frozen forecasts, measured orbits and failed numerical controls. The [subsequent certificate](docs/chua-infinity-certificate.md) bounds the exact invariant manifolds and proves the organizing tangency. It checks entire flight intervals, curvature, dissipativity and parameter unfolding. A classical theorem then supplies the infinite cascade; the finite sequence is not extrapolated into a proof.

## Where layers overlap

Different attracting motions can coexist at the same parameters. Their basins determine which one a trajectory reaches. Following only one initial state can therefore miss a family even on a dense parameter grid.

![Symmetric coexistence islands connected through a three-cycle region](docs/assets/chua-overlap-story.png)

For two selected families, the symmetric coexistence region splits into islands while a region with **three attracting cycles** connects their primary branches. Nearby parameters support four cycles. The [overlap study](docs/chua-overlap-geometry.md) follows symmetric and asymmetric branches, traces 18 contact curves and checks representative motions independently. [Orbit portraits](docs/assets/chua-overlap-orbits.png) and [magnified section points](docs/assets/chua-overlap-sections.png) show the different outcomes.

## Varying both diode slopes

The [m₀–m₁ section](docs/chua-slope-plane.md) and the two volumes below vary both slopes independently. The volumes use **(m₀,m₁,α)** at β = 100/7 and **(m₀,m₁,β)** at α = 9.

![Interior views of the m0-m1-alpha and m0-m1-beta parameter volumes](docs/assets/chua-slope-volumes-interior.png)

Each volume contains 152,561 sampled parameter nodes with three initial states. These are finite-time classifications; they can miss thin windows and small basins. [Continued five- and seven-peak families](docs/chua-slope-families.md) resolve parts of this structure and establish coexistence. The [volume report](docs/chua-slope-volumes.md) provides sections, numerical limits and an offline interactive viewer. [Slow rotation](docs/assets/chua-slope-volumes.gif).

## Model

```math
\dot x=\alpha[y-x-h(x)],\qquad \dot y=x-y+z,\qquad \dot z=-\beta y,
```

```math
h(x)=m_1x+\frac{m_0-m_1}{2}\left(|x+1|-|x-1|\right).
```

m₀ and m₁ are the inner and outer diode slopes. Some conventions absorb the x term into h and shift both slopes by one. The classic starting section uses α = 9, β = 100/7, m₀ = −8/7 and m₁ = −5/7. Background: [Chua circuit](https://doi.org/10.4249/scholarpedia.1488), [Kennedy's circuit primer](https://doi.org/10.1109/81.246141), and [Muthuswamy's model discussion, Eqs. 11–14](https://www.bharathm.org/pubs/MuthuswamyWSPCSampleChapter-ChuaCircuit.pdf).

## Reproduce the research

The interval certificate and its independent checks run on a CPU:

```bash
python -m venv .venv
# PowerShell: .\.venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements-certificate.txt
python research/certify_chua_homoclinic_tangency.py --digits 100
python research/check_chua_interval_derivatives.py
python research/check_chua_global_connection.py
python research/manifest_refined_chua_atlas.py --check
```

See the [certificate reproduction guide](results/chua-layer-infinity/README.md) for the 80-digit repeat, center reconstruction, rejected controls and figures. The recorded environment used Python 3.14.3 and python-flint 0.9.0.

The original scanner requires an OpenCL device with FP64 support:

```bash
python -m pip install -r requirements.txt
python run.py devices
python run.py all --grid 21 --transient 100 --duration 200 --output outputs/quick
```

[Methods](docs/methods.md) covers the full scans, classifications, numerical controls and array layouts. The GPU surveys were run on an NVIDIA GTX 1080; scan timings exclude independent CPU checks and plotting.

| Evidence | Location |
|---|---|
| Interval proof, positive and negative controls | [Certificate files](results/chua-layer-infinity/README.md) |
| Published forecasts and 75/85-peak checks | [Prediction files](results/chua-layer-prediction/README.md) |
| Family meshes, ancestry and sample audits | [Atlas report](docs/chua-refined-atlas.md) |
| Symmetry and coexistence | [Overlap files](results/chua-overlap-geometry/README.md) |
| Wider chaos surveys and large-α scaling | [Research guide](docs/research-guide.md#broader-surveys-and-scaling) |
| Large historical arrays and offline reports | [Release v1.1.0](https://github.com/Griboedow/chua-gpu-atlas/releases/tag/v1.1.0) |

The release predates the layer prediction and interval proof; their data are in this repository. The [manifest](results/chua-layer-atlas/manifest.json) records hashes of the research files. Historical failed checks remain available and are identified in their reports. A positive finite-time Lyapunov exponent is treated as a chaos candidate, not as a proof of an asymptotic attractor.

## Scope and related work

Periodic windows, coexistence and homoclinic cascades are established phenomena. This repository contributes explicit parameter studies, reproducible forecasts and a computer-assisted certificate for the stated transect. Historical priority for these Chua structures has not been established. The certificate is inspectable mathematical and numerical work; it has not been checked by a formal proof assistant.

Investigations of other physical and biological models are maintained in [Chaos Attractor Lab](https://github.com/Griboedow/chaos-attractor-lab). Chua data and reports remain here.

## License

Original code, text, figures and datasets use the [MIT License](LICENSE). Dependencies retain their own licenses; bundled Plotly notices are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Publications are linked, not redistributed.
