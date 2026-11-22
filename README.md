# Chua GPU Atlas

This atlas maps periodic windows in Chua’s circuit. It follows individual periodic orbits as circuit parameters change and identifies settings where different attractors coexist. GPU scans locate candidate windows, numerical continuation traces the orbits, and interval calculations certify selected results.

[Browse the 34 illustrated families](docs/chua-layer-cards.md) · [Find a research report](docs/README.md) · [Run a CPU example](#reproducibility)

## Reading the atlas

Each card pairs a family’s sampled stability region (its **layer**) with a reference orbit and one period of capacitor voltage. The 34 reference orbits use different parameter settings. A04 and A01 both have 25 voltage peaks per period. Numerical knot analysis identifies A04 as a $(2,5)$ torus knot and A01 as an unknot. [Compare their orbit diagrams](docs/chua-knot-atlas.md)

![A04 family card showing its sampled stability region, periodic orbit, 25-peak voltage trace and knot diagram](docs/assets/chua-layer-cards/A04.png)

## The system and the question

The state consists of the two capacitor voltages $x,y$ and the inductor current $z$, all in dimensionless units. The classical atlas uses zero inductor loss ($\gamma=0$). The circuit studies below also treat positive loss. The diode characteristic $h$ has inner slope $m_0$, outer slope $m_1$ and breakpoints at $x=\pm1$.

```math
\dot x=\alpha(y-x-h(x)),\qquad \dot y=x-y+z,\qquad \dot z=-\beta y-\gamma z,
```

```math
h(x)=m_1x+\tfrac12(m_0-m_1)(|x+1|-|x-1|).
```

On the classical transect, $\alpha=9$, $m_0=-8/7$ and $m_1=-5/7$. Varying $\beta$ reveals windows in which periodic orbits attract. [Parameter conventions and survey methods](docs/methods.md#equations-and-parameter-conventions)

## Research questions

| Question | Where to start |
|---|---|
| Where do periodic windows accumulate? | [Infinitely many windows near a certified parameter](docs/chua-infinity-certificate.md) · [Numerical returns through 245 peaks](docs/chua-nested-returns.md) |
| Which attractors coexist, and what separates their basins? | [Hidden periodic coexistence](docs/chua-hidden-infinity-certificate.md) · [A certified local basin boundary](docs/chua-basin-and-coexistence.md) |
| Can two attracting orbits be linked? | [Proof of a non-split link between S02 and S08 at the same parameters](docs/chua-link-certificate.md) |
| When do all bounded trajectories approach the equilibrium set? | [Exact parameter criterion for all memoryless, locally Lipschitz diode laws](docs/chua-global-relaxation.md) |
| What happens at large $\alpha$? | [Quadratic scaling limit](docs/chua-scaling-structure.md) · [Finite-rate passage through a fold](docs/chua-dynamic-fold.md) |
| What can one capacitor voltage reveal? | [Mean voltage at symmetry breaking](docs/chua-symmetry-dc.md) · [Spectrum from diode crossings](docs/chua-diode-spectrum.md) · [Reconstruct the other two coordinates](docs/chua-voltage-reconstruction.md) |

### Hidden attractors in a lossy circuit

How precisely must a capacitor be charged to keep a trajectory near a saddle cycle? At one published circuit setting, bringing the initial voltage ten times closer to the threshold adds **2.899 turns** to the predicted residence. In numerical integrations the forecast gave the exact exit count for **23 of 24** starts and missed by one turn on the other. [Forecast, plot and numerical checks](docs/chua-saddle-residence.md) · [Published circuit, Figs. 9–10](https://doi.org/10.1007/s11071-022-08078-y)

The circuit reports also compare measured and predicted periods, test component tolerances and locate coexisting hidden cycles. [Circuit reports and data](docs/README.md#hidden-attractors-and-circuit-experiments)

## Reproducibility

From the repository root, install the CPU dependencies and check a point on the convergence boundary.

```bash
python -m pip install -r requirements-certificate.txt
python research/chua_global_relaxation.py --classify 9 25 0
```

At $(\alpha,\beta,\gamma)=(9,25,0)$ the command returns an inner diode slope $m_0=-5/9$ and a squared angular frequency $\omega^2=20$. With that slope the circuit admits an exact sinusoidal orbit.

To replay the 24 voltage starts, run `python research/observed_saddle_validate.py --figure`. It writes the exit counts and figure to `outputs/chua-saddle-residence/`. [Other CPU calculations](research/README.md) · [Recorded data](results/README.md) · [Run the FP64 OpenCL survey](docs/methods.md)

## Background

My university coursework involved studying the Chua system by scanning two parameters. I spent months running calculations. I did not understand the subject very well, and my computer was weak. Much of the experience was simply waiting for another calculation to finish.

I decided to return to the problem and see whether there had been anything interesting hiding in it. This time I could use a GPU, vary three parameters, and look inside the shapes instead of just waiting for them to appear.

[MIT license](LICENSE) · [Third-party notices](THIRD_PARTY_NOTICES.md) · [Citation](CITATION.cff) · [Other dynamical systems](https://github.com/Griboedow/chaos-attractor-lab)
