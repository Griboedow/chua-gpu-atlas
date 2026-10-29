# Light-mediated operators: a source and equivalence audit

**Status: the single operator is an additional physical interpretation of our earlier relay limit; the coupled system now has independent numerical chaos evidence.** Guo et al. [1] provide an experimentally motivated optical feedback network. Its single-operator heat-balance model is exactly equivalent to the normalized relay limit already studied in this repository. The separate [six-state coupled investigation](coupled-operators-research.md) now includes a GPU survey, independent CPU spectra, a finite-width control and rotation GIF. None of these establishes scientific novelty or apparatus-level chaos.

The 56-page supplement was retrieved and its model pages inspected, including visual checks of pp. 37 and 40. Downloads remain local under `outputs/`; they are not redistributed under this repository's MIT license. The supplement filename is `41467_2025_66395_MOESM1_ESM.pdf`.

## Continuous temperature matters

Use the mechanical and heat-balance equations, SI Eqs. (6) and (10). In inertial time, write the dimensionless displacement, velocity and thermally preferred displacement as `w,v,z`:

```text
w' = v
v' = -beta_bar*v + z - w
theta*z' = P_bar*H(w0_bar-w) - z
theta = tau_heat/tau_inertial
```

The thermal state remains continuous at a switching event. The source's Eqs. (14–15) display heating from zero and cooling from the fully illuminated limit. Those special exponential responses omit general initial temperatures and do not by themselves specify how temperatures are carried between events. Resetting temperature at every crossing would change the hybrid system. We have not recovered the authors' implementation, so our interpretation is explicitly anchored in their heat balance.

Set `s=t/tau_heat`, `delta=w0_bar/P_bar`, and

```text
U = delta - w/P_bar
V = -theta*v/P_bar
Y = z/P_bar
k = theta^2,  b = beta_bar*theta
```

Then exactly, away from the switching surface,

```text
U' = V
V' = k*(delta-U-Y) - b*V
Y' = H(U) - Y.
```

These are the equations in the [slider relay-limit notebook](lce-slider-relay-limit.md). Unlike the photochemical slider's singular weak-light limit, this photothermal system has a linear heat-loss law directly. It therefore realizes the same ideal mathematical model at finite light power. The algebra audit checks 128 parameter/state samples; it is not a new trajectory survey or a proof of experimental chaos.

## How close is the published reference?

The source's reference has `theta=20`, `beta_bar=0.8`, hence `k=400,b=16`. Our existing chaotic relay candidate has `k=816.0170,b=13.49168,delta=0.0309042`, corresponding to `theta≈28.5660,beta_bar≈0.47230`. Those are changed parameters. No measured apparatus has been shown to reach this regime.

For scale only, retaining the source's 25 mm beam length and deflection coefficient `0.0026/mW` at 100 mW would require a threshold near 0.201 mm for that delta. The source uses a 5 mm single-operator threshold. Supplementary Fig. 35 reports a 2 mm laser spot, and several coupled-operator captions report 2–3 mm spots. However, the main article's Methods / Light excitation describes focusing to approximately **100 micrometres near the baffle edge**. These quantities cannot automatically be equated with the optical switching width; the beam propagation geometry and the location at which a diameter is measured matter. Our earlier reading mentioned only the 2 mm figure, which was incomplete. Neither reported diameter establishes a measured displacement-to-power transition law for our model.

A finite transition width can substantially smooth the switch. Our [earlier smoothing results](lce-slider-switching.md) make this an essential physical question, not a minor correction. They concern a nearby finite-intensity photochemical model and do not establish a switch-width threshold for this apparatus.

The paper describes an onset as a Hopf bifurcation. In an ideal discontinuous switch model a classical smooth Hopf theorem cannot be inferred just from oscillation onset; the actual optical transition law and equilibrium need to be specified.

## Coupled operators: source equations

Normalize each operator's displacement and thermal displacement by its own illuminated limit. For the common mechanical and thermal scales used in SI Eqs. (17–20), time is again `t/tau_heat`:

```text
x_minus' = v_minus
v_minus' = k*(y_minus-x_minus) - b*v_minus
y_minus' = H(x_plus-delta_plus) - y_minus

x_plus' = v_plus
v_plus' = k*(y_plus-x_plus) - b*v_plus
y_plus' = H(delta_minus-x_minus) - y_plus.
```

Thus four parameters remain: `k,b,delta_minus,delta_plus`. Independent beam powers affect the two deltas through normalization. This follows the source's positive and negative optical operators; an arbitrary delay or extra state has not been inserted. The implementation is in `research/audit_light_operators.py`.

For positive `k,b`, each temperature fraction stays in `[0,1]` when initialized there. Each mechanical subsystem is a stable linear oscillator under bounded forcing, so forward trajectories are bounded. This does not prove chaos, exclude long transients, or resolve simultaneous/grazing switching events. The completed [event-driven calculation and explicitly assumed smooth-profile control](coupled-operators-research.md) keep those limits explicit.

## Reproduction and references

```powershell
.venv/Scripts/python.exe research/audit_light_operators.py
```

Output: `results/light-operators/equivalence-audit.json`. This script checks the single-operator algebra only; the separate coupled-system notebook gives its survey commands and results.

1. Guo, H., Li, K., Yang, J., Li, D., Liu, F., and Zeng, H. (2025). Light-mediated communication in responsive materials ranging from individual self-oscillators to feedback-driven network. *Nature Communications* **16**, 11315. [doi:10.1038/s41467-025-66395-3](https://doi.org/10.1038/s41467-025-66395-3). Supplementary Notes 2.1–2.3, pp. 34–47. [Open primary article and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC12722257/).
2. Ge, D., Hong, Q., Liu, X., and Liang, H. (2024). See the primary model and complete citation in the [LCE slider notebook](lce-slider-research.md). [doi:10.3390/polym16233298](https://doi.org/10.3390/polym16233298).
3. Older relay-system comparisons and their access limits are recorded in the [relay-limit audit](lce-slider-relay-limit.md). Exact equivalence here narrows the novelty claim rather than establishing priority.
