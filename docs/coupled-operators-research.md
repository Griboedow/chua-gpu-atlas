# Chaos candidates in two light-mediated material operators

**Status: independently checked numerical chaos candidates in a published-model interpretation; scientific novelty and experimental realization remain unconfirmed.** The single operator is equivalent to our previous relay limit. The coupled six-state system is a distinct investigation, and its chaotic examples survive an independent numerical ODE solver. One example also survives an explicitly assumed smooth optical transition.

![Coupled operator evidence](assets/coupled_operators_1026.png)

![Camera rotation of a physical projection](assets/coupled_operators_1026.gif)

The three-dimensional picture is a **projection of a six-dimensional flow**, not its complete phase portrait. It shows the two material displacements and the difference between their thermal states. The folded bands and multiple return-section branches are interesting numerical structure, but they do not establish a new topological type of attractor.

## What physical system is being modeled?

Guo et al. [1] connect two responsive cantilevers with constant laser beams. One baffle opens the beam reaching its partner, while the other closes the returning beam. A material's temperature and deformation therefore affect its partner's illumination. The paper demonstrates self-oscillation and communication between materials; it does not supply a chaos analysis in the sections we inspected.

Our [source audit](light-operators-audit.md) records the complete supplement and a consequential ambiguity. We use the continuous heat balance in SI Eq. (10), mechanical Eqs. (17–18), and the illumination conditions of Eqs. (19–20). Temperature remains continuous at switches. The source's displayed exponential step responses do not specify arbitrary initial temperatures; an implementation that resets heat at every switch would be different. We have not recovered the authors' implementation.

Each displacement and thermally preferred displacement is divided by that operator's own fully illuminated deflection. Time is divided by the heat-relaxation time. The model is

```text
x_minus' = v_minus
v_minus' = k*(y_minus-x_minus) - b*v_minus
y_minus' = H(x_plus-delta_plus) - y_minus

x_plus' = v_plus
v_plus' = k*(y_plus-x_plus) - b*v_plus
y_plus' = H(delta_minus-x_minus) - y_plus.
```

Here `k=theta^2`, `b=beta_bar*theta`, and `theta=tau_heat/tau_inertial`. The two deltas are the baffle thresholds divided by the individual illuminated deflections. No external time delay, noise, or extra dynamic state has been added.

## Structural checks

For `k,b>0`, the thermal fractions remain in `[0,1]`. The mechanical subsystems are stable linear filters driven by those bounded fractions. Consequently all forward trajectories from finite states are bounded. This is a model-level boundedness argument, not a statement about the validity of arbitrary deflections in a real cantilever.

With `P(D)=(D+1)(D^2+bD+k)`, the exact pair of jerk equations is

```text
P(D) x_minus = k H(x_plus-delta_plus)
P(D) x_plus  = k H(delta_minus-x_minus).
```

Every open affine branch has the same stable linear part. The flow divergence is `-2*(b+1)`. At a transverse switching event the saltation matrix has determinant one, because the changed heating derivative is tangential to the switching surface. The full Lyapunov spectrum must therefore sum to that divergence, providing a useful numerical control.

For interior thresholds there is no equilibrium strictly inside an affine branch: the binary steady states would require simultaneously `x_minus=x_plus` and `x_plus=1-x_minus`. **This is not an equilibrium-free-system claim.** A Filippov convexification permits an equilibrium on the intersection of the switching surfaces. Our event calculation follows transverse crossings and explicitly rejects unresolved grazing and simultaneous events.

The transformation

```text
T(xm,vm,ym,xp,vp,yp) = (xp,vp,yp,1-xm,-vm,1-ym)
(delta_minus,delta_plus) -> (delta_plus,1-delta_minus)
```

maps parameter-related systems exactly. `T^4=identity`. It is a symmetry of one fixed system when both thresholds equal `1/2`, giving a fourfold symmetry in six-dimensional state space. The two plotted candidates do not have those symmetric parameters. Random algebra checks of the transformation and jerk reduction are saved in `coupled-structure.json`.

## Numerical search and independent checks

The FP64 OpenCL code uses exact affine segments and locates the first crossing on a finite bracketing mesh. It then applies the saltation matrix to the flow tangent. This avoids integrating across a discontinuity with an ordinary unsplit RK step. Mesh refinement is still necessary: sign bracketing alone cannot rigorously exclude two roots inside one bracket.

The scout contains 2,048 parameter tuples with two independent starts each. Its ranges are `k=10..10,000`, `beta_bar=b/sqrt(k)=10^-1.3..10^0.2`, and each delta between `0.005` and `0.995`, using a scrambled Sobol design. The source's Fig. 37 reference replaces the first tuple. This broad search changes physical ratios substantially; it is not a scan confined to measured apparatus parameters.

| Stage | Settings | Result |
| --- | --- | --- |
| Scout, 4,096 trajectories | 1,000 discarded + 4,000 observed switching events; 64 brackets per mechanical period | All completed; 436 preliminary positive candidates; 34.8 s on GTX 1080 |
| Every preliminary candidate | Additional 4,000 discarded + 16,000 observed events; 256 brackets per period | All 436 completed; 422 retained, representing 213 distinct parameter tuples; 233.5 s |
| Reference A, scout index 1026 | Independent DOP853, 200 discarded + 800 observed thermal-time units | Positive leading exponent in all four windows; full six-exponent spectrum |
| Reference B, scout index 3054 | Independent DOP853, 100 discarded + 400 observed thermal-time units | Positive leading exponent in all four windows; repeated with tenfold tighter tolerance and half the QR interval |

The screening criterion is mean leading exponent above `0.005` and all four window means above `0.002`, after successful integration. Retained fractions depend on the chosen sampling distribution and integration lengths; they are not basin probabilities or confidence levels.

An additional 256-trajectory local check uses 32 new random starts per reference and perturbation size, with the same long event counts and 256-point bracketing. All trajectories complete. At reference A, retained counts are `32/32,31/32,23/32,14/32` for log-parameter perturbations `0,0.001,0.01,0.03`; reference B retains `32/32` in all four groups. Every parameter is multiplied independently by `exp(uniform(-j,j))`. This shows local numerical persistence and parameter sensitivity, not a proof of structural stability.

Implementation validation gives an affine-segment discrepancy below `6e-14` against DOP853, fixed-time saltation-Jacobian discrepancy below `9e-10` against finite differences, and a CPU/GPU twelve-event state discrepancy below `4e-14`. The sum of each independent CPU spectrum agrees with the divergence to about `3e-12` or better.

At the source Fig. 37 parameters, the same solver instead returns a strongly stable four-switch periodic orbit: period `1.8608154588` in thermal-time units and largest Poincare multiplier magnitude `0.0101660`. Independent DOP853 closes that orbit within `8.3e-14`. Tiny remaining multipliers have limited relative precision. This is a useful regular control and does not claim exact reproduction of the unavailable author implementation.

| Quantity | Reference A | Reference B |
| --- | --- | --- |
| `k` | 74.1567393414624 | 8608.828902752228 |
| `b` | 1.784180421282182 | 4.7177230976207 |
| `delta_minus` | 0.4139517339970917 | 0.8733635635208339 |
| `delta_plus` | 0.25683679452165964 | 0.8404887881875038 |
| Refined GPU leading exponent, first start | 0.22182 | 4.42050 |
| Independent CPU spectrum | `(0.22650,-0.00373,-0.56073,-1.27571,-1.80234,-2.15235)` | `(4.29812,-0.00402,-0.59409,-4.20581,-4.73502,-6.19464)` |
| Tighter CPU leading exponent | See smooth control below | 4.30254 |

The nearly zero second exponent is consistent with the neutral direction of an autonomous flow. It is not evidence of hyperchaos. Positive finite-time exponents and long bounded trajectories remain numerical evidence, not a proof that the asymptotic object is a strange attractor rather than an extremely long transient.

![Higher-frequency physical projection and section](assets/coupled_operators_3054.png)

Reference A has a structured projected return section with several curved branches. Reference B has a broader projected section and much faster expansion. Neither appearance alone establishes a new attractor class. Plotting is densified with exact branch flow between CPU-verified event/QR knots; the maximum branch-end mismatch is below `1.1e-11`, so fast mechanical motion is not represented by sparsely connected QR samples.

## Finite optical width: a stronger physical test

Replace `H(z)` with `S(z)=1/(1+exp(-z/w))` in both illumination laws. This is an **explicit modified model**, not a measured beam profile. The published source gives both approximately 100 micrometres at the baffle edge in Methods and 2–3 mm spot sizes in several supplemental captions. Those reports must not be converted directly to a switching width without the optical geometry.

At reference A, `w=0.005` retains a positive leading exponent under independent DOP853 integration:

```text
(0.22207, 0.00228, -0.46622, -1.38822, -1.79217, -2.14610).
```

The four leading window estimates are `0.14054,0.26572,0.23281,0.24921`. Thus this particular observation does not depend on an infinitely sharp relay. For a logistic profile, the 10–90% transition spans `2*ln(9)*w≈0.02197` of the individual illuminated deflection. This is a dimensionless sensitivity test, not apparatus calibration.

A second run discards 400 thermal-time units and observes 1,600, using `rtol=1e-12` and half the QR interval. It gives leading exponent `0.217994`, with window values `0.21321,0.22287,0.19936,0.23653`; its second exponent is `-0.001236`. The spectrum sum agrees with the analytic divergence within `1.3e-12`.

The smooth model has a unique equilibrium: composing an increasing illumination response with the decreasing response gives a strictly decreasing fixed-point equation. Let `a_minus,a_plus` be the positive magnitudes of their slopes there. Its characteristic equation is

```text
[(lambda+1)(lambda^2+b*lambda+k)]^2 + k^2*a_minus*a_plus = 0.
```

An imaginary pair occurs when

```text
omega^2 = k/(b+1)
a_minus*a_plus = b^2*(k+b+1)^2 / [k*(b+1)^3].
```

For reference A this gives a linear Hopf candidate at `w≈0.0632608463`, `omega≈5.16091404`. It does not locate the chaos boundary or prove the nonlinear Hopf bifurcation type. The algebra and direct eigenvalue check are reproducible.

## Prior art and what remains open

General relay chaos is old. Our [relay literature audit](lce-slider-relay-limit.md) already covers several established examples and explicitly lists inaccessible originals. The coupled jerk representation above makes further control-theory equivalence searches necessary. A targeted search for the source title/DOI plus chaos did not identify the same result, which is insufficient to establish novelty.

An OpenAlex snapshot for source work `W4416430930` on 29 October 2026 returned five citing works. The access levels are recorded rather than treating an unexamined citation as cleared:

| Citing work | Inspection and relevance |
| --- | --- |
| Peng et al., [10.1016/j.chaos.2026.118687](https://doi.org/10.1016/j.chaos.2026.118687) | Publisher highlights: delayed reaction–diffusion, hydrogel swelling and periodic pulsation. Full equations not retrieved; separate candidate, not cleared for absence of chaos. |
| [10.1016/j.ccr.2026.217958](https://doi.org/10.1016/j.ccr.2026.217958) | Citation metadata only: review on fueling dissipative self-assembly. Full text not audited. |
| Cheng et al., [10.1007/s10118-026-3646-1](https://doi.org/10.1007/s10118-026-3646-1) | Publisher abstract: pre-twisted hydrogel rings, figure-of-eight reconfiguration and irregular tumbling. No faithful dynamic model retrieved yet; irregular tumbling is not itself a chaos diagnosis. |
| Boron et al., [10.1002/advs.76810](https://doi.org/10.1002/advs.76810) | Primary main text: tendon-driven origami tentacles and contact-mediated entanglement, with Cosserat-rod simulation. Its finite grasping operation is not evidence of an autonomous sustained chaotic attractor. |
| Yuan et al., [10.1016/j.chaos.2026.118988](https://doi.org/10.1016/j.chaos.2026.118988) | Indexed bibliographic record/abstract only; publisher full text inaccessible in this pass. Light-powered dominoes and impact-mediated cycling remain a possible hybrid-model target. |

Other retrieved relay-search results concerned delayed neural relays or a third intermediary oscillator between already chaotic systems. The shared word “relay” does not establish equation equivalence. A comprehensive classification of two stable third-order filters with these cross-coupled switching laws remains unfinished.

## Reproduction

From the repository root, with NumPy, SciPy, PyOpenCL, Matplotlib and Pillow installed:

```powershell
python research/coupled_operators.py
python research/coupled_operators.py --scout --transient 1000 --events 4000 --subdivisions 64
python research/coupled_operators.py --refine results/light-operators/coupled-scout.npz --transient 4000 --events 16000 --subdivisions 256
python research/check_coupled_operators.py --index 1026 --transient 200 --observation 800 --output results/light-operators/coupled-cpu-1026
python research/check_coupled_operators.py --index 3054 --transient 100 --observation 400 --output results/light-operators/coupled-cpu-3054
python research/check_coupled_operators.py --index 3054 --rtol 1e-12 --qr-interval .025 --output results/light-operators/coupled-cpu-3054-tight
python research/check_coupled_operators.py --index 1026 --width .005 --transient 200 --observation 800 --output results/light-operators/coupled-smooth-1026-w005
python research/check_coupled_operators.py --index 1026 --width .005 --transient 400 --observation 1600 --rtol 1e-12 --qr-interval .025 --output results/light-operators/coupled-smooth-1026-w005-tight
python research/audit_coupled_operators.py
python research/check_coupled_operators.py --periodic-control --rtol 1e-12 --qr-interval .02 --output results/light-operators/coupled-periodic-control-cpu
python research/audit_coupled_operators.py --robustness
python research/render_coupled_operators.py --gif
```

The saved NPZ files include all surveyed parameters, starts, statuses, exponents/windows and final states; independent CPU traces include physical coordinates and Poincare sections. CPU verification requires no GPU. The analytic GPU implementation currently supports underdamped mechanics only.

## References

1. Guo, H., Li, K., Yang, J., Li, D., Liu, F., and Zeng, H. (2025). Light-mediated communication in responsive materials ranging from individual self-oscillators to feedback-driven network. *Nature Communications* **16**, 11315. [doi:10.1038/s41467-025-66395-3](https://doi.org/10.1038/s41467-025-66395-3). Supplementary Notes 2.1–2.3, pp. 34–47; main-text Methods / Light excitation.
2. di Bernardo, M., Johansson, K. H., and Vasca, F. (2001). Self-oscillations and sliding in relay feedback systems: symmetry and bifurcations. *International Journal of Bifurcation and Chaos* **11**, 1121–1140. [doi:10.1142/S0218127401002584](https://doi.org/10.1142/S0218127401002584). A primary relay comparison, not the same two-operator transfer structure.
3. Mallet-Paret, J., and Smith, H. L. (1990). The Poincare–Bendixson theorem for monotone cyclic feedback systems. *Journal of Dynamics and Differential Equations* **2**, 367–421. [doi:10.1007/BF01054041](https://doi.org/10.1007/BF01054041). Author/institutional abstract inspected. A possible route for an overdamped smooth-model analysis; no theorem is applied to the discontinuous underdamped candidates here.
4. Letellier, C., and Aguirre, L. A. (2012). Required criteria for recognizing new types of chaos: Application to the “cord” attractor. *Physical Review E* **85**, 036204. [doi:10.1103/PhysRevE.85.036204](https://doi.org/10.1103/PhysRevE.85.036204). Publisher abstract inspected; our figures do not constitute that topological analysis.
