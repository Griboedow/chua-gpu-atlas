# A concrete thermofluidic candidate: the NIFTE heat engine

Screening date: 30 October 2026, 00:20 UTC. This initial note provides the published equations and an independent boundedness argument. No trajectory integration was performed in this screening block. **The subsequent [implementation and numerical study](nifte-research.md) now reproduces the source cycles and records a large search without a validated chaos candidate.** The original transcription and source-access evidence below remain available.

The non-inertive-feedback thermofluidic engine (NIFTE) converts a constant hot/cold temperature difference into fluid motion through evaporation and condensation. Markides et al. [1] introduced a smooth saturating temperature profile. Wang et al. [2] give explicit equations, thermal-loss corrections and coexisting stable/unstable cycles. The latter study is the precise implementation target. Its calibrated thermal-loss model has **nine states**, not six: four additional RC states supplement the five original states.

## Exact equations and sign conventions

Use the dimensionless state `(p,d,q,u,v)` for the source's adiabatic vapor pressure, displacer hydrostatic pressure, power-cylinder hydrostatic pressure, feedback flow and power-cylinder flow. Hats are omitted below; time is `t_hat=t/tau`. These pressures and flows are signed deviations, not positive absolute pressures. In particular, reversing the signs of both flow variables without changing the pressure equations would change the model.

This is a compact transcription of [2], Eqs. 7a–7e and 8. Define

```math
a=L_d+L_f,\qquad D=L_l a+L_p(a+L_l),\qquad
c=\frac{\tau U^0}{P^0},\qquad b=\frac{\tau P^0}{U^0},
\qquad f(d)=F\tanh(Kd).
```

Then the five-state NTP model is

```math
\dot p=\frac{\tau[f(d)-p]}{R_{th}C_{ad}}+
\frac{c(u+v)}{C_{ad}},\qquad
\dot d=\frac{cu}{C_d},\qquad
\dot q=\frac{cv}{C_p},
```

```math
\dot u=\frac{b[L_lq-L_pp-(L_p+L_l)d]
-\tau[L_lR_f+L_p(R_f+R_l)]u-\tau L_pR_lv}{D},
```

```math
\dot v=\frac{b[L_ld-ap-(a+L_l)q]
+\tau(L_lR_f-aR_l)u-\tau aR_lv}{D}.
```

The dots in these three blocks mean derivatives with respect to dimensionless time. In particular, the last term in `u_dot` is **negative**. The `u` coefficient in `v_dot` can have either sign; its two physical contributions must not be combined as a sum of magnitudes.

For the NTP-TL model, introduce `z_i`, the dimensionless pressure of loss capacitor `i`. Retain the other four equations and replace/add

```math
\dot p=\frac{\tau[f(d)-p]}{R_{th}C_{ad}}+
\frac{c(u+v)}{C_{ad}}
-\frac{\tau(p-\sum_{i=1}^n z_i)}{R_{tl,0}C_{ad}},
```

```math
\dot z_i=\frac{\tau(p-\sum_{j=1}^n z_j)}{R_{tl,0}C_{tl,i}}
-\frac{\tau z_i}{R_{tl,i}C_{tl,i}},\qquad i=1,\ldots,n.
```

The published fit uses `n=4`. A one-capacitor, six-state model would require a new fit and validation; simply discarding three capacitors is not a reproduction. The main closure is smooth `tanh`, with no events, impacts or switching. Wang's later thesis [3] discusses a separate switched-linear approximation; that is not the model specified here.

## Parameters and published controls

The following numerical data are from [2], Table 2, visually checked in the publisher PDF. Resistances have units `kg m^-4 s^-1`, inertances `kg m^-4`, and capacitances `m^4 s^2 kg^-1`.

| Parameter | Nominal value | Unit |
| --- | ---: | --- |
| `F` (`P_th_max` in source) | 1.37 | dimensionless |
| `K` (`K_th`) | 330 | dimensionless |
| `R_f` | 2.13e6 | resistance |
| `R_l` | 4.08e6 | resistance |
| `R_th` | 5.02e8 | resistance |
| `L_d` | 1.80e5 | inertance |
| `L_f` | 4.74e6 | inertance |
| `L_l` | 1.27e7 | inertance |
| `L_p` | 3.77e5 | inertance |
| `C_ad` | 1.76e-9 | capacitance |
| `C_d` | 7.35e-8 | capacitance |
| `C_p` | 7.43e-8 | capacitance |
| `R_tl,0` | 3.15e7 | resistance |
| `R_tl,1 ... R_tl,4` | 8.82e6, 1.72e7, 3.83e7, 2.81e8 | resistance |
| `C_tl,1 ... C_tl,4` | 8.61e-10, 5.95e-9, 1.79e-8, 2.66e-8 | capacitance |
| `P^0` | 1.013e5 | Pa |
| `U^0` | 8.00e-4 | m^3 s^-1 |
| `tau` | 5.00 | s |

For the control in [2], Table 3, change only `L_f` to `4.74e5`. The five-state stable cycle has period **1.31 s in the table** (1.30 s in the surrounding prose) and multipliers `1`, `0.1028`, `0.3093`, `-0.4704 +/- 0.2959i`. The unstable cycle has period `0.88 s` and multipliers `1`, `0.0348`, `0.7958`, `-0.3848 +/- 1.2078i`. With the four thermal-loss capacitors, the reported stable cycle has period `1.69 s`. These are published reference values, not our calculations. No complete cycle initial condition was recovered from the inspected tables.

Two transcription traps deserve attention. Table 3 and its preceding paragraph attach an extra `s^-1` to `L_f`, inconsistent with Table 2 and the dimensional flow equation; use the inertance units above. Also, the printed RC impedance sum in Eq. 6 starts at zero, while the state equations, fit and table distinguish a series `R_tl,0` from four parallel pairs indexed 1–4. The implementation target is the explicit Eqs. 7a' and 7f with those four pairs.

## Our analytical controls

The following calculations are direct consequences of the equations above; they are not claimed as new theorems or unpublished results.

### Passive part and global mathematical boundedness

Return temporarily to dimensional deviations `(P,D_p,Q,U,V)` and physical time. Write

```math
M=\begin{pmatrix}a+L_l&L_l\\L_l&L_p+L_l\end{pmatrix},
\qquad \mathcal F(D_p)=P^0 F\tanh(KD_p/P^0).
```

The two flow equations are exactly equivalent to

```math
M\begin{pmatrix}\dot U\\\dot V\end{pmatrix}
=-\begin{pmatrix}P+D_p\\P+Q\end{pmatrix}
-\begin{pmatrix}R_f+R_l&R_l\\R_l&R_l\end{pmatrix}
\begin{pmatrix}U\\V\end{pmatrix}.
```

The mass matrix is positive definite for positive inertances. Define

```math
E=\frac12\left(C_{ad}P^2+C_dD_p^2+C_pQ^2
+aU^2+L_pV^2+L_l(U+V)^2\right).
```

Differentiating and cancelling the pressure-flow terms gives

```math
\dot E=\frac{P\mathcal F(D_p)-P^2}{R_{th}}
-R_fU^2-R_l(U+V)^2.
```

For the linear subsystem with input `mathcal F=0`, this is nonpositive. Its zero-dissipation set has `P=U=V=0`. A trajectory remaining in that set must also have `D_p=Q=0`, by the two flow equations. The positive quadratic energy and this invariant-set argument imply asymptotic stability of the passive linear subsystem; finite-dimensional linearity then gives exponential stability.

The nonlinear thermal input is uniformly bounded by `P^0 F`. Variation of constants for an exponentially stable linear system therefore bounds every solution and gives an ultimate bound independent of its starting state. Smoothness and at-most-linear growth ensure global existence. This establishes mathematical boundedness for positive component parameters, without claiming a chaotic attractor.

For the loss model, add `sum(C_tl,i Z_i^2)/2` to the energy. The extra dissipation is

```math
-\frac{(P-\sum_i Z_i)^2}{R_{tl,0}}
-\sum_i\frac{Z_i^2}{R_{tl,i}}.
```

The same argument applies for any finite positive RC network, including the published four-pair fit.

An algebraic implementation check at 64 seeded random states gives a maximum scaled energy-identity residual of `7.3e-16`; the explicit and mass-matrix flow derivatives agree within `3.3e-19 m^3 s^-2`. This is not a trajectory or attractor calculation.

### Equilibrium, symmetry and a spectrum check

At an NTP equilibrium, the flow equations give `U=V=0` and `D_p=Q=-P`. The remaining equation becomes `mathcal F(D_p)+D_p=0`, so the origin is the unique equilibrium for positive `F,K`. The loss network adds a positive linear coefficient to `D_p` and preserves uniqueness. All terms are odd, giving the symmetry `X -> -X`.

In dimensional time the five-state divergence is the constant

```math
\nabla\cdot f=-\frac{1}{R_{th}C_{ad}}
-\frac{L_lR_f+L_p(R_f+R_l)+aR_l}{D}<0.
```

For `n` loss states add

```math
-\frac{1}{R_{tl,0}C_{ad}}
-\sum_{i=1}^n\frac{1/R_{tl,0}+1/R_{tl,i}}{C_{tl,i}}.
```

Multiply by `tau` when comparing against exponents per dimensionless time. This identity is useful for checking a full variational implementation. Volume contraction and boundedness do not imply chaos.

No planar reduction or elementary monotone-cyclic-feedback exclusion was established. The given equations have several bidirectional couplings rather than one incoming link per state. In particular, the `p,u` reciprocal derivatives have opposite signs, whose product survives diagonal sign changes. A direct sign flip cannot turn this into a system with all off-diagonal derivatives of the same sign. This does not rule out a more sophisticated coordinate change, cone argument or exclusion theorem.

## Physical domain and interpretation

The equations define a smooth mathematical flow on all of real state space. Physical interpretation is narrower: the reconstructed liquid interfaces must remain inside their cylinders, the vapor region must retain positive volume, and absolute pressure must remain positive. The relation `D_p=rho_wf*g*y` reconstructs the displacer level. Linear inertance/compliance/resistance values and fixed thermophysical properties are modeling assumptions. A bounded mathematical solution need not satisfy these operating assumptions.

The RC fit describes thermal-loss impedance over `0.1–100 rad s^-1`, not unrestricted microscopic dynamics. Do not label behavior produced at frequencies far outside its calibration range as a prediction for the apparatus. Likewise, changing circuit-analogue parameters independently can represent a mathematically valid model but may not correspond to an achievable geometry. A physical parameter map and explicit interface limits should precede claims about an experimental attractor.

## Literature check and priority limits

The full primary texts of [1] and [2] were inspected, including the explicit equations and conclusions. The 2018 primary precursor [4] already reports stable/unstable cyclic-state multiplicity. These features are prior work. Exact-name, exact-identifier and chaos/attractor/Lyapunov searches did not establish an exact-model chaos report in the inspected material.

The OpenAlex snapshot returned **3 citing works for [2] and 40 for [1]**. All returned titles/identifiers were screened; this is not a full-text audit of forty papers. The direct 2020 experimental descendant [5] concerns zeotropic working-fluid performance; its publisher preview was inspected. Two other direct citations are reviews, retained as search leads only. A related 2024 experiment [6] explicitly reports quasiperiodicity and mode transitions in a regenerator-equipped thermofluidic oscillator. Its abstract is not a demonstration of chaos in the NTP equations, but it prevents treating complex thermofluidic oscillations as a new general idea.

Wang's 2020 thesis [3] was checked through its institutional record. Its full PDF returned HTTP 429, so the thesis equations and all later developments remain incompletely checked. A single citation index and search silence cannot establish priority. No equivalence classification against general saturating-feedback electrical circuits has been completed.

The additional field-emission nanowire lead [7] was not promoted to a second implementation target: its full articulated model requires geometry-dependent electrostatic quantities obtained from a field solve. Replacing them with convenient guessed functions would change the model. This decision does not exclude chaos in that device class.

## Next implementation block

1. Implement the literal five-state equations and the equivalent dimensional mass-matrix form independently; verify their equality and energy identity.
2. Reproduce the stable cycle at the published reduced-`L_f` setting. Use shooting/continuation for the unstable cycle; a generic forward trajectory will not reproduce it. State any chosen initial conditions as our choices.
3. Check the full Lyapunov sum against the constant divergence and the trivial Floquet multiplier against one.
4. Add all four published thermal-loss states and reproduce the `1.69 s` control before a broad GPU survey. Explore the already published coexistence/instability neighborhoods first, while tracking physical-domain and approximation limits.

The [machine-readable screening record](../results/candidate-audits/physical-candidate-screening.json) contains the equations, numerical parameters, reference spectra, source access record and citation-index snapshot. Downloaded source PDFs and page renders are kept only in ignored local outputs.

## References

1. Markides, C. N., Osuolale, A., Solanki, R., and Stan, G.-B. V. (2013). Nonlinear heat transfer processes in a two-phase thermofluidic oscillator. *Applied Energy* **104**, 958–977. [DOI](https://doi.org/10.1016/j.apenergy.2012.11.056); [author-hosted full paper](https://gstan.bg-research.cc.ic.ac.uk/Markides_Stan_App_Energy_2013.pdf).
2. Wang, Y., Markides, C. N., and Chachuat, B. (2019). Optimization-based investigations of a two-phase thermofluidic oscillator for low-grade heat conversion. *BMC Chemical Engineering* **1**, 20. [Full primary article](https://doi.org/10.1186/s42480-019-0019-3); [publisher PDF](https://link.springer.com/content/pdf/10.1186/s42480-019-0019-3.pdf). Equations 7–8 and Tables 2–3 checked against PDF pages 6, 8 and 10; notation above is rearranged. The article is [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
3. Wang, Y. (2020). *Two-phase thermofluidic engines for low-grade heat recovery: system analysis and supporting algorithms*. Imperial College London doctoral thesis. [Institutional record](https://doi.org/10.25560/86023). Record/abstract inspected; full PDF unavailable in this pass.
4. Wang, Y., Markides, C. N., and Chachuat, B. (2018). Optimization-based investigations of a thermofluidic engine for low-grade heat recovery. *IFAC-PapersOnLine* **51**(18), 690–695. [DOI](https://doi.org/10.1016/j.ifacol.2018.09.343); [primary proceedings preprint](https://skoge.folk.ntnu.no/prost/proceedings/ADCHEM-2018/media/files/0136.pdf).
5. Tan, J., Luo, J., Huang, J., Wei, J., and Jin, T. (2020). A closed two-phase thermofluidic oscillator with zeotropic mixtures for low-grade heat recovery. *Energy* **211**, 118691. [DOI](https://doi.org/10.1016/j.energy.2020.118691). Publisher preview inspected, not full equations.
6. Zhou, Q., Luo, J., Chen, Y., Sheng, S., and Jin, T. (2024). Characterizing the mode transition features of a two-phase thermofluidic oscillator. *Applied Thermal Engineering* **240**, 122225. [DOI](https://doi.org/10.1016/j.applthermaleng.2023.122225). Publisher abstract inspected.
7. Lazarus, A., de Langre, E., Manneville, P., Vincent, P., Perisanu, S., Ayari, A., and Purcell, S. (2010). Statics and dynamics of a nanowire in field emission. *International Journal of Mechanical Sciences* **52**, 1396–1406. [DOI](https://doi.org/10.1016/j.ijmecsci.2010.01.011); [author-hosted paper](https://ilm.univ-lyon1.fr/images/ILM/Equipe_17/articles/articlePNEC/lazarusmechanics.pdf). Model and electrostatic-closure sections inspected; comprehensive priority search not completed.
