# NIFTE: physical operating limits and prior-model audit

Audit date: 30 October 2026. **The published NTP equations are a concrete mathematical target, but the inspected sources do not supply a complete numerical operating envelope for a realizable apparatus.** The checks below provide state conversions and a geometry map. They do not certify an experimental chaotic regime or scientific priority. No trajectories were integrated for this audit.

The equations and numerical controls are recorded in [the implementation screening note](physical-candidate-screening.md). This follow-up uses the primary 2011, 2013 and 2019 papers [1–3], with a separate check of related literature. The [machine-readable record](../results/nifte/physical-priority-audit.json) separates published values, our derived quantities and unresolved assumptions.

## Pressures, interface positions and missing bounds

Keep the dimensionless state convention `(p,d,q,u,v)` from the screening note. These pressure coordinates are **deviations**, not absolute pressures. With the published reference pressure and liquid density, the reconstruction is

```math
P_{\mathrm{abs}}=P^0(1+p),\qquad
y_d=\frac{P^0d}{\rho_{\mathrm{wf}}g},\qquad
y_p=\frac{P^0q}{\rho_{\mathrm{wf}}g}.
```

For `P0=101300 Pa`, `rho_wf=621 kg/m^3` and `g=9.81 m/s^2`, one unit of `d` or `q` means **16.6283 m** of interface displacement. Thus a seemingly small coordinate such as `d=0.01` represents about 166 mm. The nominal fluid is n-pentane; water is the pumped load fluid. The nominal equilibrium temperature is 309 K [2].

`p > -1` is a necessary absolute-pressure condition under this reference convention. It is a weak test: pressure positivity alone does not validate constant density, linear compliance or phase-equilibrium approximations. The additional thermal-loss coordinates are auxiliary RC-network pressures, not separate absolute chamber pressures; imposing `z_i > -1` has no established physical justification.

The inspected texts give the following geometric facts [2]:

| Quantity | Reported value | What it does **not** establish |
| --- | ---: | --- |
| Hot exchanger height | 90 mm | Upper mechanical stroke limit |
| Cold exchanger height | 60 mm | Lower mechanical stroke limit |
| Exchanger separation | 3–4 mm | Total cylinder travel |
| Nominal displacer/power lengths in Table 1 | 0.26 / 0.55 m | The two signed clearances from equilibrium |
| Nominal vapor volume | 1.94e-4 m^3 | A separate stop/event law |

Appendix A of [2] calibrates the smooth temperature profile using approximately 95% saturation over a 90 mm hot exchanger. This is a heat-transfer scale, **not a hard state boundary**. The ODE continues smoothly beyond it. Cylinder end locations, liquid-inventory constraints, overflow conditions, tank head and gas ingestion were not recovered as a closed set of state inequalities from [1–3]. Adding impacts or clipping the interface at the exchanger edge would change the published model.

If both cylindrical free surfaces account for the varying vapor volume, constant cross-sections imply the geometric reconstruction

```math
V_{\mathrm{geom}}=V^0-A_dy_d-A_py_p
=V^0-P^0(C_dd+C_pq).
```

This is our **conditional geometry diagnostic**, not an additional equation published in the five-state model. It needs confirmation of the chamber layout and equilibrium fill levels. Likewise, a bottom clearance would require `H_d0+y_d>0`, with the actual equilibrium liquid height `H_d0`. None of these missing offsets should silently be replaced by an arbitrary symmetric interval. Positive volume is necessary for a consistent reconstruction, but this formula alone does not certify the apparatus geometry.

## Which component variations correspond to geometry?

Equation 1 of [3] maps fixed-fluid cylindrical components to the electrical analogy:

```math
R=\frac{128\mu l}{\pi d_{\mathrm{tube}}^4}
=\frac{8\pi\mu l}{A^2},\qquad
L=\frac{\rho l}{A},\qquad A=\frac{\pi d_{\mathrm{tube}}^2}{4}.
```

Here the feedback branch uses working-fluid properties, while the load uses pumped-fluid properties. Inverting these relations gives

```math
A=\frac{8\pi\mu L}{\rho R},\qquad
l=\frac{8\pi\mu L^2}{\rho^2R}.
```

Consequently, `R` and `L` are not independent knobs of one *fixed* pipe. They can be selected jointly by changing two geometric variables. If their ratios to nominal are `r` and `ell`, the corresponding diameter and length ratios are

```math
\frac{d_{\mathrm{tube}}}{d_{\mathrm{tube},0}}
=\sqrt{\frac{\ell}{r}},\qquad
\frac{l}{l_0}=\frac{\ell^2}{r}.
```

For example, the published control with feedback inertance reduced tenfold and resistance unchanged [3, Table 3] maps to a **1.581 mm diameter, 1.499 mm long** uniform pipe, using the stated fluid properties. The nominal pair maps to 5.000 mm and 149.862 mm. The reduced length lies outside the `0.05–0.5 m` feedback-length range examined in [2]. This does not invalidate the article's mathematical cycle control. It does mean that changing `Lf` alone should not be described as merely shortening the nominal tube tenfold. A separate valve could alter the resistance map, but its required calibration has not been supplied here.

For either hydrostatic cylinder, the inverse map is

```math
A=\rho g C,\qquad l_{\mathrm{eff}}=gLC.
```

At nominal component values, the inferred displacer/power areas are `4.4776e-4 / 4.5264e-4 m^2`; the effective liquid-column lengths are **0.12979 / 0.27479 m**. Inserting the Table 1 lengths `0.26 / 0.55 m` from [2] directly into its inertance formula instead gives `3.604e5 / 7.556e5 kg m^-4`, approximately twice the nominal inertances `1.80e5 / 3.77e5`. The earlier analysis [1] explicitly interprets the inertial lengths as time-averaged liquid-column lengths. A half-filled-cylinder interpretation is plausible, but the factor-of-two convention was not explicitly reconciled in the inspected tables. Preserve the published RLC values when reproducing cycles; keep this geometry ambiguity visible.

Changing cylinder diameter also changes the vapor-volume contribution unless the top chamber or fill is adjusted. Thus keeping `Cad` fixed while varying `Cp` or `Cd` is a compensating design assumption, not an automatic property of a modified tube. At fixed pressure and adiabatic index,

```math
V^0=\gamma P^0 C_{ad}.
```

The nominal values imply `V0=1.94334e-4 m^3`, consistent with the reported rounded volume.

The saturator parameters have physical interpretations [2,3]. In our notation, with `s_T=(dT/dP)_sat`,

```math
\Delta T=2P^0s_TF,\qquad
\beta=\frac{K\rho g}{P^0}.
```

At `F=1.37`, `K=330`, the implied difference is 80.215 K and the inverse spatial scale is `19.846 m^-1`. Independent changes of `F` and `K` change the imposed temperature profile. Changing fluid or mean pressure also changes the density, viscosity, saturation slope and thermal resistance, so those properties cannot all remain physically fixed by default.

The four-pair thermal-loss fit in [3] is calibrated over `0.1–100 rad/s`. A redesign affecting glass dimensions or heat transfer generally needs a new loss-impedance fit. The source uses linear laminar drag and lumped inertia [1,3]. Report Reynolds numbers and displacement amplitudes as approximation diagnostics; do not turn a conventional flow threshold into an exact theorem about this oscillatory apparatus. For the nominal feedback pipe, `Re_f ≈ 580317 |u|` with the adopted flow scale.

## An exact circuit-family comparison

The model already has the standard **Lur'e feedback structure**: a linear dynamical network driven by one static scalar saturation. The use of `tanh` is therefore not a new chaos-generating mechanism. A useful exact comparison is the transfer from an independent dimensional thermal-pressure input `w` to the displacer-pressure deviation `D`. Define, with physical-time Laplace variable `s` and `a=Ld+Lf`,

```math
A_p(s)=1+L_pC_ps^2,\qquad
A_d(s)=1+R_fC_ds+aC_ds^2,
```

```math
B(s)=C_ds(L_ls+R_l)A_p(s)
+[1+C_ps((L_p+L_l)s+R_l)]A_d(s).
```

Eliminating the two flows and the other pressures from the published five-state equations gives our identity

```math
G(s)=\frac{D(s)}{w(s)}
=-\frac{A_p(s)}{(1+R_{th}C_{ad}s)B(s)
+R_{th}s[C_dA_p(s)+C_pA_d(s)]}.
```

Closing `w=P0 F tanh(K D/P0)` recovers NTP. The denominator is degree five and the relative degree is three. For nominal components, the numerator zeros are at `s=±5.97496 i s^-1`. The formula agrees with an independent matrix-resolvent calculation at 32 complex frequencies to a maximum relative discrepancy of `1.01e-14`. This is an algebraic check, not a chaos test. The exact transfer function, nonlinear gain and state reconstruction are more informative comparison targets than a projected attractor's shape.

Sprott [4, Eq. 12 and Table I] already published chaotic third-order equations with `tanh` saturation. Those particular equations have three equilibria, while positive-parameter NTP has one. This rules out a global bijective conjugacy of the entire flows to those particular examples; it does not exclude a local attractor equivalence, an invariant subsystem, or another known higher-order circuit. A complete classification of saturating-feedback realizations has not been performed.

## Priority and access record

- The 2013 and 2019 full primary papers describe periodic regimes and their bifurcations. Their occurrences of the word *chaos* are bibliographic: Strogatz's textbook in [2] and a Pyragas control reference in [3]. The latter is cited as a possible way to stabilize an unstable cycle, not as a demonstrated chaotic NIFTE regime.
- Exact-model queries combining NIFTE, its expanded name, both article titles, chaos, chaotic and Lyapunov did not locate a verified exact-model chaos demonstration in this pass. This is a bounded search result, not evidence that none exists. The earlier citation snapshot contains three citations of [3] and forty of [2]; not all forty full texts have been read.
- A related 2024 TPTO experiment explicitly reports quasiperiodic mode transitions [5]. A further external-load experiment [6] examines mode transitions with a van der Pol model. Publisher-indexed abstracts/previews were available, but full equations were not inspected; direct publisher fetches returned access/rate errors. These are related apparatus studies, not confirmed equivalences to the five-state or nine-state NTP equations.
- Autonomous thermoacoustic chaos was already reported experimentally by Yazaki in 1993 [7]. The broad claim “first chaos in a thermally driven fluid oscillator” is therefore untenable regardless of what a future NIFTE scan finds.
- Wang's thesis record [8] describes the nonlinear model and a separate switched-linear approximation. An ordinary retry of its linked PDF returned HTTP 429 again. No bypass or author contact was attempted. Its full content remains a priority gap.

For subsequent numerical work, report mathematical boundedness, numerical chaos evidence, physical reconstruction and novelty review as separate outcomes. A surviving chaotic numerical candidate would justify a focused continuation and literature audit; it would not yet establish a new experimental attractor.

## References

1. Markides, C. N., and Smith, T. C. B. (2011). A dynamic model for the efficiency optimization of an oscillatory low-grade heat engine. *Energy* **36**, 6967–6980. [DOI](https://doi.org/10.1016/j.energy.2011.08.051); [institutional full text](https://spiral.imperial.ac.uk/server/api/core/bitstreams/c7601b88-278c-4f37-8bd1-df7eb9851579/content). Sections 2.2 and 3.1–3.4 inspected.
2. Markides, C. N., Osuolale, A., Solanki, R., and Stan, G.-B. V. (2013). Nonlinear heat transfer processes in a two-phase thermofluidic oscillator. *Applied Energy* **104**, 958–977. [DOI](https://doi.org/10.1016/j.apenergy.2012.11.056); [author-hosted paper](https://gstan.bg-research.cc.ic.ac.uk/Markides_Stan_App_Energy_2013.pdf). Tables 1–3 and Appendix A visually inspected.
3. Wang, Y., Markides, C. N., and Chachuat, B. (2019). Optimization-based investigations of a two-phase thermofluidic oscillator for low-grade heat conversion. *BMC Chemical Engineering* **1**, 20. [Primary article](https://doi.org/10.1186/s42480-019-0019-3). Equations 1–8 and the design discussion inspected; [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
4. Sprott, J. C. (2000). Simple chaotic systems and circuits. *American Journal of Physics* **68**, 758–763. [Author-hosted primary paper](https://sprott.physics.wisc.edu/pubs/paper245.pdf). Equation 12 and its discussion inspected.
5. Zhou, Q., Luo, J., Chen, Y., Sheng, S., and Jin, T. (2024). Characterizing the mode transition features of a two-phase thermofluidic oscillator. *Applied Thermal Engineering* **240**, 122225. [DOI](https://doi.org/10.1016/j.applthermaleng.2023.122225). Publisher-indexed abstract, not full equations.
6. Zhou, Q., Chen, Y., Sheng, S., He, T., and Jin, T. (2024). Experimental investigation on the mode transition characteristics of a phase-change thermofluidic oscillator with external loads. *Applied Thermal Engineering* **257**, 124192. [Publisher page](https://www.sciencedirect.com/science/article/abs/pii/S135943112401860X). Publisher-indexed preview only.
7. Yazaki, T. (1993). Experimental observation of thermoacoustic turbulence and universal properties at the quasiperiodic transition to chaos. *Physical Review E* **48**, 1806–1818. [Primary abstract](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.48.1806). Spontaneous and forced regimes are distinguished in the abstract.
8. Wang, Y. (2020). *Two-phase thermofluidic engines for low-grade heat recovery: system analysis and supporting algorithms*. Imperial College London doctoral thesis. [Institutional record](https://spiral.imperial.ac.uk/entities/publication/ff79c230-8eff-4695-8e58-e07d63edca1e). Linked PDF returned HTTP 429 on this retry.
