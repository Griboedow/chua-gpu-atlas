# Solar-array thermal flutter: source gate

Audit: 30 October 2026 UTC. **Excluded from the current numerical search: no unambiguous, reproducible autonomous finite-state model was recovered.** This is a model-definition decision, not a proof that thermal flutter cannot be chaotic. No simulation, fitted forcing, extra thermal state or replacement mode truncation was introduced.

## Source inspected

The final Motaharifard–Daneshjou–Bakhtiari paper [1] was read through its complete Europe PMC XML and publisher PDF text. It reports limit cycles, internal resonance, quasiperiodicity and period-doubling-related instability, but a positive Lyapunov-spectrum calculation was not found. The displayed thermal equations include prescribed harmonic heating; the eclipse section supplies a time-dependent illumination-related speed ratio. Later structural equations use apparently fixed thermal coefficients. These stages must not be silently identified with one autonomous system.

Equation (40) contains three second-order mechanical coordinate families and acceleration-dependent nonlinear terms. Thus an explicit implementation requires inversion of a state-dependent mass matrix. One mode per family would give six mechanical states, but this audit did not verify that truncation as the numerical implementation. The paper also discusses dropping in-plane motion, yielding a bending–torsion reduction. Neither choice supplies a verified full parameter/initial-state tuple.

Coefficient definitions or numerical values sufficient to reconstruct all $C_{uwj},C_{wwj},C_{awj}$ terms were not located. The XML contains no supplementary-material attachment; data are offered on request. Tables give useful physical constants, but not a self-contained executable closure. [1]

## Checks that matter before implementation

The following are our consistency checks of the displayed equations, not proposed corrections.

* **Forcing versus autonomy.** Equations (9), (15), (18) retain $\cos(\omega t)$. Adding a phase variable would represent an externally driven system, not establish spontaneous oscillation under constant heating. Freezing the eclipse schedule or deleting the harmonic term would change the model. Conversely, treating Eq. (40) as autonomous requires a justified derivation of its fixed coefficients.
* **Thermal normalization.** In Eq. (15), let $K$ be the coefficient multiplying the perturbation temperature. Equation (19) defines $k_n^m=K/(\rho c)$, while Eq. (18) retains $\rho c$ on the acceleration and damping terms and uses $k_n^m$ for stiffness. As displayed, these do not follow from (15) by substituting a steady mean temperature alone. No correcting factor is assumed here.
* **Boundary consistency.** Equation (10) specifies zero values at the rectangular boundary. An individual cosine product in Eq. (11), plus a nonzero mean, does not automatically satisfy those conditions. Cancellation constraints or a different boundary convention would need to be established before truncation.
* **Time dependence.** Equation (20) has a decaying transient plus an undamped harmonic particular solution. Equation (22), as encoded, places the harmonic terms inside the decay factor; later load expressions contain $e^{-1/(2\tau_0)}$ without $t$. A source-confirmed choice is necessary. PDF image rendering failed, so these specific typography concerns remain based on the authoritative XML and text, not visual equation verification.
* **Physical parameter meaning.** The text defines $U=\sqrt{k/(\rho c\tau_0)}$ as thermal-wave speed and describes ratios to light speed, including values above one. Using Table 2 and literal $U=0.95c_{\mathrm{light}}$ gives $\tau_0\approx5.82\times10^{-23}$ s. This arithmetic does not validate that interpretation. A dimensional/reference-speed explanation is needed before treating the plotted ratio as a physical control.

The strong conclusion is limited: an exact-model sweep would currently require undocumented choices. No inference about the correctness of every result in the paper follows from this gate.

## A useful mathematical exclusion within the thermal subsystem

For constant positive coefficients, the printed mean-temperature equation (14), with $A=T_n^{(0)}$, is

```math
\tau_0\rho c\ddot A+\rho c\dot A+
\frac{\varepsilon\sigma}{t_s}A^4=
\frac{\alpha_s q}{\pi t_s}.
```

Its elementary energy function satisfies

```math
E=\frac{\tau_0\rho c}{2}\dot A^2+
\frac{\varepsilon\sigma}{5t_s}A^5-
\frac{\alpha_s q}{\pi t_s}A,
\qquad \dot E=-\rho c\dot A^2.
```

On any compact invariant subset of the physical domain $A\ge0$, the largest invariant subset of $\dot E=0$ consists of the unique positive equilibrium. Therefore this **isolated mean-temperature equation** cannot supply a chaotic attractor. With the mean frozen at equilibrium, the displayed perturbation equation is a stable linear, harmonically forced oscillator for positive coefficients. This does not exclude chaos in the coupled nonlinear mechanics; it prevents misidentifying thermal transient complexity as an independently chaotic thermal oscillator.

## Prior work and forward check

The direct 2025 predecessor [2] already studies nonlinear modal interaction during penumbra transitions. Its full equations were not retrieved, so exact equivalence was not established. Flexible-panel spacecraft chaos also predates both papers: Aslanov [3] studies aerodynamic/gravitational attitude dynamics, a different mechanism. This is broad prior art, not an asserted conjugacy.

The OpenAlex snapshot returned two citing works. The accessible primary sections of Li et al. [4] concern spin-modulated thermal bending and model-admissibility limits; they do not establish chaos in [1]'s equations. Zhang et al. [5] is a free-vibration/modal study according to its primary abstract; full text was not inspected. Index coverage and this bounded reading do not establish absence of an exact-model chaos publication.

**Disposition:** preserve the audit; do not manufacture a small ODE from the incomplete coefficient set. Reopening requires a fully specified reduction, forcing convention, parameter normalization and source-control trajectory. A three-dimensional projection of several mechanical/thermal modes would not itself establish unusual attractor topology.

## References

1. O. Motaharifard, K. Daneshjou and M. Bakhtiari, “Nonlinear stability and vibration of flexible spacecraft solar arrays under thermally induced flutter during the penumbra phase,” *Scientific Reports* **16**, 9856 (2026). [DOI](https://doi.org/10.1038/s41598-026-38274-4). [Complete XML read](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13018175/fullTextXML). [Publisher PDF](https://www.nature.com/articles/s41598-026-38274-4.pdf).
2. “Investigation of nonlinear vibration in spacecraft with flexible solar arrays during penumbra phase,” *Aerospace Science and Technology* **167**, 110615 (2025). [Primary publisher excerpt](https://www.sciencedirect.com/science/article/abs/pii/S1270963825006868). [DOI](https://doi.org/10.1016/j.ast.2025.110615).
3. V. S. Aslanov, “Chaotic attitude dynamics of a LEO satellite with flexible panels,” *Acta Astronautica* **180**, 538–544 (2021). [Author-hosted primary manuscript](https://aslanov.ssau.ru/papers/DAC_Chaotic_attitude_dynamics_of_a_LEO_satellite_with_flexible_panels.pdf). [DOI](https://doi.org/10.1016/j.actaastro.2020.12.055). Abstract inspected; no equation-equivalence claim.
4. O. Li et al., “Spin-Modulated Thermoelastic Response of a Flexible Spacecraft Appendage Under Attitude-Dependent Solar Radiation,” *Aerospace* **13**, 645 (2026). [Primary article](https://www.mdpi.com/2226-4310/13/7/645). DOI 10.3390/aerospace13070645.
5. S. Zhang et al., “Semi-analytical modeling and free-vibration analysis of composite thin-walled lenticular tubes,” *International Journal of Mechanical Sciences* **327**, 111933 (2026). [Primary abstract](https://www.sciencedirect.com/science/article/abs/pii/S0020740326007861). DOI 10.1016/j.ijmecsci.2026.111933.

Source files remain in the ignored cache; their publisher licenses are not replaced by this repository's MIT license. [Machine-readable audit](../results/candidate-audits/thermal-flutter-source-audit.json).
