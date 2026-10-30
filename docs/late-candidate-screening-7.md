# Late source screening: collective motors and soft fluidic valves

Status: 30 October 2026. The usable outcome is an **old, fully specified molecular-motor model whose authors explicitly searched for chaos without reporting a positive example**. It is not evidence of an unpublished attractor.

## Selected exact target

Alouges et al. (2025), *Some mathematical models for flagellar activation mechanisms*, provides an exact five-state first-harmonic subsystem, published periodic controls, and a local Hopf analysis. Addition and subtraction of its motor-density coefficients close a three-state mechanical core. We then recovered an exact equivalence to Oriola et al. (2014), *Subharmonic oscillations of collective molecular motors*, Eq. (10).

The [separate source audit](motor-pair-source-audit.md) gives the coordinate transformation, two essential parameters, physical rate restrictions, an invariant bounded region, Jacobian, control settings and the independent Python implementation. The 2014 paper already describes intricate periodic orbits and performs a Lyapunov calculation. A new drawing of such a cycle would not be a new chaotic attractor, and the 2025 two-row interpretation does not create a new dynamical mechanism.

The remaining question is narrow: does this exact existing core have bounded chaotic regimes outside its previously reported periodic examples? No positive answer or complete impossibility theorem was established by this screening. The separate CPU control reproduces and checks a reached periodic branch before a possible survey by the parent task.

## Rejected alternatives from this block

**Soft hysteretic valve, van Laake, Comoretto and Overvelde (2024).** The baseline Eqs. (14–15) evolve pressure difference and overdamped valve displacement, giving two continuous states. For a regular autonomous branch this is not a three-state chaos target. The modified piecewise spring law needs its own switching/well-posedness analysis; a blanket smooth planar theorem should not be applied across discontinuities. No additional inertia or memory was invented.

**Networks of fluidic relaxation oscillators, Frankevoort (2023).** The institutional thesis is from **2023**, not 2025 as a search-index date initially suggested. Its single-oscillator Eq. (2.9) adds a low-pass-filter flow variable to the two pressure states. The text explicitly introduces that third state as a numerical regularization of an algebraic flow relation; the continuation study changes its time constant to reduce stiffness. It is a complete published mathematical model, but chaos in it could depend on artificial relaxation. A faithful network study should separately audit the underlying hybrid valve model, network dimension, and the limit of the artificial time constant. This block did not obtain a new physical chaos candidate by promoting the regularization to measured valve inertia.

## Scope of source retrieval

Primary full text was read for the 2024 valve paper and 2023 thesis, for the 2025 motor preprint, and for the 2014 motor paper through an author-uploaded copy and its detailed thesis treatment. For the motors, forward discovery and closest retrieved successors are recorded in the [machine-readable audit](../results/motor-pair/source-audit.json). Search counts and missing results are not proofs of priority. No source code was executed or imported, and no large simulation was performed by this source-screening agent.

## References

1. F. Alouges et al. (2025), *Mathematical Models and Methods in Applied Sciences* **35**, 2395–2424. [DOI 10.1142/S0218202525500423](https://doi.org/10.1142/S0218202525500423); [complete preprint](https://arxiv.org/html/2409.03506v2).
2. D. Oriola et al. (2014), *EPL* **107**, 18002. [DOI 10.1209/0295-5075/107/18002](https://doi.org/10.1209/0295-5075/107/18002).
3. L. C. van Laake, A. Comoretto, J. T. B. Overvelde (2024), “On the coexistence of pressure regulation and oscillation modes in soft hysteretic valves,” *Journal of Fluids and Structures* **126**, 104090. [DOI](https://doi.org/10.1016/j.jfluidstructs.2024.104090); [institutional full paper](https://ir.amolf.nl/pub/10955/16970publishedVersion.pdf).
4. M. G. H. Frankevoort (2023), *Modeling networks of fluidic relaxation oscillators*, bachelor thesis, Eindhoven University of Technology, 25 January 2023. [Institutional full thesis](https://pure.tue.nl/ws/portalfiles/portal/293050364/1366157_NetworksOfFluidicRelaxationOscillators.pdf).
5. L. C. van Laake et al. (2022), “A fluidic relaxation oscillator for reprogrammable sequential actuation in soft robots,” *Matter* **5**, 2898–2917. [DOI](https://doi.org/10.1016/j.matt.2022.06.002); [author data and code deposit, not executed](https://zenodo.org/records/6576063).
