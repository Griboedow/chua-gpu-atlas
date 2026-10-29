# One-transistor organic electrochemical oscillator: count the states first

Research notebook, 29 October 2026. **The published autonomous two-state model is not a candidate for a classical chaotic strange attractor on a compact regular domain.** This is a dimensional exclusion, not a negative numerical survey or a statement about every physical organic transistor.

Bisquert and Tessler [1] obtain sustained oscillations using one organic electrochemical transistor and passive RC components. The two independent states are the external capacitor voltage `u` and the internal electrochemical potential `v`. With constant applied voltage, their Eqs. 31–32 can be written

$$\dot v=G(u,v)=\frac{-v+[R_1u+R_t(v)V_{app}]/[R_1+R_t(v)]}{R_3C_\mu},$$
$$\dot u=\frac{-u/R_0+(V_{app}-u)/[R_1+R_t(v)]-C_\mu G(u,v)}{C_0}.$$

All parameters are fixed here. On a region where the capacitances and resistances are positive and the constitutive function is smooth, this is an autonomous planar vector field. Transconductance changes its stability, not its state dimension. The Poincare–Bendixson restriction rules out a classical recurrent chaotic strange attractor for bounded regular motion; a positive short-time stretching rate near an unstable equilibrium would not evade that restriction.

The source's stationary S-shaped current-voltage relation does not itself introduce a switching memory state. These ODEs remain single valued as long as their explicit constitutive law remains regular. One should not replace a smooth branch by arbitrary hysteretic jumps and then label the resulting chaos as a property of the original model.

For an actual device, additional ionic relaxation, spatial transport, external periodic forcing or a network could increase the effective dynamical dimension. Each would require a specified physical model and a separate literature audit. The paper's suggestion of richer internal dynamics is a possible research direction, not a supplied third-state equation. No GPU sweep was run for the two-state equations because it cannot override the dimensional restriction.

The full primary XML and MathML equations were retrieved from the official Europe PMC endpoint; source provenance is saved in [oect-source-audit.json](../results/candidate-audits/oect-source-audit.json). No source PDF or XML is redistributed with our code.

## Reference

1. Bisquert, J., and Tessler, N. (2025). A one-transistor organic electrochemical self-sustained oscillator model for neuromorphic networks. *Newton* **1**(8), 100207. [DOI](https://doi.org/10.1016/j.newton.2025.100207); [open primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12501206/). In particular Eqs. 9, 17–19 and 31–32. Our exclusion applies to this autonomous two-state model, not to arbitrary higher-dimensional OECT circuits.
