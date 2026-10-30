# Mineral protocells: planar oscillations and irreversible membrane growth

Audit: 30 October 2026 UTC. **Excluded from the present search for a sustained strange attractor.** The exact 2024 model has two dynamic states. Its three-state growing-membrane predecessor has a strictly increasing thickness during active precipitation and was designed to describe eventual rupture. The older authors already considered transient chaos and reported a negative sensitivity check. No simulations were run for this audit.

## Exact 2024 model and physical scope

Ding, Cardoso, and Cartwright [1] couple precipitate concentration to internal pressure. The two fluid speeds are algebraic functions, not extra dynamic states. The primary article's equations (1)–(3) were checked against its rendered PDF, especially page 3. Dropping hats on dimensionless variables, the complete reduced system is

```math
\dot c=[M-cp+c^2p]H(1-c),\qquad
\dot p=N+[-p+cp]H(1-c).
```

Here $M,N>0$ are constants and $H$ is the Heaviside function. The active physical interior is $0<c<1$, $p>0$. The membrane thickness and flow areas are fixed. The concentration and osmotic-pressure gradients supply energy; the model assumes them maintained over the interval considered. It also uses the approximation that osmotic pressure is much larger than internal pressure. A pressure excursion invalidating that approximation is not a new physical attractor.

For reference, the dimensional definitions are

```math
u_{\mathrm{in}}=\frac{k_{\mathrm{in}}}{\mu}\frac{p_o-p}{L_m},\qquad
u_{\mathrm{out}}=\frac{k_{\mathrm{out}}}{\mu}\frac{p}{L_m}
\left(1-\frac{c}{c^*}\right)H\left(1-\frac{c}{c^*}\right),
```

```math
\dot c=\left(\frac{u_{\mathrm{in}}c_o}{L_r}
-\frac{u_{\mathrm{out}}c}{R}\right)H(u_{\mathrm{out}}),\qquad
\dot p=\gamma(u_{\mathrm{in}}A_{\mathrm{in}}-u_{\mathrm{out}}A_{\mathrm{out}}),\qquad
\gamma=\frac{2L_mE}{R_pV_p}.
```

The nondimensional parameters and scales are

```math
N=\frac{p_ok_{\mathrm{in}}A_{\mathrm{in}}V_pR_p}
{2EL_mk_{\mathrm{out}}A_{\mathrm{out}}^2R},\qquad
M=N\frac{c_oA_{\mathrm{out}}R}{c^*A_{\mathrm{in}}L_r},
```

```math
t_s=\frac{\mu R_pV_p}{2EA_{\mathrm{out}}k_{\mathrm{out}}},\qquad
p_s=\frac{2EL_mA_{\mathrm{out}}R}{V_pR_p},\qquad c_s=c^*.
```

There is no hidden third state in these definitions. Assigning an independent relaxation equation to permeability, osmotic pressure or flow would create a different model and is not authorized by these equations.

## Why the oscillations are not a chaos lead

Inside $c<1$, this is a smooth autonomous planar polynomial flow. The planar limit-set restriction excludes a regular bounded strange attractor in a compact subset of that interior. At $c>1$, the printed cutoff gives $\dot c=0$ and $\dot p=N>0$, not recurrent bounded dynamics. The exact convention at the switch must be specified for numerical integration; it does not provide a reset or another state. We have not replaced it with a smoothing or invented rupture-and-reformation rule.

An algebraic check also makes the known stability mechanism transparent. For $0<M<N$, the interior equilibrium is

```math
c_e=\frac{M}{N},\qquad p_e=\frac{N^2}{N-M}.
```

Our differentiation gives

```math
\det J_e=N-M>0,\qquad
\mathrm{tr}\,J_e=\frac{N(2M-N)}{N-M}-\frac{N-M}{N}.
```

The trace crosses zero at

```math
M=N\left[1+N-\sqrt{N^2+N}\right].
```

This is the linear oscillatory-instability condition, not by itself a proof of the type or stability of a nonlinear cycle. The source's Fig. 3 already maps stable, oscillatory and explosive regimes. A fixed-point/periodic control would be useful for a tutorial, but a GPU chaos survey of this exact planar model has no justified target.

## The three-state predecessor does not repair the problem

The full author manuscript of Ding et al. [2] has concentration $c$, pressure $p$, and membrane thickness $L_m$ as states. Its equations (2a)–(2c), visually checked on manuscript page 6, retain the concentration/pressure balances above and add

```math
\dot L_m=K\left[-\dot c+\frac{u_{\mathrm{in}}c_{\mathrm{Si}}}{L_r}\right]
H(u_{\mathrm{out}}),\qquad K=\frac{MR}{\rho_s}>0.
```

Here $K$ abbreviates the positive dimensional material prefactor printed in that paper; its $M$ is not the dimensionless control parameter above. In the active interior, where $u_{\mathrm{out}}>0$ and $c>0$, substitution of the concentration balance yields the exact identity

```math
\dot L_m=\frac{K}{R}u_{\mathrm{out}}c>0.
```

Consequently there is no periodic orbit or recurrent compact invariant set wholly within an active physical region on which this derivative stays positive. More specifically, any compact invariant set lying strictly inside that region would have a positive lower bound on $\dot L_m$, contradicting bounded thickness. This is a scoped structural exclusion; it is not a theorem about arbitrary boundary extensions, dissolution, repair, or a newly forced model.

The primary paper describes loss of permeability and rupture following all three tested regimes. On manuscript page 12 it explicitly raises possible **transient chaos**, then reports no change in explosion time when initial dimensionless pressure was varied from 0 to at least 10 and initial concentration from 0 to 0.4. The authors interpret the experimental spread as likely stochastic rather than produced by their model. This is already a published chaos question with a limited negative test, not evidence of a sustained strange attractor. It also does not rigorously exclude all transient chaotic scattering.

The 2024 paper explicitly identifies its frozen thickness as a short-time approximation to this growing system. A finite sequence of oscillations before membrane rupture cannot be promoted to an asymptotic attractor by discarding the thickness coordinate in a plot.

## Source access, earlier work, and later citations

The complete 2024 primary PDF and the 2019 growing-membrane author manuscript were obtained from Cambridge's institutional repository and retained only in the ignored source cache. The latter's publisher lists a supplement, but its download returned HTTP 403; the institutional manuscript bundle contains only the main manuscript. Thus the complete growth-dependent permeability law and numerical parameter table from that supplement are **not independently verified here**, and the growing model is not presented as a simulation-ready package.

The cited 2019 *Interface Focus* concentration-cycle paper [3] and the 2020 chemical-nanopump paper [4] establish the earlier research lineage. Their publisher/index records were checked, but the older full-text XML retrieval failed; they are not substitutes for the two full primary texts read here.

OpenAlex forward metadata returned five citing records for [1] and ten for [3]. These were screened at title/DOI level, not all read in full. They include later work on membrane growth and fracture; the search did not verify a higher-state, maintained-reservoir oscillator with a closed set of equations that evades the restrictions above. This is a bounded screening result, not a claim that chaos is absent throughout chemical-garden physics or that an unreturned paper does not exist.

## References

1. Y. Ding, S. S. S. Cardoso, and J. H. E. Cartwright, “Dynamics of the osmotic lysis of mineral protocells and its avoidance at the origins of life,” *Geobiology* **22**, e12611 (2024). [DOI](https://doi.org/10.1111/gbi.12611); [full institutional primary PDF](https://api.repository.cam.ac.uk/server/api/core/bitstreams/1fb23cec-d7ec-411c-93a6-4e0bfe2010e3/content).
2. Y. Ding, C. M. Gutiérrez-Ariza, C. I. Sainz-Díaz, J. H. E. Cartwright, and S. S. S. Cardoso, “Exploding Chemical Gardens: A Phase-Change Clock Reaction,” *Angewandte Chemie International Edition* **58**, 6207–6213 (2019). [DOI and supplement listing](https://doi.org/10.1002/anie.201812331); [full institutional author manuscript](https://api.repository.cam.ac.uk/server/api/core/bitstreams/ec8ecd2c-0956-4d94-8832-d49a9254860d/content).
3. Y. Ding, J. H. E. Cartwright, and S. S. S. Cardoso, “Intrinsic concentration cycles and high ion fluxes in self-assembled precipitate membranes,” *Interface Focus* **9**, 20190064 (2019). [DOI](https://doi.org/10.1098/rsfs.2019.0064); [indexed primary record](https://pubmed.ncbi.nlm.nih.gov/31641435/).
4. Y. Ding, J. H. E. Cartwright, and S. S. S. Cardoso, “Convective flow driven by a chemical nanopump,” *Physical Review Fluids* **5**, 082201(R) (2020). [Primary publisher record](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.082201).
