# Microfluidic source audit: the final supplement retains the discrepancy

Audit date: 30 October 2026 UTC. **The final published supplement agrees with the arXiv supplement on the coefficients that fail our periodic source control.** This rules out a simple preprint-versus-publication explanation. It does not identify a unique correction or exclude a coexisting cycle with an untested basin.

## What was checked

The primary model is Case, Angilella, and Motter's pressure-driven five-channel network [1]. Two flow constraints reduce its six listed variables to four independent states. The [earlier extraction](late-candidate-screening-4.md) records all equations and dimensional conversions; the [numerical report](microfluidic-research.md) records independently implemented source controls. This source-audit block performed no trajectory integration.

The final published `aay6761_SM.pdf` was obtained through the public [Europe PMC supplementary-files endpoint](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7220308/supplementaryFiles), rather than inferred from extracted arXiv text. Its 11 pages include the paper's cover sheet, supplement, and figures. Pages 5, 9, and 10 were rendered and visually checked: S1–S2, S5–S13, and Fig. S5 respectively. The file is 2,367,829 bytes; its SHA256 is `885b93467f60324f762da7bbe7556ff71ebd26199f0c7964f2360689a5b22357`. It remains in the ignored local source cache, not among the project's MIT-licensed results. The retrieved archive contains the supplement and supporting media, but no author implementation or raw initial-condition table.

The relevant formulas, signs and constants agree with [arXiv v1](https://arxiv.org/pdf/2006.09400). Crossref metadata returned no correction relation for the DOI; this bounded metadata check is not proof that no correction exists elsewhere.

## The printed coefficient is negative on the plotted flow scale

Source S1 and S2 require

```math
\eta(Q_1)=p+qQ_1,\qquad
\xi(Q_1)=\frac{(r^{\ast})^3\eta(Q_1)-Q_3^{\ast}}{r^{\ast}}.
```

Both coefficients therefore vary with $Q_1$. Freezing $\xi$ while retaining variable $\eta$ would change the source model. S13 prints

```math
r^{\ast}=5.5,\quad Q_3^{\ast}=25,\quad
p=0.12,\quad q=0.8\times10^{-5},\quad \epsilon=10.9.
```

Although the main text describes $\eta$ and $\xi$ as positive, these constants imply

```math
\xi(Q_1)>0\quad\Longleftrightarrow\quad
Q_1>\frac{Q_3^{\ast}/(r^{\ast})^3-p}{q}
=3782.87002254.
```

Fig. S5 plots total flow roughly between 0.4 and 0.9 mL per second per millimeter, equivalent to dimensionless $Q_1+Q_2$ between 400 and 900. With forward flow, $Q_1$ is below that total. Thus the printed $\xi$ is negative throughout that plotted range. The dimensional conversion follows the primary supplement, not a guessed plotting rescaling.

At the independently calculated equilibrium for equal 250 Pa inlet pressures, $Q_1=694.316102856$, the literal formulas give

```math
\eta=0.125554528823,\qquad \xi=-0.747430048563.
```

The [numerical controls](microfluidic-research.md) show convergence to this stable equilibrium from three stated starts, despite the periodic example displayed at 250 Pa in Fig. S5. They do not exhaust the basin structure. No factor-of-ten change to $q$, sign repair, frozen coefficient, or refitted relaxation rate has been adopted. Several such modifications could change the dynamics, and the primary evidence does not select one.

## Which source comparisons are legitimate

The constant $\epsilon$ is a fitted relaxation parameter. The flow dependence of $\eta$ and consequent dependence of $\xi$ are explicit. The source keeps the base resistance $L_b=42$ constant after noting that making it flow-dependent did not markedly improve the fit. None of these facts licenses silently changing the closure to recover a desired oscillation.

Fig. S5 is the reduced dynamical model's bifurcation diagram. Fig. S2's oscillation periods and Fig. S3's traces are Navier–Stokes results. Agreement with those periods is not a source control for the four-state closure. The isolated-channel model and its fitted constants in the main text are another, lower-dimensional model and should not be mixed into the network equations.

The absence of author initial conditions or an implementation leaves an unresolved reproducibility gap. The exact code, units, signs and algebra have been cross-checked in our [numerical](microfluidic-research.md) and [structural](microfluidic-structure.md) reports. The present evidence supports pausing the large chaos scan until the periodic control or a justified source correction is established.

## Structural limits and priority

For positive $\eta$ and negative $\xi$, the vortex equation has

```math
\frac{\partial\dot r}{\partial r}
=\epsilon\left[\xi-3\eta(r-r^{\ast})^2\right]<0.
```

This is damping in one state variable, not a theorem excluding oscillation or chaos in the coupled system. In the eliminated flow equations, vortex-to-flow derivatives contain $R'(r)=2\gamma(r-r_b)$ and change sign across $r_b$. A standard fixed-sign cyclic-feedback theorem therefore cannot be applied directly over a domain spanning that point. The separate [structural audit](microfluidic-structure.md) also supplies a counterexample to global forward-flow invariance. Physical interpretation must be checked along trajectories rather than assumed from positive inlet pressures alone.

Shao, Angilella, and Motter [2] already demonstrate chaotic Eulerian flow-rate dynamics in a different serialized rigid network using Navier–Stokes simulations. This excludes a broad first-chaos claim for rigid fluidic feedback. Their model is not the exact four-state S5–S13 closure audited here. The bounded exact-title, equation and forward-citation checks did not establish whether chaos in this particular closure has already been reported. Failure to retrieve such a result is not evidence of scientific novelty.

## References

1. D. J. Case, J.-R. Angilella, and A. E. Motter, “Spontaneous oscillations and negative-conductance transitions in microfluidic networks,” *Science Advances* **6**, eaay6761 (2020). [DOI](https://doi.org/10.1126/sciadv.aay6761); [full primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7220308/); [official supplement archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7220308/supplementaryFiles); [author manuscript with supplement](https://arxiv.org/pdf/2006.09400).
2. Y. Shao, J.-R. Angilella, and A. E. Motter, “Emergent oscillations and chaos in noncompliant microfluidic networks,” *Physical Review Fluids* **10**, 054401 (2025). [DOI](https://doi.org/10.1103/PhysRevFluids.10.054401); [full author manuscript](https://arxiv.org/html/2505.00068v1).
