# Material-memory screening: a three-state model that does not need a chaos scan

Audit date: 2026-10-30 UTC. This is a literature and analytical screening result; no trajectories were integrated in this block. Neither screened system is promoted as a new attractor candidate.

## 1. Halide-perovskite memristor: an exact three-state negative-feedback loop

Bisquert and Guerrero [1] study a mixed ionic/electronic memristor and report simulated self-sustained spiking under current drive. The model follows Berruet et al. [2]. With constant applied current, its three states are voltage $u$, delayed internal current $i$, and occupation $f$. Equations (1), (11), and (12) are:

```math
C_m\dot u=I-\frac{u}{R_b}-i,
\qquad \tau_d\dot i=i_0f-i,
\qquad \tau_k\dot f=1-f-E(u)f,
\qquad E(u)=\exp\!\left(-\frac{u-V_T}{V_m}\right).
```

All capacitance, resistance, time, current-scale, and voltage-scale parameters are positive; $I$ is the imposed constant drive. In particular, the last equation is **not** the constant-time-constant relaxation $\tau_k\dot f=f_{ss}(u)-f$.

Figure 5b of [1] uses $R_b=1$, $i_0=10$, $V_T=1$, $V_m=0.05$, $(\tau_m,\tau_d,\tau_k)=(0.1,0.1,0.1)$, $I=6$, with $\tau_m=R_bC_m$. Figure 5c changes $\tau_k$ to 1. Equations and figure labels were checked visually in the institutional manuscript, pages 7 and 10. The paper predicts these spiking regimes; its conclusion does not claim their experimental observation.

### Physical domain and boundedness: our derivation

The strip $0\le f\le1$, $0\le i\le i_0$ is forward invariant. At the two occupation boundaries the vector field points inward, while the current equation preserves its interval. Voltage satisfies

```math
u_-\le u(t)\le u_+,
\qquad
u_- = \min\{u(0),R_b(I-i_0)\},
\qquad
u_+ = \max\{u(0),R_b I\}.
```

Let $E_-=E(u_+)>0$ and $E_+=E(u_-)\lt \infty$. Comparison with scalar linear equations gives

```math
\liminf_{t\to\infty}f(t)\ge\frac{1}{1+E_+}>0,
\qquad
\limsup_{t\to\infty}f(t)\le\frac{1}{1+E_-}\lt 1.
```

Thus every physical trajectory has a compact omega-limit set strictly inside the occupation interval. The current also eventually stays away from its boundary values. There is exactly one equilibrium: its voltage solves

```math
I=\frac{u_*}{R_b}+\frac{i_0}{1+E(u_*)},
```

whose right-hand side is strictly increasing from minus infinity to plus infinity.

### The theorem mapping, including equilibrium-containing limit sets

Set $g=\log[f/(1-f)]$, write $\sigma(g)=(1+e^{-g})^{-1}$, and order the coordinates as $(u,g,i)$. The transformed vector field is smooth on all of $\mathbb R^3$:

```math
\dot u=\frac{I-u/R_b-i}{C_m},
\qquad
\dot g=\frac{1+e^{-g}-E(u)(1+e^g)}{\tau_k},
\qquad
\dot i=\frac{i_0\sigma(g)-i}{\tau_d}.
```

Its three cyclic predecessor derivatives have fixed, strictly nonzero signs:

```math
\frac{\partial\dot u}{\partial i}=-\frac1{C_m}\lt 0,
\qquad
\frac{\partial\dot g}{\partial u}=\frac{E(u)(1+e^g)}{\tau_kV_m}>0,
\qquad
\frac{\partial\dot i}{\partial g}=\frac{i_0\sigma(g)(1-\sigma(g))}{\tau_d}>0.
```

Each equation depends only on its own coordinate and its cyclic predecessor. This is an irreducible monotone cyclic system with negative loop sign. The comparison bounds above ensure compactness after the logit transformation; a boundary occupation at the initial instant enters the interior immediately.

The relevant result is the Mallet-Paret–Smith theorem [3]. The original 1990 proof was not accessible in this bounded audit. We checked two full primary research papers that explicitly state the result: Wang–Zhou [4], introduction equations (1.1)–(1.2), state the planar embedding of bounded omega-limit sets; Czechowski–Piliouras [5], Definition 1 and Theorem 2, pages 321–322, state the alternatives including equilibria and connecting orbits, on an open convex domain with strictly nonzero cyclic derivatives.

For this model the hypotheses map to the equations above and the equilibrium is unique. The allowed limit-set alternatives are an equilibrium, a periodic orbit, or an equilibrium together with connecting homoclinic orbits. **This is not a claim that every solution converges to a periodic orbit.** In particular, an equilibrium-containing limit set is not discarded by invoking only the weaker equilibrium-free periodic-orbit statement.

On this theorem application, the exact autonomous, constant-current three-state model is excluded from a search for a recurrent strange attractor. This is an analytical screening conclusion, not a numerical failure to find chaos. It does not cover extra circuit states, networks, nonmonotone reaction laws, stochastic models, or periodic external forcing. Those would require their own published equations and separate analysis.

### A coefficient inconsistency in the accessible manuscript

Independently expanding the characteristic determinant of the equilibrium Jacobian gives

```math
\det(\lambda I-J)=(\lambda+a)(\lambda+b)(\lambda+c)+d,
\quad
a=\frac1{\tau_m},\quad b=\frac1{\tau_d},
\quad c=\frac1{f_*\tau_k},
\quad d=\frac{i_0(1-f_*)}{C_m\tau_d\tau_kV_m}.
```

Therefore $c_0=abc+d$. Equation (20) on page 9 of the accessible institutional manuscript displays only $d$, omitting the positive term $abc$. Its displayed Jacobian agrees with the original ODE. This mismatch was checked visually, rather than inferred from PDF text extraction. The final publisher equation image was not independently recovered, so this is specifically a **manuscript consistency observation**, not an assertion about an uninspected final typeset equation or an erratum. No numerical Hopf boundary should be copied from that coefficient without recomputation. The cyclic-feedback conclusion does not depend on this discrepancy.

### Citation and equivalence screen

An OpenAlex forward-citation query for [1] returned 44 records; 38 included abstracts. Titles and available abstracts were screened for chaos, chaotic, Lyapunov, and Poincare terms, with no matches. These metadata are not a full-text audit of all citing papers. Related primary work [2] supplies the earlier physical model; the connection to monotone cyclic oscillators is structural, not a claim that this memristor is globally conjugate to a particular Goodwin oscillator.

The search did not establish priority for applying [3] to this exact device model. The useful outcome is a reproducible reason to avoid a GPU scan, independent of whether someone has already written down the same application.

## 2. LCO/LCMO electrothermal neuronal oscillators: two states, not three

Jardali et al. [6] use temperature plus capacitor voltage for LaCoO3, and temperature plus inductor current for La0.7Ca0.3MnO3. Their Methods equations (6a–6b) and (7a–7b), and section 4.4 explicitly identify two ODEs. The apparent material inductance/capacitance does not add a third independent state. With the AC probe removed and a constant bias, each compact oscillator is planar. They are excluded from this autonomous three-or-more-state search; an AC-driven experiment would be a different nonautonomous question. We did not transcribe the material-fit coefficients or run the authors' code, because the dimensional exclusion already settles the present screening criterion.

## Reproducibility and limits

The machine-readable companion is [material-memory-screening.json](../results/candidate-audits/material-memory-screening.json). Source PDFs/XML are cached locally under ignored `outputs/candidate_research/`; no author code was executed, copied into the project, or modified. No access barriers were bypassed. No new viable chaos candidate or priority claim emerges from this block.

## References

1. J. Bisquert and A. Guerrero, “Dynamic Instability and Time Domain Response of a Model Halide Perovskite Memristor for Artificial Neurons,” *J. Phys. Chem. Lett.* **13**, 3789–3795 (2022). [DOI](https://doi.org/10.1021/acs.jpclett.2c00790); [institutional manuscript](https://repositori.uji.es/bitstreams/7cfd1818-1271-4bb8-92db-5a91d7ac1666/download); [public full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9974066/fullTextXML).
2. M. Berruet et al., “Physical Model for the Current–Voltage Hysteresis and Impedance of Halide Perovskite Memristors,” *ACS Energy Lett.* **7**, 1214–1222 (2022). [DOI](https://doi.org/10.1021/acsenergylett.2c00121).
3. J. Mallet-Paret and H. L. Smith, “The Poincare-Bendixson theorem for monotone cyclic feedback systems,” *J. Dyn. Differ. Equ.* **2**, 367–421 (1990). [DOI and publisher abstract](https://doi.org/10.1007/BF01054041).
4. Y. Wang and D. Zhou, “Transversality for Cyclic Negative Feedback Systems,” arXiv:1412.1540 (2014). [Full primary manuscript](https://arxiv.org/html/1412.1540).
5. A. Czechowski and G. Piliouras, “Poincaré-Bendixson Limit Sets in Multi-Agent Learning,” *AAMAS 2022*, pp. 318–326. [Proceedings PDF](https://aamas.csc.liv.ac.uk/Proceedings/aamas2022/pdfs/p318.pdf).
6. F. Jardali, G. Svatek, R. S. Williams, T. D. Brown, and P. J. Shamberger, “Material-Driven Neuronal Oscillators and Filters via Active Reactance in CC-NDR and VC-NDR Electro-Thermal Memristors,” *Advanced Physics Research* **5**, e70138 (2026). [Primary full text](https://doi.org/10.1002/apxr.70138).
