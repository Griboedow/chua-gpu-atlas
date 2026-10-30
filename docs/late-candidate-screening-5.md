# Late screening V: photothermal Marangoni oscillators

Audit: 30 October 2026 UTC. **A physically distinct three-state oscillator and six-state pair are identified, but this block does not deliver a simulation-ready source control.** The retrieved 2021 main paper and supplement omit numerical values for the fitted absorption, effective mass, drag, and coupling coefficients. No trajectories were integrated and no chaos or priority claim is made.

## Physical target and exact scope

Kim et al. [1] place photothermal hydrogel particles at an air–water interface in stationary patterned illumination. Unequal heating of two particle halves changes surface tension; finite thermal relaxation delays the force relative to displacement. The source has three independent states per translating oscillator: position, velocity, and temperature difference. A pair has six states. Maintained illumination supplies energy and the surrounding water removes heat.

The experimental rings are particularly interesting for geometry: three particles exhibit frustration and switching phase relationships; four particles form changing phase clusters. These observations are not a diagnosis of deterministic chaos. The article gives the single-particle and pair equations, but not a complete, calibrated vector model for its frustrated ring experiments. Extending the pair formula to a nine- or twelve-state network would need to be labeled an additional modeling assumption. Rotating an existing attractor or presenting a three-dimensional projection would not establish geometric novelty.

## Source equations

Use the signs in the source's numerical equations 5–8. Write its temperature difference as $T$ and disk thickness as $t_d$ to avoid confusing thickness with time:

```math
\dot x=v,\qquad
\dot v=\frac{\gamma_T aT-\eta_{\rm eff}v}{m_{\rm eff}},\qquad
\dot T=\frac{\Delta Q(x)-T/R_{\rm th}(v)}{C_{\rm th}(v)}.
```

Source equations 2–4 define

```math
h(v)=\frac{2k_w}{L}+0.664\,\mathrm{Pr}^{1/3}k_w
\sqrt{\frac{\rho_w|v|}{a\mu}},\qquad
\mathrm{Pr}=\frac{c_p\mu}{k_w},
```

```math
R_{\rm th}(v)=\frac{1}{h(v)A},\qquad
C_{\rm th}(v)=c_p\rho_m A t_d
\left(1+\frac{h(v)t_d}{3k_w}\right),\qquad A=\frac{a^2}{2},
```

```math
m_{\rm eff}=\Lambda_m m,\qquad
\eta_{\rm eff}=\Lambda_{\rm drag}\,6\pi\mu\left(\frac{0.566a}{2}\right).
```

Here the source approximates a disk by a square of width $a$, with each lumped element covering half its area. It specifies $L\sim a$, rather than a numerical choice of $L$. The heat input $\Delta Q(x)$ is obtained from overlap of the illuminated pattern with the two halves. The local approximation $\Delta Q=2I_0\alpha W x$ is stated only near the trap center; extending it to arbitrary displacement would change the model's confinement.

For the pair, source equation 8 replaces the acceleration by

```math
\dot v_i=\frac{\gamma_T aT_i-\eta_{{\rm eff},i}v_i
\mathbin{\pm}\Lambda_{\rm th}/r^2}{m_{{\rm eff},i}}.
```

The source calls $r$ the separation between particles and assigns plus to the left particle and minus to the right. A common-coordinate implementation must explicitly resolve its axis and temperature conventions before interpreting these as repulsive forces. Trap-center spacing and instantaneous particle separation must not be silently interchanged. The separate small-amplitude derivation writes equation 10 with the opposite algebraic sign of the thermal force from equation 6; the numerical equations and analytical threshold calculation therefore require a sign-convention audit. This is an implementation ambiguity, not evidence that the experimental mechanism is wrong.

## Parameters actually supplied

Table S1 of [1], visually checked in the author-hosted PDF, supplies the following physical constants:

| Quantity | Value |
|---|---:|
| Specific heat $c_p$ | 4180 J/(kg K) |
| Dynamic viscosity $\mu$ | $0.890\times10^{-3}$ N s/m² |
| Signed temperature derivative $\gamma_T$ | $-0.1379\times10^{-3}$ N/(m K) |
| Water and gel thermal conductivity | 0.6 W/(m K) each |
| Water and gel density | 997 kg/m³ each |
| Model thickness $t_d$ | 10 µm |
| Model width $a$ | 300 µm |

The same table names $\alpha$, $\Lambda_{\rm drag}$, $\Lambda_m$, and $\Lambda_{\rm th}$ **without values**. The main text explicitly says these are fitted. The value $\alpha=0.5$ occurs in a simplified threshold example and is not identified as the absorption used for the numerical controls. Neither a documented global heat-overlap function nor a complete fitted control tuple was recovered.

Useful future controls are Fig. 1E (2.8 W/cm² and a 140 µm trap), Fig. 1G (position–temperature cycle), and Fig. 2E (two particles at 510 µm trap spacing). The source adjusts the two effective masses to produce approximately 9% uncoupled frequency mismatch. Matching only the quoted light intensity and dimensions would therefore not reproduce its calculation. The paper also acknowledges a mismatch in the onset intensity and omitted short-range interactions.

## Low-cost mathematical checks

For the isolated oscillator, suppose the supplied illumination closure is bounded, $|\Delta Q(x)|\le Q_{\max}$, and all dissipative parameters are positive. Since $h(v)\ge h_0=2k_w/L$, comparison at the temperature boundary gives

```math
|T(t)|\le T_{\max}=\max\left\{|T(0)|,\frac{Q_{\max}}{h_0A}\right\}.
```

The linear drag then gives

```math
|v(t)|\le\max\left\{|v(0)|,
\frac{|\gamma_T|aT_{\max}}{\eta_{\rm eff}}\right\}.
```

These are our conditional bounds, not source claims. They do **not** prove bounded position, persistence inside the trap, collision avoidance for the pair, or existence of an attractor. Pairwise inverse-square forcing additionally requires a positive separation bound. Escape and collision regions cannot be counted as bounded chaotic candidates.

The velocity dependence through $\sqrt{|v|}$ is generally not differentiable at zero velocity. Any variational/Lyapunov implementation must handle this explicitly; smoothing it changes the literal model. Its extra velocity-to-temperature coupling and the global illumination geometry also prevent a direct application of a simple monotone cyclic-feedback exclusion. Three states alone neither imply chaos nor rule it out.

## Exact nondimensionalization and relation to known oscillator classes

The absence of a calibrated source tuple does not forbid a clearly labeled mathematical study of the equation family. It does prevent claiming experimental prediction or source-control reproduction. The remaining global heat-input closure and coordinate conventions must still be specified explicitly rather than guessed.

For any chosen positive heat scale $Q_*$, define

```math
h_0=\frac{2k_w}{L},\quad b=\frac{h_0t_d}{3k_w},\quad
C_0=c_p\rho_mAt_d(1+b),\quad
\tau_0=\frac{C_0}{h_0A},\quad T_* =\frac{Q_*}{h_0A}.
```

Use $s=t/\tau_0$, $X=x/a$, $V=\tau_0v/a$, $\Theta=T/T_{\ast}$, and $q(X)=\Delta Q(aX)/Q_{\ast}$. With primes denoting differentiation in $s$, equations 5–7 are exactly

```math
X'=V,\qquad V'=\Gamma\Theta-\delta V,\qquad
\Theta'=F(V)\big[q(X)-H(V)\Theta\big],
```

```math
H(V)=1+\kappa\sqrt{|V|},\qquad
F(V)=\frac{1+b}{1+bH(V)},\qquad
\kappa=\frac{0.664\,\mathrm{Pr}^{1/3}k_w}{h_0}
\sqrt{\frac{\rho_w}{\mu\tau_0}},
```

```math
\delta=\frac{\eta_{\rm eff}\tau_0}{m_{\rm eff}},\qquad
\Gamma=\frac{\gamma_T T_*\tau_0^2}{m_{\rm eff}}.
```

The sign of $\Gamma$ retains the signed $\gamma_T$ in the numerical source equations. In particular, equation 6 with Table S1 gives $\Gamma\lt0$. For a pair with common geometric/thermal scales, the additional acceleration is $\mathbin{\pm}K_i/R^2$, where $R=r/a$ and $K_i=\Lambda_{\rm th}\tau_0^2/(m_{{\rm eff},i}a^3)$.

This separates parameter roles. Geometry and fluid properties determine $b$ and $\kappa$ jointly; they are not independent experimental knobs. Effective mass changes $\delta$, $\Gamma$, and $K$ simultaneously; effective drag changes $\delta$; illumination and absorption determine the heat scale and may also affect pair coupling. A mathematical family may vary these dimensionless groups independently, but that larger parameter domain must not be presented as independently realizable controls. The source fitted coefficients do not establish a physical parameter range for such a survey. The shape of $q$ also depends on the fixed illumination and trap geometry, not just one amplitude.

For $\Gamma\ne0$, elimination of $\Theta$ is an exact linear change from $(X,V,\Theta)$ to $(X,X',X'')$. It gives the scalar jerk equation

```math
X'''+[\delta+\lambda(X')]X''+
\delta\lambda(X')X'-\Gamma F(X')q(X)=0,\qquad
\lambda(V)=F(V)H(V).
```

Thus the model belongs to the broad, well-established autonomous third-order feedback/jerk class, where chaotic examples have long been known [6]. This derivation neither identifies it with a particular known chaotic jerk model nor establishes novelty within that class. In the **changed-model limit** of velocity-independent heat exchange ($\kappa=0$), it reduces to $X'''+(\delta+1)X''+\delta X'-\Gamma q(X)=0$. Dropping the thermal state instead would be another approximation. Neither operation is an exact reduction of the published velocity-dependent system to a planar Rayleigh, van der Pol, or unforced Duffing oscillator. The thermal-memory mechanism and this jerk representation are more precise comparisons than an unsupported visual-equivalence claim.

## Closest literature and priority limits

The direct predecessor is Hauser, Sundaram, and Hayward's 2018 *Photothermocapillary Oscillators* [2]. Its three-state thermal-mechanical model includes a capillary restoring force from a curved droplet and a Gaussian illumination profile. The main paper already discusses irregular motion near a possible transition to chaos, without presenting a Lyapunov confirmation. Thus, even a later positive result would not be the first suggestion of chaos in this physical family. Its cooling law lacks the 2021 zero-velocity diffusive term, and its fitted constants and geometry cannot be transplanted as the 2021 parameter tuple. A resemblance to a Duffing oscillator is not an established coordinate equivalence.

A forward-citation metadata query returned 57 records for the 2021 paper. The closest later primary article examined was Tanjeem et al. (2024) [3], which uses larger traps and lower illumination to keep particles static for logic operations. It also studies particle escape caused by collective repulsion. This reinforces the need to test confinement, but does not establish chaos in the six-state oscillator model. Its supplementary material may contain relevant dynamical calibration; full official supplement retrieval was unsuccessful (HTTP 403), so that route remains unaudited. Search snippets are not treated as a substitute for inspecting its equations and values.

Searches for the exact title and Marangoni optical traps combined with chaos/chaotic/Lyapunov did not identify a verified chaos calculation for the literal 2021 model. This limited retrieval is **not proof of novelty**. The 57 citing works were screened by metadata, not all read in full. Other thermocapillary, active-droplet, and optothermal oscillators cannot be assumed equivalent or irrelevant without comparing their equations.

## Additional exclusions in this search block

- Marchegiani, Braggio, and Giazotto (2020) [4] provide an autonomous three-state thermoelectric SQUID: junction voltage, load current, and Josephson phase. Its primary paper already labels an intermediate-load regime chaotic in Fig. 3. The abstract emphasizes oscillations, so an abstract-only audit would have missed this prior result. It is excluded as a first-chaos target; the figure alone does not certify a permanent strange attractor.
- Sharp and Kalogirou (2025) [5] report regular motion of fluidically levitated bubbles, drops, and spheres. The retrieved primary treatment supplies force-balance and frequency-scaling relations, not a closed dynamic wake/particle state model. Appending a wake oscillator would be a new closure, so it is not ready for exact-model screening here.

## Decision

Retain the Marangoni pair as a **source-gap lead**, not a GPU-ready candidate. Its exact nondimensional equation family is available for a later mathematical study once the global heat-input closure and conventions are resolved; missing fit values alone are not a mathematical prohibition. The next useful step is recovering a complete, consistent control specification from public materials, including the 2024 supplement if accessible, then reproducing the isolated cycle and two-particle synchronization. No parameter fitting, arbitrary ring extension, regularization, simulation, or new-attractor claim was performed in this block. Source PDFs remain in the ignored research cache; repository MIT licensing does not relicense those papers.

## References

1. H. Kim et al. *Coupled oscillation and spinning of photothermal particles in Marangoni optical traps*. PNAS **118**, e2024581118 (2021). [DOI](https://doi.org/10.1073/pnas.2024581118), [author-hosted main PDF](https://subras.github.io/papers/PNAS_21.pdf), [author-hosted supplement](https://subras.github.io/papers/PNAS_SI_21.pdf).
2. A. W. Hauser, S. Sundaram, R. C. Hayward. *Photothermocapillary Oscillators*. Phys. Rev. Lett. **121**, 158001 (2018). [DOI](https://doi.org/10.1103/PhysRevLett.121.158001), [author-hosted main PDF](https://subras.github.io/papers/PRL_18.pdf), [supplement](https://subras.github.io/papers/PRL_PTCO_SI.pdf).
3. N. Tanjeem, K. M. Kreienbrink, R. C. Hayward. *Modulating photothermocapillary interactions for logic operations at the air–water interface*. Soft Matter **20**, 1689–1693 (2024). [Primary full text](https://pubs.rsc.org/en/content/articlehtml/2024/sm/d3sm01487h).
4. G. Marchegiani, A. Braggio, F. Giazotto. *Phase-tunable thermoelectricity in a Josephson junction*. Phys. Rev. Research **2**, 043091 (2020). [DOI](https://doi.org/10.1103/PhysRevResearch.2.043091), [author manuscript](https://arxiv.org/abs/2006.10631).
5. J. S. Sharp, A. Kalogirou. *Fluidic levitation of bubbles, drops and solid spheres*. J. Fluid Mech. **1023**, A17 (2025). [DOI](https://doi.org/10.1017/jfm.2025.10837).

6. J. C. Sprott. *Some Simple Chaotic Jerk Functions*. Am. J. Phys. **65**, 537–543 (1997). [DOI](https://doi.org/10.1119/1.18585), [author-hosted primary](https://sprott.physics.wisc.edu/pubs/PAPER229.HTM).
