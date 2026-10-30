# A thermal spring pair and a valve lead

Audit date: 30 October 2026, UTC. This pass recovered **one implementable four-state published reduced model**, plus one valve model that is not ready for faithful implementation. No trajectories were integrated and no new chaotic attractor was found in this audit. The purpose is to narrow the next calculation, not to assert that the model has never been studied for chaos.

Machine-readable equations, source access limits, search scope, and algebra checks are in [the audit record](../results/candidate-audits/late-physical-candidates-2.json).

## 1. Thermally powered coupled LCE springs

Li et al. study two liquid crystal elastomer spring oscillators coupled by another spring in a maintained temperature gradient. Their published Equation 14 is a pair of second-order autonomous equations. The thermal field has already been eliminated: the implemented reduced model has **four independent states**, not six. The reported examples include synchronization and decay. [Primary article and full MathML source](https://www.mdpi.com/2073-4360/15/16/3349), [repository XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10458843/fullTextXML).

### Exact reduction of the printed equations

Below, bars on the paper's dimensionless variables are omitted; the symbol `k` denotes its dimensionless coupling stiffness. Define

```math
s=\alpha\beta,\qquad
A=\frac{Ks}{e^s-1},\qquad
B=\frac{Ks\tau(1-e^s+se^s)}{(e^s-1)^2},\qquad
d=a_0+B.
```

At zero `s`, use the continuous limits `A=K` and `B=0`. Numerically, `expm1(s)` and a small-argument series avoid cancellation. The source's damping law contains an absolute value: `(a0+a1|v|)v`. This was checked in the MathML of Equations 4 and 14; flattening XML text alone loses its vertical bars.

Using the paper's Equation 2 to eliminate the second prestretch gives the common constant term and equilibrium:

```math
C=(A-K)\lambda_1-K\alpha Q,\qquad
w_e=-\frac{C}{A+2k}.
```

With `qi=wi-we`, Equation 14 becomes exactly

```math
\dot q_i=v_i,\qquad
\dot v_i=-Aq_i-k(q_1+q_2)-dv_i-a_1|v_i|v_i,
\qquad i=1,2.
```

This algebraic simplification is our analysis of the published equation; priority for the simplification has not been established. The mathematical assumptions used below are `A>0`, `k>=0`, and `a1>0`.

For the active case `d<0`, the further rescaling

```math
\sigma=\sqrt A\,t,\qquad
q_i=\frac{-d}{a_1\sqrt A}x_i,\qquad
\mu=\frac{-d}{\sqrt A},\qquad \kappa=\frac{k}{A}
```

leaves only **two independent dynamical parameters**:

```math
x_i''=-x_i-\kappa(x_1+x_2)+\mu(1-|x_i'|)x_i',
\qquad i=1,2.
```

Here primes mean differentiation with respect to `sigma`. Reflecting the second displacement and velocity changes the sum coupling to the usual elastic difference coupling. Thus this is a generalized Rayleigh oscillator pair with quadratic velocity drag. Changing `a1` changes the displacement scale; changing the prestretch or `Q` translates the equilibrium. A grid over all those parameters would repeatedly sample equivalent dimensionless dynamics, although their physical admissibility can differ.

### A source-grounded starting point

The source Figure 2 uses the following dynamic parameters; its velocity seeds include `(0.1,0.5)` and `(0.1,-0.5)`.

| Parameter | Value |
|---|---:|
| `K` | 14 |
| `k` | 8 |
| `alpha` | -0.5 |
| `beta` | 0.1 |
| `tau` | 0.15 |
| `a0` | 0.02 |
| `a1` | 0.2 |

For this point, the derived values are `A=14.3529165451`, `d=-0.03337492709`, **`mu=0.00880947943`**, and **`kappa=0.557378006`**. The displacement scale is `0.04404739715`. These are coefficient conversions, not simulated findings. The inspected Figure 2 description does not fully specify the positions, prestretch, and thermal offset needed to reconstruct an absolute physical trajectory; these must not be invented and labelled as the author's control.

### Analytic exclusions and useful controls

For the shifted system define

```math
E=\frac{v_1^2+v_2^2}{2}
 +\frac{A(q_1^2+q_2^2)}{2}
 +\frac{k(q_1+q_2)^2}{2}.
```

Direct differentiation gives

```math
\dot E=-d(v_1^2+v_2^2)-a_1(|v_1|^3+|v_2|^3).
```

When `d>=0` and `a1>0`, this radially unbounded energy is nonincreasing. The only invariant trajectory with both velocities zero is the sole equilibrium. LaSalle's principle therefore gives convergence to that equilibrium: this parameter region cannot contain a bounded nonstationary attractor.

When `d<0`, the bound

```math
\dot E\leq\frac{8(-d)^3}{27a_1^2}
```

rules out finite-time energy explosion. It does **not** prove a uniformly bounded absorbing region or a chaotic attractor. A subsequent [structural audit](rayleigh-pair-structure.md) supplies an explicit absorbing-set proof and seven independent modal Floquet controls; it still does not establish chaos or a global active-regime exclusion.

The planes `(q2,v2)=(q1,v1)` and `(q2,v2)=(-q1,-v1)` are invariant. Each carries a planar autonomous oscillator. A bounded strange attractor cannot be wholly confined to either plane; a collective-chaos candidate must persist off both. These planes are useful implementation controls and starting points for transverse-stability analysis.

The canonical-vector-field, energy-derivative, and equilibrium-shift identities were checked at 256 seeded random states. Maximum residuals were respectively `8.9e-16`, `1.7e-16`, and `1.4e-15`. These are algebra checks only; no integration or GPU calculation occurred.

### Physical and novelty limits

The temperature gradient supplies continuing external energy, so this is not a closed finite-fuel oscillator. Physical use still requires positive fiber lengths, admissible prestretch, and tensile forces in the unshifted coordinates. A compressed or slack fiber cannot be accepted merely because the mathematical ODE continues smoothly.

The force relation traces to a small-thermal-time asymptotic constitutive approximation in Zhou et al. The earlier paper distinguishes that approximation from its fuller constitutive model. At the source point, `tau*sqrt(A)=0.568`; our algebra alone supplies no small-error certificate. Large thermal times or fast motion need validation against the full thermal model. Also, the pair paper's Table 2 lists a positive `alpha` range while its examples use `alpha=-0.5`; a scan must retain the explicitly intended contraction sign. [Primary constitutive-model preview](https://www.sciencedirect.com/science/article/abs/pii/S0020740322002971).

The bounded literature check found these substantial limitations on novelty:

- A 2024 author follow-up studies three coupled thermal LCE oscillators and reports full, partial, and asynchronous behavior. Its complete equations were not accessible in this pass. [Wu et al., Physical Review E 109, 024701](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.109.024701).
- A September 2026 paper studies nonidentical light-fueled coupled LCE oscillators. Metadata describes phase locking, Arnold tongues, and bistability; full equations were not recovered. Exact equivalence and its treatment of irregular states remain open checks. [Dai et al., DOI 10.1631/jzus.a2600177](https://doi.org/10.1631/jzus.a2600177).
- Chaos in coupled Rayleigh oscillators is established. Sprott and Hoover's Equation 18 uses cubic drag, nonidentical frequencies, and multiplicative coupling, so it is not an established equivalence to this quadratic-drag elastically coupled pair. It nevertheless rules out a broad claim of discovering the mechanism. [Primary author-hosted paper, Section 4.4](https://sprott.physics.wisc.edu/pubs/paper492.pdf).

The audit screened 30 forward-citation records for the three-oscillator follow-up and 55 for the constitutive paper, mostly at title/identifier level. Several concern explicitly chaotic LCE pendulums or other actuators. The pair paper's zero-citation index result was demonstrably insufficient because independent author/title searches found follow-ups. **No exact Equation 14 chaos demonstration was verified; absence of such a result has not been established.** This is a modest-priority, precise mathematical target, not a new-attractor claim.

## 2. Pressure servo valve: hold

Huang et al. provide a detailed pressure-servo-valve model and discuss self-excited oscillations. The relief subsystem alone appears to require pressures, valve displacement, and velocity; the full pilot/cylinder/feedback system adds dynamics. This pass did not recover and audit all port, contact, and parameter closures, so it does not certify the complete independent-state count or supply an executable surrogate. [Primary article](https://www.mdpi.com/2076-3417/12/17/8477).

A source consistency issue also needs resolution. Equations 15–16 give `ph_dot=(beta/Vh)*(Qin-kj*ph)`. For constant inflow and zero initial pressure its exponential step response is a formula for `ph`. The source's Equation 17 labels that response as `ph_dot`; its right side has pressure units rather than pressure/time. This is a reason to verify the original model carefully, not permission to silently repair all its conventions.

Hydraulic-valve chaos is already documented by Hayashi et al. (1997) and the pressure-relief-valve literature reviewed by Hős and Champneys (2012). That excludes a broad first-valve-chaos claim, while leaving exact topology distinctions unresolved. Positive cavity volumes, admissible pressure drops, spool limits, and seat-contact behavior must all be specified before simulations could test the Huang system. **Hold this lead; do not create missing physical closures.** [Hayashi et al.](https://doi.org/10.1006/jfls.1997.0096), [Hős and Champneys primary preview](https://www.sciencedirect.com/science/article/abs/pii/S0167278911001229).

## Follow-up: narrower equivalences and the 2026 supplementary equations

This follow-up supersedes the initial metadata-only access statement for the 2026 paper. Its primary abstract and references were checked, and its **official, publicly linked seven-page supplementary PDF** was recovered. The main article remains subscription-only; it was not purchased or bypassed. [Publisher landing page](https://link.springer.com/article/10.1631/jzus.A2600177), [official supplement](https://media.springernature.com/original/springer-static/esm/art:10.1631%2Fjzus.A2600177/MediaObjects/11582_2026_177_MOESM1_ESM.pdf).

### The coupling is conventional after reflection

Put `z1=x1`, `z2=-x2`. Oddness of the velocity force gives, with no approximation or time change,

```math
z_1''=-z_1-\kappa(z_1-z_2)+\mu(1-|z_1'|)z_1',
\qquad
z_2''=-z_2-\kappa(z_2-z_1)+\mu(1-|z_2'|)z_2'.
```

The potential becomes `(z1^2+z2^2)/2+kappa*(z1-z2)^2/2`. The original sum coupling therefore cannot be used to claim an unconventional interaction mechanism.

There is also an exact state-space equivalence to **elastically coupled nonsmooth van der Pol oscillators**. Let `M` have diagonal entries `1+kappa` and off-diagonal entries `-kappa`, and let `f(v)=mu*(1-|v|)*v` act componentwise. The reflected pair is `z''=-Mz+f(z')`. Setting `u=2z'` and differentiating yields

```math
u_i''+\mu(|u_i|-1)u_i'+u_i+\kappa(u_i-u_j)=0,
\qquad j\ne i.
```

This does not add two spurious states. From `(u,b)`, where `b=u'`, the original state is recovered by

```math
z'=u/2,\qquad z=M^{-1}\bigl(f(u/2)-b/2\bigr).
```

For `kappa>=0`, `det(M)=1+2*kappa>0`. The forward map `(z,z')->(u,u')` and its inverse are continuously differentiable, including at zero velocity; its Jacobian determinant is `16*det(M)`. It is therefore an invertible coordinate description of the same flow. A striking change in projection under this transformation would not establish a new dynamical mechanism.

The scalar equation `u''+mu*(|u|-1)*u'+u=0` is already an established nonsmooth van der Pol model. Buică, Llibre, and Makarenkov explicitly give it in Section 3 and credit Hogan's 2003 limit-cycle analysis; their own work treats periodic forcing and periodic-solution stability. This is exact prior art for the **uncoupled oscillator law**, not a verified chaos result for the present autonomous coupled pair. [Primary author-hosted paper, Section 3](https://www.math.ubbcluj.ro/~abuica/BuLlMa-SIAM-2009.pdf).

These identities were checked for 256 seeded random states, including zero velocities, over `mu` from `0.001` to `10` and `kappa` from `0.001` to `100`. Maximum reflection and inverse-map residuals were `5.7e-14` and `1.9e-14`; the transformed-vector-field identity agreed to the arithmetic precision of this check. No trajectories were integrated.

### The 2026 paper is a much closer precedent than its title alone suggested

The supplement explicitly expands the thermal field for small `tau` in S1–S8. Visually inspected Equation S13 contains the first-order perturbation terms

```math
\delta_i\bigl(D_0y_{i0}-c_i|D_0y_{i0}|D_0y_{i0}\bigr)
   +k(y_{j0}-y_{i0}).
```

Thus its reduction uses the same quadratic velocity drag and elastic difference coupling. S19–S21 derive amplitude and phase equations, and S4/Figure S1 discuss coexistence of in-phase and anti-phase locking and loss of synchronization. No chaos claim was found in these seven pages. This directly identifies a close nonidentical generalization of the same dynamical family. However, S13 is a perturbation equation: without main-text Equation 10 and its parameter definitions, it does not certify every coefficient or the entire model's equivalence. “Asynchronous” in the abstract or supplement does not mean chaotic.

The remaining exact-model priority gap is consequently narrower: any candidate must be compared against this nonidentical-pair analysis, including its identical limit, before describing it as previously unreported.

### Other older matches that do not establish exact equivalence

| Primary source inspected | Relevant distinction |
|---|---|
| Maeda and Yokomori (1999), Equation 7 | A pair of Rayleigh oscillators with elastic and viscous coupling, cubic velocity friction, a moving-belt offset, and potentially different frequencies. Period-doubling evidence is reported, but this is not the exact quadratic-drag identical pair. |
| Buică, Llibre, and Makarenkov (2009), Section 3 | Exact scalar absolute-value damping law after the transformation above; the analyzed forcing and periodic solutions do not supply an autonomous pair chaos result. |
| Li et al. (2023), *Light-Fueled Synchronization of Two Coupled Liquid Crystal Elastomer Self-Oscillators* | The primary description specifies rotating LCE/passive bars linked by a torsion spring and a bounded illuminated angular region. Its similar title does not make it the same translational linear-gradient pair. |
| Pandey (2026), *Damping-induced amplitude gain in nonlinearly damped coupled oscillators* | The primary abstract specifies Duffing oscillators with external forcing on one oscillator. It is not evidence for the target autonomous pair. |

An additional 35 forward-citation records for the 2009 nonsmooth-oscillator paper were screened at title/identifier level. This did not verify a chaos result for the exact elastically coupled pair; it is not exhaustive full-text clearance. The source mapping supports a precise implementation target but lowers its priority as a candidate for an unusual new mechanism.

## References

1. K. Li, H. Wu, B. Zhang, Y. Dai, and Y. Yu (2023). Heat-Driven Synchronization in Coupled Liquid Crystal Elastomer Spring Self-Oscillators. *Polymers* **15**, 3349. [doi:10.3390/polym15163349](https://doi.org/10.3390/polym15163349).
2. L. Zhou, C. Du, W. Wang, and K. Li (2022). A thermally-responsive fiber engine in a linear temperature field. *International Journal of Mechanical Sciences* **225**, 107391. [doi:10.1016/j.ijmecsci.2022.107391](https://doi.org/10.1016/j.ijmecsci.2022.107391).
3. H. Wu, B. Zhang, and K. Li (2024). Synchronous behaviors of three coupled liquid crystal elastomer-based spring oscillators under linear temperature fields. *Physical Review E* **109**, 024701. [doi:10.1103/PhysRevE.109.024701](https://doi.org/10.1103/PhysRevE.109.024701).
4. Y. Dai, K. Wang, X. Jiang, and K. Li (2026). Phase-locking in nonidentical coupled light-fueled liquid crystal elastomer self-oscillators. *Journal of Zhejiang University-SCIENCE A*. [doi:10.1631/jzus.a2600177](https://doi.org/10.1631/jzus.a2600177). Metadata-only audit.
5. J. C. Sprott and W. G. Hoover (2017). Harmonic Oscillators with Nonlinear Damping. *International Journal of Bifurcation and Chaos* **27**, 1730037. [doi:10.1142/S0218127417300373](https://doi.org/10.1142/S0218127417300373).
6. J. Huang, Q. Zhang, F. Zhao, X. Liu, and T. Wang (2022). Analysis and Suppression of Self-Excited Oscillations in Pressure Servo Valve System. *Applied Sciences* **12**, 8477. [doi:10.3390/app12178477](https://doi.org/10.3390/app12178477).
7. S. Hayashi, T. Hayase, and T. Kurahashi (1997). Chaos in a hydraulic control valve. *Journal of Fluids and Structures*. [doi:10.1006/jfls.1997.0096](https://doi.org/10.1006/jfls.1997.0096).
8. C. Hős and A. R. Champneys (2012). Grazing bifurcations and chatter in a pressure relief valve model. *Physica D*. [doi:10.1016/j.physd.2011.05.013](https://doi.org/10.1016/j.physd.2011.05.013).
9. A. Buică, J. Llibre, and O. Makarenkov (2009). Asymptotic Stability of Periodic Solutions for Nonsmooth Differential Equations with Application to the Nonsmooth van der Pol Oscillator. *SIAM Journal on Mathematical Analysis* **40**, 2478–2495. [doi:10.1137/070701091](https://doi.org/10.1137/070701091).
10. S. J. Hogan (2003). Relaxation oscillations in a system with a piecewise smooth drag coefficient. *Journal of Sound and Vibration* **263**, 467–471. [doi:10.1016/S0022-460X(02)01431-1](https://doi.org/10.1016/S0022-460X(02)01431-1). Identified and its equation verified through the primary 2009 paper; not an independently completed full-text audit.
11. I. Maeda and M. Yokomori (1999). Systems of Multiply Connected Nonlinear Oscillators and Their Application to Seismic Phenomena. *Journal of the Faculty of Science, Hokkaido University, Series 7, Geophysics* **11**, 633–655. [Institutional primary manuscript](https://hdl.handle.net/2115/8853).
12. K. Li, B. Zhang, Q. Cheng, Y. Dai, and Y. Yu (2023). Light-Fueled Synchronization of Two Coupled Liquid Crystal Elastomer Self-Oscillators. *Polymers* **15**, 2886. [doi:10.3390/polym15132886](https://doi.org/10.3390/polym15132886).
13. A. Pandey (2026). Damping-induced amplitude gain in nonlinearly damped coupled oscillators. *Nonlinear Dynamics* **114**, article 1. [Primary abstract](https://link.springer.com/article/10.1007/s11071-025-11869-8). Online publication: 3 December 2025.
