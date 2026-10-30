# NIFTE: gain redundancy, linear modes, and local competition

Audit date: 30 October 2026 UTC. These are independent algebraic controls for the five-state NTP and nine-state NTP-TL equations transcribed in [the physical candidate screening](physical-candidate-screening.md). The primary equations on page 6 of Wang et al. [1] were checked again against the publisher PDF. No sign mismatch was found. This block contains no trajectory integration or newly established attractor.

The most useful consequences are an exact reduction of the two saturation parameters to one dynamical gain, independent characteristic-polynomial checks, and a concrete point where two oscillatory modes simultaneously become marginal. The local cubic calculation near that point predicts competition between periodic modes. Periodic-state multiplicity and even coexistence of two stable cycles were already reported in [1]; we do not claim those phenomena as new.

Reproduce the calculations with:

```bash
python research/audit_nifte_structure.py
```

The script needs NumPy and SciPy and reads the source parameter record in `results/candidate-audits/physical-candidate-screening.json`. Its output is [structural-audit.json](../results/nifte/structural-audit.json). Every eigenvalue, angular frequency, and period in this note uses **physical seconds**, not the source's dimensionless time `t/tau`.

## One dynamical gain, two amplitude parameters

Both exact models have the form

```math
\dot X=A X+bF\tanh(Kc^TX),\qquad c^TX=d,
```

where `b` has only its pressure entry nonzero and `A` includes all passive components, including all four RC pairs when present. For `F>0`, set `X=F Y`. Then

```math
\dot Y=A Y+b\tanh(Gc^TY),\qquad G=FK.
```

Thus, with component parameters fixed, two positive pairs `(F1,K1)` and `(F2,K2)` having the same product are exactly conjugate by `X2=(F2/F1) X1`. Initial states must also be scaled. Periods, Floquet multipliers, and asymptotic Lyapunov exponents are unchanged. A scan over both `F` and `K` duplicates dynamical regimes along each curve of constant `FK`.

Physical interpretation is more restrictive: amplitudes and reconstructed interface positions do change under this scaling. Changing the physical thermal drive may also change component properties in a more complete apparatus model. The conjugacy is exact for the stated equations with those components held fixed.

An alternative useful normalization is `xi=K X`, which gives `xi'=A xi+G b tanh(xi_d)`. The nominal gain is `G=1.37*330=452.1`.

## Independent transfer function and characteristic polynomial

Use physical time, let `a=Ld+Lf`, and eliminate the two flow variables from the dimensional pressure equations. Define the hydraulic polynomial matrix

```math
H(s)=\begin{pmatrix}
1+C_d(R_f+R_l)s+C_d(a+L_l)s^2 & C_pR_ls+C_pL_ls^2\\
C_dR_ls+C_dL_ls^2 & 1+C_pR_ls+C_p(L_p+L_l)s^2
\end{pmatrix},
\qquad J(s)=\det H(s).
```

Let

```math
B(s)=1+L_pC_ps^2,\qquad
C(s)=1+C_dR_fs+C_da s^2,
```

```math
Q(s)=(1+R_{th}C_{ad}s)J(s)
+R_{th}s[C_dB(s)+C_pC(s)].
```

The passive transfer from the imposed thermal-pressure input to the displacer pressure is exactly

```math
T_5(s)=-\frac{B(s)}{Q(s)}.
```

The origin has thermal feedback slope `G`, so its characteristic polynomial, up to a positive constant, is

```math
Q(s)+G B(s).
```

This is a separate calculation from assembling the five component ODEs. In particular, the numerator is quadratic rather than a constant. The passive model is Hurwitz by the energy argument below; since the zeros of `B` are purely imaginary, no pole-zero cancellation reduces the five-state transfer to a lower order.

For all four thermal-loss pairs define

```math
W(s)=\prod_{i=1}^4(1+R_iC_i s),\qquad
V(s)=R_0W(s)+\sum_{i=1}^4R_i\prod_{j\ne i}(1+R_jC_j s).
```

The loss admittance is `W/V`. The full nine-state transfer and origin characteristic polynomial become

```math
T_9(s)=-\frac{V(s)B(s)}{V(s)Q(s)+R_{th}W(s)J(s)},
```

```math
V(s)Q(s)+R_{th}W(s)J(s)+G V(s)B(s).
```

This retains all four RC states. At the nominal and reduced-`Lf` controls the hydraulic roots are nonreal, whereas the four roots of `V` are negative real numbers, so they do not cancel through the numerator. A six-state replacement is not an exact reduction of this fit.

## Linear imaginary-pair controls

Write `Q(s)=q0+q1*s+...+q5*s^5`. Since `B(i*omega)` is real, a nonzero imaginary pair for the five-state model requires

```math
q_1-q_3\omega^2+q_5\omega^4=0,\qquad
G=-\frac{q_0-q_2\omega^2+q_4\omega^4}
{1-L_pC_p\omega^2}.
```

Only positive `G` and `omega` are retained, and a zero denominator is not a valid finite-gain candidate. For the nine-state model the script solves the analogous real polynomial condition from the imaginary part of `(VQ+Rth WJ)*conjugate(V)` on the imaginary axis. Both constructions are checked against eigenvalues of the independently assembled component matrix.

| Model and `Lf` | Gain `G` at imaginary pair | Angular frequency, rad/s |
|---|---:|---:|
| 5D, nominal `4.74e6` | 32.23150 | 2.23600 |
| 5D, nominal `4.74e6` | 3191.27226 | 6.74378 |
| 5D, source control `4.74e5` | 210.30738 | 4.80736 |
| 5D, source control `4.74e5` | 370.30653 | 7.14020 |
| 9D, nominal `4.74e6` | 37.11052 | 2.15245 |
| 9D, source control `4.74e5` | 226.77732 | 3.71972 |
| 9D, source control `4.74e5` | 2560.57861 | 7.72235 |

The corresponding pairs cross toward positive real part as `G` increases. These are linear crossing controls. A nonlinear Hopf statement additionally needs its nondegeneracy conditions. At the second crossing another pair is already unstable where recorded in the JSON; it is not a second loss of stability of a previously stable equilibrium.

At the five-state source control, `G=452.1`, the origin has the spectrum

```math
-4.82296180,\quad
0.46667154\pm5.28802699i,\quad
0.17724691\pm7.19451774i.
```

There are two unstable oscillatory modes. This helps interpret the already published two-cycle behavior; it does not imply chaos. At the nine-state source control only one pair is unstable.

## A focused two-mode point and its cubic local dynamics

For a fifth-order real polynomial, two distinct imaginary pairs plus one negative real eigenvalue have the factorization

```math
Q(s)+G B(s)=q_5(s+\gamma)(s^2+\omega_1^2)(s^2+\omega_2^2).
```

Coefficient matching gives

```math
\gamma=\frac{q_4}{q_5},\qquad
G=\frac{q_4q_1}{q_5}-q_0,\qquad
q_2+G L_pC_p=\frac{q_4q_3}{q_5}.
```

Solving the last scalar equation for `Lf`, with all other source components nominal, gives

```math
L_f\simeq286255.183178673,\qquad G\simeq302.871676391795.
```

For `F=1.37`, this means `K≈221.074216344376`. The eigenvalues are

```math
-4.00681552,\qquad \pm5.17564246i,\qquad \pm7.34206644i.
```

The script checks the polynomial factorization and a nonsingular two-parameter unfolding in `(G,log Lf)`, with two finite-difference step sizes. The frequency ratio is approximately `1.41858`; the modes have neither a cubic 1:1 nor a cubic 1:3 resonance. This is a numerical double-Hopf candidate with additional local algebra, not an interval-certified bifurcation or a claim of priority.

The cubic calculation is especially transparent because the nonlinearity is odd. In `xi=K X`,

```math
\dot\xi=A_G\xi-\frac{G}{3}b\xi_d^3+O(\xi_d^5).
```

Normalize each critical right eigenvector to have displacer component one, and let its dual left eigenvector be normalized biorthogonally. The complex gain derivative is `alpha_j=d lambda_j/dG`. Projecting the resonant cubic monomials gives the radial amplitude equations

```math
\dot r_j=r_j\left(\mu_j+a_{jj}r_j^2+a_{jk}r_k^2\right)
+\text{higher-order terms},\qquad j\ne k,
```

```math
a_{jj}=-G\,\mathrm{Re}(\alpha_j),\qquad
a_{jk}=-2G\,\mathrm{Re}(\alpha_j).
```

There are no quadratic center-manifold terms contributing to this cubic order. Direct eigenvector projection and differentiation of the characteristic polynomial agree. Numerically, in the specified amplitude normalization,

```math
(a_{jk})\simeq
\begin{pmatrix}
-0.61960318&-1.23920636\\
-1.55007159&-0.77503579
\end{pmatrix}.
```

The determinant is negative. Therefore any mixed positive-amplitude equilibrium of this **leading cubic radial system** is a saddle: its radial Jacobian has determinant `4*r1^2*r2^2*det(a)<0`. The local approximation predicts competing periodic modes and permits bistability of pure-mode cycles; it does not predict an attracting mixed torus. This is consistent with the known coexistence of stable cycles in [1]. Higher-order effects and dynamics away from this small-amplitude neighborhood remain open. No global chaos exclusion follows.

## Symmetry, equilibrium, boundedness, and low-gain controls

All terms are odd, so `X -> -X` is an exact symmetry. An attractor may be centrally symmetric or have a distinct partner. If a single periodic orbit is invariant as a set under this symmetry, the antipodal operation corresponds to half a period. Two different symmetry-related invariant sets should not be counted as a new mechanism merely because their projections differ.

At equilibrium, both flows vanish and `d=q=-p`. For the five-state model the remaining equation is `F*tanh(Kd)+d=0`; it has only the root zero. For the loss model, the total series resistance is `Rsum=R0+sum(Ri)` and the equation becomes `F*tanh(Kd)+(1+Rth/Rsum)*d=0`, again with a unique root for positive components and positive `F,K`.

The passive boundedness argument in the screening note is sound. Its positive quadratic stored energy has nonpositive derivative when the thermal input is zero. The zero-dissipation set requires vapor pressure and both flows to vanish; invariance then forces both hydrostatic pressures to vanish. With thermal losses it also forces every loss-capacitor pressure to vanish. LaSalle's argument gives asymptotic stability, and finite-dimensional linearity makes the passive matrix Hurwitz. The actual thermal input is bounded by `F`, so variation of constants gives a global ultimate mathematical bound. This does not enforce cylinder-interface limits or all experimental operating assumptions.

There is also a reproducible conservative global low-gain control. For the passive dimensionless-state matrix in physical time, solve

```math
A^TP+PA=-I,\qquad P>0,
```

and define

```math
S=(Pb)c^T+c(Pb)^T.
```

For two trajectories the secant slope of `F*tanh(Kd)` lies in `[0,G]`. Their quadratic separation therefore decays if `I-G*S` is positive definite; a sufficient strict bound is `G<1/lambda_max(S)`. The script stores the numerical matrix `P`, equation residual, and a control at half that gain bound. The resulting control gains are approximately `0.03577, 0.05423, 0.74610, 0.77445` for 5D nominal, 9D nominal, 5D reduced-`Lf`, and 9D reduced-`Lf`, respectively. These are conservative numerically constructed certificates, not sharp thresholds or interval-certified proofs. The analytical implication of a valid certificate is global incremental contraction to the unique equilibrium.

The rank-one smooth nonlinearity places both models in the general Lur'e feedback class. It does not reduce their state dimension or identify them with a particular previously known chaotic oscillator. The reciprocal pressure/flow derivatives have opposite signs, so the elementary orthant sign-flip argument for a cooperative or competitive system is unavailable. No stronger monotonicity theorem or global no-chaos result was established here.

## Verification and source precision

For all four parameter/model controls, the independent transfer polynomials and component matrices agree to approximately `4e-15` in the stored scaled coefficient and frequency checks. Odd symmetry is exact in the numerical check; amplitude-scaling residuals are at floating-point roundoff. The passive quadratic certificates have positive margins and small equation residuals. These are algebraic transcription checks, not attractor validation.

The source computes cycles using CRONOS single shooting, IPOPT and CVODES/SUNDIALS, and obtains monodromy matrices through CVODES sensitivities. It does not state an integration tolerance, shooting tolerance, or higher-precision replacement for the rounded Table 2 parameters in the inspected text. The explicit `1e-3` tolerance applies to a separate outer Nelder–Mead design optimization. Consequently, a small mismatch with the printed Floquet multipliers should remain documented as unresolved; parameter rounding or source solver tolerance are possibilities, not demonstrated explanations. Table 3 itself gives nine multipliers for the thermal-loss control.

Strongly contracting nine-dimensional modes also require care: a determinant formed from a double-precision monodromy matrix may be meaningless when its smallest singular values fall below numerical resolution. The constant trace and Liouville identity are the appropriate independent reference; an inaccurate direct determinant must not be reported as a successful spectrum check.

## References

1. Y. Wang, C. N. Markides, and B. Chachuat, “Optimization-based investigations of a two-phase thermofluidic oscillator for low-grade heat conversion,” *BMC Chemical Engineering* **1**, 20 (2019). [DOI and primary article](https://doi.org/10.1186/s42480-019-0019-3); [publisher PDF](https://link.springer.com/content/pdf/10.1186/s42480-019-0019-3.pdf). Equations 7–8, Tables 2–3, computational methods and the published cycle-multiplicity section were inspected. The article is CC BY 4.0.
2. C. N. Markides, A. Osuolale, R. Solanki, and G.-B. V. Stan, “Nonlinear heat transfer processes in a two-phase thermofluidic oscillator,” *Applied Energy* **104**, 958–977 (2013). [DOI](https://doi.org/10.1016/j.apenergy.2012.11.056). Original NTP-model lineage; the precise implementation target and parameter source for this audit are [1].
