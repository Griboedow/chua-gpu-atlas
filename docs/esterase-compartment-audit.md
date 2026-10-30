# Esterase compartment: source and structural audit

Audit date: 30 October 2026. Candidate: Leathard, Beales and Taylor, *Chaos* **33**, 123128 (2023), [doi:10.1063/5.0180256](https://doi.org/10.1063/5.0180256).

**Status: hold numerical screening until the supplementary model is recovered.** The published article specifies the reaction network and a constant external reservoir, but does not fully enumerate the transport term and baseline parameter of every species. Supplement S1 is explicitly identified as the authoritative XPPAUT model. We did not obtain it in this audit. No author code was executed, and no trajectory or parameter scan was run.

This is a plausible narrowly defined literature lead, not a new-chaos result. There is also an important mathematical gate: if either total product is retained without removal, its monotone accumulation rules out a bounded recurrent attractor with sustained reaction. The main diagram alone is insufficient to decide whether that condition holds in S1.

## What was recovered

The [institutional published PDF](https://eprints.whiterose.ac.uk/207599/1/123128_1_5.0180256.pdf) was read and its model pages visually checked. The four reactions are hydrolysis of ethyl acetate into acetic acid and ethanol, dissociation of each product, and water dissociation. The article reports pH oscillations and bistability, using a pH-dependent Michaelis-Menten rate and an external solution whose concentrations are constant. It does not demonstrate a chaotic regime.

Use the following unambiguous species names rather than interpreting the diagram's product indices:

| Symbol | Concentration |
| --- | --- |
| S | Ethyl acetate |
| A, B | Acetic acid and acetate |
| E, F | Ethanol and ethoxide |
| h, o | Hydrogen and hydroxide ions |

The reaction rates are

```math
R=\frac{V_{\max}S}{(K_M+S)(1+K_{es2}/h+h/K_{es1})},
\qquad r_A=k_2A-k_{2r}hB,
```

```math
r_E=k_3E-k_{3r}hF,\qquad r_W=k_w-k_{wr}ho.
```

Here is our stoichiometric expansion, with unresolved transport explicitly left as J. It is **not** a claim that the complete published simulation has been recovered:

```math
\begin{aligned}
\dot S&=-R+J_S, & \dot A&=R-r_A+J_A, & \dot B&=r_A+J_B,\\
\dot E&=R-r_E+J_E, & \dot F&=r_E+J_F,\\
\dot h&=r_A+r_E+r_W+J_h, & \dot o&=r_W+J_o.
\end{aligned}
```

The article's generic transport form is `J_X = k_X (X_external - X)`. Substrate and hydroxide exchange are explicit; proton transport is discussed, but its baseline is not given in the main text. Product-specific exchange, reservoir values and initial values must be recovered from S1. Assigning them zero merely because a diagram omits arrows would be an unjustified closure.

| Directly specified quantity | Value |
| --- | --- |
| k2 | 7.3e5 s^-1 |
| k2r | 4.5e10 M^-1 s^-1 |
| k3 | 1.25e-6 s^-1 |
| k3r | 1e10 M^-1 s^-1 |
| kw | 1e-4 M s^-1 |
| kwr | 1e10 M^-1 s^-1 |
| Figure 3 example | S_external = 2.56e-3 M, Vmax = 1.02e-3 M s^-1 |
| Figure 3 transport/feed | kOH = 0.045 s^-1, kS = 0.0014 s^-1, OH_external = 1e-4 M |

Missing baseline values include KM, the individual enzyme protonation constants, proton transport, all product transport choices, and initial concentrations. A pH optimum of 7 constrains a combination of the protonation constants, not their separate values or the curve width. The paper's enzyme-concentration annotations also contain arithmetic inconsistencies; reproduce the explicit Vmax values rather than infer a calibrated enzyme concentration from those annotations.

## What the reaction network already tells us

The following deductions are our algebraic checks of the recovered stoichiometry. They do not fill in missing physical parameters.

### Open reservoir and finite fuel are different models

The published reservoir is maintained at constant composition. Introducing a finite external volume and depleting its substrate would add states and change the model. A finite-fuel batch simulation cannot establish a sustained attractor of the maintained-reservoir system, and the converse does not predict indefinite oscillations in a sealed experiment.

For nonnegative concentrations and rates, the chemical reaction terms preserve nonnegativity. With the stated substrate exchange and positive kS,

```math
0\leq S(t)\leq\max(S(0),S_{\rm external}).
```

At h = 0 the continuous form of R is obtained by multiplying numerator and denominator by h. Water dissociation supplies a positive contribution to the proton derivative. These are useful implementation checks. They do not provide an absorbing bound for every species.

### A product-removal test is essential

Define total acetate `T_A = A+B`, total ethanol `T_E = E+F`, and the charge combination of the modeled ions `q = h-o-B-F`. Direct cancellation gives

```math
\dot T_A=R+J_A+J_B,\qquad
\dot T_E=R+J_E+J_F,
```

```math
\dot q=J_h-J_o-J_B-J_F.
```

If `J_A+J_B` is identically zero, then `T_A` is nondecreasing, and strictly increases wherever R is positive. Any putative periodic orbit with sustained reaction would give a positive integral of R over its period, contradicting return of T_A. The same argument applies to a compact recurrent invariant set in the positive-concentration interior, where R has a positive minimum. It also applies independently to trapped total ethanol.

This is a **conditional exclusion**, not a criticism or refutation of the published oscillations: S1 may contain removal terms omitted from the schematic. A pH trace alone cannot settle boundedness of all seven concentrations. If a product is deliberately omitted as an accumulating spectator in a reduced model, any attractor claim must be restricted to that stated reduced system and its justified feedback assumptions.

### Seven concentrations do not imply seven slow degrees of freedom

The reaction stoichiometric matrix has rank four. With all exchange turned off, three independent linear quantities are conserved: `S+A+B`, `S+E+F`, and `h-o-B-F`. The three acid-base reaction columns alone have rank three. Both rank calculations and the cancellation of the conserved quantities were checked directly; the conservation identities use integer arithmetic.

In the instantaneous acid-base-equilibrium limit, define

```math
K_A=k_2/k_{2r},\quad K_E=k_3/k_{3r},\quad K_W=k_w/k_{wr}.
```

Their values are `1.622222e-5 M`, `1.25e-16 M`, and `1e-14 M^2`. The fast equilibrium relations give

```math
A=\frac{T_Ah}{h+K_A},\quad B=\frac{T_AK_A}{h+K_A},\quad
E=\frac{T_Eh}{h+K_E},\quad F=\frac{T_EK_E}{h+K_E},\quad o=K_W/h.
```

Consequently h is determined by

```math
q=h-\frac{K_W}{h}-\frac{T_AK_A}{h+K_A}-\frac{T_EK_E}{h+K_E}.
```

The right-hand side is strictly increasing for h positive and nonnegative totals:

```math
1+\frac{K_W}{h^2}+\frac{T_AK_A}{(h+K_A)^2}+\frac{T_EK_E}{(h+K_E)^2}>0.
```

It ranges from negative to positive infinity. Thus there is a unique positive h for each finite q and nonnegative pair of totals. This produces a candidate four-variable slow description `(S,T_A,T_E,q)` once the actual transport law is known. It is a singular-limit derivation, not an exact reduction of finite rate constants. Do not impose `q=0`: unspecified counterions and exchange cannot be replaced by an invented electroneutrality constraint on this subset of ions.

If both products are trapped, `T_A-T_E` is additionally constant, but the monotonicity obstruction still applies. If removal is present, it can change both the effective dimension and recurrence. The present evidence supports neither an unconditional no-chaos theorem nor a positive candidate.

### Numerical consequences

The acetic-acid dissociation coefficient alone is 7.3e5 s^-1, while substrate exchange in the illustrated case is 0.0014 s^-1. A naive explicit RK4 GPU scan of the complete chemistry would have severe stability restrictions. Reconstruct the exact model and reproduce its baseline with an implicit CPU solver first. A GPU implementation would need a suitable stiff method or a separately validated slow reduction. Artificially slowing acid-base reactions could create misleading dynamics and would be a new model experiment.

## Priority and later-source checks

The 2023 paper's classification concerns steady and oscillating pH, not Lyapunov spectra. Searching the journal name *Chaos* is therefore not evidence that this particular system was shown to be chaotic.

OpenAlex work [W4390277561](https://openalex.org/W4390277561) returned four citing records on this audit date, representing three distinct works:

| Citing work | Checked scope and relevance |
| --- | --- |
| Le Cacheux et al., *Experimentally Guided Iterative Parameter Estimation for Predictive Chemical Oscillator Models* (2025), [doi:10.1002/anie.202511413](https://doi.org/10.1002/anie.202511413) | Primary article text describes an Fmoc/amine chemical network, not the esterase compartment. Its German-edition DOI is the fourth duplicate record. |
| Shumba et al., *Tunable photochemical properties of silica gel Belousov-Zhabotinsky micro-reactors* (2025), [doi:10.1063/5.0266767](https://doi.org/10.1063/5.0266767) | Primary abstract concerns BZ microparticles. Full equivalence/chaos analysis was not performed. |
| *Chemical Feedback and Collective pH-Dynamics in Enzyme-Powered Artificial Cells* (2026), [doi:10.1021/acs.jpcb.6c03926](https://doi.org/10.1021/acs.jpcb.6c03926) | Retrieved article passages concern collective pH feedback and bistability in other enzyme systems. Not an exhaustive full-model comparison. |

Older papain models are relevant mechanism-level prior art. Ipsen and Schreiber's [primary paper](https://doi.org/10.1063/1.1311980), Eq. 21, describes a two-variable substrate/proton reduction. The accessible primary excerpt was inspected; the full PDF download was blocked. An autonomous smooth planar reduction cannot supply an ordinary strange attractor, but this does not exclude the incompletely recovered seven-species esterase model. Matching the full transport and elimination assumptions matters more than matching the word “enzyme.”

Focused title, esterase/papain, Lyapunov and chaotic-oscillation searches did not identify a primary demonstration of chaos in the exact 2023 compartment equations. This is a retrieval result with a serious source gap, not evidence of novelty.

## Source gap and next concrete action

Both White Rose deposits and the Southampton deposit expose the main paper only. The publisher's accessible minimal article body also refers to S1 without exposing its file. A normal browser visit reached a security-verification page; no challenge was solved. The author's [2025 thesis record](https://etheses.whiterose.ac.uk/id/eprint/38169/) states an embargo until **9 February 2027**, so its full text was not accessed or requested.

The next useful step is recovery of publicly accessible S1, with the exact product fluxes and baseline parameters checked before implementation. If both product totals lack removal, prioritize the monotonicity exclusion over a chaos search. If removal permits bounded recurrence, begin from the published oscillatory parameter set and validate the slow dimension before allocating GPU time.

Machine-readable extraction, unknowns, conditional deductions and citation coverage: [esterase-compartment.json](../results/candidate-audits/esterase-compartment.json). The cached source PDF and page images remain outside tracked results; no third-party source code was copied into the project.
