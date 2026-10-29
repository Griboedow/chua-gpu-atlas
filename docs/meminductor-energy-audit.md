# Moving-coil meminductor: reciprocity before a chaos search

Research notebook, 29 October 2026.

**Finite electrical response alone does not produce recurrent chaos in a reversible moving-coil model driven by a fixed DC voltage through a fixed positive resistance.** A shifted energy decreases exactly, without eliminating the electrical state. However, the force and flux equations printed in the Dinavahi–Harris article do not together satisfy the reciprocity condition of that model. The conditional exclusion below therefore does not settle the experimental apparatus, especially its optical feedback circuit.

This is an application of standard electromechanical energy methods [3], not a claim of a new general theorem or a refutation of the reported oscillations. It extends the narrower [constant-current mechanical argument](candidate-research.md#physical-meminductor-a-useful-exclusion-before-simulation).

## What the source actually specifies

The 2025 article [1], Eqs. (2), (3), (5), (9) and (13), gives the following pair, with the mechanical coordinate increasing upward:

$$
\Phi=L(x)i,\qquad L(x)=c\mu(x),\quad c>0,
$$
$$
\mu(x)=8-\frac7{1+\exp[-250(x-0.03)]},\qquad
F(x,i)=a(x)i,\quad a(x)=K_1\mu(x)/x^4>0,
$$
$$m\ddot x=F-k_a\dot x-mg.$$

The force is **linear in current** because the moving coil interacts with a permanent magnet. Replacing it by a pure reluctance force proportional to current squared would change the source model. Likewise, the experimental optical circuit has a position-dependent photoresistor branch and is not a fixed-voltage, fixed-resistance drive.

For a reversible magnetic transducer with a differentiable coenergy function $C(x,i)$, the same sign convention requires

$$\Phi=\partial_i C,\qquad F=\partial_x C,
\qquad \boxed{\partial_i F=\partial_x\Phi}.$$

Taking the printed pair literally gives

$$\partial_iF=a(x)>0,\qquad \partial_x\Phi=L'(x)i<0\quad(i>0).$$

The cross derivatives have opposite signs. This conclusion does not depend on the positive scale chosen for inductance. Thus these two expressions cannot, together, be the complete constitutive law of a reversible magnetic transducer with these two ports and coordinates. They can still be separate empirical approximations; omitted magnetic bias, additional internal states, or a different meaning of the measured inductance must be resolved before constructing a closed energy-consistent circuit model.

The earlier experimental paper [2] distinguishes motion-dependent flux from the quasi-static estimate and constructs an effective inductance from flux divided by current. Such an effective ratio is not automatically the differential self-inductance needed in a terminal state equation. We have not reconstructed a calibrated magnetic energy from those measurements.

**Supplement recovered, 29 October 2026:** the official Europe PMC supplementary archive contains `41598_2025_12529_MOESM4_ESM.docx`. Its complete text and equations, together with the embedded contact diagram, were inspected. It gives illustrative sinusoidal current/inductance signals, a measured winding resistance of **628 ohms**, and limitations of subtracting the resistive voltage before integrating flux. It describes collisions during bipolar-current fingerprint measurements, including the vertical device hitting its lower magnet. This does not imply that every DC optical-feedback trajectory has such a collision. The supplement supplies neither a calibrated dynamical photoresponse nor a reciprocal magnetic energy or impact law. The missing closure is therefore still unresolved, now after inspecting the supplement rather than leaving its contents unknown. Retrieval metadata and source hashes are in [the supplement audit](../results/meminductor-energy/supplement-audit.json); the source document is not redistributed under our MIT license.

## A reciprocal model including a permanent magnet

One possible reversible constitutive form is

$$\Phi=L(x)i+\psi(x),\qquad L(x)>0,$$
$$H(x,\Phi)=\frac{[\Phi-\psi(x)]^2}{2L(x)}.$$

Here $\psi$ is the permanent-magnet flux linkage, and $H$ is stored magnetic energy. An additional zero-current magnetic potential can be included in the mechanical potential $U(x)$. Differentiation at fixed flux gives

$$i=\partial_\Phi H,\qquad
F=-\partial_xH=\tfrac12L'(x)i^2+\psi'(x)i.$$

Consequently the voltage equation contains both motional terms:

$$L(x)\dot i=V_s-Ri-[L'(x)i+\psi'(x)]\dot x.$$

Choosing $\psi'=a$ reproduces the source's linear force term, but also requires the reluctance contribution $L'i^2/2$ and the associated back EMF. This is a **changed, uncalibrated completion**, not an assertion that these are the fitted equations of the published device.

## Exact DC exclusion with the electrical state retained

More generally, let $H(x,\Phi)$ be a smooth reversible magnetic energy. Assume a regular autonomous model

$$
\dot x=v,\qquad m\dot v=-U'(x)-b v-\partial_xH,
\qquad \dot\Phi=V_s-Ri,\quad i=\partial_\Phi H,
$$

with $m,b,R>0$ and constant $V_s$. Set $i_0=V_s/R$ and define

$$W=\tfrac12mv^2+U(x)+H(x,\Phi)-i_0\Phi.$$

Direct differentiation cancels mechanical and magnetic power transfer:

$$
\dot W=-b v^2+(i-i_0)(V_s-Ri)
=\boxed{-b v^2-R(i-i_0)^2}.
$$

There is no assumption that $L/R$ is small. On any compact forward-invariant region where the equations remain regular, the invariant part of $\dot W=0$ consists of equilibria: $v=0$, $i=i_0$, $\dot\Phi=0$, and mechanical force balance. LaSalle's invariance argument excludes a non-equilibrium recurrent attractor there, including periodic, quasiperiodic and chaotic recurrent motion.

This identity alone does **not** establish global boundedness, uniqueness of equilibrium, or convergence to a single point when the equilibrium set is non-isolated. An invariant set can contain connecting transients; the exclusion concerns recurrent non-equilibrium dynamics. Singular constitutive laws and additional hysteretic states require a separate analysis.

The same argument works for a static strictly increasing series voltage law $\rho(i)$, with $V_s=\rho(i_0)$:

$$\dot W=-b v^2-(i-i_0)[\rho(i)-\rho(i_0)]\le0.$$

A constant current source shunted by a **constant** resistor also reduces to a fixed Thevenin drive. A photoresistor that changes with position or its own state does not: its effective source and resistance vary together, and the constant-$i_0$ cancellation need not hold. This is why the optical experiment remains a separate candidate.

### Passive impacts do not automatically remove the exclusion

The recovered supplement makes contact worth considering explicitly. At a fixed wall, suppose position and flux are continuous and the velocity resets as $v^+=-e v^-$, with restitution $0\le e\le1$. Under these additional assumptions the same shifted energy jumps by

$$W^+-W^-=-\tfrac12m(1-e^2)(v^-)^2\le0.$$

Thus a passive mechanical bounce cannot replenish energy lost between impacts in the reciprocal fixed-DC model. A regular recurrent orbit with nonzero integrated dissipation cannot close its energy balance. This observation is not a global hybrid convergence theorem: sticking, accumulating impacts, impulsive electrical changes and additional internal states need their own laws. In particular, the measured voltage spikes in the source are not a calibrated proof of flux continuity. The identity is a conditional audit, not a fitted simulation of those collisions.

## Why a missing back EMF can mislead a numerical search

For the literal pair $\Phi=Li$, $F=ai$, the natural trial energy $E=mv^2/2+U+Li^2/2$ obeys

$$\dot E=i\dot\Phi-bv^2+[ai-\tfrac12L'i^2]v.$$

The last term is unaccounted power with no fixed sign. A numerically accurate integration of this closure would not establish energy consistency.

Locally, a quasi-static expansion of its electrical equation gives

$$\delta i\simeq-\frac{L'i_0}{R}v.$$

With $L'<0$ and $F_i=a>0$, the force correction points along velocity: apparent negative damping. In a reciprocal model, both the force derivative and motional-voltage coefficient are instead the same quantity $h=L'i_0+\psi'$, giving

$$\delta i\simeq-hv/R,\qquad \delta F\simeq-h^2v/R.$$

That correction is dissipative. These are local fast-response expansions, not proofs of oscillation in the inconsistent closure; the exact energy exclusion above does not rely on the expansions.

## Reproducible algebra and numerical control

Run from the repository root:

```bash
python research/audit_meminductor_energy.py
python research/audit_meminductor_energy.py --impacts
```

The [script](../research/audit_meminductor_energy.py) audits the source cross derivatives and independently checks the shifted-energy gradient and its derivative at 128 random states. It then integrates a deliberately illustrative reciprocal model with

$$L=1+0.3\tanh x,\quad \psi=0.2\tanh x,\quad U=x^2/2,$$
$$m=1,\quad b=0.3,\quad i_0=0.8,$$

for $R=0.1,1,10$, starting at $(x,v,i)=(-1.5,1,-0.5)$. This model has positive bounded inductance and a confining shifted energy. These dimensionless choices are **not a fit to the apparatus**.

- Maximum direct energy-derivative identity error: $1.43\times10^{-14}$.
- Maximum gradient finite-difference discrepancy: $7.22\times10^{-10}$.
- Maximum integrated energy-balance residual over 200 time units: $8.12\times10^{-13}$.
- All three illustrative trajectories approach the same equilibrium, $(x,v,i)\approx(0.241618,0,0.8)$.

These CPU controls check the implementation and algebra; the exclusion follows from the analytical identity. A GPU parameter sweep cannot overturn that identity inside its assumptions, so no bulk search was spent on this class. Compact data, solver tolerances and the explicit calibration disclaimer are in [results/meminductor-energy](../results/meminductor-energy/).

## What remains open

For the actual device, the next useful inputs are a consistent magnetic flux/force calibration, the optical transfer function, the photoresistor's response dynamics, and friction/contact laws. The [assumed optical-feedback model](photofeedback-research.md) already demonstrates an instability mechanism, but its negative chaos survey applies only to that assumed model. Neither its extra state nor the permanent-magnet completion above should be mistaken for measured device physics.

A forward-citation snapshot returned two works: a five-MOSFET meminductor emulator explicitly studied for chaos [4], and a hardware-security review [5]. The emulator uses a different active circuit. The review distinguishes emulator demonstrations from physical devices and does not supply the missing moving-coil equations. Neither is evidence that this particular physical apparatus has a newly established chaotic attractor; nor does this small, incomplete citation index establish a publication gap.

## References

1. Dinavahi, A., and Harris, H. R. (2025). Edge-of-chaos operation and persistent dynamics for neuromorphic meminductor computing. *Scientific Reports*, **15**, 40672. [DOI](https://doi.org/10.1038/s41598-025-12529-y). [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12630749/).
2. Dinavahi, A., Yamamoto, A., and Harris, H. R. (2023). Physical evidence of meminductance in a passive, two-terminal circuit element. *Scientific Reports*, **13**, 1817. [DOI](https://doi.org/10.1038/s41598-022-24914-y). [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9892601/).
3. Haus, H. A., and Melcher, J. R. (1989). *Electromagnetic Fields and Energy*. Prentice Hall, Section 11.7, Macroscopic Magnetic Forces. [MIT open text](https://web.mit.edu/6.013_book/www/chapter11/11.7.html).
4. Kalra, A., Pandey, N., and Rai, S. K. (2026). Compact design of MOSFET-only meminductor emulator and its application in chaos and neural spike generator. *Integration* **109**, 102740. [DOI](https://doi.org/10.1016/j.vlsi.2026.102740). Publisher preview inspected; distinct active emulator.
5. Romero, F. J., et al. (2026). Memory circuit elements for hardware security: a review from memristors to memcapacitors and meminductors. *Frontiers in Electronics* **7**, 1919384. [DOI](https://doi.org/10.3389/felec.2026.1919384). Full primary review inspected; it is a literature map, not a calibrated device model or a general proof that memory elements generate chaos.
