# Why the two-oscillator search has an extra constraint

Research note, 29 October 2026. This is a direct algebraic control for the fixed-center Fourier model of Leyva et al. [1], not a claim of a new integrability theorem or of a conserved quantity in the full experimental apparatus.

**Two oscillators have a first integral after common phase is removed.** Their three apparent angular variables are therefore not three freely evolving degrees of freedom on a regular level set. This gives a structural explanation for the negative two-oscillator GPU controls in the [main notebook](synchronematic-research.md).

Let $\delta=\theta_2-\theta_1$, let $a,b$ be the two orientations measured from the same line of centers, and let $w=d^{-3}>0$. Define the odd parts of the phase and orientation interaction functions as

$$S_\theta(\delta)=\frac{F_\theta(\delta)-F_\theta(-\delta)}2,
\qquad S_\chi(\delta)=\frac{F_\chi(\delta)-F_\chi(-\delta)}2.$$

The reduced pair equations are

$$\dot\delta=-2wS_\theta(\delta)\cos a\cos b,$$
$$\dot a=wF_\chi(\delta)\sin a\cos b,\qquad
\dot b=wF_\chi(-\delta)\sin b\cos a.$$

On a chamber where $\sin a\sin b\ne0$ and $S_\theta(\delta)\ne0$,

$$\frac{d}{dt}\log\left|\frac{\sin a}{\sin b}\right|
=2wS_\chi(\delta)\cos a\cos b.$$

Hence

```math
\boxed{I=\log\left|\frac{\sin a}{\sin b}\right|
+\int^\delta\frac{S_\chi(u)}{S_\theta(u)}\,du,\qquad \dot I=0.}
```

The even Fourier coefficients cancel from this identity. No assumption that the entire field is a gradient is required.

For the printed coefficients,

$$S_\theta=\sin\delta(1.758+1.248\cos^2\delta),$$
$$S_\chi=\sin\delta(1.582-5.436\cos^2\delta).$$

Thus the ratio has no singularity inside a phase chamber. An elementary antiderivative is available: if the ratio is $(A+B\cos^2u)/(C+D\cos^2u)$ with $C,C+D>0$ and $D\ne0$, use

```math
\frac BD u+
\frac{A-BC/D}{\sqrt{C(C+D)}}\,
\arctan_{\rm continuous}
\left(\sqrt{\frac C{C+D}}\tan u\right).
```

The continuous branch matters; a raw principal arctangent would create false jumps. The implementation uses `atan2` and a phase-chamber index.

## What this does and does not exclude

The planes $\sin a=0$, $\sin b=0$ and $\sin\delta=0$ are invariant for the printed coefficients. Inside each chamber, a regular level of $I$ is a two-dimensional invariant surface. Ordinary autonomous strange-attractor dynamics cannot occur as a smooth recurrent flow on a compact regular part of such a surface. The boundary planes themselves have at most two independent states. Interior critical points of this integral have $\cos a=\cos b=0$ and are equilibria of the pair equations.

This is a useful analytical obstruction to interpreting a finite positive exponent in this pair model as a three-dimensional chaotic attractor. It is not a complete classification of all boundary connections or singular limits. The integral is singular at some chamber boundaries, so its numerical value should not be monitored there without a different coordinate chart.

The cancellation does **not** carry over directly to three or more oscillators: each orientation and phase difference then receives other pair contributions. Consequently this result does not prove that four oscillators are minimal for collective chaos. It also does not extend automatically to evolving positions, noise, or the untruncated physical interaction law.

## Reproduction

```bash
python research/audit_synchronematic_pair.py
```

The [independent audit](../research/audit_synchronematic_pair.py) checks both the printed coefficient set and the alternative NumPy-source sign. For each, 256 random interior states with varied pair direction and distance satisfy $\nabla I\cdot F=0$ to at most $1.78\times10^{-15}$. Finite-difference gradients agree within $1.20\times10^{-9}$. Two independent DOP853 control trajectories conserve $I$ within $4.89\times10^{-15}$ over their recorded regular intervals. These computations check the algebra and implementation; they do not replace the derivation.

Results and settings: [n2-first-integral/audit.json](../results/synchronematic/n2-first-integral/audit.json). The article supplement and the earlier dissertation [2] already analyze pair synchronization and orientation dynamics; this note does not claim that the identity has never appeared elsewhere.

## References

1. Leyva, S. G., Zhang, Z., Olvera de la Cruz, M., and Bishop, K. J. M. (2026). Self-oscillating synchronematic colloids. *Nature Communications*, **17**, 1841. [DOI](https://doi.org/10.1038/s41467-026-68552-8). The stated identity applies to the specified separable Fourier approximation.
2. Zhang, Z. (2023). *Quincke Oscillators: Dynamics, synchronization, and assembly of self-oscillating colloids*. PhD dissertation, Columbia University. [DOI](https://doi.org/10.7916/0pgk-pw45).
