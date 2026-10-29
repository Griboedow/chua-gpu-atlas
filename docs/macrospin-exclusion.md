# Why the single-macrospin candidate does not need a chaos scan

Arun et al. [1] report self-oscillations driven by self-induced torque (SIT) in a magnetic double tunnel junction. We inspected their explicit macrospin equations. Their noiseless constant-current system has no independent circuit-memory state. The three Cartesian magnetization components obey `|m|=1`; the effective field and pinned magnetization are axial. A large parameter scan would therefore be misplaced.

The following reduction is our own algebra applied to their Eq. 1. Write `z=m_z`, `p=e_z`, `w=m cross p`, `q=m cross (m cross p)`, and `kappa=gamma G`. Every scalar coefficient, including `alpha`, `G`, and `Hs`, depends at most on `z` under these assumptions. With `H_eff=Hz(z) p`, the implicit equation is

$$\dot m=\gamma(\beta H_s-H_z)w+\alpha m\times\dot m+\gamma H_s q+\kappa[p\cdot(m\times\dot m)]q.$$

For `dot m=a w+b q`, use `m cross w=q`, `m cross q=-w` and `p dot (m cross dot m)=-a(1-z^2)`. Define

$$c(z)=\alpha-\kappa(1-z^2),\qquad d(z)=1+\alpha c(z).$$

On a regular branch with `d(z)!=0`, solving the two linear equations gives

$$a=\frac{\gamma[(\beta-\alpha)H_s-H_z]}{d},\qquad
b=\frac{\gamma[(1+\beta c)H_s-cH_z]}{d}.$$

Consequently,

$$\dot z=\frac{\gamma(1-z^2)[cH_z-(1+\beta c)H_s]}{d},\qquad
\dot\varphi=\frac{\gamma[H_z-(\beta-\alpha)H_s]}{d}.$$

The tilt is governed by a **scalar autonomous ODE**, with azimuth slaved to it. A scalar trajectory cannot repeatedly cross its equilibria on a smooth unique-solution branch. Constant tilt permits steady precession, which explains periodic Cartesian components without providing a chaotic mechanism. A singular denominator is a failure of the regular ODE description, not evidence for chaos.

`python research/audit_macrospin.py` compares this reduction against an independent Cartesian 3-by-3 linear solve at 128 random algebra samples. It also checks tangency to the unit sphere. These samples are identity checks, not a fit to the device. Results are in `results/macrospin-exclusion/algebra.json`.

The exclusion applies to the deterministic axial single-spin model. Thermal noise, another dynamical magnetization, circuit memory, or external forcing changes the question. It does not exclude chaos in spintronics generally.

## Reference

1. **Arun, R., Gopal, R., Chandrasekar, V. K., and Lakshmanan, M.** (2026). Self-oscillations induced by self-induced torque in magnetic double tunnel junction. [arXiv:2601.07300v1](https://arxiv.org/html/2601.07300v1), Eqs. 1–3. This audit specifies the inspected preprint version; bibliographic aggregators return inconsistent journal labels.
