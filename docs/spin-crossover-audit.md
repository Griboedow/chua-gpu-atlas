# Spin crossover: a structural screening result

Research note, 30 October 2026. **The two interaction limits illustrated in Orlov et al. are poor targets for a new autonomous strange attractor.** One is explicitly planar; the other has an exact first integral. The model with both interactions remains a separate question, but its population domain and a source parameter discrepancy need resolution before a large numerical search.

The audit concerns the literal equations in [arXiv:2606.25506v1](https://arxiv.org/abs/2606.25506), not an independently verified implementation of the later journal article. The arXiv HTML and PDF both give Eqs. 16 and 26 below. The PDF was visually checked on pp. 7–8. Its linearization is numbered 31–32; the HTML numbers it 35–36. The published [Physical Review B article](https://doi.org/10.1103/n5cq-qxdq) is accessible at abstract level, but its full text returned an authorization error. No correction was inferred from it.

## An exact reduction of the magnetic interaction limit

Write $n=n_{HS}=\tau+1/2$, and let $m$ denote the source's magnetic variable. Define

```math
p=\frac{\theta}{\exp(D_\tau/T)+\theta},\qquad
D_\tau=\Delta_s-T\log g_\tau-J_\tau(n-1/2),\qquad D_s=J_snm.
```

Energies are expressed in kelvin, so $k_B=1$. The source's Eq. 26 becomes

```math
\dot m=\Gamma J_s n\left[pS B_S(SD_s/T)-nm\right],\qquad
\dot n=\Gamma J_\tau(p-n)+\frac{m}{n}\dot m.
```

For $J_s>0$, its thermal equation uses $\kappa=m$:

```math
\dot T=-\alpha(T-T_R)+I[1+(\rho-1)n]-H\dot m,
\qquad H=\Delta H/C_p.
```

Our direct differentiation gives

```math
\frac{d}{dt}(n^2-m^2)=2\Gamma J_\tau n(p-n).
```

Therefore, at **$J_\tau=0$**, the quantity

```math
\boxed{C=n^2-m^2}
```

is exactly conserved. On the regular positive-population branch, $n=\sqrt{C+m^2}$, so substitution leaves an autonomous two-state system in $(m,T)$. A classical autonomous strange attractor is excluded on compact regular regions of this planar branch. This applies regardless of constant illumination or latent-heat coefficient. It does not classify singular limits at $n=0$.

The source's Figure 11 uses precisely $J_\tau=0$. Plotting temperature, magnetization and population does not create a third independent degree of freedom. Conversely, setting $J_s=0$ and $m=0$ gives the source's explicit two-state $(\tau,T)$ photothermal model. Neither restriction excludes chaos in a spatial model or in a different kinetic closure.

The identity is checked at 512 random states to $8.9\times10^{-16}$. Independent DOP853 integration of the three equations and Radau integration on the exact level surface agree within $1.7\times10^{-10}$ over a regular 20-unit control; the invariant drifts by less than $10^{-15}$. This control relaxes and is not offered as a reproduction of Figure 11's periodic orbit. Numerical checks support the algebra, which supplies the exclusion.

## Why the combined-interaction model is still unresolved

When both $J_s$ and $J_\tau$ are positive, the displayed derivative of $C$ generally does not vanish. The first-integral argument must not be extended to this case.

There is a separate physical-domain issue. A population requires $0\leq n\leq1$. At the admissible boundary point $(n,m,T)=(1,1.9,20)$, with the reference magnetic parameters $J_s=112$, $\Delta_s=309.5$, $g_\tau=3$, and $\Gamma=0.05$, the literal vector field gives:

| Elastic interaction | Population derivative |
| --- | ---: |
| $J_\tau=0$ | $\dot n=1.042399$ |
| $J_\tau=152$ | $\dot n=1.063097$ |

Both point outside the population interval. This is a counterexample to **global invariance of the nominal population domain**, not proof that every trajectory is unphysical. A useful search would first identify a justified invariant operating region or establish which revised closure was intended. Clipping $n$ or replacing Eq. 26 with a probability-preserving kinetic model would change the equations and is not done here.

## A parameter convention that must be resolved

At the chosen photothermal equilibrium $(\tau,T)=(0,T_{eq})$, the literal source linearization gives

```math
\mathrm{tr}\,J=\Gamma J_\tau\left(\frac{J_\tau}{4T_{eq}}-1\right)-\alpha.
```

Substituting the listed $J_\tau=152$, $T_{eq}=112.6$, and $\Gamma=1/19$ yields

```math
\mathrm{tr}\,J=-5.300178-\alpha.
```

It cannot vanish for positive $\alpha$, whereas the source's Figure 3 places a center at $\alpha\simeq2.8$. Its printed critical-reservoir formula similarly evaluates to $176.543$ K rather than the stated $78.8$ K. Finite differences of the literal vector field reproduce our Jacobian within $7.5\times10^{-10}$.

Removing the factor $1/4$ from the trace expression gives an intercept of $2.799290$, which suggests a normalization discrepancy. **This is only a diagnostic; it does not identify a justified corrected model.** Replacing one coefficient while retaining the rest would be arbitrary. The separately rounded $\Delta_s=788$ versus $T_{eq}\log g=788.2$ cannot explain a factor of four. The journal version or original simulation conventions must be checked before claiming reproduction or error in the final publication.

## Reproduction and search limits

```bash
python research/audit_spin_crossover.py
```

The [independent implementation](../research/audit_spin_crossover.py) evaluates the $S=2$ Brillouin terms through their equivalent finite spin sum, avoiding a removable singularity at zero magnetization. [Audit results](../results/spin-crossover/audit.json) retain all parameters, residuals and source-PDF SHA-256; [control data](../results/spin-crossover/regular-control.npz) retain both solver trajectories. No GPU survey is warranted for the excluded limits.

Exact-identifier, title, author and first-integral searches did not establish an earlier occurrence of our displayed identity. That limited search does not establish mathematical novelty. Older photothermal oscillations are already published. The 2014 original was inspected at abstract level; the 2024 precursor was likewise verified bibliographically, but its full equations were not retrieved. The [access record](../results/spin-crossover/literature-audit.json) states these limits.

## References

1. Orlov, Yu. S., Paklin, N. N., Nikolaev, S. V., Shneyder, E. I., and Dudnikov, V. A. (2026). *Dynamic dissipative structures in bistable magnetically ordered spin crossover systems: Self-oscillations of magnetization*. Physical Review B **114**, 134311. [DOI](https://doi.org/10.1103/n5cq-qxdq). Audited equations: [arXiv v1](https://arxiv.org/abs/2606.25506v1).
2. Boukheddaden, K., Paez-Espejo, M., Varret, F., and Sy, M. (2014). *Autocatalytic spin-crossover transition: Nonlinear dynamics induced by a photothermal instability*. Physical Review B **89**, 224303. [DOI](https://doi.org/10.1103/PhysRevB.89.224303).
3. Orlov, Yu. S., Nikolaev, S. V., and Paklin, N. N. (2024). *Photoinduced nonlinear dynamics of strongly correlated systems with spin crossover: Autocatalytic spin transition*. JETP Letters **119**, 227–232. [DOI](https://doi.org/10.1134/S0021364023603962).
