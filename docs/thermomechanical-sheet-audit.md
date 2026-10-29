# Thermomechanical sheet: model audit before searching

**Status: 29 October 2026. No chaos survey or positive result yet.** Shen's 2025 article [1] studies a heated elastic sheet using displacement, velocity and temperature. It reports regular oscillations in its examples. The full publisher PDF was retrieved and its equation page visually checked. A cubic velocity damping term is introduced in the numerical section, although absent from initial Eq. 1; that term must be included when investigating the simulated model.

Write `theta=T-T0` and `Delta=TH-T0`. The equations to investigate are

$$m\ddot A+c\dot A+(k_0+\alpha\theta)A+\beta\dot A^3=\gamma\theta,
\qquad \tau\dot\theta=\Delta-\theta-\mu A.$$

## An equilibrium inconsistency

Our substitution of `theta*=Delta-mu A*` gives

$$-\alpha\mu A_*^2+(k_0+\alpha\Delta+\gamma\mu)A_*-\gamma\Delta=0.$$

The paper's Eq. 3 omits the quadratic term. Its displayed expression is therefore not an exact equilibrium of the printed equations in general. It can be a small-deflection approximation only when that term is negligible; we must not use it uncritically as the center of a numerical scan.

For an illustrative mathematical example `m=c=k0=gamma=mu=tau=Delta=1`, `alpha=-0.1`, the exact positive displacement is `0.5124921973`. The displayed expression gives `0.5263157895`, with nonzero force-balance residual `0.0277008310`. These are deliberately simple audit values, **not** an experimental fit or a reproduction of a source figure.

## Exact dimensionless model

For positive `Delta,mu,tau`, put

$$s=t/\tau,\quad X=\mu A/\Delta,\quad V=\tau\mu\dot A/\Delta,\quad Y=\theta/\Delta.$$

Our transformation yields

$$X'=V,\qquad V'=GY-K(1-SY)X-BV-DV^3,\qquad Y'=1-Y-X,$$

where

$$K=\frac{k_0\tau^2}{m},\quad G=\frac{\gamma\mu\tau^2}{m},\quad
B=\frac{c\tau}{m},\quad S=-\frac{\alpha\Delta}{k_0},\quad
D=\frac{\beta\Delta^2}{m\tau\mu^2}.$$

The physical equilibrium satisfies `X*=1-Y*` and `K(1-SY*)(1-Y*)=G Y*`. For positive `K,G` and nonnegative `S`, there is a unique root with `0<Y*<min(1,1/S)` (interpret `1/S` as infinity for `S=0`). The other algebraic branch is outside this interval.

Set `Ke=K(1-SY*)` and `F=G+KSX*`. The characteristic polynomial is

$$p(\lambda)=(\lambda+1)(\lambda^2+B\lambda+K_e)+F.$$

At this physical equilibrium all coefficients are positive. Local stability is equivalent to `B(Ke+B+1)>F`; equality is a linear Hopf candidate with squared frequency `Ke+B`. The cubic damping has zero derivative at equilibrium, but affects nonlinear saturation. This local criterion follows from the literal equations and retains the displacement-dependent feedback omitted by substituting the plate temperature directly for the sheet temperature.

The divergence is `-1-B-3DV^2`. Negative divergence does **not** prove boundedness. The constraints `X>=0`, `0<=Y<=1`, and positive effective stiffness `1-SY>0` must be monitored separately: the smooth equations do not automatically enforce contact with the plate or all thermal constraints. Excursions outside those conditions cannot be advertised as a physically valid sheet attractor.

## Reproducibility and next checks

`python research/thermomechanical_sheet.py` independently checks the coordinate transformation on 32 samples, the analytic Jacobian by finite differences, and the equilibrium characteristic polynomial. Recorded errors are `1.8e-15` for coordinate conversion and `8.9e-10` for the Jacobian. Results are in `results/thermomechanical-sheet/algebra-audit.json`.

Before a broad search: establish the source's cubic damping scale and units, audit forward citations and equivalent thermal/jerk models, then use GPU ensembles with domain-exit flags and independent CPU controls. The paper labels `mu` dimensionless while using displacement in the dimensional thermal balance; a dimensional calibration needs clarification. Exact-title/DOI chaos queries found no matching result in the checked material, which is not proof of novelty.

An OpenAlex snapshot returned five forward records, saved in `results/thermomechanical-sheet/literature-forward.json`. These include two LCE actuator papers, two apparently less direct material/thermal papers and a thermomechanical metamaterial review. Their exact citation context has not yet been checked; index metadata alone does not settle coverage or relevance.

## Reference

1. **Shen, X.** (2025). Heat Conduction and Its Related Interdisciplinary Areas: Self-Excited Oscillations in a Thermomechanical Elastic Sheet. *Chinese Physics Letters*, **42**, 090801. [doi:10.1088/0256-307X/42/9/090801](https://doi.org/10.1088/0256-307X/42/9/090801). [Publisher article](https://cpl.iphy.ac.cn/article/doi/10.1088/0256-307X/42/9/090801). [Publisher PDF](https://cpl.iphy.ac.cn/en/article/pdf/preview/10.1088/0256-307X/42/9/090801.pdf).
