# Exact N-row molecular-motor model: source and control audit

Status: 30 October 2026. The published N-row equations have been independently implemented and their eight-row resting and oscillatory controls reproduced. This is **not a new model or a new chaotic-attractor claim**.

The source is Alouges et al. (2025), *Some mathematical models for flagellar activation mechanisms*, specifically arXiv v2 Eq. (15), the constant-rate assumption, the cosine interaction potential, Table 1 and Fig. 9. Its motor ensemble is deterministic at the probability-density level, with fixed chemical activity. A depleting ATP reservoir is not an additional state of this model.

## Exact closure and state count

Let the absolute filament shifts obey $x_0=x_N=0$, and define the row displacements $\Delta_j=x_j-x_{j-1}$. Then $\sum_j\Delta_j=0$. Each row has cosine and sine coefficients $p_j^c,p_j^s$ in its motor probability density. For a cosine potential, only these coefficients enter the mechanical force. The other Fourier modes cannot feed back into this first-harmonic subsystem under the source assumptions.

The force-difference equations and $\sum_jv_j=0$ give, exactly,

```math
v_j=-\frac{k}{\eta}\Delta_j-
\frac{\pi U}{\eta}\left(p_j^s-\frac1N\sum_{i=1}^{N}p_i^s\right).
```

This mean coupling is a consequence of the published mechanical constraint, not an added network interaction. With $a=a_0>0$, introduce

```math
\tau=at,\quad X_j=\frac{2\pi\Delta_j}{\ell},\quad
C_j=\ell p_j^c,\quad S_j=\ell p_j^s,
\quad A=\frac{a_1}{a},\quad B=\frac{b_1}{a}.
```

For the source active sign $A<0$, put $c_j=C_j/A$, $s_j=S_j/A$ and

```math
\rho=\frac{k}{\eta a},\qquad
\gamma=\frac{2\pi^2U}{\eta\ell^2a},\qquad
r=-\gamma A,\qquad b=B/A.
```

The autonomous normalized equations are

```math
\begin{aligned}
V_j&=-\rho X_j+r(s_j-\bar s),\qquad
\bar s=\frac1N\sum_{i=1}^{N}s_i,\\
c_j'&=-c_j-V_js_j+1,\\
s_j'&=-s_j+V_jc_j+b,\\
X_j'&=V_j,\qquad \sum_{j=1}^{N}X_j=0.
\end{aligned}
```

[motor_rows.py](../research/motor_rows.py) stores $(c_1,\ldots,c_N,s_1,\ldots,s_N,X_1,\ldots,X_{N-1})$, reconstructing $X_N=-\sum_{j<N}X_j$. The dimension is **$3N-1$**: 8 for three rows, 11 for four rows, and 23 for the source's eight rows. A three-dimensional figure of the full eight-row motion would be a projection. No translation mode, clipping or numerical constraint projection is added.

## Jacobian and structural checks

Let $P=I-\mathbf1\mathbf1^T/N$, and let $D$ stack $I_{N-1}$ above a row of minus ones. In the state ordering above the Jacobian is

```math
J=\begin{pmatrix}
-I&-\mathrm{diag}(V)-r\mathrm{diag}(s)P&\rho\mathrm{diag}(s)D\\
\mathrm{diag}(V)&-I+r\mathrm{diag}(c)P&-\rho\mathrm{diag}(c)D\\
0&rP_{1:N-1,:}&-\rho I_{N-1}
\end{pmatrix}.
```

The code checks this Jacobian by complex differentiation, reconstructs velocities independently by solving the original force-difference linear system, and checks the exact two-row correspondence. For $N=2$,

```math
Y=-\frac{s_1-s_2}{2},\quad
Z=1-\frac{c_1+c_2}{2},\quad
h=\rho X_1+rY
```

obey the [three-state Oriola core](motor-pair-source-audit.md) exactly, for any fixed $b$. This equivalence check is especially useful for detecting factors of two or a reversed transport sign.

For each row,

```math
\frac12(c_j^2+s_j^2)'=-c_j^2-s_j^2+c_j+bs_j.
```

Hence the disk $c_j^2+s_j^2\leq1+b^2$ is invariant from the source's equilibrium-density starts. Since $X_j'=-\rho X_j+r(s_j-\bar s)$, these disks also bound every relative displacement; for instance $|X_j|\leq\max(|X_j(0)|,2r\sqrt{1+b^2}/\rho)$. This proves boundedness on the stated invariant region, not chaos.

The first-harmonic source rates require $|A|\sqrt{1+b^2}\leq1/2$ to stay nonnegative. With equilibrium means, $\ell Q_j$ is bounded by $1/2\pm |A|\sqrt{c_j^2+s_j^2}$, so physical density bounds can be checked directly. Varying $b$ changes the prescribed transition-rate asymmetry. It is not part of the fixed-$b=1$ numerical controls below.

The unique equilibrium has $c_j=1,s_j=b,X_j=0$. Its eigenvalues consist of $N+1$ copies of $-1$ and $N-1$ pairs satisfying

```math
\lambda^2+(\rho+1-r)\lambda+\rho=0.
```

At $r=\rho+1$, the full system has $N-1$ simultaneous critical pairs. For $N>2$, it cannot be identified near this point with one three-state core plus strictly stable passive variables. This does not establish a new chaotic mechanism.

The identical-row equations have simultaneous permutation symmetry. Subspaces where groups of rows share the same $(c,s,X)$ are invariant. In particular, two equal-sized clusters reduce to the two-row problem. The published eight-row examples and our reproduced oscillation have two groups of four. The single-core global Lyapunov argument must not simply be asserted for the entire $3N-1$ system.

## Source control tuples and initial conditions

The source uses $a_0=c\sqrt{\Omega_{\rm src}}$ and $a_1=b_1=d\Omega_{\rm src}$, with Table 1's mechanical constants. Evaluating chemical activity numerically in nN nm, so $\Omega_0=15k_BT=0.064002$, reproduces its printed $d=-56.4588$. The [two-row audit](motor-pair-source-audit.md) records this unit convention and the original constants; no coefficient is fitted here.

| Fig. 9 control | $\delta$ | $a_0$ (s⁻¹) | $\rho$ | $r$ | $b$ | $A=B$ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Below onset | −0.1 | 240.00374997 | 3.95827149 | 4.75514605 | 1 | −0.01355033 |
| Above onset | +0.1 | 265.33412898 | 3.58039127 | 4.75514605 | 1 | −0.01498045 |

Here $\Omega_{\rm src}=(1+\delta)\Omega_0$. Along this particular source activity path, $r$ stays fixed while $\rho$ changes. An independent $(\rho,r)$ survey would therefore explore a broader mathematical family than that one source path, and would need a separate physical interpretation.

The source sets all motor densities initially to equilibrium and draws the **absolute** shifts $x_1,\ldots,x_7$ uniformly from $[-0.01,0.01]$ nm, with endpoints fixed. Our seed is 20251030; it is our reproducibility choice, not a recovered source seed. Differences of these absolute positions are converted to $X_j$. Independently drawing relative displacements would change the initial prescription and could violate the sum constraint.

## Independent numerical controls

Reproduction commands:

```bash
python research/motor_rows.py --audit-output outputs/motor-rows-algebra.json
python research/motor_rows.py --control-output outputs/motor-rows-controls.json
```

The complete parameters, initial state and numerical diagnostics are saved in [source-audit.json](../results/motor-rows/source-audit.json). No author code is executed.

The algebra audit uses 30 random samples at each of $N=2,3,8$. Maximum scaled errors are $6.1\times10^{-16}$ for the analytic Jacobian, $1.7\times10^{-15}$ for the independent force solve, and $6.5\times10^{-16}$ for the two-row reduction.

The independent CPU runs used DOP853 and Radau, requested relative tolerance $2\times10^{-10}$ and absolute tolerance $2\times10^{-12}$. Both were compared over $\tau\leq20$ and over the shared longer horizon $\tau\leq300$; DOP853 continued to $\tau=1000$.

| Control | Maximum shared-horizon state difference | Observed outcome |
| --- | ---: | --- |
| Below onset | $1.17\times10^{-11}$ | Equilibrium reached; final DOP853 derivative below $7\times10^{-16}$ |
| Above onset | $6.19\times10^{-11}$ | Sustained two-cluster oscillation, with peak-to-peak $X_j\simeq0.688242$ |

For the positive control, sampled return intervals are approximately 3.285893 in $\tau$. They were obtained by interpolating sampled crossings, not by periodic-orbit shooting. Our seed selects clusters $\{1,2,5,7\}$ and $\{3,4,6,8\}$; a different partition from the figure is expected because the source does not specify its random draws. The constraint residual stays below $1.2\times10^{-16}$, and physical density bounds remain within $0.4788\leq\ell Q_j\leq0.5212$. The source regimes are qualitatively reproduced. No Lyapunov spectrum, chaotic regime or transverse cluster-stability claim is inferred from these controls.

All six integrations completed; total control wall time was about 53 seconds. There were no hidden solver failures or filtered trajectories. The implicit and explicit methods agreeing at the reported resolution does not prove all higher-dimensional regimes are correctly resolved.

## Priority and interpretation

The two-row core is already present in Oriola et al. (2014), which explicitly searched for chaos and reported periodic behavior. The published N-row extension supplies a genuine higher-dimensional target, but not an unstudied family by definition. Related flagellar, coupled-motor and active-filament models already have extensive nonlinear-dynamics literature.

The bounded primary-source audit and forward-discovery record are inherited from the [two-row audit](motor-pair-source-audit.md). Supplementary exact-title, author and N-row queries found no additionally verified exact-model positive-chaos result. This is not proof of absence or priority. The citing JFM paper on bending–compression explicitly prescribes deformation and studies its kinematic response; it is not a chaos result for these N-row equations. Its full primary text was checked. No new broad literature search or GPU survey forms part of this source-control report.

## References

1. F. Alouges, I. Anello, A. DeSimone, A. Lefebvre-Lepot, J. Levillain, “Some mathematical models for flagellar activation mechanisms,” *Mathematical Models and Methods in Applied Sciences* **35**, 2395–2424 (2025). [DOI](https://doi.org/10.1142/S0218202525500423); [complete arXiv v2](https://arxiv.org/html/2409.03506v2); [SISSA record](https://iris.sissa.it/handle/20.500.11767/147990).
2. D. Oriola, H. Gadêlha, C. Blanch-Mercader, J. Casademunt, “Subharmonic oscillations of collective molecular motors,” *EPL* **107**, 18002 (2014). [DOI](https://doi.org/10.1209/0295-5075/107/18002).
3. K. Ishimoto, J. Herault, C. Moreau, “Bending–compression coupling in extensible slender microswimmers,” *Journal of Fluid Mechanics* **1020**, A1 (2025). [Complete primary article](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/bendingcompression-coupling-in-extensible-slender-microswimmers/EF56B665D16B5603D3712959E0FA4911).

The implementation is independent. Source papers retain their original copyrights and licenses; cached source files are excluded from the repository.
