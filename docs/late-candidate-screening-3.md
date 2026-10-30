# Late screening III: mechanochemical models and a consequential version check

Audit: 2026-10-30 UTC. No trajectories were simulated. This block recovers one exact five-state follow-up problem and excludes a second model because its chaos is already explicitly published. It does **not** establish a new chaos discovery or a clean first-chaos opportunity.

## 1. A single mechanochemical Brusselator: ready to implement, with prior-art qualifications

The early preprint by Dewan, Mondal, and Sarkar [1] is titled *Oscillation death by mechanochemical feedback*. Its latest version [2] adds an author, changes the title to *Symmetry-protected phases in a 1D active solid with mechanochemical feedback*, and corresponds to PRL **137**, 078401 (2026), DOI [10.1103/p1bw-v87g](https://doi.org/10.1103/p1bw-v87g). Checking only the initial version would miss material prior art: the latest main text explicitly reports **spatiotemporal chaos**, and supplementary section X discusses Lyapunov spectra.

The precise follow-up retained here is the **single-cell**, five-state equation S2 of [2], with the source's explicit choice $f(L,t)=0$. It is distinct from the spatial model and retains no externally varying drive:

```math
\dot M=a-(b+1)M+cM^2E,
\qquad \dot E=bM-cM^2E-DE,
```

```math
\tau_D\dot D=D_0-D-\beta D(L-1),
\qquad
\tau_L\dot L_0=1-L_0-\alpha(E-E_r),
\qquad
\tau_R\dot L=L_0-L.
```

We use $E_r$ to distinguish the **fixed reference concentration** called $E_0$ in the source from an initial condition. The source parameters are $a=c=1$, $b=3$, $D_0=0.1$, $\tau_R=1$, with $\tau_D=1.1$, $\tau_L=1.2$ as defaults. The source initializes chemical variables on the uncoupled limit cycle, with $D=D_0$ and $L=L_0=1$.

### Exact normalization and implementation checks: our calculations

The uncoupled reference is determined without integrating:

```math
M_r^3-M_r^2+0.4M_r-0.1=0,
\qquad
M_r=0.6137372763943942,
\qquad
E_r=\frac{3M_r}{M_r^2+0.1}=3.8626272360560576.
```

The cubic is strictly increasing, so this real root is unique. The two chemical residuals are below $3\times10^{-16}$ in floating-point arithmetic. Initializing every state exactly at this equilibrium is not a useful oscillator test.

For $\alpha>0$, put $x=(1-L)/\alpha$, $y=(1-L_0)/\alpha$, and $\mu=\alpha\beta$. This is an invertible affine transformation, not a new constitutive law. The last three equations become

```math
\tau_D\dot D=D_0-D+\mu Dx,
\qquad
\tau_L\dot y=E-E_r-y,
\qquad
\tau_R\dot x=y-x.
```

Thus the single-cell dynamical search can use $(\mu,\tau_D,\tau_L)$, fixing $\tau_R=1$. The separate value of $\alpha$ is still needed to reconstruct physical lengths; $\beta=\mu/\alpha$ must be reported. This reduction **does not** survive inserting the spatial model's cubic regularizers without further changes.

For state order $(M,E,D,y,x)$, the full trajectory-dependent Jacobian is

```math
J=\begin{pmatrix}
-(b+1)+2cME & cM^2 & 0 & 0 & 0\\
b-2cME & -cM^2-D & -E & 0 & 0\\
0 & 0 & (-1+\mu x)/\tau_D & 0 & \mu D/\tau_D\\
0 & 1/\tau_L & 0 & -1/\tau_L & 0\\
0 & 0 & 0 & 1/\tau_R & -1/\tau_R
\end{pmatrix}.
```

In particular, an equilibrium Jacobian evaluated at $D_0,L=1$ must not be used throughout a nonlinear orbit. For a future spectrum calculation, the sum should match the time average of

```math
\mathrm{div}F=-(b+1)+2cME-cM^2-D
+\frac{-1+\mu x}{\tau_D}-\frac1{\tau_L}-\frac1{\tau_R}.
```

An exact pipeline control is available at $\mu=0$. On the stable planar chemical limit cycle with $D=D_0$, the complete five-state spectrum consists of the chemical pair $\{0,\lambda_{\mathrm{rad}}\}$, with $\lambda_{\mathrm{rad}}\lt0$, and the three exponents $-1/\tau_D,-1/\tau_L,-1/\tau_R$. The mechanical variables are then driven linear filters. This follows from the block-triangular variational equations; it is not a spectrum measured here.

### Physical interpretation and guards

The nonnegative chemical octant is forward invariant: at its respective boundaries, $\dot M=a$, $\dot E=bM$, and $\tau_D\dot D=D_0$. There is no finite fuel inventory in this maintained-feed phenomenological model. However, positivity of chemical variables does not establish boundedness or physical length positivity.

The $L_0=0$ boundary can point outward if $E>E_r+1/\alpha$. Therefore, $L,L_0>0$ must be checked throughout any proposed physical orbit; do not clip negative lengths. In normalized coordinates these conditions are $\alpha x\lt1$ and $\alpha y\lt1$. Linear-elastic validity additionally limits deformation. A bounded trace may allow a small enough $\alpha$ for positive lengths at the same $\mu$, but that is a parameter choice to disclose, not a proof of long-time boundedness or an experimental calibration.

The continuum model in [2] contains extra cubic saturation terms, and its mechanical equation conserves total length. A spatially uniform limit of that PDE is not obtained by simply replacing it with the single-cell relaxation equation. We retain the source's S2 literally and do not transfer the PDE's boundedness or chaos claims to it.

### What the priority check actually establishes

The broad family is already reported to have spatiotemporal chaos. The exact unsaturated single-cell S2 is a narrower possible verification project; this audit does not establish that its chaos is unpublished.

The latest SI section X describes QR factorization, associates five reported exponents with component variables, and reports four positive values throughout its surveyed parameter space. Its written description omits the usual logged, accumulated, time-normalized QR factors. The text alone does not specify sufficiently clearly whether these are a full spatial discretization spectrum or a local approximation. These are reasons to request an independent diagnostic calculation, **not** proof that the authors' unseen implementation is wrong. For a bounded, nonstationary recurrent autonomous orbit the flow-direction exponent is neutral, and a full five-state spectrum cannot consist of five strictly positive asymptotic exponents. The control above makes this distinction testable.

The OpenAlex record returned zero forward citations for the merged preprint/journal record at audit time, while a June 2026 primary epithelial-tissue article [3] explicitly cites the earlier preprint. This mismatch demonstrates incomplete indexing; it cannot support an absence claim. Targeted queries covered both titles, identifier, authors, mechanochemical Brusselators, and chaos. Known chemical/spatial chaos and the latest paper's own claims substantially narrow the possible contribution.

**Follow-up:** the [literal five-state study](hhs-research.md) now implements this target, verifies the planar/full-spectrum control, and completes 3,267 scouts. No retained chaos candidate emerges from 186 completed longer checks; one stiff GPU interval is incomplete. The [source audit](hhs-source-audit.md) preserves a figure-reproduction mismatch, and the [structural note](hhs-structure.md) establishes an initial-condition-dependent global bound. None of these results establishes first chaos in this already studied family.

## 2. Inertial mechanochemical solid: an already published chaotic three-state model

Sarkar et al. [4] derive a single-mode Galerkin model of a chemically driven elastic solid. In units where $\tau_{\mathrm{el}}=1$, their equation (5) is

```math
\dot U=V,
\qquad
\dot V=-U-\mathsf{Vi}V-X,
\qquad
\mathsf{Ch}\dot X=X(1-X^2)-A_uU-A_vV.
```

Figure 4e explicitly reports a chaotic strange attractor for $\mathsf{Vi}=0.16$, $\mathsf{Ch}=100$, $A_u=50$, $A_v=84$. The neighboring illustrated values $A_v=68$, $78.5$, and $106$ give a cycle, period doubling, and a larger cycle, respectively. This is decisive exact-model prior art, not an inference from a visually irregular figure. Both the initial and revised primary manuscripts were checked; [4] is PRL **135**, 258301 (2025), DOI [10.1103/19rh-3whq](https://doi.org/10.1103/19rh-3whq).

The variables are modal amplitudes, and the ODE is a published spatial approximation. Its physical derivation assumes small strain and small reaction chemical potentials. A mathematical ODE solution alone does not validate those approximations. Our direct checks show inversion symmetry $(U,V,X)\mapsto(-U,-V,-X)$ and divergence $-\mathsf{Vi}+(1-3X^2)/\mathsf{Ch}$, which is strictly negative for the stated example. Neither property is a proof of chaos.

The source explicitly cites related control-system and hearing models; no equivalence to Lorenz or Rössler is established here. The indexed forward-citation query returned three records, saved in the companion JSON. **Decision:** exclude from first-chaos searching; optionally useful as a known physical-model benchmark, without a novelty claim.

## Sources and provenance

1. P. Dewan, S. Mondal, and S. Sarkar, *Oscillation death by mechanochemical feedback*, [arXiv:2504.19655v1](https://arxiv.org/html/2504.19655v1) (28 April 2025).
2. S. Mondal, P. Dewan, L. S. Kumar, and S. Sarkar, *Symmetry-protected phases in a 1D active solid with mechanochemical feedback*, [full v3 with SI](https://arxiv.org/html/2504.19655v3) (8 July 2026); [published record](https://journals.aps.org/prl/abstract/10.1103/p1bw-v87g) (11 August 2026). The journal full text requires a subscription; the public author manuscript was used, without bypassing that restriction.
3. S. Muthukrishnan et al., *Glassy dynamics in active epithelia emerge from an interplay of mechanochemical feedback and crowding*, [Nature Communications primary article](https://www.nature.com/articles/s41467-026-74163-0) (2026). This is a related tissue model and a forward citation, not a demonstrated equivalence to S2.
4. S. Sarkar et al., *Mechanochemical feedback drives complex inertial dynamics in active solids*, [revised primary manuscript](https://arxiv.org/html/2505.18272v2) (25 December 2025); [DOI](https://doi.org/10.1103/19rh-3whq).

Companion: [machine-readable audit](../results/candidate-audits/late-candidate-screening-3.json). Public HTML and metadata snapshots are cached under ignored `outputs/candidate_research/`. No author code was copied or executed, and no new simulation, GPU job, commit, or external message was made in this screening block.
