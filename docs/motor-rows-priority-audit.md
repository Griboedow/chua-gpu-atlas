# Priority audit: chaos in the exact molecular-motor row model

Audit date: 30 October 2026, with a bounded search completed before 07:15 UTC. This is a literature assessment, not another numerical experiment.

**The numerical result remains a candidate in an existing model, with scientific priority unresolved.** The inspected full source and closest retrieved follow-ups did not supply an explicit demonstration of chaos in the exact three-row equations. However, coupled molecular-motor chaos, sarcomere chaos and chaotic cluster dynamics already have published precedents. The full 2025 thesis of a source author remains an important unread document. A distinctive projection does not settle either model equivalence or attractor topology.

The [numerical report](motor-rows-gpu-survey.md) contains the independent calculations. The [source audit](motor-rows-source-audit.md) establishes the equations and their physical interpretation. This note concerns what can responsibly be claimed about prior work.

## Exact target and equivalence standard

The target is the constant-total-rate, cosine-potential N-row system of Alouges et al. (2025), with the fixed-endpoint mechanical constraint. The independent first-harmonic subsystem has dimension $3N-1$, hence eight states for $N=3$. Its normalized parameters are $(\rho,r,b)$; the reported three-row candidate has $b=1$.

For comparison, write $z_j=c_j+i s_j$. Our algebraic rewriting of the already verified field is

```math
z_j'=(-1+iV_j)z_j+1+ib,\qquad
V_j=-\rho X_j+r\left(s_j-\frac1N\sum_k s_k\right),\qquad
X_j'=V_j,\qquad \sum_jX_j=0.
```

An exact prior result would need this field, an established invertible coordinate/time transformation, or a justified restriction containing the same candidate. A common Hopf normal form, similar biological motivation, permutation symmetry, or a similar-looking projection is insufficient. Conversely, different printed equations alone do not prove inequivalence under every possible transformation.

The two-row case already has an exact closed core equivalent to Oriola et al. (2014); that result is preserved in the [pair audit](motor-pair-source-audit.md). Equal-cluster embeddings and other invariant cluster subspaces must be distinguished from the full N-row dynamics. The present three-row trajectory does not approach a fixed two-plus-one cluster in the saved observations. This is a numerical distinction, not a proof of a new mechanism.

## What the exact source and its follow-ups establish

| Work and access level | Relevant evidence | Implication for priority |
| --- | --- | --- |
| **Alouges et al. (2025)**, complete arXiv v2, Eq. (15), Fig. 9 and conclusion [1] | The eight-row examples settle into two synchronized groups. General N-row bifurcations and initial-condition selection are left open. The current arXiv record still lists v2, dated 16 April 2025. | No explicit exact N=3 chaos demonstration was identified in the inspected source. This is a statement about that text. |
| **Anello, Alouges and DeSimone (2025)**, complete preprint [2] | Spatial elastic and hydrodynamic flagellar dynamics are added to motor activation. Its comparison of motor descriptions explicitly distinguishes attachment fractions from the Fourier coefficients of the microscopic closure; displayed long-time controls are periodic. | Related physical mechanism, but its spatial model is not automatically the isolated eight-state N=3 target. |
| **Agostinelli et al. (2026)**, primary publisher abstract/preview [3] | Cilia and photo-responsive rods are compared through reduced oscillatory models and Hopf behavior. A repository PDF request returned HTTP 403. | Full coefficient and result audit remains incomplete. No exact-model exclusion or absence claim is made from the abstract. |
| **Ishimoto, Herault and Moreau (2025)**, complete primary text checked in the earlier source audit [4] | Deformation is prescribed for a kinematic swimming calculation. | Its citation of Alouges does not itself report chaos of the autonomous N-row motor system. |

A fresh OpenAlex cited-by request for the 2025 DOI returned **four records representing three distinct works**: the three later journal papers above and a duplicate SSRN version of the flagellar-beating paper. This is discovery metadata, not a complete census. It may miss recently posted, differently indexed, unpublished or non-citing equivalent results. The exact request and records are in the [machine-readable audit](../results/motor-rows/priority-audit.json).

## Existing chaos results that limit broad claims

**Günther and Kruse (2010) already report chaotic sarcomere dynamics.** Their complete primary paper, Section III.C and Fig. 7, describes apparently chaotic solutions after secondary oscillatory bifurcations [5]. The model has four variables: two motor-binding probabilities and two linker extensions. It uses load-dependent detachment, changing filament overlap and mechanical coupling. Those are different closures from the constant-total-rate Fourier transport above; no exact transformation to the N=3 model was established here. This precedent rules out wording such as “first chaos of coupled molecular motors” or “first chaotic few-unit sarcomere model.” The paper's numerical route-to-chaos interpretation is prior published evidence, not a rigorously established mechanism for our candidate.

**A July 2026 preprint makes the broader sarcomere claim still less tenable.** Shintani studies five complex Stuart–Landau units with a phase-lagged nearest-neighbor coupling, heterogeneous detunings and mean damping [6]. Its autonomous internal system has ten real states; a common phase is applied in the observation map, rather than forcing that internal field. The full paper reports positive Lyapunov exponents under repeated integrations. It is a phenomenological amplitude model, not the exact motor probability closure. Its preprint status and the difference in equations both matter; neither makes it irrelevant as prior art on collective sarcomere chaos.

**Permutation-symmetric cluster chaos is also established.** Haugland, Tosolini and Krischer (2021) study globally coupled Stuart–Landau oscillators, cluster splitting and attractor collisions [7]. Their Eq. (1) has cubic local and mean terms, while Eq. (2) constrains the complex mean to harmonic motion. The target motor system instead has mechanical displacement states and no such fixed harmonic-mean identity. No exact equivalence was derived. The shared permutation symmetry makes cluster bifurcations and symmetry-related attractor collisions relevant hypotheses, not new mechanisms by themselves. A triangular three-coordinate cloud cannot demonstrate a topology absent from this literature.

The source-cited motor chain of Costantini and Puglisi (2024) is another useful comparison [8]. Its nearest-neighbor phase equations contain thermal noise, and its main objective is precision and thermodynamic uncertainty. It does not supply an exact deterministic N-row Fourier-field chaos result. This comparison also prevents treating every use of “motor synchronization” as the same model.

## Unresolved gaps and defensible wording

The closest unread item is **Irene Anello's 2025 doctoral thesis, _Mathematical modeling of spontaneous oscillations in flagellar and bio-inspired systems_** [9]. The official institutional listing and indexed abstract were accessible. Normal requests for the official PDF returned HTTP 403. Its full chapters, appendices and any unpublished N-row examples were not inspected. No access restriction was bypassed. A complete thesis check is required before making a strong priority claim.

Other remaining gaps are the complete 2026 cilia/rod follow-up, unindexed later material, and a systematic equivalence analysis beyond the explicit pair reduction. The bounded exact-title, author, motor-chaos and cluster-chaos queries are not an exhaustive mathematical classification. Search-result absence supplies no positive evidence that a result is new.

The supportable description is:

> We obtained sustained numerical chaos candidates in the published N-row molecular-motor model. The inspected source examples were periodic, and the bounded audit did not verify an earlier positive result for these exact equations. Related motor, sarcomere and cluster models already exhibit chaos. Scientific priority and the global bifurcation mechanism remain unresolved.

The same distinction applies to later calculations at the source's eight-row setting: retaining the published equations or activity path strengthens physical relevance, but it does not establish priority. A new model, a previously unreported parameter regime, a distinct dynamical mechanism and an unusual projection are four separate claims.

## References

1. F. Alouges et al., “Some mathematical models for flagellar activation mechanisms,” *Mathematical Models and Methods in Applied Sciences* **35**, 2395–2424 (2025). [DOI](https://doi.org/10.1142/S0218202525500423); [complete v2](https://arxiv.org/html/2409.03506v2); [version record](https://arxiv.org/abs/2409.03506).
2. I. Anello, F. Alouges, A. DeSimone, “Beating of eukaryotic flagella via Hopf bifurcation of a system of stalled molecular motors,” *European Journal of Mechanics A/Solids* **114**, 105729 (2025). [DOI](https://doi.org/10.1016/j.euromechsol.2025.105729); [complete preprint](https://arxiv.org/html/2412.06067v1).
3. D. Agostinelli et al., “Spontaneous oscillations in eukaryotic cilia and photo-responsive rods,” *Journal of the Mechanics and Physics of Solids* **210**, 106529 (2026). [Primary publisher preview](https://www.sciencedirect.com/science/article/pii/S0022509626000293).
4. K. Ishimoto, J. Herault, C. Moreau, “Bending–compression coupling in extensible slender microswimmers,” *Journal of Fluid Mechanics* **1020**, A1 (2025). [Complete primary article](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/bendingcompression-coupling-in-extensible-slender-microswimmers/EF56B665D16B5603D3712959E0FA4911).
5. S. Günther, K. Kruse, “Spontaneous sarcomere dynamics,” *Chaos* **20**, 045122 (2010). [DOI](https://doi.org/10.1063/1.3523283); [complete primary manuscript](https://arxiv.org/pdf/1202.0755).
6. S. A. Shintani, “Chaotic internal dynamics coexist with a stable temporal scaffold in a mesoscale sarcomere model informed by high-resolution recordings,” bioRxiv preprint (20 July 2026). [DOI](https://doi.org/10.64898/2026.07.20.739477); [full primary-paper copy used](https://www.researchgate.net/publication/410649922_Chaotic_internal_dynamics_coexist_with_a_stable_temporal_scaffold_in_a_mesoscale_sarcomere_model_informed_by_high-resolution_recordings).
7. S. W. Haugland, A. Tosolini, K. Krischer, “Between synchrony and turbulence: intricate hierarchies of coexistence patterns,” *Nature Communications* **12**, 5634 (2021). [Complete primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC8463560/); [DOI](https://doi.org/10.1038/s41467-021-25907-7).
8. G. Costantini, A. Puglisi, “Thermodynamic precision of a chain of motors: the difference between phase and noise correlation,” *Journal of Statistical Mechanics* **2024**, 024003. [Complete primary preprint](https://arxiv.org/html/2401.09952v1); [DOI](https://doi.org/10.1088/1742-5468/ad244b).
9. I. Anello, *Mathematical modeling of spontaneous oscillations in flagellar and bio-inspired systems*, doctoral thesis, SISSA (2025). [Official thesis listing](https://www.math.sissa.it/user/72/former_phd_students/); [official PDF, inaccessible during this audit](https://iris.sissa.it/retrieve/b97cdabf-fec2-4383-bfd1-fc36bbdd4e2f/PhD_Thesis_IreneAnello_2025%20%284%29.pdf).

Original papers and theses retain their own copyrights and licenses. Only our audit and independently derived comparisons are part of this repository; source caches are ignored.
