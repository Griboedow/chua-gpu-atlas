# Searching for an attractor with a different geometry

The added research objective is to find a robust chaotic-attractor candidate that looks substantially different from familiar examples. This complements the search for overlooked regimes in concrete physical models. The objective remains open: the figures below describe our existing candidates, not a newly established attractor type.

![Existing-candidate geometry baseline](assets/attractor_geometry_baseline.png)

## Three separate questions

1. **Does the geometry look different?** Inspect several views, original state variables and physically meaningful observables. Record the scaling and projection. Do not generate apparent novelty by changing the camera, wrapping angles, adding redundant copies, or applying a nonlinear coordinate distortion.
2. **Is the dynamical organization different?** Examine Poincare return maps, branches, symmetries, bifurcations and selected periodic orbits. Full Lyapunov spectra help characterize a candidate but do not classify its topology. In appropriate low-dimensional cases, periodic-orbit organization and linking can be more discriminating than a picture.
3. **Is that precise result new in the checked literature?** Compare the actual equations and their reductions with primary papers, not just names or thumbnails. Report the limits of the search; an absent result in a query is not a proof of novelty.

Letellier and Aguirre [1] explicitly distinguish a new system from a new attractor and analyze the unusual-looking cord attractor. Its existence is a useful reminder that the comparison set extends beyond Lorenz, Rossler and Chua. Their publisher abstract was inspected in this block; we have not reproduced the cord model or completed its full topology analysis. Gilmore [2] describes the role and dimensional scope of topological methods. These sources motivate the investigation; they are not an automated certificate of novelty.

## What the existing trajectories already tell us

The [four-unit synchronematic candidate](synchronematic-research.md) has seven effective angular states after the common phase is removed. Its saved long CPU spectrum gives a Kaplan–Yorke estimate of **3.11054**, with estimates **3.0970–3.1173** across the four observation windows. This is an estimate from finite-time exponents, not a rigorous Hausdorff dimension or a measured embedding dimension. It cautions against interpreting a three-observable picture as a faithful three-dimensional reconstruction.

The top-left panel uses collective observables:

```text
Q = mean(exp(2*i*chi))
U = mean(exp(i*(theta+chi)))
V = mean(exp(i*(theta-chi)))
C = sqrt((abs(U)**2 + abs(V)**2)/2)
plotted coordinates = (Re Q, Im Q, C).
```

These observables are invariant under the physical equivalence `(theta_i,chi_i) -> (theta_i+pi,chi_i+pi)` and under a common phase shift. An earlier version used `R=abs(mean(exp(i*theta)))`, which changes under equivalent local phase/orientation choices. The figures and GIF now use `C`; no trajectory was changed. The [priority and symmetry audit](synchronematic-priority-audit.md) explains this correction and verifies all 16 local choices along 64,000 saved states.

The bottom-left panel shows the same orbit through different invariant observables. Neither is a Poincare section. No angular discontinuities are joined by artificial straight lines, and no additional orbit is generated. Four-oscillator phase chaos is itself established prior art: Bick et al. [3] reported it in a symmetric phase network in 2011. The later audit compares its full primary equations and symmetries with our phase-plus-orientation system. We make no first-four-oscillator claim or proof of inequivalence to all known systems.

The [slider candidate](lce-slider-research.md) has three states. Its full-state view here only rescales the existing physical coordinates by `q=I/(1+I)` and the mechanical frequency scale. The visually funnel-like structure and curved two-dimensional event section do not establish a new topological class. Inserting the autonomous flow's zero exponent between the two measured transverse exponents gives a Kaplan–Yorke estimate of **2.11292**. The existing relay-equivalence audit already provides a concrete reason to compare it closely with older relay attractors.

The two dimension estimates are descriptive evidence, not a ranking of novelty. Higher dimension alone does not meet the objective either.

## Rotation preview

![Camera rotation of a three-observable projection](assets/synchronematic_projection_rotation.gif)

This preview rotates the camera around the saved synchronematic data, changing azimuth and elevation. It shows a projection of a seven-dimensional flow. It is not a new simulation, a time evolution animation, or evidence that its topology differs from known examples.

## Next comparisons

Prioritize candidates with reproducible bounded dynamics and a recognizable change in their return structure. The newly recovered [six-state coupled optical operators](light-operators-audit.md) offer a concrete next system. For a positive candidate, first validate the dynamics and switching treatment, then produce a full-state or clearly labeled observable comparison. In three dimensions, investigate return branches and periodic-orbit organization before assigning a new name. In higher dimensions, avoid inferring topology from a three-dimensional crossing pattern.

The repository retains separate labels for **visually distinctive**, **dynamically characterized**, and **novelty unconfirmed**. A successful bonus result would need an explicit account of what is different, what comparisons were made, and what remains unproved.

## Reproduction and references

```powershell
.venv/Scripts/python.exe research/attractor_geometry.py
```

The script reuses saved trajectories; no mass calculation is repeated. Numerical estimates and exact observable definitions are saved in `results/attractor-geometry/baseline.json`.

1. Letellier, C., and Aguirre, L. A. (2012). Required criteria for recognizing new types of chaos: Application to the “cord” attractor. *Physical Review E* **85**, 036204. [doi:10.1103/PhysRevE.85.036204](https://doi.org/10.1103/PhysRevE.85.036204). Publisher abstract inspected; full comparison pending.
2. Gilmore, R. (1998). Topological analysis of chaotic dynamical systems. *Reviews of Modern Physics* **70**, 1455. [doi:10.1103/RevModPhys.70.1455](https://doi.org/10.1103/RevModPhys.70.1455). Publisher abstract inspected for the method's scope.
3. Bick, C., Timme, M., Paulikat, D., Rathlev, D., and Ashwin, P. (2011). Chaos in Symmetric Phase Oscillator Networks. *Physical Review Letters* **107**, 244101. [doi:10.1103/PhysRevLett.107.244101](https://doi.org/10.1103/PhysRevLett.107.244101); [publisher issue cover](https://journals.aps.org/prl/issues/107/24). Full primary equations inspected in the later priority audit; global equivalence questions remain open.
