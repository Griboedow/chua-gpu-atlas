# Terminology and data conventions

The reports use the following terms to distinguish objects in state space, parameter space and the finite numerical catalogue. Each report specifies its equations, parameter values and return map. The definitions of hidden and self-excited attractors follow [Stankevich et al. (2017)](https://arxiv.org/abs/1710.02677) and [Kuznetsov et al. (2023)](https://doi.org/10.1007/s11071-022-08078-y).

## Periodic motion and stability

| Term | Meaning in this repository |
|---|---|
| Periodic orbit; cycle | A nonconstant trajectory that repeats after a positive time. Its period is the least such time unless stated otherwise. |
| Limit cycle | An isolated periodic orbit. This term is not used for the nonisolated periodic orbits at a linear centre. |
| Attracting periodic orbit | A periodic orbit that attracts nearby initial states. A numerical approximation is labelled as such unless existence and attraction have been certified. |
| Basin of attraction | The set of initial states whose forward trajectories approach the specified attractor. It belongs to state space. |
| Hidden attractor | An attractor whose basin avoids some open neighborhood of every equilibrium. Avoidance by a finite set of initial states does not establish hiddenness. |
| Self-excited attractor | An attractor whose basin intersects every neighborhood of at least one equilibrium. |
| Poincaré return map | The map from a transverse section to its next specified return. The section and return convention must be stated. |
| Floquet multipliers | Eigenvalues of the linearized full-period evolution. For an autonomous flow the tangential multiplier is 1; the nontrivial multipliers govern transverse stability. Eigenvalues of a symmetry-reduced half-return map are called half-return multipliers, not full-period Floquet multipliers. |
| Saddle-node bifurcation of periodic orbits; fold | A local collision of two periodic-orbit branches. "Fold" is used as its short name after the object has been specified. |
| Symmetry-breaking bifurcation | A bifurcation from a symmetric branch to branches lacking that symmetry. The reports identify the relevant map and evidence before calling it a pitchfork. |
| Homoclinic tangency | A nontransverse intersection of the stable and unstable manifolds of the same saddle. A finite-iterate tangent-line estimate is distinguished from a certificate for the exact manifolds. |

## Parameter regions and the catalogue

A **stability region** is a set of parameter values where a specified periodic-orbit branch is attracting. A **periodic window** is an interval of periodic attraction on a specified one-dimensional parameter section. The **stability interval of the symmetric cycle** ends when that cycle loses stability; attracting asymmetric or period-doubled branches can continue beyond it.

**Layer** is the atlas's visual shorthand for the stability region associated with a continued periodic-orbit family. It is not an additional dynamical invariant. A displayed **sheet** is a surface interpolating selected orbit representatives or bifurcation points, as specified in its caption. Its graphical thickness is not the stability width, and an interpolated mesh does not certify the unsampled interior or global topology.

A **catalogue entry** identifies a reference orbit under the stated normalization, symmetry and continuation convention. Numerically traced branch relationships can associate several descendants with one entry. Two inversion-related cycles count as distinct coexisting attractors, even when the catalogue groups them in one family. The 32 original mesh entries and 34 extended reference entries therefore do not constitute a proof of the number of globally distinct families.

An **itinerary** records successive visits to the affine regions. A **flight** in the implementation is a trajectory segment between switching-surface crossings. Equal itineraries or equal peak counts do not imply equal periodic orbits.

The slope prediction study uses a **bias correction**: it adds each target family's fold-location residual at the central parameter node to the uncorrected prediction throughout the small tested slope patch. The correction is held fixed before the off-centre tests; its success there does not justify global extrapolation.

For the distinct map conventions and verification dependencies in the hidden-coexistence proof, see the [review guide](chua-review-guide.md). Numerical tolerances, interval bounds and formal proof obligations remain separate in that guide.
