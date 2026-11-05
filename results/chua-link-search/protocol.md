# Searching for linked coexisting attracting cycles

This search follows the negative four-cycle experiment in `../chua-linked-attractors/`. Computation and caches remain on D:. Earlier roots, classifications, forecasts and certificates are inputs, not overwritten outputs.

## Selection before joint classification

1. Test all nine asymmetric references together with their inversion partners at the identical reference parameters. A04 and A05 have individually knotted components. Inversion guarantees coexistence but forces their pairwise linking number to zero, so this test must examine splitness rather than rely on Gauss linking alone.
2. Intersect the recorded stability intervals at shared slope/alpha nodes and test other representatives falling inside them. Prioritize configurations containing a knotted family. Select well-interior witnesses for each distinct family combination, rather than count repeated samples of the same continuation as independent discoveries.
3. Where the catalogue lacks bounds, continue additional symmetric representatives to both stability boundaries. Use the resulting overlaps to propose new common-parameter configurations. Retain failures and root identity checks.

The joint topology may be split even when its individual components are knots. A nontrivial polynomial relative to the unlink is therefore not sufficient. A diagrammatic decomposition into one-component pieces establishes a split union conditional on extraction. A nonzero Alexander minor for a multicomponent link, or a Jones polynomial different from the split-union product, gives sufficient evidence of non-splitness of that diagram. Zero/equal invariants without a separating diagram remain unresolved.

## Verification and reporting

Check controls containing both split knotted components and non-split links. Reconstruct candidate roots at fixed common parameters, reject invalid itineraries, nonpositive flights, unstable cycles, duplicate components and multiple traversals. Repeat any non-split candidate at higher precision, higher polygon density and different projections; independently integrate its full vector field. Preserve screening failures and changes under refinement. A Brunnian conclusion additionally requires all proper sublinks to be unlinks; a linked pair alone is not such a conclusion.

Publish the search ledger, final numerical witnesses or negative results, exact diagram calculations, visual comparisons and scope limitations. Polygon extraction and a certified isotopy of exact ODE curves remain separate obligations. No claim of historical novelty follows from the search alone.
