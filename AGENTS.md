# Repository conventions

- Write repository content in English; use the user's preferred language in conversation.
- Prefer established terminology from dynamical systems and numerical analysis in prose, figure labels and new identifiers. Use the definitions in `docs/terminology.md` consistently.
- Introduce a project-specific term only when an established term would be inadequate. Define it at first use and distinguish it from the underlying mathematical object.
- In particular, a plotted "layer" represents sampled stability regions of a continued periodic-orbit family. A finite catalogue is not a classification of global connected components.
- Distinguish numerical approximations, interval enclosures, analytic deductions and formal verification. State the map, parameter section and scope of every stability or counting claim.
- Preserve numerical forecast values and certificate bounds. When renaming identifiers, update producers, consumers and links together, and verify the renamed records against their original Git publication. Keep historical source-file digests tied to their producing revision.
- Keep work, environments and caches on drive D:. Run computations locally; do not introduce GitHub Actions.
- Preserve the author's coursework motivation in the main README.
