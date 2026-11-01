# Partial formal verification: tangency and storage algebra

This Lean project checks a finite-dimensional step in the [Chua proof audit](../docs/chua-coexistence-audit.md). It does not formalize the Chua flow or prove infinitely many coexisting attractors.

The separate [Storage.lean](ChuaProof/Storage.lean) file adds three real-algebra implications used by the [global-convergence argument](../docs/chua-hidden-boundary.md): `storage_balance`, `storage_dissipation` and `storage_nonincrease`. These assume the output relation and the coefficient/sector signs. Coordinate construction, properness, the Chua flow and LaSalle's theorem are not formalized. Run `python research/check_chua_storage_formal.py` from the repository root to check these declarations, audit their dependencies and reject an altered cross term. Its [record](../results/chua-hidden-boundary/storage-formal-check.json) is separate from the five-declaration tangency audit below; `Audit.lean` still reports only that earlier scope.

## Statements

- `ChuaProof.tangency_box`: continuity on a rectangle, actual state derivatives, and four uniform boundary signs imply an interior point with g = 0 and gₐ = 0. The point minimizes g along its state interval.
- `ChuaProof.quadratic_contact`: positive second state derivative and parameter derivative throughout the box give the same properties at that witness.
- `ChuaProof.graph_guards`: the ten published graph-transform guard inequalities hold as exact rational inequalities over the reals.
- Two endpoint lemmas exclude boundary minimizers using one-sided limits of difference quotients.

The derivative hypotheses use `HasDerivAt`, including at the closed state boundaries. Joint continuity is required only on the closed rectangle. Existence of a contact is a conclusion, not an input. Neither strict convexity nor uniqueness is asserted by the basic theorem.

## Reproduce

Install [Lean through elan](https://lean-lang.org/install/), then run:

```bash
cd formal
lake exe cache get Mathlib.Analysis.Calculus.LocalExtr.Rolle Mathlib.Analysis.Calculus.Deriv.Slope Mathlib.Topology.Order.Compact Mathlib.Tactic.Linarith Mathlib.Tactic.NormNum Mathlib.Data.Real.Basic
lake build
lake env lean Audit.lean
cd ..
python research/check_chua_formal.py
```

The pinned toolchain is **Lean 4.24.0**. Mathlib is pinned to **f897ebcf72cd16f89ab4577d0c826cd14afaafc7** (v4.24.0); [lake-manifest.json](lake-manifest.json) fixes its transitive dependencies. Retain this manifest when reproducing the check. The cache supplies compiled Mathlib dependencies; it does not supply proofs of the new declarations. Lean elaborates and kernel-checks those declarations during the build.

The Python checker requires only the standard library and `lake` on PATH. It checks source identity if an optional `--project` build copy is used, builds the project, and inspects the transitive dependencies printed by Lean. Only `propext`, `Classical.choice` and `Quot.sound` are accepted. It also verifies that reversing the upper-face sign breaks the proof and that a deliberately introduced extra axiom is detected. Generated control files live under the ignored `.lake/build/` directory.

Checks run locally on the research workstation. With the existing Windows installation on D:, run this command from the repository root:

```powershell
.\research\check_chua_formal_local.ps1
```

The [local launcher](../research/check_chua_formal_local.ps1) uses the repository's Python environment and the Lean installation under `D:\_TMP\chua-proof-tools`. Its temporary files and caches stay on D:. It restores the calling shell's environment afterwards.

The [recorded local result](../results/chua-proof-audit/formal-check.json) includes source hashes and axiom lists. Hashes normalize CRLF to LF.

## Trust boundary

This project depends on Lean's kernel and the imported Mathlib environment. The declarations use proof-producing tactics and depend only on the standard axioms listed in the check record. Source hashes identify the files checked; hashes are not a proof of their mathematical interpretation.

The numerical bounds on the actual map are still supplied by the separate Python/FLINT computation and the analytical manifold argument. No code here proves that those outputs satisfy the hypotheses of `tangency_box`. The Newhouse–Robinson theorem and its localization to the Chua return domain are outside this formalization. See the audit's [scope table](../docs/chua-coexistence-audit.md#4-scope-of-formal-verification).
