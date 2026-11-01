# Coexistence audit and partial formal proof

The [report](../../docs/chua-coexistence-audit.md) separates reproduced numerical hypotheses, Lean-checked implications, and the outstanding localization step for infinite simultaneous coexistence.

| File | Contents |
|---|---|
| [numerical-audit.json](numerical-audit.json) | Fresh 80/100-digit certificates and two expected shifted-box failures; source hashes and software versions |
| [formal-check.json](formal-check.json) | Lean build, five transitive axiom reports, source hashes and two rejection controls |
| [validation-summary.json](validation-summary.json) | Local verification record, preserved forecast hash, browser review and historical Linux check |
| [browser-check.json](browser-check.json) | Live Chrome review of the report, formulas, tables, figure, formal guide and main page |
| [chua-proof-audit.pdf](chua-proof-audit.pdf) | Vector version of the argument map |

From the repository root:

```bash
python -m pip install -r requirements-certificate.txt
python research/audit_chua_coexistence.py
python research/check_chua_formal.py
python research/plot_chua_proof_audit.py
python research/manifest_refined_chua_atlas.py --check
```

Install the pinned Lean environment and Mathlib cache as described in [formal/README.md](../../formal/README.md) before running the formal checker. The numerical audit uses only CPU arithmetic. Existing forecast files and earlier certificates are not overwritten by the audit.

The numerical record labels coexistence `conditional_on_localization_argument`. A successful formal check proves the statements in the Lean source under their stated hypotheses; it does not upgrade that coexistence status.

Proof checks run locally. On the research workstation, `research/check_chua_formal_local.ps1` selects the installed tools and temporary directories on D:. The earlier completed Linux run remains in `validation-summary.json` as historical evidence; its hosted workflow has been removed.
