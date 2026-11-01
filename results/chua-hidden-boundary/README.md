# Global sector and hidden-cycle separator: evidence

The [report](../../docs/chua-hidden-boundary.md) states the analytical argument, published questions, numerical findings and unresolved obligations. Parameters in this study use γ = 0.

## Reproduce

From the repository root, with the CPU environment installed:

```bash
python -m pip install -r requirements-certificate.txt
python research/derive_chua_global_sector.py
python research/study_chua_hidden_boundary.py
python research/certify_chua_basin_saddle.py
python research/check_chua_hidden_boundary.py
python research/check_chua_storage_formal.py
python research/plot_chua_hidden_boundary.py
python research/manifest_refined_chua_atlas.py --check
```

The formal checker additionally requires the pinned [Lean environment](../../formal/README.md), with `lake` on PATH and its Mathlib dependencies available. On the research workstation, keep tools and caches on D: and use:

```powershell
$env:ELAN_HOME = 'D:\_TMP\chua-proof-tools\elan'
$env:PATH = 'D:\_TMP\chua-proof-tools\elan\bin;' + $env:PATH
$env:MATHLIB_CACHE_DIR = 'D:\_TMP\chua-proof-tools\mathlib-cache'
$env:TEMP = 'D:\_TMP\Матеша\.agent-tmp'
$env:TMP = $env:TEMP
.\.venv\Scripts\python.exe -X utf8 research\check_chua_storage_formal.py
```

The numerical study includes sixteen integrations to time 700 and may take several minutes. It uses earlier accessibility data to reconstruct the outer-cycle seed. The prediction file is written before the fine-distance sweep; an existing prediction is checked for identical coefficients and is not overwritten. It was fixed after exploratory continuation, not published externally before the search.

## Records

| Files | Role and scope |
|---|---|
| `global-sector-algebra.json` | Nine exact rational identities, the coordinate matrix, a rejected wrong feedthrough sign and SymPy version |
| `global-sector-80.json`, `global-sector-100.json` | Positive storage minors and coefficient enclosures; 1,000 secondary derivative controls |
| `cycle-seeds.json` | Reconstructed symmetric cycles and complete first-return itineraries |
| `saddle-certificate-80.json`, `saddle-certificate-100.json` | Saddle and sink inclusions, full-flight validity, multiplier classification and a rejected shifted root |
| `cycle-independent-controls.json` | DOP853 and Radau checks of both cycles and their variational equations |
| `branches.json` | 180 numerical roots along each cycle branch |
| `basin-separator.json` | 71 local section points at four/six return depths and sixteen two-sided numerical trajectory outcomes |
| `linear-contact-prediction.json` | Analytically chosen coefficients and exponents, with timing/protocol stated |
| `linear-contact-results.json` | 31 actual piecewise-linear roots per contact and eight independent integration controls |
| `numerical-audit.json` | Repeated-precision overlap and numerical acceptance checks |
| `storage-formal-check.json` | Three Lean declarations, standard axiom dependencies, input hashes and rejection of an altered balance |
| `audit-corrections.json` | Export and plotting errors found and corrected before publication |
| `browser-check.json` | Live GitHub checks: 17 rendered report formulas, three figures/tables and the README; one notation incompatibility corrected |

The preliminary files `exploratory-roots.json`, `inner-lower.json`, `inner-upper.json` and `outer-to-stability.json` preserve the exploratory continuation that preceded the prediction. They do not constitute a second independent experiment or the interval evidence. The main study program reconstructs the roots and repeats the relevant branches.

Arb's joint `ball` strings include decimal rounding; the global-sector files also provide high-precision `enclosure` strings. Separately rounded `midpoint` and internal `radius` fields are descriptive and must not be combined as if they were a lossless interval serialization. The repeated-precision audit uses the joint strings.

The global theorem applies to every initial state under its stated sector assumptions. The separator curve, finite basin outcomes and asymptotic coefficient tests have numerical scope. The formal verification covers three algebraic implications, not the entire global-convergence proof.

## Figures

- [Global sector, PDF](chua-global-sector.pdf)
- [Saddle and local separator, PDF](chua-basin-saddle.pdf)
- [Linear contacts and prediction comparison, PDF](chua-linear-contact.pdf)

The corresponding PNGs are displayed in the report and generated from the same records. Computation and proof checking are local; no GitHub Actions workflow is used.
