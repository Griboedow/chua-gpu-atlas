# Linked attracting-cycle witnesses

This directory contains a finite search for non-split unions of coexisting attracting periodic orbits in the piecewise-linear Chua system. Three principal pairs pass exact joint-Jones versus split-product tests at 20 precision/resolution/projection combinations: six for S02/S08 and seven for each daughter pair. Interval calculations certify all six cycles' existence, distinctness, attraction and primitive periods. Link classification concerns numerical polygons; an ambient-isotopy certificate for the exact trajectories is absent.

[Report](../../docs/chua-link-search.md) · [Offline 3D viewer](linked-pairs.html) · [Review record](review.json) · [Frozen protocol](protocol.md)

## Principal witnesses

| Pair | Exact parameters (γ = 0) | Refined root records |
|---|---|---|
| S02/S08 | α = 8.94; β = 14.26071457889439; m₀ = −1.1473571428571427; m₁ = −0.7142857142857143 | [S02](overlaps/S02-S08/S02-d90-n384.json), [S08](overlaps/S02-S08/S08-d90-n384.json) |
| S01 daughter/inverse | α = 9; β = 14.25984090319409; m₀ = −8/7; m₁ = −5/7 | [Daughter](daughters/S01-a1e-03-d90-n384.json) |
| S02 daughter/inverse | α = 9; β = 14.290802548649257; m₀ = −8/7; m₁ = −5/7 | [Daughter](daughters/S02-a1e-03-d90-n384.json) |

Decimal strings specify exact chosen values. The first row's slopes are not replaced by nearby rational classical slopes. The daughter suffix `a1e-03` identifies the initial flight-time asymmetry constraint, not a geometric separation.

Every `*-jones-*.json` or `jones-*.json` record stores the joint polynomial, individual component polynomials, hypothetical split-union product and computation diagnostics. Polynomial dictionaries map integer exponents of A to integer coefficients. The convention is t = A⁻⁴, J(unknot) = 1 and J(split union) = (−A² − A⁻²)J₁J₂. For all three principal pairs, J differs from that split product, while the oriented linking matrix and the sampled Alexander-minor evaluations vanish.

## Contents and interpretation

| File or directory | Content |
|---|---|
| [protocol.md](protocol.md) | Search plan frozen in commit `36ed4c5`; the report discloses the subsequent daughter-search adaptation |
| [catalogue-candidates.json](catalogue-candidates.json) | 137 candidates from recorded common-node stability bounds; five selected family combinations |
| `symmetry/` | Nine original asymmetric references paired with their inverses, all split on the retained screening calculation |
| `overlaps/` | Common-parameter roots, polygons and joint diagrams for five combinations; S02/S08 is the positive witness |
| `daughters/` | S01/S02 branch probes and fixed-parameter refinements; the S02 amplitude-0.003 root is retained as unstable |
| `*.npz` | Time samples `t` and common-coordinate curve samples `xyz`, including a repeated closing endpoint; daughter records also retain `partner_xyz` |
| `joint-*.json`, `*-joint-*.json` | Extracted and reduced planar diagrams, component counts, linking matrices and screening outcomes |
| `jones-*.json`, `*-jones-*.json` | Exact joint and component Jones polynomials and split-product comparison |
| [exhaustive-jones.json](overlaps/S02-S08/exhaustive-jones.json) | Independent enumeration of all 262,144 smoothing states for an 18-crossing S02/S08 diagram |
| [invariant-controls.json](invariant-controls.json) | Eight tests, including split knotted components and known non-split links |
| [bracket-controls.json](bracket-controls.json) | Fourteen agreements between grouped and exhaustive Kauffman state sums |
| [independent-audit.json](independent-audit.json) | Fourteen additional exact polynomial controls and six outward-rounded interval separation checks, with source hashes |
| [certificates-d140.json](certificates-d140.json), [certificates-d160.json](certificates-d160.json) | All six interval-certified cycles, section images, primitive-period checks, pair separation and shifted-root negative control |
| `certificate-*-d*.json` | The same interval witness records separated by pair |
| [certificate-development.json](certificate-development.json) | Failed first-return assumption and wider-box enclosure attempt, with corrections |
| [review.json](review.json) | Replayed polynomials, numerical integration comparisons, polygon gaps and agreement of interval runs |

An early screen marked `unresolved` has not been rewritten after a later Jones calculation. Read the companion polynomial record for the final numerical conclusion. Likewise, initial daughter probes retain their original numerical scope; their later interval certificates are separate evidence.

Text source hashes normalize CRLF to LF before SHA-256, matching the repository manifest. Binary polygon hashes use the original bytes. This convention makes the input chain independent of checkout line endings.

## Reproduction on the local computer

Run from `D:/_TMP/Матеша/chua-attractor`, using the existing environment and the dependencies in [requirements-knots.txt](../../requirements-knots.txt), [requirements-precision.txt](../../requirements-precision.txt) and [requirements-certificate.txt](../../requirements-certificate.txt). Work, temporary files and caches remain on D:. No GitHub Actions are required.

To replay the published evidence without repeating the exploratory continuation:

```powershell
.\.venv\Scripts\python.exe -X utf8 research/review_chua_link_search.py
.\.venv\Scripts\python.exe -X utf8 research/chua_link_invariants.py
.\.venv\Scripts\python.exe -X utf8 research/chua_link_bracket.py
.\.venv\Scripts\python.exe -X utf8 research/certify_chua_link_witnesses.py --digits 140 --radius 1e-60
.\.venv\Scripts\python.exe -X utf8 research/certify_chua_link_witnesses.py --digits 160 --radius 1e-60
.\.venv\Scripts\python.exe -X utf8 research/audit_chua_link_algorithms.py
.\.venv\Scripts\python.exe -X utf8 research/plot_chua_link_search.py
```

The review replays 20 joint polynomials and performs 12 independent ODE integrations. It checks the frozen protocol against Git commit `36ed4c5`, so that commit must be available locally. The certificate checks use exact parameter strings and refine numerical starting guesses anew. The plot program regenerates both PNG/PDF figures and the offline viewer from the retained records.

The exploratory selection was generated with:

```powershell
.\.venv\Scripts\python.exe -X utf8 research/search_chua_links.py symmetry --digits 60 --subdivisions 96 --seeds 0
.\.venv\Scripts\python.exe -X utf8 research/search_chua_links.py catalogue
.\.venv\Scripts\python.exe -X utf8 research/search_chua_links.py overlaps --digits 60 --subdivisions 96 --seeds 0
```

The recorded catalogue makes S02/S08 case 0. Its refined comparisons are obtained with:

```powershell
.\.venv\Scripts\python.exe -X utf8 research/search_chua_links.py overlaps --cases 0 --digits 80 --subdivisions 192 --seeds 0 2
.\.venv\Scripts\python.exe -X utf8 research/search_chua_links.py overlaps --cases 0 --digits 90 --subdivisions 384 --seeds 0 2 5
.\.venv\Scripts\python.exe -X utf8 research/compare_chua_joint_jones.py --digits 60 --subdivisions 96 --seed 0
foreach ($projection in 0,2) {
    .\.venv\Scripts\python.exe -X utf8 research/compare_chua_joint_jones.py --digits 80 --subdivisions 192 --seed $projection
}
foreach ($projection in 0,2,5) {
    .\.venv\Scripts\python.exe -X utf8 research/compare_chua_joint_jones.py --digits 90 --subdivisions 384 --seed $projection
}
```

The daughter search uses the central atlas node's recorded symmetry boundary. The exploratory probe command below regenerates the source records and should be used in a separate reproduction checkout on D: if their published hashes are to remain unchanged. Existing refinement programs reuse files already present.

```powershell
.\.venv\Scripts\python.exe -X utf8 research/probe_chua_knotted_daughters.py --stage probe --families S01 S02 --amplitudes .001 .003
.\.venv\Scripts\python.exe -X utf8 research/probe_chua_knotted_daughters.py --stage analyze --families S01 S02 --amplitudes .001 --digits 60 --subdivisions 96 --seeds 0 2
.\.venv\Scripts\python.exe -X utf8 research/probe_chua_knotted_daughters.py --stage analyze --families S01 S02 --amplitudes .001 --digits 80 --subdivisions 192 --seeds 0 2
.\.venv\Scripts\python.exe -X utf8 research/probe_chua_knotted_daughters.py --stage analyze --families S01 S02 --amplitudes .001 --digits 90 --subdivisions 384 --seeds 0 2
.\.venv\Scripts\python.exe -X utf8 research/probe_chua_knotted_daughters.py --stage analyze --families S01 S02 --amplitudes .001 --digits 90 --subdivisions 768 --seeds 0
```

The small decimal β offsets and the amplitude constraint select numerical branch points. The later interval calculation proves cycles at those fixed β values; it does not certify the entire bifurcation continuation. The protocol's broad extension of missing stability boundaries was not performed.

## Scope

These data support three non-split two-component polygonal links made from coexisting attracting-cycle approximations. They do not establish an exact link type for the ODE curves, a Borromean or Brunnian triple, a complete inventory of attracting links, or historical novelty. The earlier negative four-cycle experiment remains valid at its own parameter point.
