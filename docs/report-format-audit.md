# GitHub report rendering audit

Initial published-page check: 29 October 2026, 23:49 UTC, at [`04e00a3`](https://github.com/Griboedow/chua-gpu-atlas/commit/04e00a3c9c0a6856be225d08f593e93169966340). Subsequent passes on 30 October completed at 00:00 UTC ([`e2c70b1`](https://github.com/Griboedow/chua-gpu-atlas/commit/e2c70b18a2e00d523a47c7eb54e68ee02f101b60)), 00:34 UTC ([`d73464c`](https://github.com/Griboedow/chua-gpu-atlas/commit/d73464c9e7d1dbfec35cc9ae377e9a1d42f29e87)), 01:16 UTC ([`1dc034a`](https://github.com/Griboedow/chua-gpu-atlas/commit/1dc034ad51f993509f70c06229473a3dbfd81c70)) and 02:07 UTC ([`ed9517e`](https://github.com/Griboedow/chua-gpu-atlas/commit/ed9517e4eb68d9336e958c553aa669a5f15c13e7)).

The HHS and microfluidic follow-up completed at 02:56 UTC on 30 October, using [`ecfd8a5`](https://github.com/Griboedow/chua-gpu-atlas/commit/ecfd8a5672442c9b884e839b2826f2777d68207e) and the two-page repair revision [`a99a267`](https://github.com/Griboedow/chua-gpu-atlas/commit/a99a2675f53135bc2eb112a3d7ded34716d862ea).

**Twenty-three distinct expressions in ten reports needed formatting repairs.** This includes a single-symbol formula in a bold heading, replaced by the same symbol as plain text. The mathematical expressions and research conclusions were retained. The problems were visible in GitHub's actual renderer, not merely suspected from the Markdown source.

## What was repaired

| Report | Rendering problem | Repair |
| --- | --- | --- |
| [First research pass](candidate-research.md) | Eight displays used delimiters that GitHub treated as ordinary text | Fenced math blocks |
| [LCE slider](lce-slider-research.md) | TeX spacing commands became visible punctuation | Fenced math block |
| [Photoresponse feedback](photofeedback-research.md) | A 3 × 3 Jacobian became one row with red unknown-command tokens | Fenced math block preserving row separators |
| [Rheostat oscillator](rheostat-research.md) | Matrix row separators and a spacing command were consumed by Markdown | Fenced math blocks |
| [Synchronematic pair integral](synchronematic-pair-integral.md) | Spacing commands were consumed; the renderer rejected a custom operator macro | Fenced math blocks and the standard arctangent command |
| [Synchronematic model](synchronematic-research.md) | Norm delimiters lost their escapes; a summation subscript caused a brace error | Fenced math blocks and an explicit less-than command |
| [Thermomechanical sheet](thermomechanical-sheet-audit.md) | Markdown emphasis altered an equilibrium subscript; spacing appeared as punctuation | Fenced math blocks |
| [Two-strain model](two-strain-research.md) | The renderer rejected a custom divergence operator macro | Upright divergence notation with explicit spacing |
| [HHS structural audit](hhs-structure.md) | One bold heading and two inline inequalities retained literal dollar/TeX text | Plain D in the heading; braced asterisk commands in the two inline subscripts |
| [Microfluidic screening](late-candidate-screening-4.md) | Markdown emphasis consumed two asterisk superscripts in one table-cell expression | Braced asterisk commands preserving the superscripts |

The initial nineteen repairs are in commits `c9daf5b` and `04e00a3`; the four later repairs are in `a99a267`. The pair-integral antiderivative needed both an escaping repair and an operator repair; it is counted once.

## Verification

The check loaded **23 public GitHub pages**—the README and the 22 reports that existed at the start of this audit—in Microsoft Edge through headless Playwright. GitHub's own JavaScript and math renderer executed in the browser. This was a real browser check, using the published revision, rather than a local Markdown renderer.

Across those pages, **177 math elements** produced MathML. No unresolved formulas, renderer error panels, red unknown-command tokens, or math-container overflow were found at the checked 1320 × 1000 viewport. Representative formula screenshots were also inspected, including the Jacobian before and after repair. Merely counting math elements would have missed several failures: GitHub can place an error panel or red command text inside a math element.

All **130 local relative links** checked before this audit note was added resolved to existing paths. No loaded article image reported a failed load during the initial browser pass.

### Published follow-up

The second pass loaded ten actual GitHub pages at `e2c70b1`: the four new spin-crossover, enzyme-feedback, core-shell-gel and synchronematic-priority reports; the README, candidate queue, synchronematic research, attractor geometry and synchronematic bifurcation reports; and this audit note. Its **111 math elements** rendered without errors, unresolved formulas, red command tokens or math-container overflow. Visible article images were scrolled into view and checked for successful loading. GitHub also inserts hidden GIF preview elements; these were distinguished from the visible figures rather than treated as failed images.

Screenshots from the published pages were inspected for the trace formulas, the aligned four-species equations, the enzyme-gel survey figure, the corrected collective projection and the bifurcation figure. Seven spin-crossover displays had been prepared with math fences before publication, including equivalent upright trace notation. These preventive edits are separate from the nineteen previously observed broken formulas. Local in-memory previews were used during preparation but are not counted as published-page checks.

After the second pass, the record covered 28 distinct pages, with 250 math elements across their latest checked versions. Its 33 successful page checks included five repeated pages. All 164 relative link occurrences in those documents resolved against the published Git tree at `e2c70b1`.

### Urease, esterase and thermofluidic follow-up

The third pass loaded seven published pages at `d73464c`: the new urease research, urease structural audit, esterase compartment audit and physical candidate screening reports, together with the changed README, candidate queue and enzyme-feedback report. All **68 math elements** rendered without errors, unresolved formulas, red command tokens or math-container overflow. Visible images loaded successfully. No formatting repairs were needed in this pass.

Screenshots were inspected for the NIFTE two-by-two matrices and column vectors, the urease reduction and history integral, and the aligned esterase equations. A separate screenshot of the complete published urease survey figure confirmed that all four panels, labels and legends were displayed. This additional screenshot is supporting evidence for the same page check, not an extra page in the totals.

After the third pass, the record covered 32 distinct pages, with 301 math elements across their latest checked versions. It contained 40 successful page checks and 356 formula observations including repeats. All 189 relative link occurrences in those documents resolved against the published Git tree at `d73464c`.

### Published NIFTE and Rayleigh-pair reports

The fourth pass loaded nine published pages at `1dc034a`: the three new NIFTE reports, the second physical-candidate screening and Rayleigh-pair structural report, `results/README.md`, and the changed main README, candidate queue and first physical-candidate screening. All **74 math elements** rendered without errors, unresolved formulas, red command tokens or math-container overflow. Every visible article image loaded successfully; no formatting repairs were needed.

The inspected screenshots include both complete NIFTE figures, the local cubic-mode coefficient matrix and the Rayleigh-pair Lyapunov inequality with its surrounding derivation. Supporting screenshot visits are not counted as additional page checks. Earlier in-memory previews of the NIFTE drafts remain excluded from the published totals.

After four passes, the record covered **38 distinct pages**, with **362 math elements across their latest checked versions**. It contained **49 successful page checks** (23 + 10 + 7 + 9) and 430 formula observations (177 + 111 + 68 + 74). At that point, eighteen pages retained their `04e00a3` check, seven their `e2c70b1` check, four their `d73464c` check and nine their `1dc034a` check. The exact inspected revision is retained for each page; later edits are not implicitly certified.

At `1dc034a`, all **227 relative link occurrences** across the 38 documents, including `results/README.md`, resolve to files or directories in its published Git tree. This is a path check, not a certification of external destinations or fragment anchors. The browser checks do not certify every browser or mobile layout, or the scientific content of the figures. Future changes need their own pass.

### Rayleigh numerical results and material-model screening

The fifth pass loaded eight public pages at `ed9517e`: the new Rayleigh numerical report, material-memory screening, meminductor follow-up and third candidate-screening reports, together with the main README, results README, candidate queue and updated Rayleigh structural report. All **102 math elements** rendered without errors, unresolved formulas, red command tokens or math-container overflow. This count includes inline expressions as well as display blocks. Every visible article image loaded successfully. No formatting repairs were needed.

The inspected screenshots show all three complete Rayleigh figures: the survey and stable-cycle controls, the weak-coupling saddle/long-transient comparison, and the independent trajectory's approach to and departure from the saddle cycle. The five-by-five trajectory-dependent Jacobian and its divergence formula were visually checked, as were the material model's transformed equations and three signed cyclic derivatives. Earlier in-memory previews are excluded from the published totals.

After five passes, **42 distinct pages** had been checked, with **450 math elements across their latest checked versions**. The record contained **57 successful page checks** (23 + 10 + 7 + 9 + 8) and **532 formula observations** (177 + 111 + 68 + 74 + 102). Each page was counted once in the latest-version totals: eighteen retained their `04e00a3` check, seven their `e2c70b1` check, four their `d73464c` check, five their `1dc034a` check and eight their `ed9517e` check. The fifth pass contained four new pages and four repeated pages; it did not claim that all 42 were revisited at that revision.

All **255 relative link occurrences** across those 42 documents resolve against the published Git tree at `ed9517e`. External destinations, fragment anchors, other viewports and later edits remain outside this check.

### HHS and microfluidic follow-up: detecting unparsed inline math

Eight public pages were inspected at `ecfd8a5`. Six passed: the main README, results README, candidate queue, third screening report, HHS numerical report and HHS source audit. Their **57 math elements** rendered correctly and all visible images loaded. Screenshots of both complete HHS figures were inspected.

The HHS structural report and fourth screening report initially **failed** visual inspection. Three expressions in the HHS proof and one in a microfluidic parameter-table cell remained literal dollar/TeX text. The existing MathML checks alone did not detect these failures because those expressions never became math elements. A scan of the rendered article outside code and math elements confirmed raw remnants on these two pages and none on the other six.

After the four minimal repairs were published, only the affected two pages were rechecked at `a99a267`. Their **159 math elements** (102 and 57) rendered without errors, red command tokens, unresolved formulas or container overflow. Both raw-text scans were clean. The repaired proof inequalities and parameter-table row were visually inspected; the pressure matrix and displayed flow equations had also been inspected in the initial publication. The two failed initial page attempts, containing 156 successfully typeset elements but four unparsed expressions, are recorded separately and excluded from all successful-check totals. Screenshot-only visits and preparation previews are not extra page checks.

The record now covers **46 distinct pages**, with **612 math elements across their latest successfully checked versions**. There are **65 successful page checks** (23 + 10 + 7 + 9 + 8 + 6 + 2) and **748 successful formula observations** (177 + 111 + 68 + 74 + 102 + 57 + 159). The seven revision-specific groups retain eighteen pages at `04e00a3`, seven at `e2c70b1`, four at `d73464c`, five at `1dc034a`, four at `ed9517e`, six at `ecfd8a5` and two at `a99a267`. These counts do not certify every page at the newest revision. The additional raw-text scan covers the eight pages in this block, not the previous 38 pages retroactively.

All **283 relative link occurrences** across those 46 documents resolve against the published Git tree at `a99a267`. The link matcher excludes fenced code/math and inline code before interpreting Markdown links, so mathematical bracket-and-parenthesis expressions are not misidentified as paths. External destinations and fragment anchors remain outside this check.

The [machine-readable record](../results/candidate-audits/report-format-audit.json) lists every checked page, the corrections, and the expression-preservation checks. Downloaded GitHub HTML and UI screenshots remain in the ignored local `outputs/report-format-audit/` directory.

## Conventions for further reports

Use fenced `math` blocks for display equations, especially matrices, subscripts containing asterisks, and formulas with TeX spacing commands. Keep a blank line around each block. Use standard mathematical commands where available, and verify less-than signs inside subscripts and custom operators in the published renderer. For inline asterisks, use braced `\ast` commands when Markdown emphasis would consume a literal asterisk. Check bold labels and table cells visually, and scan for raw dollars or TeX outside rendered math. Do not treat successful Markdown conversion or clean math-element metrics as proof that every intended expression was typeset.

GitHub documents both dollar delimiters and fenced math blocks, and identifies MathJax as its renderer. The repairs above follow the fenced form and were checked against the deployed implementation. [GitHub mathematical-expression documentation](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
