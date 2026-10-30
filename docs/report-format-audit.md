# GitHub report rendering audit

Initial published-page check: 29 October 2026, 23:49 UTC, at [`04e00a3`](https://github.com/Griboedow/chua-gpu-atlas/commit/04e00a3c9c0a6856be225d08f593e93169966340). Second pass completed 30 October 2026, 00:00 UTC, at [`e2c70b1`](https://github.com/Griboedow/chua-gpu-atlas/commit/e2c70b18a2e00d523a47c7eb54e68ee02f101b60); third pass completed at 00:34 UTC, at [`d73464c`](https://github.com/Griboedow/chua-gpu-atlas/commit/d73464c9e7d1dbfec35cc9ae377e9a1d42f29e87).

**Nineteen distinct formulas in eight reports needed formatting repairs.** The mathematical expressions and research conclusions were retained. The problems were visible in GitHub's actual renderer, not merely suspected from the Markdown source.

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

The repairs are in commits `c9daf5b` and `04e00a3`. The pair-integral antiderivative needed both an escaping repair and an operator repair; it is counted once in the total of nineteen formulas.

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

Across all three passes, **32 distinct pages** have been checked, with **301 math elements across their latest checked versions**. The record contains **40 successful page checks** (23 + 10 + 7) and 356 formula observations (177 + 111 + 68). The latest-version totals count each page once: eighteen pages retain their `04e00a3` check, seven their `e2c70b1` check and seven their `d73464c` check. The exact inspected revision is retained for each page; later edits are not implicitly certified.

At `d73464c`, all **189 relative link occurrences** across the 32 documents resolve to files or directories in its published Git tree. This is a path check, not a certification of external destinations or fragment anchors. The browser checks do not certify every browser or mobile layout, or the scientific content of the figures. Future changes need their own pass.

The [machine-readable record](../results/candidate-audits/report-format-audit.json) lists every checked page, the corrections, and the expression-preservation checks. Downloaded GitHub HTML and UI screenshots remain in the ignored local `outputs/report-format-audit/` directory.

## Conventions for further reports

Use fenced `math` blocks for display equations, especially matrices, subscripts containing asterisks, and formulas with TeX spacing commands. Keep a blank line around each block. Use standard mathematical commands where available, and verify less-than signs inside subscripts and custom operators in the published renderer. Do not treat successful Markdown conversion as proof that the browser completed mathematical typesetting.

GitHub documents both dollar delimiters and fenced math blocks, and identifies MathJax as its renderer. The repairs above follow the fenced form and were checked against the deployed implementation. [GitHub mathematical-expression documentation](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
