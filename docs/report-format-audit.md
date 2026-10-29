# GitHub report rendering audit

Checked 29 October 2026, 23:49 UTC, at publication revision [`04e00a3`](https://github.com/Griboedow/chua-gpu-atlas/commit/04e00a3c9c0a6856be225d08f593e93169966340).

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

All **130 local relative links** checked before this audit note was added resolved to existing paths. No loaded article image reported a failed load during the browser pass. This does not certify every external link, every browser or mobile layout, or the scientific content of the figures. Newly created reports and future changes need their own pass.

The [machine-readable record](../results/candidate-audits/report-format-audit.json) lists every checked page, the corrections, and the expression-preservation checks. Downloaded GitHub HTML and UI screenshots remain in the ignored local `outputs/report-format-audit/` directory.

## Conventions for further reports

Use fenced `math` blocks for display equations, especially matrices, subscripts containing asterisks, and formulas with TeX spacing commands. Keep a blank line around each block. Use standard mathematical commands where available, and verify less-than signs inside subscripts and custom operators in the published renderer. Do not treat successful Markdown conversion as proof that the browser completed mathematical typesetting.

GitHub documents both dollar delimiters and fenced math blocks, and identifies MathJax as its renderer. The repairs above follow the fenced form and were checked against the deployed implementation. [GitHub mathematical-expression documentation](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
