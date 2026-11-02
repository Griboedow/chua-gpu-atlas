# Combined layer cards

The [animation](../../docs/assets/chua-layer-cards.gif) and [34 static cards](../../docs/chua-layer-cards.md) associate the published parameter samples with reference motions and numerical knot types. This is a new presentation of existing results, with no new orbit or mesh calculation.

## Contents and sources

| File | Contents |
|---|---|
| [cards.json](cards.json) | All displayed characteristics, exact reference β and period strings, source paths, mesh coverage, frame order and duration |
| [cards.html](cards.html) | Family selector and optional playback; open from a local repository checkout, preserving relative paths to the card images |
| [Static cards](../../docs/assets/chua-layer-cards/) | One 1680 × 1278 PNG per reference entry, plus an introductory frame |
| [GIF](../../docs/assets/chua-layer-cards.gif) | One twelve-second introduction, then 34 fixed-camera family cards at eight seconds each; A04 and A01 first, then remaining entries in identifier order |

The introductory frame and the banner on every family card distinguish the finite catalogue from the [proved infinite cascade of attracting windows](../../docs/chua-infinity-certificate.md). The statement concerns the fixed transect α = 9, m₀ = −8/7, m₁ = −5/7, γ = 0 and the certified β* interval. It uses an interval certificate and the classical sink-cascade theorem; it does not count globally disconnected layer surfaces. This presentation reuses the established proof and does not introduce a new certificate. The viewer selects family cards directly, while the GIF begins with the introduction.

The original 32 representative surfaces come from the [README mesh refinement](../chua-layer-atlas/readme-layer-view/README.md). The two additional families, P75 and P85, use the representative roots in `../chua-followup/sheets/beta_01-search.json`, `beta_11-search.json` and `beta_21-search.json`: α = 9, m₁ = −5/7 and m₀ = −8/7 plus −0.0004, 0 and +0.0004. Lines connect these three samples only; no two-dimensional mesh is inferred for them.

The [extended registry](../chua-layer-prediction/extended-family-registry.json) supplies peak and section-return counts. The [knot data](../chua-knot-atlas/README.md) supply refined parameter strings, the section-aligned 60-digit trajectories, classification review and standard embeddings. The producer checks closure, full-period duration, peak and upward x = 0 crossing counts, symmetry labels, common parameters and the central P75/P85 continuation points before rendering. It does not replay the original numerical or diagram audits.

ρ is the transverse spectral radius at the selected reference: squaring the normalized +0.5 half-return multiplier gives 0.25 for symmetric cycles; asymmetric cycles use the full-return multiplier −0.5 and hence ρ = 0.5. These values specify a reference convention, not the stability of all mesh nodes. P65's mesh retains its documented negative half-return representative; its portrait uses the normalized reference from the knot study.

## Reproduction

Use the repository environment with the dependencies documented in the [knot reproduction guide](../chua-knot-atlas/README.md#reproduction). Run locally; the script puts its Matplotlib cache under the repository's ignored `.publication` directory.

```bash
python research/render_chua_layer_cards.py
```

For a layout check of selected cards without replacing the GIF, metadata or gallery:

```bash
python research/render_chua_layer_cards.py --preview overview A04 E06 S01 P85
```

This writes only the selected PNGs. A full run is required to synchronize the GIF and gallery after a layout change. The original mesh data, trajectories, knot diagrams and earlier animations are unchanged.

The [reading guide](../../docs/chua-layer-cards.md#reading-a-card) defines the catalogue groups, section, symmetry and knot interpretation. The [manifest](../chua-layer-atlas/manifest.json) covers the producer, data associations and figures. Knot topology of the ODE orbits remains numerical; this presentation adds no interval isotopy certificate.
