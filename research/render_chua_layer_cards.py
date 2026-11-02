"""Join the published layer meshes, reference motions and numerical knot types."""
import argparse
import io
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.publication/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
from scipy.signal import find_peaks

from plot_readme_layers import geometry, LOW, HIGH
from plot_chua_knots import TYPES, NAMES, COLORS, knot_drawing, portrait

OUT = ROOT / 'results/chua-layer-cards'
ASSETS = ROOT / 'docs/assets'
CARDS = ASSETS / 'chua-layer-cards'
MESH = ROOT / 'results/chua-layer-atlas/readme-layer-view'
KNOTS = ROOT / 'results/chua-knot-atlas'
BG, INK, MUTED = '#fafaf7', '#182d39', '#53636b'
DWELL = 8000
INTRO_DWELL = 12000
CARD_HEIGHT = 10.65
GROUPS = {'S': 'Symmetric parent family', 'A': 'Asymmetric parent family',
          'E': 'Additional off-centre family', 'P': 'Period-adding family'}
KNOT_GROUPS = ['Unknot', 'Prime torus knot', 'Prime torus knot',
               'Composite knot', 'Composite knot']


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def rasterize(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120, facecolor=BG)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert('RGB')


def overview():
    """State the established infinity result without counting global surfaces."""
    fig = plt.figure(figsize=(14, CARD_HEIGHT), dpi=120, facecolor=BG)
    fig.text(.055, .94, 'Finite catalogue and infinite accumulation', fontsize=29, weight='bold')
    fig.text(.055, .895, 'The displayed reference motions are selected examples, not an exhaustive atlas.',
             fontsize=18, color=MUTED)
    fig.text(.25, .665, '34', fontsize=108, ha='center', color=INK)
    fig.text(.745, .665, '∞', fontsize=120, ha='center', color='#287f9c')
    fig.text(.25, .61, 'computed reference cycles', fontsize=23, ha='center', weight='bold')
    fig.text(.745, .61, 'attracting windows proved', fontsize=23, ha='center', weight='bold', color='#287f9c')
    fig.text(.25, .525, 'Explicit orbits, parameter samples,\nwaveforms and numerical knot types.',
             fontsize=17, ha='center', color=MUTED, linespacing=1.6)
    fig.text(.745, .525, 'Infinitely many stable periodic windows\naccumulate as β approaches β*.',
             fontsize=17, ha='center', color=MUTED, linespacing=1.6)
    for y in [.455, .135]:
        fig.add_artist(Line2D([.055, .945], [y, y], transform=fig.transFigure, color='#d3dcde', lw=1))
    fig.text(.07, .408, 'PROVED ON THIS PARAMETER TRANSECT', fontsize=14, weight='bold')
    fig.text(.07, .355, 'α = 9,   m₀ = −8/7,   m₁ = −5/7,   γ = 0', fontsize=21)
    fig.text(.07, .298, '14.30635909996140  <  β*  <  14.30635909996152', fontsize=20, color='#287f9c')
    fig.text(.07, .204, 'Interval certificate of a quadratic homoclinic tangency\n'
             '+ classical theorem on cascades of attracting periodic orbits.', fontsize=18, linespacing=1.6)
    fig.text(.07, .088, 'These attracting windows have local extensions in parameter space.', fontsize=17)
    fig.text(.07, .043, 'Global connections of the layer surfaces remain unclassified.', fontsize=17, color=MUTED)
    return rasterize(fig)


def add_scope_banner(fig):
    # Add vertical space; preserve the original card's physical dimensions.
    scale = 10 / CARD_HEIGHT
    fig.set_size_inches(14, CARD_HEIGHT)
    for ax in fig.axes:
        x, y, width, height = ax.get_position(original=True).bounds
        ax.set_position([x, y*scale, width, height*scale])
    for label in fig.texts:
        x, y = label.get_position()
        label.set_position((x, y*scale))
    for line in fig.artists:
        if isinstance(line, Line2D):
            line.set_ydata(np.asarray(line.get_ydata())*scale)
    fig.text(.035, .962, '∞ attracting windows proved  ·  34 reference cycles shown',
             fontsize=19, weight='bold', color='#287f9c')
    fig.add_artist(Line2D([.035, .965], [.942, .942], transform=fig.transFigure, color='#d3dcde', lw=1))


def load():
    registry = {r['id']: r for r in read(ROOT / 'results/chua-layer-prediction/extended-family-registry.json')['families']}
    review = read(KNOTS / 'review.json')
    records, meshes, partial = [], {}, {}
    assert read(MESH / 'audit.json')['status'] == 'passed'
    for r in review['families']:
        id = r['id']
        root = read(KNOTS / f'{id}-d70-n192.json')
        arrays = np.load(KNOTS / f'{id}-d60-n96.npz')
        xyz, t = arrays['xyz'], arrays['t']
        source = registry[id]
        assert r['peaks'] == source['peak_count'] == len(find_peaks(xyz[:, 0])[0]), id
        section_x = xyz[:, 0].copy()
        section_x[np.abs(section_x) < 1e-12] = 0
        assert int(np.sum((section_x[:-1] <= 0) & (section_x[1:] > 0))) == source['section_returns'], id
        assert r['symmetric'] == root['symmetric'] == source['symmetric'], id
        assert np.max(np.abs(xyz[0] - xyz[-1])) < 1e-10, id
        assert np.all(np.diff(t) > 0) and abs(t[-1] / float(root['period']) - 1) < 1e-12, id
        assert [root[k] for k in ['alpha', 'm0', 'm1', 'gamma']] == ['9', '-8/7', '-5/7', '0'], id
        both = bool(xyz[:, 0].min() < -1 and xyz[:, 0].max() > 1)
        target = float(root['target_multiplier'])
        assert target == (.5 if r['symmetric'] else -.5), id
        group = TYPES.index(r['knot_type'])
        rec = dict(id=id, catalogue_group=GROUPS[id[0]], knot_type=r['knot_type'],
                   knot_name=NAMES[group], knot_group=KNOT_GROUPS[group], color=COLORS[group],
                   beta=root['beta'], period=root['period'], peaks=r['peaks'],
                   section_returns=source['section_returns'], symmetric=r['symmetric'],
                   both_outer_regions=both, full_return_spectral_radius=target**2 if r['symmetric'] else abs(target),
                   x_range=[float(xyz[:, 0].min()), float(xyz[:, 0].max())],
                   reference_in_displayed_box=LOW <= float(root['beta']) <= HIGH,
                   sources=[str((KNOTS / f'{id}-d60-n96.npz').relative_to(ROOT)).replace('\\', '/'),
                            str((KNOTS / f'{id}-d70-n192.json').relative_to(ROOT)).replace('\\', '/')])
        if (MESH / f'{id}.json').exists():
            sheet = read(MESH / f'{id}.json')
            meshes[id] = geometry(sheet)
            rec.update(coverage='original refined mesh', displayed_samples=len(meshes[id]['points']),
                       mesh_target_multiplier=sheet['target_multiplier'])
            rec['sources'].append(f'results/chua-layer-atlas/readme-layer-view/{id}.json')
        else:
            assert id in ['P75', 'P85'], id
            pts = []
            for name in ['beta_01', 'beta_11', 'beta_21']:
                path = ROOT / f'results/chua-followup/sheets/{name}-search.json'
                data = read(path)
                row = next(s for s in data['records'] if s['peaks'] == r['peaks'])
                par = data['parameters']
                assert row['status'] == 'numerical_root' and par['alpha'] == 9
                assert abs(par['m1'] + 5/7) < 1e-15
                pts.append([par['alpha'], 1000*(par['m0'] + 8/7), float(row['representative']['coordinate'])])
                rec['sources'].append(path.relative_to(ROOT).as_posix())
            partial[id] = np.array(pts)
            assert np.max(np.abs(partial[id][1] - [9, 0, float(root['beta'])])) < 1e-12
            rec.update(coverage='three fixed-alpha continuation samples; no dense mesh', displayed_samples=3,
                       continuation_points=pts)
        records.append(rec)
    assert len(records) == 34 and len(meshes) == 32 and len(partial) == 2
    # Put the two equal-peak, different-knot examples first; keep the rest in ID order.
    records.sort(key=lambda r: (['A04', 'A01'].index(r['id']) if r['id'] in ['A04', 'A01'] else 2, r['id']))
    standard = np.load(KNOTS / 'standard-knot-curves.npz')
    shapes = {name: standard[f'k{i}'] for i, name in enumerate(read(KNOTS / 'controls.json')['types'])}
    return records, meshes, partial, shapes


def frame(r, index, meshes, partial, shapes):
    color, id = r['color'], r['id']
    fig = plt.figure(figsize=(14, 10), dpi=120, facecolor=BG)
    fig.text(.035, .95, f"{id}  ·  {r['catalogue_group']}", fontsize=25, weight='bold')
    motion = 'Both outer regions' if r['both_outer_regions'] else 'One outer region'
    symmetry = 'half-period inversion symmetry' if r['symmetric'] else 'asymmetric motion'
    fig.text(.035, .908, f"{motion}  ·  {symmetry}  ·  {r['knot_group'].lower()}", fontsize=16, color=MUTED)
    fig.text(.965, .955, f'{index:02d} / 34', fontsize=16, ha='right', color=MUTED)
    for y in [.881, .424, .092]:
        fig.add_artist(Line2D([.035, .965], [y, y], transform=fig.transFigure, color='#d3dcde', lw=1))

    fig.text(.035, .847, 'PARAMETER LAYER', fontsize=14, weight='bold')
    fig.text(.535, .847, 'ATTRACTING CYCLE', fontsize=14, weight='bold')
    fig.text(.79, .847, 'KNOT TYPE', fontsize=14, weight='bold')
    ax = fig.add_axes([.015, .475, .48, .36], projection='3d', computed_zorder=False)
    points = np.vstack([g['points'][::3] for g in meshes.values()])
    ax.scatter(*points.T, s=.8, c='#c4ced3', alpha=.15, depthshade=False, zorder=1)
    if id in meshes:
        g = meshes[id]
        ax.add_collection3d(Poly3DCollection(g['faces'], facecolors=color, alpha=.80,
                                            edgecolors=color, linewidths=.16, zorder=2))
        ax.scatter(*g['points'].T, s=5, c=color, depthshade=False, zorder=3)
    else:
        pts = partial[id]
        ax.plot(*pts.T, color=color, marker='o', ms=5, lw=2.5)
        ax.text2D(.06, .83, '3 continuation samples\nDense mesh unavailable', transform=ax.transAxes,
                  fontsize=13, color=INK, bbox=dict(facecolor=BG, edgecolor='none', alpha=.94))
    if r['reference_in_displayed_box']:
        ax.scatter([9], [0], [float(r['beta'])], marker='*', s=130, c=INK, depthshade=False, zorder=20)
    ax.set(xlim=(8.94, 9.06), ylim=(-6, 6), zlim=(LOW, HIGH),
           xlabel='α', ylabel='10³ Δm₀', zlabel='β')
    ax.set_xticks([8.94, 9, 9.06]); ax.set_yticks([-4, 0, 4]); ax.set_zticks([14.24, 14.32, 14.4])
    ax.tick_params(labelsize=11, pad=0)
    ax.xaxis.labelpad=5; ax.yaxis.labelpad=8; ax.zaxis.labelpad=5
    ax.view_init(23, -55); ax.set_box_aspect([1, 1, 1.10])
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_pane_color((.97, .98, .98, 1))
    fig.text(.035, .439, '★  Reference orbit' if r['reference_in_displayed_box'] else 'Reference β lies outside this displayed box',
             fontsize=12, color=MUTED)

    arrays = np.load(KNOTS / f'{id}-d60-n96.npz')
    xyz, t = arrays['xyz'], arrays['t']
    ax = fig.add_axes([.49, .49, .29, .33], projection='3d')
    portrait(ax, xyz, color)
    ax.lines[0].set_linewidth(.9)
    fig.text(.535, .475, f"β = {float(r['beta']):.14f}", fontsize=14)
    fig.text(.535, .445, 'Full cycle · axes rescaled', fontsize=12, color=MUTED)
    ax = fig.add_axes([.79, .535, .175, .25])
    knot_drawing(ax, shapes[r['knot_type']], color)
    fig.text(.877, .485, r['knot_name'], fontsize=16, color=color, weight='bold', ha='center')
    fig.text(.877, .452, 'Standard embedding', fontsize=12, color=MUTED, ha='center')

    fig.text(.035, .389, 'SLICE AT α = 9', fontsize=14, weight='bold')
    ax = fig.add_axes([.084, .185, .355, .168])
    for other, g in meshes.items():
        active = other == id
        ax.plot(g['offsets'], g['beta'][8], color=color if active else '#cfd6d8',
                lw=2.4 if active else .7, marker='.' if active else None, ms=4,
                zorder=3 if active else 1)
    if id in partial:
        ax.plot(partial[id][:, 1], partial[id][:, 2], color=color, lw=2.4, marker='o', ms=5)
    if r['reference_in_displayed_box']:
        ax.scatter([0], [float(r['beta'])], marker='*', s=110, c=INK, zorder=5)
    ax.set(xlim=(-6, 6), ylim=(LOW, HIGH), xlabel='1000 (m₀ + 8/7)', ylabel='β')
    ax.set_xticks([-6, -3, 0, 3, 6]); ax.set_yticks([14.24, 14.32, 14.4])
    ax.tick_params(labelsize=11); ax.grid(alpha=.15)
    count = r['displayed_samples']
    fig.text(.035, .108, f'{count} plotted samples · ' + ('refined mesh' if id in meshes else 'partial coverage'),
             fontsize=13, color=MUTED)

    for x, number, label in [(.535, str(r['peaks']), 'peaks of x'),
                              (.69, str(r['section_returns']), 'section returns'),
                              (.85, f"{float(r['period']):.4f}", 'period T')]:
        fig.text(x, .374, number, fontsize=23, weight='bold', color=color)
        fig.text(x, .346, label, fontsize=12, color=MUTED)
    ax = fig.add_axes([.57, .188, .395, .117])
    ax.plot(t/t[-1], xyz[:, 0], color=color, lw=1.05)
    for level in [-1, 1]:
        ax.axhline(level, color='#a5b2b9', lw=.7, ls='--', zorder=0)
    ax.set(xlim=(0, 1), xlabel='t / T', ylabel='x')
    ax.set_xticks([0, .5, 1]); ax.tick_params(labelsize=11)
    fig.text(.535, .108, f"Attracting reference · ρ = {r['full_return_spectral_radius']:g} (full return)", fontsize=13, color=MUTED)
    fig.text(.035, .059, 'Reference: α = 9, m₀ = −8/7, m₁ = −5/7, γ = 0.  Each family has its own β.', fontsize=13)
    fig.text(.035, .028, 'Facets guide the eye; thickness is not stability width. Knot type is numerical, not interval-certified.', fontsize=12, color=MUTED)
    add_scope_banner(fig)
    return rasterize(fig)


def write_gallery(records):
    ids = sorted(r['id'] for r in records)
    intro = '''# Layer, periodic orbit and knot type

Each card joins a sampled parameter layer to one reference attracting cycle, its waveform and its numerically identified knot type. The **34 reference cycles use different β values**; they are not 34 simultaneous attractors. [Play the slow animation](assets/chua-layer-cards.gif), or choose a static card below.

**34 is the number of displayed examples, not the total number of periodic motions.** On the transect α = 9, m₀ = −8/7, m₁ = −5/7, γ = 0, [infinitely many attracting periodic windows are proved](chua-infinity-certificate.md) to accumulate at β* in (14.30635909996140, 14.30635909996152). The result combines an interval certificate with the classical sink-cascade theorem. These windows have local extensions in parameter space, but their global connections have not been classified: the result is not a count of disconnected layer surfaces. The animation starts with this distinction and retains an infinity banner on every family card. Its [introductory frame](assets/chua-layer-cards/overview.png) is also available separately.

The 32 original families have refined meshes. P75 and P85 show only three existing continuation samples at fixed α and m₁; no dense sheet is inferred for them. Grey points and lines retain the other families as context. The star marks the reference parameters when they lie in the displayed box. The orbit and waveform use the reference parameters printed on the card, not every point of the highlighted surface.

## Reading a card

- **Catalogue group:** S denotes the original symmetric parent families; A the asymmetric parents; E the additional families found away from the central transect; P the later period-adding entries. These labels record catalogue provenance. The symmetry line describes the actual reference motion. S04–S07 also belong to the period-adding sequence.
- **Motion:** the orbit visits one or both outer affine regions, x < −1 and x > 1. Half-period inversion symmetry means u(t + T/2) = −u(t). An asymmetric cycle has an inversion-related partner; the partner has the mirror knot type.
- **Counts and period:** peaks are local maxima of x over one full cycle. Section returns count upward crossings of x = 0. T is the dimensionless full period. The waveform shows one period, with dashed lines at x = ±1. Orbit coordinates are independently rescaled to reveal shape. The layer's slope offset is Δm₀ = m₀ + 8/7.
- **Stability:** ρ is the spectral radius of the transverse full-period return map. Its values 0.25 (symmetric) and 0.5 (asymmetric) follow the selected reference normalization; they do not rank whole families by stability. P65's displayed mesh retains a different reference multiplier, as documented in the [mesh audit](../results/chua-layer-atlas/readme-layer-view/README.md).
- **Knot:** the drawing is a standard embedding of the identified knot, not another circuit trajectory. Orbit topology remains numerical. Exact diagram invariants and the sampling checks are in the [knot report](chua-knot-atlas.md). Equal knot types, peak counts or catalogue groups do not identify a family uniquely.

Facets and straight segments indicate sample adjacency. They neither certify the unsampled interior nor measure a stability interval's width. Reference parameters outside the plotted β range are labelled explicitly. Values printed on a card are rounded; [source records and reproduction](../results/chua-layer-cards/README.md) retain the full precision.

## Choose a family

'''
    intro += ' · '.join(f'[{id}](#{id.lower()})' for id in ids) + '\n'
    lookup = {r['id']: r for r in records}
    for id in ids:
        r = lookup[id]
        intro += f"\n## {id}\n\n![{id}: {r['peaks']} peaks, {r['knot_name']}; layer, orbit and waveform](assets/chua-layer-cards/{id}.png)\n"
    (ROOT / 'docs/chua-layer-cards.md').write_text(intro, encoding='utf-8', newline='\n')

    html = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Chua layer cards</title>
<style>body{margin:0;background:#fafaf7;color:#182d39;font:17px system-ui}main{max-width:1400px;margin:auto;padding:20px}h1{font-size:28px;margin:0 0 8px}p{line-height:1.5}select,button{font:inherit;padding:8px;margin:4px}img{width:100%;height:auto}nav{display:flex;gap:8px;align-items:center;flex-wrap:wrap}footer{font-size:15px}</style>
<main><h1>Layer, periodic orbit and knot type</h1><p><strong>Infinitely many attracting periodic windows are proved; these 34 reference motions are selected examples.</strong> Global connections of the layer surfaces remain unclassified. <a href="https://github.com/Griboedow/chua-gpu-atlas/blob/main/docs/chua-infinity-certificate.md">Proof and parameter transect</a>.</p><p>Select a family, or play one card every eight seconds. The camera remains fixed; each reference uses its own β.</p><nav><label for="family">Family</label><select id="family"></select><button id="previous">Previous</button><button id="next">Next</button><button id="play">Play</button></nav><img id="card" width="1680" height="1278" alt=""><footer><p>32 refined meshes; P75 and P85 have only three continuation samples each. Knot types are numerical. <a href="../../docs/chua-layer-cards.md">Reading guide and static cards</a> · <a href="../../docs/assets/chua-layer-cards.gif">GIF</a></p><p>This viewer uses the card images in the repository. Keep the relative folder structure when opening it offline.</p></footer></main>
<script>const cards=__DATA__;const select=document.getElementById('family'),picture=document.getElementById('card'),play=document.getElementById('play');let timer=null;
cards.forEach((r,i)=>{const o=document.createElement('option');o.value=i;o.textContent=r.id+' · '+r.peaks+' peaks · '+r.knot_name;select.appendChild(o)});
function draw(){const r=cards[+select.value];picture.src='../../docs/assets/chua-layer-cards/'+r.id+'.png';picture.alt=r.id+': '+r.catalogue_group+', '+r.peaks+' peaks, '+r.knot_name;}
function move(delta){select.value=(+select.value+delta+cards.length)%cards.length;draw()}
function stop(){clearInterval(timer);timer=null;play.textContent='Play'}
select.addEventListener('change',()=>{stop();draw()});document.getElementById('previous').onclick=()=>{stop();move(-1)};document.getElementById('next').onclick=()=>{stop();move(1)};play.onclick=()=>{if(timer)stop();else{timer=setInterval(()=>move(1),8000);play.textContent='Pause'}};draw();</script></html>'''
    html = html.replace('__DATA__', json.dumps(records, separators=(',', ':')))
    (OUT / 'cards.html').write_text(html, encoding='utf-8', newline='\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview', nargs='+', help='Render selected cards only, without publishing the GIF or metadata.')
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True); CARDS.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 13, 'text.color': INK,
                         'axes.labelcolor': INK, 'axes.edgecolor': MUTED, 'axes.facecolor': BG,
                         'axes.spines.top': False, 'axes.spines.right': False})
    records, meshes, partial, shapes = load()
    frames = []
    if not args.preview or 'overview' in args.preview:
        intro = overview()
        intro.save(CARDS / 'overview.png')
        frames.append(intro.quantize(colors=256))
    for index, r in enumerate(records, 1):
        if args.preview and r['id'] not in args.preview:
            continue
        im = frame(r, index, meshes, partial, shapes)
        im.save(CARDS / f"{r['id']}.png")
        frames.append(im.quantize(colors=256))
        print('Rendered', index, r['id'], flush=True)
    if args.preview:
        return
    frames[0].save(ASSETS / 'chua-layer-cards.gif', save_all=True, append_images=frames[1:],
                   duration=[INTRO_DWELL] + [DWELL]*len(records), loop=0, optimize=False, disposal=2)
    metadata = dict(families=records, frame_order=['overview'] + [r['id'] for r in records],
                    frame_durations_ms=[INTRO_DWELL] + [DWELL]*len(records),
                    reference_count=len(records), frame_count=len(frames), image_size=list(frames[0].size),
                    infinity_statement='Infinitely many attracting periodic windows accumulate at the certified tangency '
                        'on alpha=9, m0=-8/7, m1=-5/7, gamma=0. The catalogue is finite; global layer connections remain unclassified.',
                    infinity_source='docs/chua-infinity-certificate.md',
                    scope='Presentation of existing numerical results; no new orbit or mesh computation. '
                          '32 refined meshes and two three-point sections. Knot types remain numerical.')
    (OUT / 'cards.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8', newline='\n')
    write_gallery(records)
    print('Published', len(frames), 'cards and the fixed-camera GIF', flush=True)


if __name__ == '__main__':
    main()
