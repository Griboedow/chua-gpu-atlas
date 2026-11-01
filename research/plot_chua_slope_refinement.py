"""Publication figures for the refined slope volumes and grazing experiment."""
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.publication/matplotlib'))
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from refine_chua_slope_atlases import OUT, flights, parameters

ASSETS = ROOT/'docs/assets'
COLORS = {'P5': '#197da1', 'P7': '#c1678a', 'new': '#713e73', 'contact': '#ca8c25'}
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'axes.titleweight': 'semibold',
    'axes.labelcolor': '#253747', 'text.color': '#253747', 'axes.spines.top': False,
    'axes.spines.right': False, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})


def read(name):
    return json.loads((OUT/name).read_text())


def save(fig, name):
    fig.savefig(ASSETS/(name+'.png'), dpi=180)
    fig.savefig(OUT/(name+'.pdf'))
    plt.close(fig)


def volumes():
    fig = plt.figure(figsize=(15, 8))
    fig.subplots_adjust(left=.015, right=.96, bottom=.19, top=.82, wspace=.035)
    grazing = read('grazing-results.json')['contact_curves']
    for col, kind in enumerate(['alpha', 'beta'], 1):
        ax = fig.add_subplot(1, 2, col, projection='3d')
        polygons, colors, points = [], [], []
        for family in ['P5', 'P7']:
            d = read(f'{family}-{kind}.json')
            rows = {tuple(r['index']): r for r in d['nodes']}
            for i in range(40):
                for j in range(40):
                    keys = [(i, j), (i+1, j), (i+1, j+1), (i, j+1)]
                    if not all(k in rows for k in keys):
                        continue
                    vertices = [[rows[k]['core'][p] for p in ['m0', 'm1', kind]] for k in keys]
                    changed = family == 'P7' and sum(len(rows[k]['core']['targets']) > 4 for k in keys) >= 2
                    color = COLORS['new' if changed else family]
                    polygons.extend([[vertices[0], vertices[1], vertices[2]], [vertices[0], vertices[2], vertices[3]]])
                    colors.extend([color]*2)
            points.extend([[r['core'][p] for p in ['m0', 'm1', kind]] for r in d['nodes']])
        ax.add_collection3d(Poly3DCollection(polygons, facecolors=colors, alpha=.82,
            edgecolors=(1, 1, 1, .25), linewidths=.12))
        contact = next(c for c in grazing if c['kind'] == kind and c['target'] == -.43)
        curve = np.array([[r['contact']['orbit'][p] for p in ['m0', 'm1', kind]] for r in contact['nodes']])
        ax.plot(*curve.T, color=COLORS['contact'], lw=2.7, zorder=20)
        values = np.array(points)
        ax.auto_scale_xyz(*values.T)
        ax.set(xlabel=r'$m_0$', ylabel=r'$m_1$', zlabel='α' if kind == 'alpha' else 'β',
            title=r'$(m_0,m_1,\alpha)$; $\beta=100/7$' if kind == 'alpha' else r'$(m_0,m_1,\beta)$; $\alpha=9$')
        ax.set_box_aspect((1.12, 1, 1))
        ax.set_proj_type('ortho')
        ax.view_init(elev=24, azim=-57)
        ax.tick_params(labelsize=9, pad=1)
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            axis.set_major_locator(MaxNLocator(5))
            axis.set_pane_color((.97, .98, .985, 1))
            axis._axinfo['grid']['color'] = (.4, .5, .6, .17)
    fig.suptitle('Continued families in the two slope volumes', x=.5, y=.96, fontsize=23, weight='semibold')
    fig.text(.5, .89, '6,713 checked cross-sections  /  6,724 requested     ·     fold, symmetry breaking and interior attracting cycle', ha='center', fontsize=12)
    fig.legend(handles=[Patch(color=COLORS['P5'], label='P5: five peaks'),
        Patch(color=COLORS['P7'], label='P7: seven peaks, 8 crossings'),
        Patch(color=COLORS['new'], label='Same P7 family, 12 crossings'),
        Line2D([0], [0], color=COLORS['contact'], lw=2.7, label='P7 contact curve')],
        loc='lower center', bbox_to_anchor=(.5, .105), ncol=2, frameon=False, columnspacing=3)
    fig.text(.5, .035, 'Surfaces connect checked stability-interval midpoints at their actual coordinates; stability widths are shown separately.\n'
        'All P7 nodes are recovered. Eleven P5 nodes remain unresolved. Facets are interpolation, not an existence proof.',
        ha='center', fontsize=10.5, color='#536573')
    save(fig, 'chua-slope-refined-volumes')


def sections():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.4))
    fig.subplots_adjust(left=.085, right=.93, bottom=.15, top=.85, hspace=.48, wspace=.30)
    for family in ['P5', 'P7']:
        d = read(f'{family}-alpha.json')
        rows = sorted([r for r in d['nodes'] if r['index'][0] == 20], key=lambda r: r['core']['m1'])
        x = [r['core']['m1'] for r in rows]
        lower = [r['symmetry']['m0'] for r in rows]
        upper = [r['fold']['m0'] for r in rows]
        axes[0, 0].fill_between(x, lower, upper, color=COLORS[family], alpha=.7, label=family)
        axes[0, 0].plot(x, lower, color=COLORS[family], lw=.65)
        axes[0, 0].plot(x, upper, color=COLORS[family], lw=.65)
        axes[0, 1].semilogy(x, [r['width_m0'] for r in rows], color=COLORS[family], lw=2.3, label=family)
    axes[0, 0].set(xlabel=r'$m_1$', ylabel=r'$m_0$', title=r'Actual stability bands: $\alpha=9$, $\beta=100/7$')
    axes[0, 0].legend(frameon=False)
    axes[0, 1].axhline(.01, color='#75818a', ls='--', label='Original volume step in m0')
    axes[0, 1].set(xlabel=r'$m_1$', ylabel='Stability width in m0', title='Widths require continuation, not grid occupancy')
    axes[0, 1].legend(frameon=True, facecolor='white', framealpha=1, edgecolor='none', fontsize=10, loc='upper right')
    grazing = read('grazing-results.json')['contact_curves']
    for col, kind in enumerate(['alpha', 'beta']):
        ax = axes[1, col]
        d = read(f'P7-{kind}.json')
        widths = np.full((41, 41), np.nan)
        for r in d['nodes']:
            widths[tuple(r['index'])] = np.log10(r['width_m0'])
        im = ax.pcolormesh(d['m1'], d['controls'], widths, shading='nearest', cmap='viridis', vmin=-4.15, vmax=-2.85)
        for target, style, label in [(1., '-', 'Fold contact'), (-.43, '--', 'Interior-cycle contact'), (-1., ':', 'Symmetry-edge contact')]:
            c = next(x for x in grazing if x['kind'] == kind and x['target'] == target)
            ax.plot([r['contact']['orbit']['m1'] for r in c['nodes']], [r['control'] for r in c['nodes']],
                    style, color='white', lw=1.7, label=label)
        ax.set(xlabel=r'$m_1$', ylabel='α' if kind == 'alpha' else 'β', title='P7 width and three distinct contact curves')
        ax.legend(frameon=True, facecolor='#183447', edgecolor='none', labelcolor='white', fontsize=9, loc='lower left')
        bar = fig.colorbar(im, ax=ax, pad=.025, fraction=.045)
        bar.set_label(r'$\log_{10}$(width in $m_0$)', fontsize=10)
    fig.suptitle('Layer width and changes of diode itinerary', y=.965, fontsize=22, weight='semibold')
    fig.text(.5, .91, 'The grazing curve crosses the interior of a stability band; it does not mark the birth of a new family.', ha='center', fontsize=11.5)
    fig.text(.5, .035, 'Boundaries are corrected periodic roots. At 32 additional cell centres, finer interpolation has a median 2.97-fold error reduction.\n'
             'The largest remaining interpolation error is 20.4% of the local width; use corrected roots for narrow-window searches.',
             ha='center', fontsize=10.5, color='#536573')
    save(fig, 'chua-slope-refined-sections')


def sample_orbit(r, times):
    f = flights(parameters(r), r)
    mats = f.matrices(r['beta'])
    durations = np.array(r['q'][2:])
    limits = np.r_[0., np.cumsum(durations)]
    starts = [np.r_[0., r['q'][:2]]]
    for duration, region in zip(durations, r['regions']):
        A, c, eig, V, W = mats[region]
        starts.append(c+((V*np.exp(eig*duration))@W).real@(starts[-1]-c))
    result = np.empty((len(times), 3))
    for i, time in enumerate(times):
        time = time % r['period']
        sign = 1 if time < r['period']/2 else -1
        time = time % (r['period']/2)
        k = min(np.searchsorted(limits, time, side='right')-1, len(durations)-1)
        A, c, eig, V, W = mats[r['regions'][k]]
        result[i] = sign*(c+((V*np.exp(eig*(time-limits[k])))@W).real@(starts[k]-c))
    return result


def grazing():
    d = read('grazing-results.json')
    prediction = read('grazing-linear-prediction.json')
    before = d['witnesses'][0]['orbit']
    after = d['witnesses'][2]['orbit']
    contact = prediction['contact']
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.subplots_adjust(left=.08, right=.95, top=.84, bottom=.15, hspace=.43, wspace=.28)
    for r, color, label in [(before, COLORS['P7'], 'Before contact'), (after, COLORS['new'], 'After contact')]:
        times = np.linspace(0, r['period'], 3500, endpoint=False)
        state = sample_orbit(r, times)
        axes[0, 0].plot(state[:, 0], state[:, 2], color=color, lw=1.2, alpha=.9, label=label)
        centre = r['q'][2]+contact['contact_time_within_outer_flight']
        local_times = centre+np.linspace(-.045, .045, 1000)
        local = sample_orbit(r, local_times)
        offset = local_times[np.argmin(local[:, 0])]
        axes[0, 1].plot(local_times-offset, (local[:, 0]-1)*1e3, color=color, lw=2, label=label)
    axes[0, 0].set(xlabel='x', ylabel='z', title='The same seven-peak motion persists')
    axes[0, 0].legend(frameon=False, fontsize=10)
    axes[0, 1].axhline(0, color='#79848c', ls='--', lw=1)
    axes[0, 1].axhspan(-.6, 0, color='#e5c997', alpha=.3)
    axes[0, 1].set(xlabel='Time relative to the local minimum', ylabel=r'$10^3(x-1)$',
                   title='A short central flight appears at the contact', xlim=(-.045, .045), ylim=(-.6, 4.5))
    positive = sorted([w for w in d['witnesses'] if w['delta'] > 0], key=lambda w: w['delta'])
    delta = np.array([w['delta'] for w in positive])
    actual = np.array([w['new_central_flight'] for w in positive])
    predicted = prediction['coefficient']*np.sqrt(delta)
    axes[1, 0].loglog(delta, predicted, '--', color='#7d8b92', lw=2, label=r'Taylor prediction: $0.957956\sqrt{\delta}$')
    axes[1, 0].loglog(delta, actual, 'o', color=COLORS['new'], ms=6, label='Corrected periodic orbit')
    axes[1, 0].set(xlabel=r'$\delta=m_1-m_{1,*}$', ylabel='Duration of the new central flight', title='Predicted and measured square-root scaling')
    axes[1, 0].legend(frameon=False, fontsize=10)
    axes[1, 0].grid(alpha=.15, which='both')
    axes[1, 1].semilogx(delta, 100*(actual/predicted-1), 'o-', color=COLORS['new'], lw=1.6)
    axes[1, 1].set(xlabel=r'$\delta=m_1-m_{1,*}$', ylabel='Relative discrepancy (%)',
                   title='Discrepancy falls to 0.0694% near contact')
    axes[1, 1].axhline(0, color='#75818a', ls='--', lw=.8)
    axes[1, 1].grid(alpha=.15)
    fig.suptitle('New crossings do not imply a new attractor', fontsize=23, y=.965, weight='semibold')
    fig.text(.5, .905, 'P7: 8 → 12 breakpoint crossings; 7 peaks throughout; dominant multiplier magnitude 0.181–0.1856',
             ha='center', fontsize=12)
    fig.text(.5, .035, 'A fixed smooth parameter line crosses the contact at m1 = −0.655374411697. DOP853 and Radau check the actual field.\n'
        'The analytic comparison follows preliminary exploration; it is not a blind forecast. The unsuccessful first protocol is retained.',
        ha='center', fontsize=10.5, color='#536573')
    save(fig, 'chua-slope-grazing-law')


if __name__ == '__main__':
    volumes()
    sections()
    grazing()
