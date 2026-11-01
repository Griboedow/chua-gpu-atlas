"""Render family-labelled continuation samples, narrow stability intervals, and overlap evidence."""
import argparse
import colorsys
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import plotly.graph_objects as go

ROOT = Path(__file__).resolve().parents[1]


def save(fig, name, assets, out):
    fig.savefig(assets / (name + '.png'), dpi=190, facecolor='white')
    fig.savefig(out / (name + '.pdf'), facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    args = parser.parse_args(); out = args.output; assets = args.assets
    assets.mkdir(parents=True, exist_ok=True)
    families = sorted(json.loads((out / 'family-registry.json').read_text())['families'], key=lambda r: r['id'])
    sheets = {r['id']: json.loads((out / (r['id'] + '-sheet.json')).read_text()) for r in families}
    colors = {r['id']: colorsys.hsv_to_rgb((k * .61803398875) % 1, .73, .72) for k, r in enumerate(families)}
    beta = np.full((len(families), 9, 9), np.nan); periods = beta.copy(); widths = beta.copy()
    for k, family in enumerate(families):
        for n in sheets[family['id']]['nodes']:
            i, j = n['index']; beta[k, i, j] = n['core']['beta']; periods[k, i, j] = n['core']['period']
            if 'fold' in n:
                widths[k, i, j] = n['fold']['beta'] - n['symmetry']['beta']
    aa = np.array(sheets[families[0]['id']]['alpha']); mm = np.array(sheets[families[0]['id']]['m0']); offsets = (mm + 8/7) * 1000
    inside = (beta >= 14.24) & (beta <= 14.4)
    count = inside.sum(axis=0); complete = np.isfinite(beta).all(axis=0)
    np.savez_compressed(out / 'display-samples.npz', family=np.array([r['id'] for r in families]), alpha=aa, m0=mm,
        beta=beta, period=periods, symmetric_core_width=widths, representative_count=count, all_families_sampled=complete)

    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 7.2), gridspec_kw={'width_ratios': [1.65, 1]}, constrained_layout=True)
    old = np.load(ROOT / 'results/chua-layer-census/atlas-comparison.npz'); bi, mi = np.where(old['regular'][20])
    axes[0].scatter((old['m0'][mi]+8/7)*1000, old['beta'][bi], c='#d2d7dc', s=14, marker='s', label='Old regular grid')
    for k, family in enumerate(families):
        id = family['id']; line = beta[k, 4]; valid = np.isfinite(line)
        axes[0].plot(offsets, line, '.-', ms=4, lw=.85, color=colors[id], label=id)
    axes[0].set(xlim=(-6.1, 6.1), ylim=(14.24, 14.4), xlabel=r'$10^3(m_0+8/7)$', ylabel=r'$\beta$', title=r'Labelled representative branches at $\alpha=9$')
    axes[0].grid(alpha=.12)
    im = axes[1].imshow(count, origin='lower', extent=[-6.75, 6.75, 8.9325, 9.0675], aspect='auto', cmap='viridis', interpolation='nearest')
    for i in range(9):
        for j in range(9):
            axes[1].text(offsets[j], aa[i], str(count[i,j]) + ('*' if not complete[i,j] else ''), ha='center', va='center', fontsize=8,
                color='white' if count[i,j] < (count.min()+count.max())/2 else '#182124')
    axes[1].set(xlabel=r'$10^3(m_0+8/7)$', ylabel=r'$\alpha$', title='Catalogued representatives inside the beta box')
    fig.colorbar(im, ax=axes[1], fraction=.045, pad=.03, label='Count at each sampled (alpha, m0)')
    fig.suptitle('The refined atlas: 32 catalogued window families', fontsize=17)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=11,fontsize=8,frameon=False)
    axes[1].set_title('Representatives inside the beta box\n* = some family samples unresolved',fontsize=11)
    save(fig, 'chua-refined-atlas', assets, out)

    interactive = go.Figure()
    for k, family in enumerate(families):
        id = family['id']; rgb = 'rgb(%d,%d,%d)' % tuple(round(c*255) for c in colors[id]); x=[]; y=[]; z=[]; custom=[]
        # Separate row/column polylines leave holes at unresolved mesh nodes.
        for axis in [0, 1]:
            for j in range(9):
                for i in range(9):
                    u,v = (i,j) if axis == 0 else (j,i)
                    if inside[k,u,v]:
                        x.append(aa[u]); y.append(offsets[v]); z.append(beta[k,u,v]); custom.append([periods[k,u,v], widths[k,u,v]])
                    else:
                        x.append(None); y.append(None); z.append(None); custom.append([None,None])
                x.append(None); y.append(None); z.append(None); custom.append([None,None])
        interactive.add_trace(go.Scatter3d(x=x, y=y, z=z, customdata=custom, mode='lines+markers', name=f'{id} | parent: {family["peak_count"]} maxima',
            line=dict(color=rgb,width=3), marker=dict(color=rgb,size=2.4), connectgaps=False,
            hovertemplate=f'<b>{id}</b><br>alpha=%{{x:.5f}}<br>m0 offset x1000=%{{y:.3f}}<br>beta=%{{z:.11f}}<br>T=%{{customdata[0]:.6f}}<br>symmetric-cycle stability width=%{{customdata[1]:.3e}}<extra></extra>'))
    interactive.update_layout(height=850,
        margin=dict(l=20,r=20,t=55,b=85), scene=dict(xaxis_title='alpha',yaxis_title='1000 (m0 + 8/7)',zaxis_title='beta',
        xaxis=dict(range=[8.94,9.06]),yaxis=dict(range=[-6,6]),zaxis=dict(range=[14.24,14.4]),aspectmode='cube'),
        legend=dict(itemsizing='constant',title='Click to hide; double-click to isolate'),
        updatemenus=[dict(type='buttons',direction='right',x=0,y=1.07,xanchor='left',yanchor='top',buttons=[dict(label=label,method='update',args=[{'visible': mask}]) for label,mask in [
            ('All', [True]*len(families)),('Off-center discoveries',[r['id'].startswith('E') for r in families]),
            ('Adding sequence',[r['id'] in ['S04','S05','S06','S07','P55','P65'] for r in families])]])])
    fragment = interactive.to_html(full_html=False, include_plotlyjs=True, config={'displaylogo':False})
    html = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>Refined Chua layer atlas</title>
<style>body{font:16px system-ui;margin:24px;color:#172333}p{max-width:1100px;line-height:1.5}h1{margin-bottom:8px}</style>
<h1>32 catalogued window families</h1><p>Rotate, zoom, and select families in the legend. Each dot is a numerically closed attracting cycle.
Lines show mesh adjacency, not verified attraction at every interpolated point. Their visual width is unrelated to physical layer thickness.
The atlas is clipped to alpha 8.94–9.06, beta 14.24–14.40 and m0 = −8/7 ± 0.006, with m1 = −5/7.</p>
<p>This is a finite family catalogue, not an exhaustive count of windows or a proof of global connectedness.
Missing samples remain unresolved. Symmetry partners and period-doubled descendants are grouped with their traced parent family.
Hover shows the full period and, where continued, the stability width of the symmetric cycle.</p>''' + fragment + '</html>\n'
    # The bundled Plotly shader sources contain incidental trailing spaces.
    html='\n'.join(line.rstrip() for line in html.splitlines())+'\n'
    (out/'atlas.html').write_text(html,encoding='utf-8',newline='\n')

    oldadding = json.loads((ROOT/'results/chua-layer-census/census-analysis.json').read_text())['adding_sequence']
    added = json.loads((out/'extended-adding.json').read_text())
    peaks = [r['maxima'] for r in oldadding] + [r['peak_count'] for r in added]
    corewidth = [r['symmetric_core_width'] for r in oldadding] + [r['symmetric_core_width'] for r in added]
    duration = [r['period'] for r in oldadding] + [r['core']['period'] for r in added]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    axes[0].semilogy(peaks, corewidth, 'o-', color='#9653a2', lw=2)
    axes[0].axhline(.004,color='#b26030',ls='--',label='Old beta spacing: 0.004')
    for x,y in zip(peaks,corewidth): axes[0].annotate(f'{y:.2e}',(x,y),xytext=(0,-17),textcoords='offset points',ha='center',fontsize=9)
    axes[0].set(xlabel='Maxima per full cycle',ylabel='Fold-to-symmetry width in beta',xlim=(10,70),ylim=(1e-8,1e-2),title='Six shrinking stability intervals');axes[0].legend(fontsize=9);axes[0].grid(alpha=.15)
    axes[1].plot(peaks,duration,'o-',color='#197b89',lw=2)
    axes[1].set(xlabel='Maxima per full cycle',ylabel='Full period',title='Approximately 21.6 extra time units per member');axes[1].grid(alpha=.15)
    fig.suptitle('15 → 25 → 35 → 45 → 55 → 65',fontsize=17)
    save(fig,'chua-refined-adding',assets,out)

    orbit = np.load(out/'four-cycles.npz'); baseline = np.load(out/'overlap-basins.npz'); fine = np.load(out/'basin-capture-full.npz')
    palette=['#177e89','#82b5ba','#b04c7b','#e5a5bd']; names=['F1a','F1b','F2a','F2b']
    fig, axes = plt.subplots(2,2,figsize=(11.5,8.7),constrained_layout=True)
    for pair in range(2):
        ax=axes[0,pair]
        for j in range(2*pair,2*pair+2):
            xyz=orbit[names[j]+'_xyz'];ax.plot(xyz[:,0],xyz[:,2],lw=1,color=palette[j],label=names[j])
        ax.set(xlabel='x',ylabel='z',title=f'Family {pair+1}: two inversion partners');ax.legend();ax.grid(alpha=.1)
    ax=axes[1,0];ax.imshow(fine['labels'].reshape(201,201),origin='lower',extent=[.19,.30,-.54,-.30],aspect='auto',cmap=ListedColormap(palette),vmin=.5,vmax=4.5,interpolation='nearest')
    ax.set(xlabel='Initial y (x = 0)',ylabel='Initial z',title='Computed outcomes at step 0.000625')
    x=np.arange(4); coarsecounts=np.bincount(baseline['labels'],minlength=5)[1:]/40401*100; finecounts=np.bincount(fine['labels'],minlength=5)[1:]/40401*100
    ax=axes[1,1];ax.bar(x-.18,coarsecounts,.36,color=palette,alpha=.5,label='Step 0.00125');ax.bar(x+.18,finecounts,.36,color=palette,label='Step 0.000625')
    ax.set(xticks=x,xticklabels=names,ylabel='Fraction of the sampled rectangle (%)',ylim=(0,38),title='Similar fractions; pointwise labels disagree')
    ax.legend(fontsize=9)
    ax.set_title(f'Similar fractions; only {np.mean(baseline["labels"]==fine["labels"])*100:.2f}% pointwise agreement',fontsize=10)
    fig.suptitle('Four attracting cycles at the same parameters',fontsize=17)
    fig.supxlabel('alpha = 9.06, beta = 14.3458, m0 = −8/7 + 0.006.  The computed basin image is not pointwise converged.',fontsize=10)
    save(fig,'chua-four-cycle-overlap',assets,out)

    behavior=json.loads((out/'family-behavior.json').read_text()); features={r['id']:r for r in behavior}
    with (out/'family-catalog.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(behavior[0]),lineterminator='\n');writer.writeheader()
        writer.writerows(sorted(behavior,key=lambda r:r['reference_beta']))
    report=dict(families=len(families),samples=int(np.isfinite(beta).sum()),inside_samples=int(inside.sum()),
        representative_count_range=[int(count.min()),int(count.max())],fully_sampled_parameter_nodes=int(complete.sum()),
        missing_samples={r['id']:int(np.isnan(beta[k]).sum()) for k,r in enumerate(families) if np.isnan(beta[k]).any()},
        adding=dict(maxima=peaks,widths=corewidth,periods=duration,old_spacing_over_thinnest=float(.004/corewidth[-1])),
        basin_fractions=dict(baseline=coarsecounts.tolist(),half_step=finecounts.tolist(),max_change_percentage_points=float(max(abs(coarsecounts-finecounts)))),
        single_outer_region_families=[r['id'] for r in behavior if not r['both_outer_regions']])
    (out/'plot-summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
