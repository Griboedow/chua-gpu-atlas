"""Shared-coordinate portraits, exact-polynomial evidence and offline viewer."""
from search_chua_links import ROOT,OUT,read
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plotly.offline import get_plotlyjs
from plot_chua_linked_attractors import scene,save,BG,INK,MUTED,COLORS


def load():
    certificates=read(OUT/'certificates-d160.json')['witnesses'];datasets={}
    for cert in certificates:
        name=cert['name']
        if name=='S02-S08':
            paths=[OUT/'overlaps'/name/f'{id}-d90-n384.json' for id in ['S02','S08']]
            roots=[read(p) for p in paths];curves=[np.load(p.with_suffix('.npz'))['xyz'] for p in paths]
            poly=read(OUT/'overlaps'/name/'jones-d90-n384-s0.json')
            names=['S02','S08'];title='Two different attracting families'
        else:
            family=name[:3];p=OUT/'daughters'/f'{family}-a1e-03-d90-n384.json';root=read(p)
            roots=[root,root];c=np.load(p.with_suffix('.npz'))['xyz'];curves=[c,-c]
            poly=read(p.with_name(p.stem+'-jones-s0.json'))
            names=[family+' daughter', 'Inversion partner'];title=family+': a pair after symmetry breaking'
        datasets[name]=dict(names=names,curves=curves,title=title,parameters=cert['parameters'],
            periods=[float(r['period']) for r in roots],radii=[r['spectral_radius'] for r in roots],
            section=[np.array([[float(v['midpoint']) for v in point] for point in r['upward_section_points']]) for r in cert['records']],
            polynomial=poly,note='Non-split numerical link. Gauss linking number: 0. The joint Jones polynomial differs from the split union of these two components. Existence, distinctness and attraction of both cycles are interval-certified; the link type is not.')
    negative=read(OUT/'symmetry/A04-d60-n96-s0.json');path=ROOT/negative['source'];c=np.load(path)['xyz'];root=read(path.with_suffix('.json'))
    datasets['A04-control']=dict(names=['A04','Inversion partner'],curves=[c,-c],title='A04: two knots in a split union',
        parameters={k:root[k] for k in ['alpha','beta','m0','m1','gamma']},periods=[float(root['period'])]*2,
        radii=[float(max(abs(np.roots([1,-float(root['trace']),float(root['determinant'])]))))]*2,section=[],
        note='Negative control. These components are individually knotted, but the joint numerical diagram separates into two one-component pieces. Individual knotting does not imply mutual linking.')
    return datasets


def overview(data):
    fig=plt.figure(figsize=(14,9.5));gs=fig.add_gridspec(2,2,left=.025,right=.975,top=.83,bottom=.18,height_ratios=[3.8,1.35],wspace=.15,hspace=.08)
    fig.text(.045,.956,'Coexisting attracting motions can form a non-split link',fontsize=22,weight='bold',color=INK)
    fig.text(.045,.907,'Three numerical witnesses. Two discovery routes: overlapping families and symmetry breaking.',fontsize=12,color=MUTED)
    for col,name in enumerate(['S02-S08','S01-daughters']):
        d=data[name];ax=fig.add_subplot(gs[0,col],projection='3d');scene(ax,d['curves'],d['names'])
        ax.set_title('S02 + S08' if col==0 else 'S01 daughter + inversion partner',fontsize=14,pad=10,weight='bold')
    ax=fig.add_subplot(gs[1,0]);ax.set_axis_off()
    ax.text(.09,.84,'Different motions, one set of parameters',fontsize=13,weight='bold')
    ax.text(.09,.54,'Periods: 28.22965 and 10.84575\nBoth attracting; both individually inversion-symmetric.',fontsize=11,color=MUTED,linespacing=1.6)
    ax.text(.09,.06,'Gauss linking number = 0\nExact joint Jones polynomial ≠ split-union polynomial',fontsize=11,color=INK,linespacing=1.6)
    ax=fig.add_subplot(gs[1,1]);pts=[x[0] for x in data['S01-daughters']['section']];center=np.mean(pts,axis=0)
    for i,p in enumerate(pts):
        ax.scatter((p[0]-center[0])*1e6,(p[1]-center[1])*1e6,s=85,color=COLORS[i],zorder=3)
        ax.annotate('Daughter' if i==0 else 'Partner',((p[0]-center[0])*1e6,(p[1]-center[1])*1e6),xytext=(8,-10 if i==0 else 4),textcoords='offset points',color=COLORS[i],fontsize=10)
    ax.set(xlim=(-4,8),ylim=(-4.5,4.5),xlabel='(y − common center) × 10⁶',ylabel='(z − common center) × 10⁶')
    ax.grid(alpha=.2);ax.set_title('Nearly coincident at full scale; distinct on x = 0',fontsize=11)
    fig.text(.045,.085,'Every portrait preserves the shared state-space coordinates. The section close-up uses one common translation and scale.\nLeft: α = 8.94, β = 14.26071457889439. Right: α = 9, β = 14.25984090319409. Full slopes and certificates are in the report.',fontsize=10,color=MUTED,linespacing=1.5)
    fig.text(.045,.025,'Exact interval checks certify the attracting cycles. Linking is supported by refined numerical curves and exact diagram invariants;\nan ambient-isotopy certificate connecting each polygon to its exact ODE orbit remains open.',fontsize=10,color=MUTED,linespacing=1.4)
    save(fig,'chua-link-search-overview')


def evidence(data):
    fig,axs=plt.subplots(3,1,figsize=(14,9));fig.subplots_adjust(left=.09,right=.97,top=.81,bottom=.21,hspace=.65)
    fig.text(.06,.948,'Zero linking number does not mean a split pair',fontsize=24,weight='bold',color=INK)
    fig.text(.06,.892,'For a split union: V(K₁ ⊔ K₂) = (−A² − A⁻²) V(K₁) V(K₂).\nThe exact coefficients below violate this identity in every witness.',fontsize=12,color=MUTED,linespacing=1.6)
    for ax,(name,d) in zip(axs,list(data.items())[:3]):
        pol=d['polynomial'];actual={int(k):v for k,v in pol['joint_polynomial'].items()};expected={int(k):v for k,v in pol['split_union_polynomial'].items()}
        powers=np.array(sorted(set(actual)|set(expected)));ax.axhline(0,color=MUTED,lw=.6)
        ax.bar(powers-.6,[actual.get(int(p),0) for p in powers],width=1.2,color=COLORS[0],label='Measured joint diagram')
        ax.bar(powers+.6,[expected.get(int(p),0) for p in powers],width=1.2,color=COLORS[1],alpha=.8,label='Hypothetical split union')
        limit=max(abs(powers));expected_limit=max(abs(p) for p in expected)
        ax.axvspan(-limit-3,-expected_limit-2,color=COLORS[0],alpha=.06);ax.axvspan(expected_limit+2,limit+3,color=COLORS[0],alpha=.06)
        ax.set(xlim=(-limit-3,limit+3),ylabel='Coefficient',xlabel='Power of A')
        ax.set_title(name.replace('-daughters',' daughters').replace('S02-S08','S02 + S08')+f'     |     joint degrees ±{limit}; split degrees ±{expected_limit}',fontsize=12,loc='left',weight='bold')
        ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    axs[0].legend(frameon=False,fontsize=10,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.43))
    fig.text(.06,.085,'Integer coefficients are computed exactly from planar diagrams. The additional outer coefficients alone exclude a split union.\nAgreement across precision, sampling density and projections is recorded alongside the complete polynomials.',fontsize=11,color=INK,linespacing=1.5)
    fig.text(.06,.027,'Each comparison uses the actual knot polynomial of each component. This avoids mistaking individually knotted but separated curves\nfor a mutually linked pair. Equal polynomials would be inconclusive; here they are unequal.',fontsize=10,color=MUTED,linespacing=1.4)
    save(fig,'chua-link-search-evidence')


def viewer(data):
    serial={}
    for name,d in data.items():
        serial[name]={k:v for k,v in d.items() if k not in ['curves','section','polynomial']}
        serial[name]['curves']=[np.round(c,13).tolist() for c in d['curves']]
        serial[name]['section']=[p.tolist() for p in d['section']]
    template=(ROOT/'research/chua_link_search_viewer.html').read_text(encoding='utf-8')
    html=template.replace('__PLOTLY__',get_plotlyjs()).replace('__DATA__',json.dumps(serial,separators=(',',':')))
    (OUT/'linked-pairs.html').write_text('\n'.join(line.rstrip() for line in html.splitlines())+'\n',encoding='utf-8',newline='\n')


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.facecolor':BG,'axes.facecolor':BG,
        'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':INK,'text.color':INK,'xtick.color':MUTED,'ytick.color':MUTED})
    data=load();overview(data);evidence(data);viewer(data);print('Wrote portraits, coefficient figure and offline viewer',flush=True)


if __name__=='__main__':main()
