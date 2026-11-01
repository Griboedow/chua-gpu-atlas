"""Plot the layer census, direct atlas comparison, and coexistence controls."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT=Path(__file__).resolve().parents[1]
COLORS=['#086b9a','#d66b18']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/chua-layer-census')
    p.add_argument('--assets',type=Path,default=ROOT/'docs/assets')
    args=p.parse_args();out=args.output;args.assets.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
    a=json.loads((out/'census-analysis.json').read_text())
    first=json.loads((ROOT/'results/chua-layers/extended.json').read_text())['surfaces']
    second=json.loads((out/'second-surfaces.json').read_text())
    atlas=np.load(out/'atlas-comparison.npz')
    def save(fig,name):
        fig.savefig(args.assets/(name+'.png'),dpi=180,bbox_inches='tight')
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained',sharey=True)
    xx=np.linspace(8.94,9.06,301);aa=np.linspace(8.94,9.06,7);mm=np.linspace(-8/7-.006,-8/7+.006,7)
    for ax,k in zip(axes,[20,33]):
        m0=atlas['m0'][k];x,y=np.meshgrid(atlas['alpha'],atlas['beta'],indexing='ij')
        mask=atlas['regular'][:,:,k]
        ax.scatter(x[mask],y[mask],s=15,c='#b44f92',alpha=.7,label='Original regular nodes',zorder=3)
        for j,rows in enumerate([first,second]):
            vals={}
            for key in ['fold','symmetry','period_doubling']:
                field=np.array([r[key]['beta'] for r in rows]).reshape(7,7)
                vals[key]=RegularGridInterpolator((aa,mm),field)(np.column_stack([xx,np.full(len(xx),m0)]))
            ax.fill_between(xx,vals['period_doubling'],vals['fold'],color=COLORS[j],alpha=.15)
            ax.plot(xx,vals['fold'],color=COLORS[j],lw=1.8,label=f'Family {j+1}: fold / first doubling')
            ax.plot(xx,vals['period_doubling'],color=COLORS[j],lw=1.8)
            ax.plot(xx,vals['symmetry'],color=COLORS[j],lw=1,ls=':')
        ax.set(xlim=(8.94,9.06),ylim=(14.24,14.4),xlabel=r'$\alpha$',title=f'$m_0$ = {m0:.7f}')
        ax.ticklabel_format(useOffset=False);ax.grid(alpha=.15)
    axes[0].set_ylabel(r'$\beta$');axes[0].legend(loc='upper left',fontsize=8)
    axes[1].text(.03,.95,'Overlap permits two stable cycles',transform=axes[1].transAxes,va='top',fontsize=9)
    fig.suptitle('Two continued five-maximum families inside the original atlas',fontsize=14)
    fig.supxlabel('Shading: fold to first doubling; dotted: symmetry breaking. Boundaries interpolated from 7 × 7 points.',fontsize=9)
    save(fig,'chua-layer-atlas-comparison')

    fig=plt.figure(figsize=(12,8),layout='constrained');gs=fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[0,:]);data=np.load(out/'tight.npz');mask=data['period']>0
    ax.scatter(data['beta'][mask],data['period'][mask],s=3,c='#b2b2b2',label='GPU peak recurrence candidates')
    anchors=a['anchors'];ax.scatter([r['beta'] for r in anchors],[r['maxima'] for r in anchors],s=35,c='#086b9a',zorder=4,label='17 independently checked window representatives')
    for r in anchors:
        if r['maxima'] in [9,11,45,60] or r['beta']>14.319:
            offset=(-8,7) if abs(r['beta']-14.34990625)<1e-9 else ((10,7) if abs(r['beta']-14.3512265625)<1e-9 else (0,7))
            ax.annotate(str(r['maxima']),(r['beta'],r['maxima']),xytext=offset,textcoords='offset points',ha='center',fontsize=8)
    ax.set(xlim=(14.24,14.4),ylim=(0,100),xlabel=r'$\beta$',ylabel='Maxima per full cycle',title=r'Fixed transect: $\alpha=9$, $m_0=-8/7$, $m_1=-5/7$')
    ax.legend(loc='upper left',fontsize=8);ax.ticklabel_format(useOffset=False)
    ax=fig.add_subplot(gs[1,0]);seq=a['adding_sequence'];n=np.array([r['maxima'] for r in seq])
    ax.semilogy(n,[r['symmetric_stability_width'] for r in seq],'o-',color=COLORS[0],label='Fold to symmetry breaking')
    ax.set(xlabel='Maxima per full cycle',ylabel=r'Symmetric-cycle stability width in $\beta$',title='Four resolved members of a shrinking sequence',xticks=n)
    ax.grid(alpha=.2);ax.legend(fontsize=8)
    for r in seq:ax.annotate(f"{r['symmetric_stability_width']:.2e}",(r['maxima'],r['symmetric_stability_width']),xytext=(5,5),textcoords='offset points',fontsize=8)
    ax=fig.add_subplot(gs[1,1])
    for r in seq:
        dist=np.array(r['distances_to_parent_section']);i=int(dist.argmin());dist=np.roll(dist,-i)
        ax.semilogy(np.arange(len(dist)),dist,'o-',ms=4,label=f"{r['maxima']} maxima")
    ax.set(xlabel='Successive upward section returns (starting at closest)',ylabel='Distance to unstable five-maximum cycle',title='Longer cycles pass closer to the unstable parent')
    ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('A finer census reveals more windows; it does not determine a global layer count',fontsize=14)
    save(fig,'chua-layer-census')

    c=json.loads((out/'coexistence.json').read_text());d=np.load(out/'coexistence.npz')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for j,r in enumerate(c['cycles']):
        xyz=d[f'family{j+1}_xyz'];t=d[f'family{j+1}_t'];T=r['cycle']['period']
        label=f'Family {j+1}: T = {T:.6f}'
        axes[0].plot(xyz[:,0],xyz[:,2],color=COLORS[j],lw=1,label=label)
        axes[1].plot(t,xyz[:,0],color=COLORS[j],lw=1.5,label=label)
    axes[0].set(xlabel='x',ylabel='z',title='Two different attracting orbits')
    axes[1].set(xlabel='Time',ylabel='x',title='Both have five maxima per full cycle')
    for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.15)
    fig.suptitle(r'Coexistence at $\alpha=9.04$, $\beta=14.352$, $m_0=-8/7+0.004$, $m_1=-5/7$',fontsize=13)
    save(fig,'chua-layer-coexistence')


if __name__=='__main__':main()
