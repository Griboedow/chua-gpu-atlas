"""Rebuild publication figures from the layer continuation and local surfaces."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from matplotlib.patches import Patch

from trace_chua_window import Chua
from chua_layer_bifurcations import Flights


def partner_y(model,r,regions):
    state=np.array([0.,*r['q'][:2]])
    for region,duration in zip(regions,r['q'][2:]):
        A=model.jacobian(r['beta'],region==0)
        B=np.zeros((4,4));B[:3,:3]=A;B[0,3]=-model.alpha*(model.m0-model.m1)*region
        u=np.r_[state,1.];end=(expm(B*duration)@u)[:3]
        if region==0 and state[0]>.9 and end[0]<-.9:
            t=brentq(lambda t:(expm(B*t)@u)[0],0,duration)
            return float(-(expm(B*t)@u)[1])
        state=end
    raise RuntimeError('No negative section found')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    root=Path(__file__).resolve().parents[1]
    p.add_argument('--output',type=Path,default=root/'results/chua-layers')
    p.add_argument('--assets',type=Path,default=root/'docs/assets')
    args=p.parse_args();out=args.output;assets=args.assets;assets.mkdir(parents=True,exist_ok=True)
    a=json.loads((out/'affine-central.json').read_text());b=json.loads((out/'extended.json').read_text())
    c=json.loads((out/'cascade.json').read_text());audit=json.loads((out/'audit.json').read_text())
    fold=a['boundaries']['1'];sym=a['boundaries']['-1'];pd=b['period_doubling'];pd2=c['second_period_doubling']
    model=Chua();family=Flights(model,a['regions'],a['targets'],True)
    saddle=max(b['fold_partners'],key=lambda r:r['spectral_radius']);q=np.array(saddle['q']);saddles=[]
    for beta in np.linspace(saddle['beta'],14.308,100):
        q=family.solve(float(beta),q);saddles.append(family.record(float(beta),q))
    central=sorted(a['central']+[fold],key=lambda r:r['beta'])
    asym=sorted(b['asymmetric'],key=lambda r:r['beta'])
    partners=[partner_y(model,r,b['full_regions']) for r in asym]
    (out/'plot-data.json').write_text(json.dumps(dict(saddle_branch=saddles,asymmetric_partner_section_y=partners,
                                                   asymmetric_partner_beta=[r['beta'] for r in asym]),indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,layout='constrained')
    blue='#2764b0';green='#14826e';orange='#bd5e12';gray='#7b8090'
    for ax in axs:
        for l,r,col in [(sym['beta'],fold['beta'],blue),(pd['beta'],sym['beta'],green),(pd2['beta'],pd['beta'],orange)]:
            ax.axvspan(l,r,color=col,alpha=.07)
        for r in [fold,sym,pd,pd2]:ax.axvline(r['beta'],color='#a1a6af',lw=.8,ls=':')
        ax.grid(alpha=.15)
        ax.ticklabel_format(useOffset=False)
    def branch(ax,rows,color,label,other=None):
        x=np.array([r['beta'] for r in rows]);y=np.array(other if other is not None else [r['q'][0] for r in rows]);stable=np.array([r['spectral_radius']<1+1e-8 for r in rows])
        ax.plot(x,np.where(stable,y,np.nan),color=color,lw=2,label=label)
        ax.plot(x,np.where(~stable,y,np.nan),color=color,lw=1.3,ls='--')
    branch(axs[0],central,blue,'Symmetric cycle')
    branch(axs[0],asym,green,'Symmetry-related pair')
    branch(axs[0],asym,green,None,partners)
    axs[0].plot([r['beta'] for r in saddles],[r['q'][0] for r in saddles],color=gray,ls='--',label='Saddle partner')
    axs[0].scatter([fold['beta'],sym['beta']],[fold['q'][0],sym['q'][0]],color='#202630',zorder=4,s=20)
    axs[0].set_ylabel('y at the upward x = 0 section')
    axs[0].set_title('One thin layer contains several periodic regimes\nSolid: attracting cycle; dashed: unstable continuation',loc='left',fontweight='bold')
    axs[0].legend(loc='upper left',frameon=False)
    for rows,color,label in [(central,blue,'Symmetric, full period'),(asym,green,'Asymmetric, full period'),(c['records'],orange,'Doubled orbit')]:
        rr=sorted(rows,key=lambda r:r['beta']);axs[1].plot([r['beta'] for r in rr],[r['full_multipliers'][0][0] for r in rr],color=color,label=label,lw=1.8)
    axs[1].axhline(1,color='black',lw=.8);axs[1].axhline(-1,color='black',lw=.8)
    axs[1].set_ylim(-1.6,2.1);axs[1].set_ylabel('Dominant transverse Floquet multiplier')
    axs[1].set_xlabel('beta decreases to the right  →')
    axs[1].set_xlim(14.3315,14.3075);axs[1].legend(loc='center left',bbox_to_anchor=(0,.28),frameon=False,fontsize=9)
    for r,label in [(fold,'Fold'),(sym,'Symmetry\nbreaking'),(pd,'5 → 10\nmaxima'),(pd2,'Second\ndoubling')]:
        axs[1].text(r['beta'],1.76,label,ha='center',va='center',fontsize=9)
    fig.suptitle('Chua atlas: alpha = 9, m0 = −8/7, m1 = −5/7',fontsize=14)
    fig.savefig(assets/'chua-layer-bifurcations.png',dpi=180)
    fig.savefig(out/'bifurcations.pdf');plt.close(fig)
    fig=plt.figure(figsize=(11,7),layout='constrained');ax=fig.add_subplot(projection='3d')
    rows=b['surfaces'];aa=np.array([r['alpha'] for r in rows]).reshape(7,7);mm=np.array([r['m0'] for r in rows]).reshape(7,7)
    for key,col in [('fold',orange),('symmetry',blue),('period_doubling',green)]:
        z=np.array([r[key]['beta'] for r in rows]).reshape(7,7)
        ax.plot_surface(aa,mm,z,color=col,alpha=.65,edgecolor='white',linewidth=.25)
    ax.set_xlabel('alpha',labelpad=8);ax.set_ylabel('m0',labelpad=10);ax.set_zlabel('beta',labelpad=10)
    ax.set_title('The layer is bounded and subdivided by inclined bifurcation surfaces\n49 continuation points on each surface; m1 = −5/7',pad=15)
    ax.view_init(elev=24,azim=-53);ax.ticklabel_format(useOffset=False)
    ax.legend(handles=[Patch(color=orange,label='Fold'),Patch(color=blue,label='Symmetry breaking'),Patch(color=green,label='First period doubling')],loc='upper left',frameon=False)
    fig.savefig(assets/'chua-layer-surfaces.png',dpi=180);fig.savefig(out/'surfaces.pdf');plt.close(fig)
    print('saved chua-layer-bifurcations.png and chua-layer-surfaces.png',flush=True)


if __name__=='__main__':main()
