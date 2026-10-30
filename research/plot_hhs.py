"""Render saved literal single-cell controls and searches; no integrations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hhs import ROOT


def main():
    plt.rcParams.update({'font.size':10,'axes.titlesize':11})
    source=np.load(ROOT/'uncoupled-cycle.npz');s=np.load(ROOT/'scout.npz');r=np.load(ROOT/'scout-refined.npz')
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained');ax=axes[0,0]
    z=source['states'];ax.plot(z[:,0],z[:,1],color='#2166ac',lw=1.5)
    ax.set(xlabel='Chemical concentration M',ylabel='Chemical concentration E',title='Exact uncoupled periodic control (μ = 0)\nPeriod = 13.01945 units of τR')
    ax=axes[0,1];x=np.arange(5)
    ax.plot(x,source['expected_triangular_exponents'],'o',ms=9,mfc='none',color='#d95f02',label='Exact block-triangular prediction')
    ax.plot(x,source['exponents'],'x',ms=7,color='#2166ac',label='Full five-state Floquet calculation')
    ax.axhline(0,color='.6',lw=.6);ax.set(xticks=x,xticklabels=np.arange(1,6),xlabel='Sorted exponent index',ylabel='Exponent per unit of τR',title='Neutral flow direction and four stable directions\nIndependent spectrum agreement better than 10⁻¹³')
    ax.legend(fontsize=8)
    ax=axes[1,0];mu=s['parameters'][:195:3,0];means=s['lle'][:195].reshape(65,3)
    for i in range(3):ax.plot(mu,means[:,i],'.-',ms=3,lw=.7,label=f'Start {i+1}')
    ax.axhline(0,color='.5',lw=.7);ax.set(xlabel='Mechanochemical feedback μ',ylabel='Finite-time largest estimate / τR',title='Source-time-constant section: three starts\n500 discarded + 1,000 observed')
    ax.legend(fontsize=8)
    ax=axes[1,1]
    selected=np.flatnonzero(s['preliminary_candidate']);missing=0
    for i in selected:
        j=np.flatnonzero(r['source_indices']==i)[0]
        if r['status'][j]!=0:missing+=1;continue
        ax.plot([0,1],[s['lle'][i],r['lle'][j]],'o-',ms=4,lw=1,alpha=.8)
    ax.axhline(0,color='.5',lw=.7);ax.axhline(.0005,color='.6',ls=':',lw=1)
    ax.set(xticks=[0,1],xticklabels=['Scout','Longer, tighter follow-up'],ylabel='Mean largest estimate / τR',
        title=f'Follow-up of {len(selected)} preliminary candidates\n{missing} incomplete cases omitted from this comparison')
    fig.suptitle('Mechanochemical single cell: controls before chaos claims',fontsize=15)
    out=Path(__file__).resolve().parents[1]/'docs/assets/hhs-controls.png';fig.savefig(out,dpi=170);plt.close(fig);print(out)
    cpu=np.load(ROOT/'cpu-row1445-spectrum.npz');cycle=np.load(ROOT/'cpu-row1445-cycle-check.npz')
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained');ax=axes[0];x=np.arange(1,5)
    ax.plot(x,s['lle_windows'][1445],'o-',color='#d95f02',label='Scout: four 250-unit windows')
    ax.plot(x,cpu['spectrum_windows'][:,0],'s-',color='#2166ac',label='Later CPU: four 1,000-unit windows')
    ax.plot(x,cpu['neutral_flow_windows'],'kx',ms=8,label='Neutral flow-direction control')
    ax.axhline(0,color='.5',lw=.7);ax.set(xticks=x,xlabel='Window within each separate observation',ylabel='Finite-time exponent per unit of τR',
        title='A positive window can measure phase speed\nThe late estimates match a neutral direction')
    ax.legend(fontsize=8,loc='center left')
    z=cycle['states'];ax=axes[1];ax.plot(z[:,0],z[:,1],color='#2166ac')
    ax.set(xlabel='Chemical concentration M',ylabel='Chemical concentration E',
        title='Closed stable five-state cycle, shown in (M,E)\nT = 12.46623; dominant transverse multiplier = 0.005297')
    fig.suptitle('Row 1445: independent full-spectrum and original-coordinate checks',fontsize=13)
    out=out.with_name('hhs-neutral-control.png');fig.savefig(out,dpi=170);plt.close(fig);print(out)


if __name__=='__main__':main()
