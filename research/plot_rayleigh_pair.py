"""Draw saved Rayleigh-pair data; no new trajectory searches or integrations."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rayleigh_pair import ROOT
from summarize_rayleigh_pair import completed_refinement


def main():
    plt.rcParams.update({'font.size':10,'axes.titlesize':11})
    scout=np.load(ROOT/'scout.npz');refined=completed_refinement()
    fig,axes=plt.subplots(2,2,figsize=(12.5,9.5),layout='constrained')
    mu=scout['mu_grid'];kappa=scout['kappa_grid'];n=len(mu)*len(kappa)*3
    value=scout['normalized_lle'][:n].reshape(len(mu),len(kappa),3).max(axis=2)
    ax=axes[0,0];im=ax.pcolormesh(mu,kappa,value.T,cmap='RdBu_r',vmin=-.02,vmax=.02,shading='nearest')
    ax.set(xscale='log',yscale='log',xlabel='Self-excitation / quadratic drag, μ',ylabel='Elastic coupling, κ',
        title='Scout: largest estimate across three starts\nPositive color alone does not establish chaos')
    fig.colorbar(im,ax=ax,label='Finite-time exponent per scaled time τ',extend='both')
    ax=axes[0,1];selected=np.flatnonzero(scout['preliminary_candidate'])
    for i in selected:
        row=np.flatnonzero(refined['source_indices']==i)[0]
        ax.plot([0,1],[scout['normalized_lle'][i],refined['normalized_lle'][row]],'o-',lw=1,ms=4,alpha=.7)
    ax.axhline(0,color='.4',lw=.7);ax.axhline(.01,color='.6',ls=':',lw=1)
    ax.set(xticks=[0,1],xticklabels=['Scout','Longer, tighter follow-up'],ylabel='Mean exponent per scaled time τ',
        title='All 13 preliminary candidates lose retention\n468 selected follow-ups completed; zero retained')
    cycle=np.load(ROOT/'cpu-row1865-cycle-check.npz');spec=np.load(ROOT/'cpu-row1865-spectrum.npz')
    ax=axes[1,0];s=cycle['states'];ax.plot(s[:,0],s[:,1],color='#2166ac',lw=1.5)
    ax.set(xlabel='Canonical displacement x₁',ylabel='Canonical velocity v₁',
        title='A stronger scout signal closes to a stable cycle\nRow 1865: x₂ = x₁, v₂ = v₁; Tσ = 6.12234')
    ax=axes[1,1];x=np.arange(1,5)
    ax.plot(x,scout['lle_windows'][1865],'o-',color='#d95f02',label='Scout windows')
    ax.plot(x,spec['spectrum_windows'][:,0],'s-',color='#2166ac',label='Later CPU leading exponent')
    ax.plot(x,spec['flow_neutral_windows'],'kx',ms=8,label='Neutral flow-direction control')
    ax.axhline(0,color='.4',lw=.7);ax.set(xticks=x,xlabel='Window within each separate observation',ylabel='Finite-time exponent per canonical σ',
        title='The late weak sign matches a neutral direction\nOther three CPU exponents are negative')
    ax.legend(fontsize=8)
    fig.suptitle('Quadratic-drag thermal spring pair: survey and periodic-orbit controls',fontsize=15)
    out=Path(__file__).resolve().parents[1]/'docs/assets/rayleigh-pair-search.png';fig.savefig(out,dpi=170);plt.close(fig);print(out)
    q=np.load(ROOT/'quadrature-cycle.npz');weak=np.load(ROOT/'weak-coupling-long.npz')
    fig=plt.figure(figsize=(12,5),layout='constrained');ax=fig.add_subplot(1,2,1,projection='3d');s=q['states']
    ax.plot(s[:,0],s[:,2],s[:,1],color='#d95f02',lw=1.7)
    ax.set(xlabel='x₁',ylabel='x₂',zlabel='v₁',title='An off-modal saddle cycle, not an attractor\n3D projection of the four-state orbit; Tσ = 6.27722')
    ax.view_init(25,-60)
    ax=fig.add_subplot(1,2,2);x=np.arange(1,5)
    for i,source in enumerate(weak['source_indices']):
        ax.plot(x,weak['normalized_lle_windows'][i],'o-',label=f'Original start {source%3+1}')
    ax.axhline(0,color='.4',lw=.7);ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax.set(xticks=x,xlabel='Window in the very long observation',ylabel='Finite-time exponent per scaled time τ',
        title='Three starts at the same weak-coupling point\n90,546 σ each; endpoints approach modal planes')
    ax.legend(fontsize=9)
    fig.suptitle('Weak coupling: slow transients and an unstable periodic control',fontsize=14)
    out=out.with_name('rayleigh-pair-weak-coupling.png');fig.savefig(out,dpi=170);plt.close(fig);print(out)
    shadow=np.load(ROOT/'saddle-shadowing.npz');fig,ax=plt.subplots(figsize=(11,4.6),layout='constrained')
    for name,label,color in [('quadrature','Off-modal saddle cycle','#d95f02'),('symmetric','Stable symmetric cycle','#2166ac'),('antisymmetric','Stable antisymmetric cycle','#1b9e77')]:
        ax.plot(shadow['times'][::5]/1000,shadow[name+'_distance'][::5],color=color,lw=1.5,label=label)
    scale=scout['time_scales'][193]
    ax.axvspan(50*scale/1000,150*scale/1000,color='.5',alpha=.12,label='GPU scout observation')
    ax.axvspan(350*scale/1000,750*scale/1000,color='#7570b3',alpha=.09,label='GPU refined observation')
    ax.axhline(.01,color='.5',ls=':',lw=1)
    ax.set(yscale='log',xlabel='Canonical time σ / 1,000',ylabel='Phase-minimized distance in four-state space',
        title='The actual weak-positive trajectory visits a periodic saddle, then leaves\nIndependent CPU trajectory; distance includes exchange and sign symmetry copies')
    ax.legend(fontsize=8,ncol=2,loc='lower right')
    out=out.with_name('rayleigh-pair-saddle-visit.png');fig.savefig(out,dpi=170);plt.close(fig);print(out)


if __name__=='__main__':main()
