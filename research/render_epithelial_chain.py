"""Render a physical-domain audit, without labeling short positive exponents as discoveries."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from epithelial_chain import threshold


def main():
    base=Path('results/epithelial-chain');s=np.load(base/'n16-scout/scan.npz');r=np.load(base/'n16-refined/scan.npz')
    meta=json.loads((base/'n16-refined/summary.json').read_text());ids=np.array(meta['parent_indices']);cpu=json.loads((base/'cpu-340/summary.json').read_text())
    p=s['parameters'];gain=p[:,2]*p[:,3]/threshold(16,p[:,0],p[:,1]);positive=s['lle']>.005
    candidate=positive&(s['lle_windows'].min(axis=1)>.002)
    physically_clear=(s['domain'][:,[0,2]]>0).all(axis=1)
    retained=(r['lle']>.005)&(r['lle_windows'].min(axis=1)>.002)
    summary=dict(scouts=len(p),completed=int((s['status']==0).sum()),preliminary_candidates=int(candidate.sum()),
                 candidate_positive_cell_lengths=int((candidate&(s['domain'][:,0]>0)).sum()),
                 positive_cell_lengths=int((s['domain'][:,0]>0).sum()),positive_lengths_and_concentration=int(physically_clear.sum()),
                 positive_all_three_domain_checks=int((s['domain'][:,:3]>0).all(axis=1).sum()),
                 selected_followups=len(ids),retained_selected=int(retained.sum()),retained_selected_positive_lengths=int((retained&(r['domain'][:,0]>0)).sum()),
                 selection=meta['selection'],cpu_control=cpu,
                 limitation='Only 4 of 343 preliminary candidates selected for long follow-up because every preliminary candidate had negative cell lengths. Other candidates remain finite-time mathematical observations.')
    (base/'survey-summary.json').write_text(json.dumps(summary,indent=2))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(14.5,4.4),layout='constrained')
    bad=~physically_clear
    ax[0].scatter(gain[bad],s['lle'][bad],s=11,c='#cf7a60',alpha=.65,label='Length or concentration fails')
    ax[0].scatter(gain[physically_clear],s['lle'][physically_clear],s=14,c='#277da8',label='Length + concentration positive')
    ax[0].set(xscale='log',xlabel='Gain / finite-chain instability threshold',ylabel='Largest finite-time exponent',title='16 cells: 1,024 GPU scouts');ax[0].legend(fontsize=8)
    ax[1].scatter(s['domain'][:,0],s['lle'],s=11,c='#7a8e9d',alpha=.6)
    ax[1].scatter(s['domain'][candidate,0],s['lle'][candidate],s=13,c='#cc593e',label='Passes preliminary LLE screen')
    ax[1].axvline(0,color='black',lw=1);ax[1].set(xlabel='Minimum cell length / mean length',ylabel='Largest finite-time exponent',title='All 343 candidates have cell crossings');ax[1].legend(fontsize=8)
    ax[2].plot(np.arange(1,5),s['lle_windows'][340],'o-',color='#bb7a32',label='Short GPU scout')
    local=int(np.flatnonzero(ids==340)[0]);ax[2].plot(np.arange(1,5),r['lle_windows'][local],'o-',color='#277da8',label='Long GPU, half step')
    ax[2].plot(np.arange(1,5),cpu['lle_windows'],'o-',color='#49854b',label='Independent DOP853')
    ax[2].axhline(0,color='black',lw=.7);ax[2].set(xlabel='Observation window (different horizons)',ylabel='Largest finite-time exponent',title='Domain-clear weak-positive control #340');ax[2].legend(fontsize=8)
    fig.suptitle('Cell-train audit: mathematical expansion must be checked against the physical domain',fontsize=14)
    fig.savefig('docs/assets/epithelial_chain_audit.png',dpi=170);plt.close(fig)
    print(summary)


if __name__=='__main__':main()
