"""Saved-data figures for the core-shell gel's negative and unresolved survey."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core_shell_gel import ROOT


def main():
    scout=np.load(ROOT/'scout.npz');ref=np.load(ROOT/'scout-refined.npz')
    control=np.load(ROOT/'source-reduced-control.npz');cycle=np.load(ROOT/'scout-cycle-1378.npz')
    meta=json.loads((ROOT/'scout-cycle-1378.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    ax[0,0].plot(control['states'][:,0],control['states'][:,1],color='#185c84',lw=1.3)
    ax[0,0].set(xlabel='Proton concentration x',ylabel='Polymer fraction phi',title='Published two-state approximation: periodic control')
    ok=scout['status']==0
    ax[0,1].scatter(scout['parameters'][ok,6],scout['lle'][ok],s=7,c='#9aa8b1',alpha=.4,rasterized=True)
    ids=np.flatnonzero(scout['preliminary_candidate'])
    ax[0,1].scatter(scout['parameters'][ids,6],scout['lle'][ids],s=45,c='#b84f31',label='Four preliminary candidates')
    ax[0,1].axhline(0,color='black',lw=.7);ax[0,1].set_xscale('log')
    ax[0,1].set(xlabel='Glucose relaxation rate lambda_g',ylabel='Finite-time leading exponent',title='4,096 full-model trajectories; 27 scout cases unresolved')
    ax[0,1].set_ylim(-.02,.012);ax[0,1].legend(fontsize=8)
    ax[0,1].text(.03,.03,'Negative values below -0.02 are outside this view',transform=ax[0,1].transAxes,fontsize=8)
    selected=[1188,1378,1575,2882]
    for sid in selected:
        row=int(np.flatnonzero(ref['source_indices']==sid)[0])
        if ref['status'][row]==0:ax[1,0].plot(np.arange(1,5),ref['lle_windows'][row],'.-',label=f'Scout {sid}')
    ax[1,0].axhline(0,color='black',lw=.7);ax[1,0].set(xticks=[1,2,3,4],xlabel='Successive 1,000-unit observation windows',ylabel='Finite-time leading exponent',title='Three signals lose persistence; fourth hits step budget')
    ax[1,0].legend(fontsize=8)
    ax[1,1].plot(cycle['states'][:,2],cycle['states'][:,3],color='#168075',lw=1.5)
    ax[1,1].set(xlabel='Proton concentration x',ylabel='Polymer fraction phi',title='Independent closed orbit at scout 1378')
    ax[1,1].text(.15,.55,f'Period {meta["period"]:.6f}\nLargest transverse multiplier {meta["spectral_radius"]:.2e}',transform=ax[1,1].transAxes,fontsize=10)
    fig.suptitle('Core-shell enzyme gel: no validated chaos candidate in this survey',fontsize=16)
    target=Path(__file__).resolve().parents[1]/'docs/assets/core_shell_gel_survey.png'
    fig.savefig(target,dpi=160);plt.close(fig);print(target)


if __name__=='__main__':main()
