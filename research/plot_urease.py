"""Render saved regular controls and the full-network GPU screening evidence."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from urease import ROOT,SCALE


def main():
    scout=np.load(ROOT/'scout.npz');ref=np.load(ROOT/'reference.npz')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    ax=axes[0,0];keep=ref['times']>=8000
    ax.plot((ref['times'][keep]-8000)/60,-np.log10(SCALE*ref['states'][keep,1]),color='#176c91',lw=1)
    ax.set(xlabel='Time within late reference segment (min)',ylabel='pH',title='A  Regular control: 17.79399 min\nPublished period: 17.79 +/- 0.01 min')
    ax=axes[0,1];p=scout['parameters'];equal=p[:,7]==p[:,5]
    for mask,color,label in [(~equal,'#176c91','Unequal product/substrate rates'),(equal,'#bc5835','Equal product/substrate rates')]:
        ax.scatter(p[mask,5],scout['lle'][mask],s=8,c=color,alpha=.45,label=label)
    ax.axhline(0,color='black',lw=.7);ax.set_xscale('log');ax.set_ylim(-.002,.0005)
    ax.set(xlabel='Substrate transport kS (1/s)',ylabel='Mean finite-time exponent (1/s)',title='B  1,024 complete GPU trajectories')
    ax.legend(fontsize=7,loc='lower left')
    ax.text(.04,.94,'None has four positive windows.\nValues below -0.002 are outside this view.',transform=ax.transAxes,fontsize=8,va='top')
    ax=axes[1,0]
    for i,color in [(127,'#bc5835'),(494,'#23856a'),(0,'#176c91')]:
        ax.plot(np.arange(1,5),scout['lle_windows'][i],'.-',color=color,label=f'Scout {i}' if i else 'Source control')
    ax.axhline(0,color='black',lw=.7);ax.set_xticks(range(1,5));ax.legend(fontsize=8)
    ax.set(xlabel='Successive 10,000-second windows',ylabel='Finite-time exponent (1/s)',title='C  Positive means without persistence')
    ax=axes[1,1]
    for file,color,label in [('source-cycle','#176c91','Source'),('scout-cycle-127','#bc5835','Scout 127'),('scout-cycle-494','#23856a','Scout 494')]:
        data=np.load(ROOT/(file+'.npz'));meta=json.loads((ROOT/(file+'.json')).read_text());state=data['states']
        ax.plot(state[:,0]*SCALE*1e6,-np.log10(state[:,1]*SCALE),color=color,lw=.9,label=f'{label}: {meta["period_minutes"]:.3f} min')
    ax.set(xlabel='Urea concentration (micromolar)',ylabel='pH',title='D  Independently closed stable cycles')
    ax.legend(fontsize=8)
    fig.suptitle('Full four-species urease network: no validated chaos candidate',fontsize=15)
    output=Path(__file__).resolve().parents[1]/'docs/assets/urease_survey.png'
    fig.savefig(output,dpi=160);plt.close(fig);print(output)


if __name__=='__main__':main()
