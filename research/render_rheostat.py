"""Aggregate completed follow-ups and draw the rheostat model audit."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    base=Path('results/rheostat');records={};counts={};timings={}
    for name in ['refined','refined-retry']:
        data=np.load(base/name/'scan.npz');meta=json.loads((base/name/'summary.json').read_text())
        timings[name]=meta['timing'];counts[name]=meta['status_counts']
        for i,index in enumerate(meta['parent_indices']):
            records[index]=dict(index=index,status=int(data['status'][i]),lle=float(data['lle'][i]),windows=data['lle_windows'][i].tolist())
    summary=dict(scout_count=4096,scout_candidates=0,distinct_followups=len(records),
                 attempt_statuses=counts,final_statuses={str(i):sum(x['status']==i for x in records.values()) for i in set(x['status'] for x in records.values())},
                 timings=timings,records=list(records.values()))
    assert len(records)==24 and all(x['status']==0 for x in records.values())
    summary['retained_followups']=sum(x['lle']>.005 and min(x['windows'])>.002 for x in records.values())
    summary['maximum_followup_mean']=max(x['lle'] for x in records.values())
    (base/'survey-summary.json').write_text(json.dumps(summary,indent=2))
    print('Completed',len(records),'retained',summary['retained_followups'],'max mean',summary['maximum_followup_mean'])
    print('Original scout777',records[777])

    control=json.loads((base/'controls/summary.json').read_text())
    boundary=np.load(base/'controls/linear-boundary.npz')
    trace=np.load(base/'controls/strongest-scout-control.npz')['state']
    fig,axes=plt.subplots(1,3,figsize=(13,4.3),layout='constrained')
    ax=axes[0];r=boundary['rate']
    ax.loglog(r,boundary['gain_critical'],label='Full 3D linear threshold',color='#276eab')
    ax.loglog(r,.2*r,'--',label='Fast-thermal approximation',color='#d48636')
    ax.axhline(1,color='black',lw=.8)
    low,high=control['linear_example']['rate_boundaries']
    ax.axvspan(low,high,alpha=.12,color='#c94242',label='Unstable equilibrium at g=1')
    ax.set(xlabel='Thermal rate r (inverse relaxation time)',ylabel='Feedback gain g',title='A  Instability requires an intermediate rate',ylim=(.01,100))
    ax.legend(fontsize=7,loc='upper center')
    ax=axes[1]
    ax.plot(trace[:,0],trace[:,1],color='#276eab',lw=1.4)
    ax.set(xlabel='Displacement X',ylabel='Velocity V',title='B  Largest scout exponent: stable cycle')
    orbit=control['orbits'][-1]
    ax.text(.03,.97,f"Period {orbit['period']:.6f}\nLeading transverse multiplier 0.678365",transform=ax.transAxes,va='top',fontsize=8,
            bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
    ax=axes[2]
    near=control['orbits'][:3];excess=np.array([x['excess_gain'] for x in near])
    ax.loglog(excess,[x['leading_amplitude_X'] for x in near],'--',color='#d48636',label='Leading linear amplitude law')
    ax.loglog(excess,[x['amplitude_X'] for x in near],'o',color='#276eab',label='Independently closed orbits')
    ax.set(xlabel='Excess gain g - g_H',ylabel='Displacement half-amplitude',title='C  Quadratic drag: amplitude scales linearly')
    ax.legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Rheostat-driven elastomer: no retained chaos in 4,096 scouts and 24 long checks\nExact full-model reduction; finite sampling is not a global exclusion theorem',fontsize=12)
    fig.savefig('docs/assets/rheostat_audit.png',dpi=170)


if __name__=='__main__':main()
