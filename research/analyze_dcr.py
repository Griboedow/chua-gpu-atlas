"""Summarize archived DCR controls without starting new integrations."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results/dcr'


def main():
    labels=['ida-source-short','ida-analytic-1','ida-analytic-1-tight',
            'ida-source-1-tight','ida-analytic-1-fine','ida-44-30','ida-44-30-tight']
    ref=np.load(DATA/'ida-analytic-1-fine.npz')
    rows=[];hashes={}
    fig,axes=plt.subplots(1,3,figsize=(13,3.7),layout='constrained')
    colors=['#2b78a3','#e6922e','#7c5295','#2a987e']
    for i,label in enumerate(labels):
        meta=json.loads((DATA/f'{label}.json').read_text())
        a=np.load(DATA/f'{label}.npz')
        for suffix in ['json','npz']:
            path=DATA/f'{label}.{suffix}'
            payload=path.read_bytes()
            if suffix=='json':payload=payload.replace(b'\r\n',b'\n')
            hashes[str(path.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(payload).hexdigest()
        row=dict(label=label,success=meta['success'],settings=meta['settings'],
                 wall_seconds=meta['wall_seconds'],callback_calls=meta['callback_calls'])
        if meta['success']:
            row.update(actual_last_output_time=float(a['time'][-1]),nfev=meta['nfev'],njev=meta['njev'],
                       displacement_range=[float(a['states'][:,0].min()),float(a['states'][:,0].max())],
                       sampled_contact_range=[int(a['contacts'].min()),int(a['contacts'].max())])
            if len(a['time'])==len(ref['time']):
                assert np.array_equal(a['time'],ref['time'])
                err=abs(a['states']-ref['states'])
                row.update(max_state_difference_from_fine=float(err.max()),
                           max_displacement_difference_from_fine=float(err[:,0].max()),
                           end_state_difference_from_fine=float(err[-1].max()))
                axes[0].plot(a['time'],a['states'][:,0],lw=1,label=label.removeprefix('ida-'),color=colors[i-1])
                if label!='ida-analytic-1-fine':
                    axes[1].semilogy(a['time'],np.maximum(err[:,0],1e-16),lw=1,color=colors[i-1],label=label.removeprefix('ida-'))
        else:
            row.update(last_callback_time=meta['last_callback_time'],message=meta['message'],
                       interpretation='Budget-unresolved. Saved callback state can be an unaccepted solver trial; not a certified trajectory endpoint.')
        rows.append(row)
    axes[0].set(xlabel='Dimensionless time',ylabel='Backbone displacement x',title='Short strict-contact controls')
    axes[0].legend(fontsize=6.5)
    axes[1].set(xlabel='Dimensionless time',ylabel='Absolute displacement difference',title='Against finest archived run',ylim=(1e-12,1e-3))
    axes[1].legend(fontsize=6.5)
    axes[2].step(ref['time'],ref['contacts'],where='post',lw=.7,color='#2a987e')
    axes[2].set(xlabel='Dimensionless time',ylabel='Number of active contacts',title='Fine run: sampled contact count')
    fig.suptitle('DCR: numerical source controls only — no attractor classification',fontsize=12)
    asset=ROOT/'docs/assets/dcr-controls.png'
    fig.savefig(asset,dpi=170);plt.close(fig)
    result=dict(reference='ida-analytic-1-fine',runs=rows,source_sha256=hashes,
                source_hash_convention='SHA256: exact NPZ bytes; JSON bytes with CRLF normalized to LF for portable Git checkouts.',
                scope='Five completed short controls and two budget-unresolved longer attempts. Agreement away from contact events does not define a continuous-time contact-selection law.')
    (DATA/'control-summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
