"""Merge equivalent transported roots, retaining their branch ancestry."""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output;existing=json.loads((out/'sheet-ancestry.json').read_text());extras=json.loads((out/'offcenter-ancestry.json').read_text())
    groups=[]
    for r in existing+extras:
        if r['status']!='identified':continue
        match=next((g for g in groups if g['symmetric']==r['symmetric'] and abs(g['core']['beta']-r['root']['beta'])<1e-7 and abs(g['core']['period']-r['root']['period'])<1e-5),None)
        source=dict(id=r['id'],chain=r['chain'])
        if match:match['sources'].append(source)
        else:groups.append(dict(id=r['id'] if not r['id'].startswith('off') else None,symmetric=r['symmetric'],core=r['root'],
                               regions=r['regions'],targets=r['targets'],target=.5 if r['symmetric'] else -.5,sources=[source]))
    new=sorted([g for g in groups if g['id'] is None],key=lambda g:g['core']['beta'])
    for k,g in enumerate(new):g['id']=f'E{k+1:02d}'
    model=Chua(rtol=2e-12,step=.02)
    for g in groups:
        r=g['core'];u=np.r_[0.,r['q'][:2]];orbit=model.integrate(r['beta'],u,r['period'])
        g['peak_count']=len(orbit.t_events[3]);g['section_returns']=len(orbit.t_events[2][orbit.t_events[2]<r['period']-1e-6])
        g['source_count']=len(g['sources'])
    report=dict(family_count=len(groups),existing_count=len(existing),additional_count=len(new),families=groups,
        identity_rule='Same inversion-symmetry class, normalized multiplier, reference-plane beta within 1e-7, and full period within 1e-5, after tracing symmetry and doubling parents.',
        scope='A numerical catalog of continued window families, not proof that all periodic windows in the box have been found or that no families connect outside the continued patches.')
    (out/'family-registry.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    (out/'additional-templates.json').write_text(json.dumps(new,indent=2)+'\n',encoding='utf-8',newline='\n')
    roots=[g for g in groups if g['id'] in ['A03','A05']]
    (out/'parent-templates.json').write_text(json.dumps(roots,indent=2)+'\n',encoding='utf-8',newline='\n')
    for g in sorted(groups,key=lambda g:g['core']['beta']):print(g['id'],g['peak_count'],g['core']['beta'],g['core']['period'],[s['id'] for s in g['sources']],flush=True)
    print('families',len(groups),'additional',len(new),flush=True)


if __name__=='__main__':main()
