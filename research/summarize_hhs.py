"""Recompute counts from saved HHS single-cell scans without integration."""
import json
import numpy as np
from hhs import ROOT


def main():
    s=np.load(ROOT/'scout.npz');r=np.load(ROOT/'scout-refined.npz')
    stages={}
    for name,a in [('scout',s),('refinement',r)]:
        ok=a['status']==0
        stages[name]=dict(trajectories=len(ok),completed=int(ok.sum()),statuses={str(int(k)):int(v) for k,v in zip(*np.unique(a['status'],return_counts=True))},
            positive_means=int(np.sum(ok&(a['lle']>0))),all_four_positive_windows=int(np.sum(ok&np.all(a['lle_windows']>0,axis=1))),
            maximum_completed_mean=float(np.max(a['lle'][ok])) if ok.any() else None,
            minimum_sampled_length_at_alpha_point1=float(np.min(a['minimum_observed_length_at_alpha_point1'][ok])) if ok.any() else None)
    rows=np.flatnonzero(s['preliminary_candidate']);records=[]
    for i in rows:
        j=int(np.flatnonzero(r['source_indices']==i)[0]);done=r['status'][j]==0
        records.append(dict(source_row=int(i),parameters=s['parameters'][i].tolist(),scout_mean=float(s['lle'][i]),
            refined_status=int(r['status'][j]),refined_mean=float(r['lle'][j]) if done else None,
            refined_windows=r['lle_windows'][j].tolist() if done else None,retained=bool(r['retained_candidate'][j])))
    result=dict(scope='Literal fixed-reference five-state S2; not a reproduction of all published figures or a novelty claim.',stages=stages,
        preliminary_count=int(len(rows)),retained_count=int(np.sum(r['retained_candidate'])),preliminary_followups=records,
        incomplete_source_indices=r['source_indices'][r['status']!=0].tolist(),
        wall_seconds={name:json.loads((ROOT/(file+'.json')).read_text())['wall_seconds'] for name,file in [('scout','scout'),('refinement','scout-refined')]},
        qualifications=['Incomplete trajectories remain unclassified.','Positive length at accepted observation states is not physical calibration or a continuous-time proof.','Small positive completed exponents below thresholds are not all classified as periodic.'])
    result['qualifications'][0]='Incomplete GPU intervals are not exponent-classified; a separate endpoint orbit check does not complete those intervals.'
    result['independent_endpoint_controls']=[]
    for row in [1445,3009]:
        path=ROOT/f'cpu-row{row}-cycle-check.json'
        if path.exists():
            d=json.loads(path.read_text())
            result['independent_endpoint_controls'].append(dict(source_row=row,record=path.name,status=d['status'],period=d.get('period'),
                maximum_transverse_multiplier=d.get('maximum_transverse_multiplier'),scope='Local numerical orbit control, not completion of a missing GPU interval or global basin classification.'))
    (ROOT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8');print(json.dumps(result['stages'],indent=2));print('retained',result['retained_count'])


if __name__=='__main__':main()
