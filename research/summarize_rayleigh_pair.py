"""Summarize saved searches, merging exact resumptions without repeating runs."""
import json
from pathlib import Path
import numpy as np
from rayleigh_pair import ROOT


def completed_refinement():
    with np.load(ROOT/'scout-refined.npz') as archive:
        result={key:archive[key].copy() for key in archive.files}
    with np.load(ROOT/'scout-refined-resumed.npz') as resumed:
        rows=resumed['resumed_rows']
        if not np.array_equal(result['source_indices'][rows],resumed['source_indices']):
            raise ValueError('Resume source rows do not match')
        for key in result:
            if result[key].ndim and result[key].shape[0]==len(result['status']) and key in resumed:
                result[key][rows]=resumed[key]
    return result


def main():
    scout=np.load(ROOT/'scout.npz');refined=completed_refinement();weak=np.load(ROOT/'weak-coupling-long.npz')
    selected=np.flatnonzero(scout['preliminary_candidate'])
    follow_rows=np.array([np.flatnonzero(refined['source_indices']==i)[0] for i in selected])
    records=[]
    for i,row in zip(selected,follow_rows):
        records.append(dict(source_index=int(i),parameters=scout['parameters'][i].tolist(),
            scout_scaled_mean=float(scout['normalized_lle'][i]),scout_scaled_windows=scout['normalized_lle_windows'][i].tolist(),
            followup_scaled_mean=float(refined['normalized_lle'][row]),followup_scaled_windows=refined['normalized_lle_windows'][row].tolist()))
    final=weak['final_states'];dplus=np.linalg.norm(final[:,:2]-final[:,2:],axis=1)/np.sqrt(2)
    dminus=np.linalg.norm(final[:,:2]+final[:,2:],axis=1)/np.sqrt(2)
    meta=[json.loads((ROOT/(name+'.json')).read_text()) for name in ['scout','scout-refined','scout-refined-resumed','weak-coupling-long']]
    result=dict(scope='Finite-time negative search with independent periodic-orbit controls; no global exclusion or novelty claim.',
        scout=dict(trajectories=len(scout['status']),completed=int(np.sum(scout['status']==0)),preliminary=int(len(selected)),
            all_four_positive_signs=int(np.sum(np.all(scout['normalized_lle_windows']>0,axis=1)))),
        refinement=dict(trajectories=len(refined['status']),completed_after_exact_resume=int(np.sum(refined['status']==0)),
            retained=int(np.sum(refined['retained_candidate'])),all_four_positive_signs=int(np.sum(np.all(refined['normalized_lle_windows']>0,axis=1))),
            maximum_scaled_mean=float(np.max(refined['normalized_lle'])),
            resume_scope='Two interrupted rows continued with unchanged state, tangent, clocks, raw window logs and statistics.'),
        preliminary_followups=records,
        weak_followup=dict(source_indices=weak['source_indices'].tolist(),completed=int(np.sum(weak['status']==0)),retained=int(np.sum(weak['retained_candidate'])),
            total_sigma=weak['elapsed_sigma'].tolist(),canonical_means=weak['lle'].tolist(),scaled_windows=weak['normalized_lle_windows'].tolist(),
            endpoint_distance_to_symmetric_plane=dplus.tolist(),endpoint_distance_to_antisymmetric_plane=dminus.tolist(),
            qualification='Distance to invariant planes at endpoints, not a proof of complete basin coverage or asymptotic orbit convergence.'),
        elapsed_gpu_wall_seconds={name:item['wall_seconds'] for name,item in zip(['scout','refinement','resume','weak_followup'],meta)},
        source_files=['scout.npz','scout-refined.npz','scout-refined-resumed.npz','weak-coupling-long.npz'])
    (ROOT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['scout','refinement','weak_followup']},indent=2))


if __name__=='__main__':main()
