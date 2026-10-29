"""Recheck the original wide survey with the corrected full-field DOP853 reference.

The v1.0 affine-branch event routine is superseded: missed grazing crossings
could leave its integration on the wrong branch. Historical event_resolved
JSON files are preserved; new results use a distinct full_field_checks name.
"""
import numpy as np
from explore_extent import OUT
from run import save_json
from independent_reference import reference_lle

if __name__=='__main__':
    data=np.load(OUT/'verified.npz')
    keep=(data['classes']==4).any(axis=1)
    indices=np.flatnonzero(keep)
    chosen=np.unique([indices[np.argmax(data['parameters'][indices,0])],
                      indices[np.argmax(data['parameters'][indices,1])],
                      indices[np.argmin(data['parameters'][indices,2])]])
    results=[]
    for idx in chosen:
        p=data['parameters'][idx]
        seed=int(np.argmax(np.where(data['classes'][idx]==4,data['lle'][idx],-np.inf)))
        print('Checking parameters',p,'seed',seed,flush=True)
        result=reference_lle(p,data['final_states'][idx,seed])
        result.update({'parameters':p.tolist(),'gpu_lle':float(data['lle'][idx,seed]),'seed_index':seed})
        results.append(result)
        save_json(OUT/'full_field_checks.json',{'method':'DOP853 evaluates the actual full piecewise field at every stage; crossing events are diagnostic and do not select the branch.',
                   'checks':results,'scope':'Selected examples only; starts from the end of the GPU observation, then discards an additional transient.'})
        print(result,flush=True)
