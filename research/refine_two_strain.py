"""Refine all screened positives and the twelve largest finite-time estimates."""
import argparse
import numpy as np
from two_strain import ROOT, save
from two_strain_gpu import TwoStrainGPU

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--step-control',action='store_true');args=parser.parse_args()
    if args.step_control:
        source=np.load(ROOT/'slice-refined.npz');indices=np.array([17,32,39])
        data=TwoStrainGPU().scan(source['parameters'][indices],source['final_log_states'][indices],
            transient=100000,duration=400000,rtol=1e-11,max_step=.25,max_attempts=10000000)
        data['refined_index']=indices;save(ROOT/'step-control',data)
        print(data['lle']);print(data['lle_windows']);raise SystemExit(0)
    source=np.load(ROOT/'slice-scout.npz')
    indices=np.union1d(np.flatnonzero(source['preliminary_candidate']),np.argsort(source['lle'])[-12:])
    data=TwoStrainGPU().scan(source['parameters'][indices],source['final_log_states'][indices],
        transient=100000,duration=400000,rtol=1e-10,max_step=.5,max_wall_seconds=900)
    data['source_index']=indices;data['start_kinds']=source['start_kinds'][indices]
    ok=data['status']==0;positive=ok&(data['lle']>1e-5)&np.all(data['lle_windows']>5e-6,axis=1)
    data['preliminary_candidate']=positive
    data.update(lanes=len(indices),completed=int(ok.sum()),retained=int(positive.sum()),
        scope='Longer and tighter continuation of all scout positives plus top twelve estimates; surviving finite-time expansion does not prove attraction')
    save(ROOT/'slice-refined',data)
    print({k:data[k] for k in ['lanes','completed','retained','wall_seconds']},flush=True)
