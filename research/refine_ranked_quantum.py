"""Follow up large mean exponents even when the four-window screen fails.

This deliberately broadens the follow-up to reduce bias against intermittency.
It still cannot exclude weak chaos or basins not sampled by the chosen starts.
"""
from pathlib import Path
import json
import numpy as np
from quantum_neuron import QuantumGPU
from search_quantum import save

if __name__=='__main__':
    root=Path('results/candidate-research');p=[];starts=[];parents=[]
    old=set(json.loads((root/'scout/summary.json').read_text())['candidate_indices'])
    for name in ['scout','geometry-scout']:
        a=np.load(root/name/'scan.npz');chosen=[]
        for i in np.argsort(a['lle'])[::-1]:
            if a['status'][i] or (name=='scout' and i in old):continue
            chosen.append(i)
            if len(chosen)==24:break
        for i in chosen:
            p.append(a['parameters'][i]);starts.append(a['starts'][i]);parents.append({'scan':name,'index':int(i),'scout_lle':float(a['lle'][i])})
    out=root/'ranked-refined';out.mkdir(parents=True,exist_ok=True)
    (out/'parents.json').write_text(json.dumps(parents,indent=2))
    gpu=QuantumGPU();result=gpu.scan(np.array(p),np.array(starts),transient=5000,duration=8000,rtol=1e-9,atol=1e-12,
                                    progress=lambda f,s:print(f'{f:.1%}, {s:.1f}s',flush=True))
    save(out,result)
