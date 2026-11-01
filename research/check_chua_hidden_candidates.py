"""Longer, smaller-step controls and equilibrium accessibility of search candidates."""
import os
import sys
import numpy as np
from chua_research_followup import ROOT,OUT,dump
from trace_chua_window import Chua
import json


def main():
    os.environ['PYOPENCL_CACHE_DIR']=str(ROOT/'.publication/opencl-cache')
    os.environ['TEMP']=str(ROOT.parent/'.agent-tmp');os.environ['TMP']=os.environ['TEMP']
    sys.path.insert(0,str(ROOT));from chua_gpu import GPU
    source=json.loads((OUT/'hidden-search/outer-screen.json').read_text())['candidates']
    parameters=np.array([[r['parameters'][k] for k in ['alpha','beta','m0','m1']] for r in source])
    initial=[]
    for row in source:
        p=row['parameters'];A=Chua(p['alpha'],p['m0'],p['m1']).jacobian(p['beta'],True)
        values,V=np.linalg.eig(A);v=V[:,np.argmax(values.real)].real;v/=np.linalg.norm(v)
        initial.append([row['initial'],row['final'],*(sgn*eps*v for eps in [1e-4,1e-7,1e-10] for sgn in [-1,1])])
    initial=np.array(initial);gpu=GPU();rows=[]
    for dt in [.0025,.00125]:
        result=gpu.scan(parameters,initial_states=initial,dt=dt,transient=1200,duration=3200)
        np.savez_compressed(OUT/'hidden-search'/f'candidate-control-{dt}.npz',parameters=parameters,initial=initial,
            **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        for i,row in enumerate(source):
            p=row['parameters'];xe=(p['m1']-p['m0'])/(1+p['m1']);eq=np.array([[xe,0,-xe],[-xe,0,xe]])
            distances=np.linalg.norm(result['final_states'][i,:,None,:]-eq[None,:,:],axis=2).min(axis=1)
            rec=dict(parameter_index=row['parameter_index'],parameters=p,dt=dt,
                classes=result['classes'][i].tolist(),lle=result['lle'][i].tolist(),
                windows=result['lle_windows'][i].tolist(),distance_to_outer_equilibria=distances.tolist())
            rows.append(rec);print('control',row['parameter_index'],dt,rec['classes'],flush=True)
        dump(OUT/'hidden-search/candidate-controls.json',dict(records=rows,
            seed_order=['original finite-amplitude start','continued last state','negative/positive unstable eigenvector at radii 1e-4,1e-7,1e-10'],
            scope='Finite-time controls. Equilibrium-neighborhood exclusion and homoclinic tangency remain separate proof obligations.'))


if __name__=='__main__':main()
