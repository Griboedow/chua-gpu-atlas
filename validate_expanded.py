"""Independently check representative extremes of the expanded parameter cloud."""
import json
import numpy as np
from adaptive_gpu import AdaptiveGPU
from explore_expanded import OUT
from run import save_json
from independent_reference import reference_lle

def main():
    data=np.load(OUT/'verified.npz');indices=np.flatnonzero((data['classes']==4).any(axis=1))
    p=data['parameters'][indices]
    selected=[indices[np.argmax(p[:,0])],indices[np.argmax(p[:,1])],indices[np.argmin(p[:,2])]]
    distant=indices[(p[:,0]>768)]
    if len(distant):
        selected.append(distant[np.argmin(data['parameters'][distant,2])])
    selected=list(dict.fromkeys(map(int,selected)));g=AdaptiveGPU();checks=[]
    for idx in selected:
        point=data['parameters'][idx]
        seed=int(np.argmax(np.where(data['classes'][idx]==4,data['lle'][idx],-np.inf)))
        initial=data['final_states'][idx,seed]
        print('Independent expanded check:',point.tolist(),flush=True)
        refined=g.scan(point[None,:],initial_states=initial[None,None,:],transient=20,duration=80,
                       rtol=1e-9,atol=1e-12)
        reference=reference_lle(point,initial,transient=20,duration=80)
        check={'parameters':point.tolist(),'seed_index':seed,'initial_state':initial.tolist(),
               'survey_lle':float(data['lle'][idx,seed]),'matched_adaptive_lle':float(refined['lle'][0,0]),
               'matched_adaptive_windows':refined['lle_windows'][0,0].tolist(),
               'independent':reference,'positive_in_all_independent_windows':bool(min(reference['windows'])>.005 and reference['lle']>.01)}
        checks.append(check)
        save_json(OUT/'independent_checks.json',{'checks':checks,
            'method':'DOP853 evaluates the full piecewise field at every stage; crossing events are diagnostics. Matched starting state, transient and observation in an additional tighter adaptive GPU run.',
            'reference_correction':'The original v1.0 affine-branch event routine failed at the maximum-alpha example. Full-field evaluation prevents a missed crossing from trapping integration on the wrong branch.',
            'scope':'Selected examples only. Survey and independent checks use different observation lengths. Positive signs do not prove asymptotic chaos.'})
        print(json.dumps(check),flush=True)

if __name__=='__main__':main()
