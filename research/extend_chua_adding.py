"""Construct longer period-adding candidates by inserting a parent-cycle return."""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_layer_bifurcations import from_orbit
from chua_fast_flights import FastFlights
from audit_chua_census import independent

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    old=json.loads((ROOT/'results/chua-layer-census/census-analysis.json').read_text())['adding_sequence']
    model=Chua(rtol=2e-12,step=.02);last=old[-1];previous=old[-2]['fold']['beta'];records=[]
    for count in [55,65,75,85]:
        edge=last['fold'];b=edge['beta'];q=np.array(edge['q']);T=edge['period']
        predicted=b+.4*(b-previous)
        parent=last['parent_five_cycle'];pf,pq=from_orbit(model,b,[*parent['q'][:2],parent['period']],False)
        parentfast=FastFlights(model,pf.regions,pf.targets,False);pq=parentfast.solve(b,pq)
        sol=model.integrate(b,[0.,*q[:2]],T/2)
        times=sol.t_events[2];states=sol.y_events[2]
        indices=np.flatnonzero((times>1e-6)&(times<T/2-1e-6))
        i=indices[np.argmin(np.linalg.norm(states[indices,1:]-pq[:2],axis=1))]
        t=times[i];ends=np.cumsum(q[2:]);k=int(np.searchsorted(ends,t));starts=np.r_[0.,ends[:-1]]
        regions=[a['region'] for a in edge['arcs']]
        assert regions[k]==0 and pf.regions[0]==pf.regions[-1]==0
        durations=np.r_[q[2:2+k],t-starts[k]+pq[2],pq[3:-1],pq[-1]+ends[k]-t,q[3+k:]]
        newregions=regions[:k]+pf.regions+regions[k+1:]
        targets=[newregions[j] or newregions[j+1] for j in range(len(newregions)-1)]
        candidate=FastFlights(model,newregions,targets,True);guess=np.r_[q[:2],durations]
        record=dict(maxima_requested=count,predicted_beta=predicted,insertion_time=float(t),status='unresolved')
        try:
            bf,qf=candidate.boundary(predicted,guess,1);fold=candidate.record(bf,qf,True)
            bs,qs=candidate.boundary(bf,qf,-1);sym=candidate.record(bs,qs,True)
            width=bf-bs
            if width<=0 or width>last['symmetric_core_width']*.7:raise RuntimeError('Unexpected adding width')
            bc=bs+.45*width;qc=candidate.solve(bc,qs);core=candidate.record(bc,qc,True)
            if max(r['region_violation'] for r in [fold,sym,core])>1e-7:raise RuntimeError('Invalid itinerary')
            orbit=model.integrate(bc,[0.,*qc[:2]],core['period']);peak_count=len(orbit.t_events[3])
            var=independent(core);radau=Chua(rtol=2e-13,step=.01,method='Radau').integrate(bc,[0.,*qc[:2]],core['period'])
            closure=float(max(abs(radau.y[:,-1]-np.r_[0.,qc[:2]])))
            record.update(fold=fold,symmetry=sym,core=core,regions=newregions,targets=targets,
                peak_count=peak_count,symmetric_core_width=width,independent_variational=var,radau_closure=closure)
            if peak_count!=count or core['spectral_radius']>=.999 or not var['success'] or abs(var['dominant'])>=.999 or var['multiplier_error']>2e-5 or var['closure']>2e-7 or closure>2e-7:
                raise RuntimeError('Independent orbit check failed')
            record['status']='verified_attracting_cycle'
            previous=b;last=dict(fold=fold,symmetric_core_width=width,parent_five_cycle=parent)
        except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:record['error']=str(e)
        records.append(record);(out/'extended-adding.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(count,record['status'],record.get('fold',{}).get('beta'),record.get('symmetric_core_width'),record.get('radau_closure'),record.get('error'),flush=True)
        if record['status']!='verified_attracting_cycle':break


if __name__=='__main__':main()
