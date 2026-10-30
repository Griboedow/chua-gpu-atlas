"""Close off-center candidates and transport their orbit identity to the reference plane."""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_layer_bifurcations import from_orbit
from chua_fast_flights import FastFlights
from audit_chua_census import independent
from refine_chua_satellites import advance

ROOT=Path(__file__).resolve().parents[1]


def verify(m0,beta,start):
    model=Chua(m0=m0,rtol=2e-12,step=.025)
    warm=model.integrate(beta,start,120.,dense=False);late=model.integrate(beta,warm.y[:,-1],500.)
    times,states=late.t_events[2],late.y_events[2];record=dict(alpha=9.,m0=m0,beta=beta,status='no_recurrence')
    for lag in range(1,min(60,len(times)//3)):
        if np.max(abs(states[-lag:]-states[-2*lag:-lag]))>2e-6:continue
        v=np.r_[states[-1,1:],times[-1]-times[-1-lag]];family,q=from_orbit(model,beta,v,False)
        r=family.record(beta,q,True);u=np.r_[0.,q[:2]];T=r['period'];orbit=model.integrate(beta,u,T)
        var=independent(r);radau=Chua(m0=m0,rtol=2e-12,step=.02,method='Radau').integrate(beta,u,T)
        error=float(max(abs(radau.y[:,-1]-u)));half=float(np.max(abs(orbit.sol(np.linspace(0,T/2,1001))+orbit.sol(np.linspace(T/2,T,1001)))))
        record.update(status='unresolved',cycle=r,peak_count=len(orbit.t_events[3]),section_period=lag,
                      half_symmetry_error=half,radau_closure=error,variational=var,regions=family.regions,targets=family.targets)
        if error>2e-7 or not var['success'] or abs(var['dominant'])>=.9999 or var['multiplier_error']>2e-5 or r['region_violation']>1e-7:return record
        record['status']='verified_attracting_cycle'
        try:
            symmetric=half<2e-5
            f,qq=from_orbit(model,beta,[*q[:2],T],symmetric);f=FastFlights(model,f.regions,f.targets,symmetric)
            target=.5 if symmetric else -.5
            b,qq=f.boundary(beta,qq,target);local=f.record(b,qq,True)
            if abs(local['period']/T-1)>.1 or local['region_violation']>1e-7:raise RuntimeError('Canonical local root changed itinerary')
            template=dict(regions=f.regions,targets=f.targets,symmetric=symmetric,target=target)
            row=dict(alpha=9.,m0=m0,core=local)
            for m in np.linspace(m0,-8/7,25)[1:]:row=advance(template,row,9.,float(m))
            record['canonical']=dict(core=row['core'],symmetric=symmetric,regions=f.regions,targets=f.targets,target=target)
        except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:record['canonical_error']=str(e)
        return record
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--resume',action='store_true');args=p.parse_args();out=args.output
    candidates=json.loads((out/'offcenter-candidates.json').read_text());raw=np.load(out/'offcenter-raw.npz');path=out/'offcenter-cycles.json'
    records=json.loads(path.read_text()) if args.resume and path.exists() else []
    for k,c in enumerate(candidates):
        if k<len(records):continue
        i=int(np.argmin(abs(raw['beta']-c['beta'])));beta=float(raw['beta'][i]);index=c['slice']*len(raw['beta'])+i
        try:r=verify(c['m0'],beta,raw['final_states'][index,0])
        except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:r=dict(status='failed',error=str(e),beta=beta,m0=c['m0'])
        r['candidate_index']=k;r['gpu_period_hint']=c['period_hint'];records.append(r)
        path.write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(k+1,'of',len(candidates),r['status'],beta,r.get('peak_count'),'canonical',r.get('canonical',{}).get('core',{}).get('beta'),r.get('canonical_error'),flush=True)


if __name__=='__main__':main()
