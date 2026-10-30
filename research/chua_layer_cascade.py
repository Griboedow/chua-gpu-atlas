"""Follow the first doubled orbit and validate independent full-field returns."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from trace_chua_window import Chua, initialize
from chua_layer_bifurcations import Flights, from_orbit

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results/chua-layers')
    OUT=parser.parse_args().output
    model=Chua(rtol=2e-12,step=.025)
    v=initialize(model,14.3097)
    family,q0=from_orbit(model,14.3097,v,False)
    result=dict(regions=family.regions,targets=family.targets,records=[])
    for beta in np.arange(14.3097,14.3079,-.0001):
        q0=family.solve(float(beta),q0)
        r=family.record(float(beta),q0)
        result['records'].append(r)
        print('doubled',beta,r['period'],r['full_multipliers'],flush=True)
        if r['full_multipliers'][0][0]<-1:
            b,q=family.boundary(beta,q0,-1)
            result['second_period_doubling']=family.record(b,q,True)
            print('PD2',b,flush=True);break
    (OUT/'cascade.json').write_text(json.dumps(result,indent=2)+'\n')
    base=json.loads((OUT/'affine-central.json').read_text())
    extended=json.loads((OUT/'extended.json').read_text())
    picks=[('symmetry_boundary',base['boundaries']['-1']),('fold',base['boundaries']['1']),
           ('asymmetric',extended['asymmetric_initial']),('first_pd',extended['period_doubling']),
           ('fold_branch_a',extended['fold_partners'][0]),('fold_branch_b',extended['fold_partners'][1]),
           ('doubled',result['records'][0]),('second_pd',result['second_period_doubling'])]
    checks=[];arrays={}
    for name,r in picks:
        u=[0.,*r['q'][:2]];T=r['period'];b=r['beta']
        one=[]
        for method,tol,step in [('DOP853',2e-12,.025),('DOP853',2e-13,.0125),('Radau',2e-12,.025)]:
            m=Chua(rtol=tol,step=step,method=method)
            sol=m.integrate(b,u,T)
            one.append(dict(method=method,rtol=tol,max_step=step,closure=float(np.max(np.abs(sol.y[:,-1]-u))),
                            breakpoint_crossings=sum(len(e) for e in sol.t_events[:2]),maxima_count=len(sol.t_events[3])))
            if method=='Radau':
                arrays[name+'_t']=np.linspace(0,T,4001);arrays[name+'_xyz']=sol.sol(arrays[name+'_t']).T
        checks.append(dict(name=name,beta=b,period=T,independent_integrations=one))
        print('verify',name,one,flush=True)
    # The opposite orbit really is distinct: shift the inverted full trajectory
    # to its positive x=0 section, then solve closure without symmetry imposed.
    r=extended['asymmetric_initial'];u=[0.,*r['q'][:2]]
    sol=model.integrate(r['beta'],u,r['period'])
    down=brentq(lambda t:sol.sol(t)[0],4.5,7.)
    mate=-sol.sol(down)
    partner,q=from_orbit(model,r['beta'],[*mate[1:],r['period']],False)
    twin=partner.record(r['beta'],q,True)
    checks.append(dict(name='symmetry_partner',record=twin,
                       section_separation=float(np.linalg.norm(q[:2]-np.array(r['q'][:2])))))
    (OUT/'verification.json').write_text(json.dumps(dict(checks=checks),indent=2)+'\n')
    np.savez_compressed(OUT/'verified-orbits.npz',**arrays)


if __name__=='__main__':main()
