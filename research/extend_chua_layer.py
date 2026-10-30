"""Trace symmetry-broken cycles, the fold partner, and local layer surfaces."""
import argparse
import json
from pathlib import Path

import numpy as np

from trace_chua_window import Chua, initialize
from chua_layer_bifurcations import Flights, from_orbit


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results/chua-layers')
    OUT=parser.parse_args().output
    model=Chua(rtol=2e-12,step=.025)
    data=json.loads((OUT/'affine-central.json').read_text())
    half=Flights(model,data['regions'],data['targets'],True)
    fold=data['boundaries']['1']; pitch=data['boundaries']['-1']
    qfold=np.array(fold['q']); bfold=fold['beta']
    U,S,V=np.linalg.svd(half.evaluate(bfold,qfold)[1])
    left,right=U[:,-1],V[-1]
    fold_diagnostics=[]
    for h in [1e-4,5e-5,2e-5]:
        r=half.evaluate(bfold,qfold)[0]
        rb=(half.evaluate(bfold+h,qfold)[0]-half.evaluate(bfold-h,qfold)[0])/(2*h)
        rqq=(half.evaluate(bfold,qfold+h*right)[0]-2*r+half.evaluate(bfold,qfold-h*right)[0])/h**2
        fold_diagnostics.append(dict(step=h,left_parameter_derivative=float(left@rb),left_quadratic_derivative=float(left@rqq)))
    partners=[]
    for sign in [-1,1]:
        q=half.solve(bfold-.0001,qfold+sign*.01*right)
        partners.append(half.record(bfold-.0001,q,True))
    result=dict(fold_diagnostics=fold_diagnostics,fold_partners=partners,asymmetric=[],surfaces=[])
    v=initialize(model,14.315)
    full,q0=from_orbit(model,14.315,v,False)
    result['full_regions']=full.regions;result['full_targets']=full.targets
    result['asymmetric_initial']=full.record(14.315,q0,True)
    print('asymmetric initial',json.dumps(result['asymmetric_initial']),flush=True)
    for betas in [np.arange(14.315,14.317251,.00025),np.arange(14.31475,14.2999,-.00025)]:
        q=q0.copy()
        for beta in betas:
            try:q=full.solve(float(beta),q)
            except RuntimeError as e:print(e,flush=True);break
            r=full.record(float(beta),q)
            result['asymmetric'].append(r)
            print('asymmetric',beta,r['period'],r['full_multipliers'],flush=True)
            if r['full_multipliers'][0][0]<-2:break
    for r in result['asymmetric']:
        if r['full_multipliers'][0][0]<-1:
            b,q=full.boundary(r['beta'],np.array(r['q']),-1)
            result['period_doubling']=full.record(b,q,True)
            print('PD',b,flush=True);break
    for alpha in np.linspace(8.94,9.06,7):
        for m0 in np.linspace(-8/7-.006,-8/7+.006,7):
            local=Chua(alpha=float(alpha),m0=float(m0))
            family=Flights(local,half.regions,half.targets,True)
            row=dict(alpha=float(alpha),m0=float(m0))
            for key,target,base in [('fold',1,fold),('symmetry',-1,pitch)]:
                b,q=family.boundary(base['beta'],np.array(base['q']),target)
                r=family.record(b,q,True)
                if r['region_violation']>1e-7:raise RuntimeError('Invalid surface itinerary')
                row[key]=r
            base=result['period_doubling']
            b,q=base['beta'],np.array(base['q'])
            # Track the asymmetric branch in small parameter increments;
            # direct corner solves can jump toward the symmetric branch.
            for weight in np.linspace(0,1,25)[1:]:
                intermediate=Chua(alpha=float(9+weight*(alpha-9)),m0=float(-8/7+weight*(m0+8/7)))
                doubled_family=Flights(intermediate,full.regions,full.targets,False)
                b,q=doubled_family.boundary(b+(row['fold']['beta']-fold['beta'])/24,q,-1)
            r=doubled_family.record(b,q,True)
            if r['region_violation']>1e-7:raise RuntimeError('Invalid period-doubling surface itinerary')
            row['period_doubling']=r
            result['surfaces'].append(row)
        print('surface alpha',alpha,flush=True)
    (OUT/'extended.json').write_text(json.dumps(result,indent=2)+'\n')
    print('saved extended',flush=True)


if __name__=='__main__':main()
