"""Independent variational integration, local sensitivities, and audit gates."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from trace_chua_window import Chua
from chua_layer_bifurcations import Flights


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results/chua-layers')
    out=p.parse_args().output
    central=json.loads((out/'affine-central.json').read_text())
    extended=json.loads((out/'extended.json').read_text())
    cascade=json.loads((out/'cascade.json').read_text())
    verification=json.loads((out/'verification.json').read_text())
    model=Chua()
    picks=[('fold',central['boundaries']['1']),('symmetry',central['boundaries']['-1']),
           ('period_doubling',extended['period_doubling']),('second_pd',cascade['second_period_doubling'])]
    independent=[]
    for name,r in picks:
        beta=r['beta'];u=np.array([0.,*r['q'][:2]])
        def rhs(t,state):
            f=model.field(beta,state[:3]);J=model.jacobian(beta,abs(state[0])<1)
            return np.r_[f,(J@state[3:].reshape(3,3)).ravel()]
        initial=np.r_[u,np.eye(3).ravel()]
        sol=solve_ivp(rhs,[0,r['period']],initial,method='DOP853',rtol=2e-12,atol=2e-13,max_step=.02)
        end=sol.y[:3,-1];M=sol.y[3:,-1].reshape(3,3);f=model.field(beta,end)
        P=((np.eye(3)-np.outer(f,[1.,0.,0.])/f[0])@M)[1:,1:]
        eig=np.linalg.eigvals(P);dominant=eig[np.argmax(abs(eig))]
        expected=r['full_multipliers'][0][0]
        rec=dict(name=name,success=sol.success,closure=float(max(abs(end-u))),
                 dominant_multiplier=float(dominant.real),expected_multiplier=expected,
                 multiplier_error=float(abs(dominant-expected)),rtol=2e-12,atol=2e-13,max_step=.02)
        independent.append(rec);print(rec,flush=True)
    gradients={}
    for name,base,target in [('fold',central['boundaries']['1'],1),('symmetry',central['boundaries']['-1'],-1),
                             ('period_doubling',extended['period_doubling'],-1)]:
        symmetric=name!='period_doubling'
        regions=central['regions'] if symmetric else extended['full_regions']
        targets=central['targets'] if symmetric else extended['full_targets']
        components={}
        for parameter,h in [('alpha',1e-4),('m0',1e-5)]:
            values=[]
            for sign in [-1,1]:
                kwargs={parameter:(9 if parameter=='alpha' else -8/7)+sign*h}
                family=Flights(Chua(**kwargs),regions,targets,symmetric)
                b,q=family.boundary(base['beta'],np.array(base['q']),target)
                values.append(b)
            components[parameter]=dict(step=h,values=values,derivative=(values[1]-values[0])/(2*h))
        gradients[name]=components
    family=Flights(model,central['regions'],central['targets'],True)
    r=min(central['central'],key=lambda r:abs(r['beta']-14.32));q=np.array(r['q']);beta=r['beta']
    analytic=family.evaluate(beta,q)[1];h=1e-6
    numeric=np.column_stack([(family.evaluate(beta,q+h*e)[0]-family.evaluate(beta,q-h*e)[0])/(2*h) for e in np.eye(len(q))])
    jacobian_error=float(np.max(abs(analytic-numeric)))
    controls=[c for row in verification['checks'] if 'independent_integrations' in row for c in row['independent_integrations']]
    surfaces=[row[key] for row in extended['surfaces'] for key in ['fold','symmetry','period_doubling']]
    gates=dict(variational_solvers_succeeded=all(r['success'] for r in independent),
               independent_multiplier_agreement=max(r['multiplier_error'] for r in independent)<2e-5,
               independent_variational_closure=max(r['closure'] for r in independent)<1e-7,
               full_field_closure=max(r['closure'] for r in controls)<1e-6,
               radau_closure=max(r['closure'] for r in controls if r['method']=='Radau')<1e-8,
               affine_jacobian=jacobian_error<1e-5,
               surface_closure=max(r['residual'] for r in surfaces)<1e-7,
               sampled_region_validity=max(r['region_violation'] for r in surfaces)<1e-7,
               transverse_breakpoint_crossings=min(r['minimum_breakpoint_speed'] for r in surfaces)>.5,
               fold_two_distinct_roots=abs(extended['fold_partners'][0]['q'][0]-extended['fold_partners'][1]['q'][0])>1e-4,
               asymmetric_symmetry_partner=verification['checks'][-1]['section_separation']>1e-3)
    summary=dict(status='passed' if all(gates.values()) else 'failed',gates=gates,
                 independent_variational=independent,local_surface_gradients=gradients,
                 jacobian_max_error=jacobian_error,surface_points=len(surfaces),
                 max_surface_residual=max(r['residual'] for r in surfaces),
                 max_full_field_closure=max(r['closure'] for r in controls),
                 max_radau_closure=max(r['closure'] for r in controls if r['method']=='Radau'),
                 scope='Numerical continuation and independent controls, not validated interval bounds or an exhaustive atlas.',
                 code_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [
                     *Path(__file__).parent.glob('*chua*layer*.py'),Path(__file__).with_name('trace_chua_window.py')]})
    (out/'audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)
    if not all(gates.values()):raise SystemExit(2)


if __name__=='__main__':main()
