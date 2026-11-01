"""Independent matrix-exponential and ODE controls for the follow-up results."""
import argparse
import hashlib
import json
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from chua_research_followup import ROOT,OUT,dump,ParameterPrecisionFlights
from continue_chua_slope_families import CoordinateFlights,extrema_check
from check_chua_slope_refinement import independent_segmented
from study_chua_trigger import setup,PARS
from trace_chua_window import Chua


def sheets():
    path=OUT/'sheets/independent-checks.json'
    checks=json.loads(path.read_text())['records'] if path.exists() else []
    for name in ['beta_00','beta_22','alpha_00','alpha_22']:
        source=OUT/'sheets'/(name+'-search.json')
        if not source.exists():continue
        data=json.loads(source.read_text())
        for row in data['records']:
            uid=name+'_'+str(row['peaks'])
            if row['status']!='numerical_root' or any(c['id']==uid for c in checks):continue
            representative=row['representative'];free=data['free_parameter'];p={k:float(v) for k,v in representative['parameters'].items()}
            solver=ParameterPrecisionFlights(row['regions'],row['targets'],p,free,digits=65,spectral=False)
            # Preserve exactly the parameter decimals used in the original solve.
            solver.parameters={k:mp.mpf(v) for k,v in representative['parameters'].items()}
            refined=solver.boundary(representative['coordinate'],representative['q'],'.5')
            delta=abs(mp.mpf(refined['coordinate'])-mp.mpf(representative['coordinate']))
            assert delta<mp.mpf('1e-23')
            fast=CoordinateFlights(p,row['regions'],row['targets'],free=free)
            r=fast.record(float(representative['coordinate']),np.array(representative['q'],float),True)
            geometry=extrema_check(fast,r)
            assert geometry['full_cycle_peaks']==row['peaks']
            ode=independent_segmented(r)
            assert ode['success'] and ode['closure']<2e-6 and ode['multiplier_error']<2e-4,ode
            # A deliberately displaced free parameter cannot satisfy the same
            # orbit equations with the unchanged shooting vector.
            w=mp.matrix([mp.mpf(v) for v in representative['q']]+[mp.mpf(representative['coordinate'])+10*mp.mpf(row['signed_width'])])
            displaced=max(abs(v) for v in solver.system(w,mp.mpf('.5')))
            assert displaced>mp.mpf('1e-12')
            checks.append(dict(id=uid,parameter_difference=mp.nstr(delta,16),refined_residual=refined['residual'],
                working_digits=65,method='mpmath matrix exponential, independent of spectral flight evaluation',
                geometry=geometry,ode=ode,displaced_parameter_residual=mp.nstr(displaced,16),
                representative_sha256=hashlib.sha256(json.dumps(representative,sort_keys=True).encode()).hexdigest(),passed=True))
            dump(path,dict(passed=all(r['passed'] for r in checks),records=checks,
                scope='Independent numerical controls, not interval certificates for these new sheets.'))
            print('independent sheet',uid,'closure',ode['closure'],flush=True)


def trigger():
    p,pc,P,level,e,V,W,seeds=setup();rows=[]
    thresholds=json.loads((OUT/'trigger/thresholds.json').read_text())['records']
    for row in thresholds:
        if row['angle'] not in [-1.2,0.,.8]:continue
        A=row['approximants'][-1]['amplitude'];direction=np.array(row['direction'])
        for side in [-1,1]:
            for method in ['DOP853','Radau']:
                model=Chua(8.4,-.95,-.05,rtol=2e-12,step=.04,method=method)
                sol=model.integrate(12,np.r_[0,(A+side*1e-6)*direction],700,dense=False)
                endnorm=float(np.linalg.norm(sol.y[:,-1]));dist=float(np.linalg.norm(sol.y_events[2][-1,1:]-pc))
                label='origin' if endnorm<1e-7 else 'cycle' if dist<1e-6 else 'unresolved'
                assert label==('origin' if side<0 else 'cycle')
                rows.append(dict(angle=row['angle'],side=side,method=method,label=label,duration=700,
                    final_norm=endnorm,last_section_cycle_distance=dist))
                dump(OUT/'trigger/independent-checks.json',dict(passed=True,records=rows,
                    scope='Separate full-field ODE controls; the continuum and uniqueness claims come from the interval certificate.'))
                print('independent pulse',row['angle'],side,method,label,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['sheets','trigger']);a=ap.parse_args()
    if a.phase=='sheets':sheets()
    else:trigger()
