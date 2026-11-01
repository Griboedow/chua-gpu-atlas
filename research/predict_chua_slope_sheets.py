"""Freeze organizer-anchored layer forecasts before solving their off-slice targets."""
import argparse
import hashlib
import json
import time
import numpy as np
import mpmath as mp
from chua_research_followup import (ROOT,OUT,BASE,dump,classical_templates,
    continue_tangency,ParameterPrecisionFlights)

DIR=OUT/'sheets'


def nodes():
    for free in ['beta','alpha']:
        for i in [-1,0,1]:
            for j in [-1,0,1]:
                p=BASE.copy();p['m0']+=i*.0004;p['m1']+=j*.0004
                yield f'{free}_{i+1}{j+1}',free,i,j,p


def initialize():
    path=DIR/'protocol.json'
    spec=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        base_parameters=BASE,slope_offsets=[-.0004,0.,.0004],
        sections=['alpha=9; predict beta','beta=14.30635909996146; predict alpha'],
        training_peaks=[35,45,55],target_peaks=[65,75,85],
        position_model='coordinate_n = organizer + A*r^(n-5) + (B+C*(n-5))*r^(2*(n-5))',
        width_model='width_n = width_55*r^(2*(n-5))',
        r_source='Inverse full-return unstable multiplier of the short saddle at the numerical organizer',
        primary_gate='Predicted midpoint lies strictly between independently solved fold and symmetry-breaking endpoints',
        secondary_gate='Relative width error below 10 percent',
        organizer_scope='Finite-iterate tangent-line approximation, with increased-iteration controls; exact tangency certified only at the original point.',
        provenance='Existing central target roots are used only as initial guesses during the later search; target endpoints do not enter forecasts. Central beta node is retrospective.',
        failure_policy='Retain forecast misses and unresolved roots without changing the frozen forecasts.')
    if path.exists():return
    dump(path,spec)


def calibrate(only=None):
    templates=classical_templates();initialize()
    for name,free,i,j,p in nodes():
        if only and name!=only:continue
        path=DIR/(name+'-calibration.json')
        if path.exists():continue
        start=time.monotonic();org=continue_tangency(p,free)
        # A second finite-iterate construction tests sensitivity of the organizer approximation.
        from chua_research_followup import approximate_tangency
        check=approximate_tangency(org['parameters'],free,
            [org['starting_coordinate']/org['expansion'],org['coordinate']],iterations=10)
        training=[]
        for peaks in [35,45,55]:
            seed=templates[peaks]
            solver=ParameterPrecisionFlights(seed['regions'],seed['targets'],p,free,digits=45)
            displacement=seed['core']['beta']-BASE['beta']
            # The alpha surface uses the local beta/alpha tangent for a starting guess only.
            guess=org['coordinate']+displacement*(1 if free=='beta' else -1/2.135)
            fold=solver.boundary(guess,seed['core']['q'],1)
            item=dict(peaks=peaks,fold=fold)
            if peaks==55:
                symmetry=solver.boundary(fold['coordinate'],fold['q'],-1)
                item['symmetry']=symmetry
                item['signed_width']=mp.nstr(mp.mpf(fold['coordinate'])-mp.mpf(symmetry['coordinate']),38)
            training.append(item)
            print(name,'training',peaks,fold['coordinate'],flush=True)
        r=mp.mpf(str(1/org['expansion']));limit=mp.mpf(str(org['coordinate']))
        ns=[(x['peaks']-5)//10-5 for x in training]
        design=mp.matrix([[r**n,r**(2*n),n*r**(2*n)] for n in ns])
        coeff=mp.lu_solve(design,mp.matrix([mp.mpf(x['fold']['coordinate'])-limit for x in training]))
        width55=mp.mpf(training[-1]['signed_width']);forecasts=[]
        for peaks in [65,75,85]:
            k=(peaks-5)//10-5
            fold=limit+coeff[0]*r**k+(coeff[1]+coeff[2]*k)*r**(2*k)
            width=width55*r**(2*k)
            forecasts.append(dict(peaks=peaks,fold=mp.nstr(fold,38),signed_width=mp.nstr(width,38),
                midpoint=mp.nstr(fold-width/2,38)))
        dump(path,dict(node=name,free_parameter=free,index=[i,j],parameters=p,
            organizer=org,organizer_control=check,
            organizer_difference=abs(check['coordinate']-org['coordinate']),training=training,
            coefficients=[mp.nstr(x,38) for x in coeff],forecasts=forecasts,
            elapsed_seconds=time.monotonic()-start))
        print(name,'calibration saved in',time.monotonic()-start,flush=True)


def freeze():
    path=DIR/'forecasts-frozen.json'
    if path.exists():raise SystemExit('Frozen forecasts already exist; refusing overwrite.')
    files=[DIR/(name+'-calibration.json') for name,*_ in nodes()]
    assert all(p.exists() for p in files),'All calibration nodes must finish before freezing'
    dump(path,dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        calibration_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        records=[{k:v for k,v in json.loads(p.read_text()).items() if k in
                  ['node','parameters','free_parameter','index','forecasts','organizer']} for p in files]))


def search(only=None):
    path=DIR/'forecasts-frozen.json';raw=path.read_bytes();frozen=json.loads(raw)
    templates=classical_templates()
    for node in frozen['records']:
        name=node['node']
        if only and name!=only:continue
        result_path=DIR/(name+'-search.json')
        if result_path.exists():continue
        records=[];free=node['free_parameter']
        for forecast in node['forecasts']:
            seed=templates[forecast['peaks']]
            solver=ParameterPrecisionFlights(seed['regions'],seed['targets'],node['parameters'],free,digits=45)
            rec=dict(peaks=forecast['peaks'],forecast=forecast,status='unresolved')
            try:
                fold=solver.boundary(forecast['fold'],seed['core']['q'],1)
                symmetry=solver.boundary(fold['coordinate'],fold['q'],-1)
                # A third root is useful for independent integration and phase portraits.
                core=solver.boundary(fold['coordinate'],fold['q'],.5)
                x=mp.mpf(forecast['midpoint']);a=mp.mpf(fold['coordinate']);b=mp.mpf(symmetry['coordinate'])
                width=a-b
                rec.update(status='numerical_root',fold=fold,symmetry=symmetry,core=core,
                    signed_width=mp.nstr(width,38),midpoint_inside=bool(min(a,b)<x<max(a,b)),
                    signed_fold_error_in_widths=float((mp.mpf(forecast['fold'])-a)/width),
                    width_relative_error=float(abs(mp.mpf(forecast['signed_width'])/width-1)),
                    regions=seed['regions'],targets=seed['targets'])
            except (ValueError,RuntimeError,ZeroDivisionError) as exc:rec['error']=str(exc)
            records.append(rec)
            dump(result_path,dict(node=name,parameters=node['parameters'],free_parameter=free,
                forecast_sha256=hashlib.sha256(raw).hexdigest(),records=records,complete=len(records)==3))
            print(name,'target',forecast['peaks'],rec['status'],rec.get('midpoint_inside'),rec.get('signed_fold_error_in_widths'),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['init','calibrate','freeze','search']);ap.add_argument('--only')
    a=ap.parse_args()
    if a.phase=='init':initialize()
    elif a.phase=='calibrate':calibrate(a.only)
    elif a.phase=='freeze':freeze()
    else:search(a.only)
