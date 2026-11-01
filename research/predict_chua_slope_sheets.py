"""Freeze tangency_estimate-anchored layer forecasts before solving their off-slice targets."""
import argparse
import hashlib
import json
import time
import numpy as np
import mpmath as mp
from chua_research_followup import (ROOT,OUT,BASE,dump,classical_templates,
    continue_tangency,ParameterPrecisionFlights)
from continue_chua_slope_families import CoordinateFlights

DIR=OUT/'sheets'


def boundary_seed(seed,p,free,guess):
    """Use a short continuation only to initialize the high-precision root."""
    q=np.array(seed['representative']['q'],float);source={k:float(seed['representative'][k]) for k in BASE}
    value=source[free];failures=[]
    for fraction in np.linspace(0,1,9):
        pars={k:source[k]+fraction*(p[k]-source[k]) for k in BASE}
        fast=CoordinateFlights(pars,seed['regions'],seed['targets'],free=free)
        try:value,q=fast.boundary(value,q,1)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            failures.append(dict(fraction=float(fraction),error=str(exc)))
            # The published forecast remains unchanged; this changes only the
            # starting guess used by the independent target solve.
            return guess,q,failures
    return value,q,failures


def continued_target_seed(peaks,p,free):
    """Implicit parameter prediction initializes shooting, not the frozen forecast."""
    calibration=json.loads((DIR/'central-residual-calibration.json').read_text())[free][str(peaks)]['fold']
    seed=classical_templates()[peaks]
    source={k:float(v) for k,v in calibration['parameters'].items()}
    if all(source[k]==p[k] for k in BASE if k!=free):
        return calibration['coordinate'],calibration['q'],[]
    w=np.r_[np.array(calibration['q'],float),float(calibration['coordinate'])]
    fraction=0.;step=.25;failures=[]
    while fraction<1:
        nxt=min(1.,fraction+step)
        oldp={k:source[k]+fraction*(p[k]-source[k]) for k in BASE}
        newp={k:source[k]+nxt*(p[k]-source[k]) for k in BASE}
        old=CoordinateFlights(oldp,seed['regions'],seed['targets'],free=free)
        new=CoordinateFlights(newp,seed['regions'],seed['targets'],free=free)
        try:
            forcing=np.zeros(len(w))
            for key in BASE:
                if key==free or newp[key]==oldp[key]:continue
                h=1e-7;values=[]
                for sign in [-1,1]:
                    pp=oldp.copy();pp[key]+=sign*h
                    ff=CoordinateFlights(pp,seed['regions'],seed['targets'],free=free)
                    values.append(ff.boundary_system(w,1))
                forcing+=(newp[key]-oldp[key])*(values[1]-values[0])/(2*h)
            trial=w-np.linalg.solve(old.boundary_system(w,1,True),forcing)
            # Double precision is only a preconditioner; stop at its practical
            # cancellation floor and let arbitrary precision establish closure.
            for _ in range(8):
                residual=new.boundary_system(trial,1);error=max(abs(residual))
                if error<2e-4:break
                change=np.linalg.solve(new.boundary_system(trial,1,True),-residual)
                for factor in [1.,.5,.25,.125,.0625]:
                    candidate=trial+factor*change
                    if min(candidate[2:-1])>0 and max(abs(new.boundary_system(candidate,1)))<error:
                        trial=candidate;break
                else:raise ValueError('Seed corrector stalled')
            if max(abs(new.boundary_system(trial,1)))>2e-3:raise ValueError('Seed residual too large')
            w=trial;fraction=nxt;step=min(.5,step*1.5)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            failures.append(dict(fraction=fraction,step=step,error=str(exc)));step/=2
            if step<1e-4:raise ValueError('Target seed continuation failed: '+str(failures[-1]))
    return w[-1],w[:-1],failures


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
        position_model='coordinate_n = estimated_tangency_parameter + A*r^(n-5) + (B+C*(n-5))*r^(2*(n-5))',
        width_model='width_n = width_55*r^(2*(n-5))',
        r_source='Inverse full-return unstable multiplier of the short saddle at the numerical tangency estimate',
        primary_gate='Predicted midpoint lies strictly between independently solved fold and symmetry-breaking endpoints',
        secondary_gate='Relative width error below 10 percent',
        tangency_estimate_scope='Finite-iterate tangent-line approximation, with increased-iteration controls; exact tangency certified only at the original point.',
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
        # A second finite-iterate construction tests sensitivity of the tangency_estimate approximation.
        from chua_research_followup import approximate_tangency
        check=approximate_tangency(org['parameters'],free,
            [org['starting_coordinate']/org['expansion'],org['coordinate']],iterations=10)
        training=[]
        for peaks in [35,45,55]:
            seed=templates[peaks]
            solver=ParameterPrecisionFlights(seed['regions'],seed['targets'],p,free,digits=45)
            displacement=seed['representative']['beta']-BASE['beta']
            # The alpha surface uses the local beta/alpha tangent for a starting guess only.
            guess=org['coordinate']+displacement*(1 if free=='beta' else -1/2.135)
            value,q,failures=boundary_seed(seed,p,free,guess)
            fold=solver.boundary(value,q,1)
            item=dict(peaks=peaks,fold=fold,seed_failures=failures)
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
            tangency_estimate=org,tangency_estimate_control=check,
            tangency_estimate_difference=abs(check['coordinate']-org['coordinate']),training=training,
            coefficients=[mp.nstr(x,38) for x in coeff],forecasts=forecasts,
            elapsed_seconds=time.monotonic()-start))
        print(name,'calibration saved in',time.monotonic()-start,flush=True)


def freeze():
    path=DIR/'forecasts-frozen.json'
    if path.exists():raise SystemExit('Frozen forecasts already exist; refusing overwrite.')
    files=[DIR/(name+'-calibration.json') for name,*_ in nodes()]
    assert all(p.exists() for p in files),'All calibration nodes must finish before freezing'
    corrections=json.loads((DIR/'central-residual-calibration.json').read_text())
    records=[{k:v for k,v in json.loads(p.read_text()).items() if k in
                  ['node','parameters','free_parameter','index','forecasts','tangency_estimate']} for p in files]
    for rec in records:
        for forecast in rec['forecasts']:
            delta=mp.mpf(corrections[rec['free_parameter']][str(forecast['peaks'])]['fold_residual'])
            forecast['corrected_fold']=mp.nstr(mp.mpf(forecast['fold'])+delta,38)
            forecast['corrected_midpoint']=mp.nstr(mp.mpf(forecast['midpoint'])+delta,38)
    dump(path,dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        calibration_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        central_calibration_sha256=hashlib.sha256((DIR/'central-residual-calibration.json').read_bytes()).hexdigest(),
        records=records))


def central_residuals():
    """A separately labelled bias-correction model uses known central target folds.

    These central points are calibration, not blind validation. Off-centre
    endpoint solves have not yet run. The original asymptotic model is retained.
    """
    path=DIR/'central-residual-calibration.json'
    if path.exists():raise SystemExit('Central calibration already exists')
    templates=classical_templates();result={}
    for free in ['beta','alpha']:
        node=json.loads((DIR/(free+'_11-calibration.json')).read_text());rows={}
        for forecast in node['forecasts']:
            seed=templates[forecast['peaks']]
            value,q,failures=boundary_seed(seed,BASE,free,float(forecast['fold']))
            solver=ParameterPrecisionFlights(seed['regions'],seed['targets'],BASE,free,digits=45)
            fold=solver.boundary(value,q,1)
            rows[str(forecast['peaks'])]=dict(fold=fold,
                fold_residual=mp.nstr(mp.mpf(fold['coordinate'])-mp.mpf(forecast['fold']),38),seed_failures=failures)
            print('central calibration',free,forecast['peaks'],fold['coordinate'],flush=True)
        result[free]=rows
    result['scope']='Both central nodes calibrate a constant fold-residual correction. Only the other 16 nodes test prediction at other parameter values. The uncorrected forecasts remain separately scored.'
    dump(path,result)


def search(only=None):
    path=DIR/'forecasts-frozen.json';raw=path.read_bytes();frozen=json.loads(raw)
    templates=classical_templates()
    for node in frozen['records']:
        name=node['node']
        if only and name!=only:continue
        result_path=DIR/(name+'-search.json')
        previous=json.loads(result_path.read_text()) if result_path.exists() else dict(records=[])
        records=[];free=node['free_parameter']
        for forecast in node['forecasts']:
            oldrec=next((r for r in previous['records'] if r['peaks']==forecast['peaks']),None)
            if oldrec and oldrec['status']=='numerical_root':
                records.append(oldrec);continue
            seed=templates[forecast['peaks']]
            solver=ParameterPrecisionFlights(seed['regions'],seed['targets'],node['parameters'],free,digits=45)
            rec=dict(peaks=forecast['peaks'],forecast=forecast,status='unresolved')
            if oldrec:rec['previous_failed_attempt']=oldrec
            try:
                value,q,failures=continued_target_seed(forecast['peaks'],node['parameters'],free)
                rec['initial_solver']='Implicit continuation from the central fold, followed by arbitrary precision'
                rec['seed_rejected_steps']=failures
                fold=solver.boundary(value,q,1)
                symmetry=solver.boundary(fold['coordinate'],fold['q'],-1)
                # A third root is useful for independent integration and phase portraits.
                representative=solver.boundary(fold['coordinate'],fold['q'],.5)
                x=mp.mpf(forecast['midpoint']);a=mp.mpf(fold['coordinate']);b=mp.mpf(symmetry['coordinate'])
                width=a-b
                corrected_midpoint=mp.mpf(forecast['corrected_midpoint'])
                rec.update(status='numerical_root',fold=fold,symmetry=symmetry,representative=representative,
                    signed_width=mp.nstr(width,38),midpoint_inside=bool(min(a,b)<x<max(a,b)),
                    corrected_midpoint_inside=bool(min(a,b)<corrected_midpoint<max(a,b)),
                    corrected_fold_error_in_widths=float((mp.mpf(forecast['corrected_fold'])-a)/width),
                    signed_fold_error_in_widths=float((mp.mpf(forecast['fold'])-a)/width),
                    width_relative_error=float(abs(mp.mpf(forecast['signed_width'])/width-1)),
                    regions=seed['regions'],targets=seed['targets'])
            except (ValueError,RuntimeError,ZeroDivisionError) as exc:rec['error']=str(exc)
            records.append(rec)
            dump(result_path,dict(node=name,parameters=node['parameters'],free_parameter=free,
                forecast_sha256=hashlib.sha256(raw).hexdigest(),records=records,complete=len(records)==3))
            print(name,'target',forecast['peaks'],rec['status'],rec.get('midpoint_inside'),rec.get('signed_fold_error_in_widths'),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['init','calibrate','central','freeze','search']);ap.add_argument('--only')
    a=ap.parse_args()
    if a.phase=='init':initialize()
    elif a.phase=='calibrate':calibrate(a.only)
    elif a.phase=='central':central_residuals()
    elif a.phase=='freeze':freeze()
    else:search(a.only)
