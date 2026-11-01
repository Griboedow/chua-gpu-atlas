"""Recompute frozen forecast scores and independently check selected cycles."""
import hashlib
import json
import subprocess
import numpy as np
import mpmath as mp
from chua_research_followup import ROOT,OUT,dump,ParameterPrecisionFlights
from continue_chua_slope_families import CoordinateFlights,extrema_check
from check_chua_slope_refinement import independent_segmented

DIR=OUT/'extended-slopes'


def main():
    path=DIR/'forecasts-frozen.json';raw=path.read_bytes()
    published=subprocess.check_output(['git','show','0c39117:results/chua-followup/extended-slopes/forecasts-frozen.json'],cwd=ROOT)
    assert json.loads(published)==json.loads(raw)
    frozen=json.loads(raw);scores=[];checks=[]
    checkpath=DIR/'independent-checks.json'
    if checkpath.exists():checks=json.loads(checkpath.read_text())['records']
    for node in frozen['records']:
        data=json.loads((DIR/(node['node']+'-search.json')).read_text())
        assert data['forecast_sha256']==hashlib.sha256(raw).hexdigest()
        for forecast,row in zip(node['forecasts'],data['records']):
            assert row['forecast']==forecast
            uid=node['node']+'_'+str(row['peaks'])
            if row['status']!='numerical_root':
                scores.append(dict(id=uid,status=row['status']));continue
            mp.mp.dps=60
            a,b=mp.mpf(row['fold']['coordinate']),mp.mpf(row['symmetry']['coordinate']);w=a-b
            mid=mp.mpf(forecast['midpoint']);corr=mp.mpf(forecast['corrected_midpoint'])
            score=dict(id=uid,status='numerical_root',raw_inside=bool(min(a,b)<mid<max(a,b)),
                corrected_inside=bool(min(a,b)<corr<max(a,b)),
                corrected_fold_error=float((mp.mpf(forecast['corrected_fold'])-a)/w),
                width_relative_error=float(abs(mp.mpf(forecast['signed_width'])/w-1)))
            scores.append(score)
            representative=row['representative'];free=node['free_parameter']
            p={k:float(v) for k,v in representative['parameters'].items()}
            fast=CoordinateFlights(p,row['regions'],row['targets'],free=free)
            r=fast.record(float(representative['coordinate']),np.array(representative['q'],float),True)
            geometry=extrema_check(fast,r)
            assert geometry['full_cycle_peaks']==row['peaks']
            assert r['region_violation']<1e-7
            score['geometry']=geometry
            if node['node'] not in ['beta_00','beta_22','alpha_00','alpha_22'] or any(c['id']==uid for c in checks):continue
            solver=ParameterPrecisionFlights(row['regions'],row['targets'],p,free,digits=65,spectral=False)
            solver.parameters={k:mp.mpf(v) for k,v in representative['parameters'].items()}
            refined=solver.boundary(representative['coordinate'],representative['q'],'.5')
            delta=abs(mp.mpf(refined['coordinate'])-mp.mpf(representative['coordinate']))
            assert delta<mp.mpf('1e-23')
            ode=independent_segmented(r)
            assert ode['success'] and ode['closure']<2e-6 and ode['multiplier_error']<2e-4,ode
            checks.append(dict(id=uid,parameter_difference=mp.nstr(delta,16),geometry=geometry,ode=ode,
                working_digits=65,method='Nonspectral mpmath matrix exponential and segmented ODE integration',passed=True))
            dump(checkpath,dict(records=checks,passed=True))
            print(uid,'independently checked',flush=True)
    dump(DIR/'scores.json',dict(forecast_commit='0c39117',records=scores,
        total=len(scores),raw_hits=sum(r.get('raw_inside',False) for r in scores),
        corrected_hits=sum(r.get('corrected_inside',False) for r in scores),
        max_width_relative_error=max(r.get('width_relative_error',0) for r in scores),
        max_corrected_fold_error=max(abs(r.get('corrected_fold_error',0)) for r in scores),
        scope='A prospective finite test at eight new slope nodes, not a uniform error bound.'))


if __name__=='__main__':main()
