"""Compare frozen forecasts to measured roots and describe their organizing motion."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.optimize import brentq
from chua_fast_flights import FastFlights
from chua_precision_flights import PrecisionFlights
from trace_chua_window import Chua
from predict_chua_layers import ROOT,OUT,old_layers,parent_at


def section_word(core,regions,targets,high=False):
    beta=float(core['beta']);q=np.array(core['q'],float)
    if high:
        f=PrecisionFlights(regions,targets,True,55)
        ev=f.evaluate(mp.mpf(core['beta']),[mp.mpf(x) for x in core['q']],False)
        arcs=[dict(region=a['region'],duration=float(a['duration']),start=np.array(list(a['start']),float)) for a in ev['arcs']]
    else:
        f=FastFlights(Chua(),regions,targets,True)
        arcs=f.record(beta,q,True)['arcs']
    times=[0.];points=[q[:2].tolist()];t0=0.
    total=float(core['period'])/2
    model=Chua()
    for arc in arcs:
        if arc['region']==0:
            A=model.jacobian(beta,True);eig,V=np.linalg.eig(A);coef=np.linalg.solve(V,arc['start'])
            def state(t):return np.real(V@(np.exp(eig*t)*coef))
            def x(t):return state(t)[0]
            grid=np.linspace(1e-10,arc['duration']-1e-10,101)
            xx=np.array([x(t) for t in grid])
            for l,r,vl,vr in zip(grid[:-1],grid[1:],xx[:-1],xx[1:]):
                if vl*vr>=0:continue
                t=brentq(x,l,r,xtol=1e-14)
                if not 1e-7<t0+t<total-1e-7:continue
                u=state(t);sgn=np.sign(model.field(beta,u)[0])
                times.append(float(t0+t));points.append((sgn*u[1:]).tolist())
        t0+=arc['duration']
    parent=parent_at(beta);delta=np.array(points)-np.array(parent['q'][:2])
    vals,V=np.linalg.eig(np.array(parent['section_matrix']));order=np.argsort(-abs(vals));V=V[:,order];vals=vals[order]
    coordinates=np.linalg.solve(V,delta.T).T.real
    ratio_checks=[]
    for i in range(len(points)):
        j=(i+2)%len(points)
        if np.linalg.norm(delta[i])<.002 and np.linalg.norm(delta[j])<.002 and abs(coordinates[i,0])>1e-7:
            ratio_checks.append(dict(index=i,ratio=float(coordinates[j,0]/coordinates[i,0]),
                predicted=parent['full_multipliers'][0][0]))
    return dict(times=times,normalized_section_points=points,
        distances_to_parent=np.linalg.norm(delta,axis=1).tolist(),
        unstable_coordinate=coordinates[:,0].tolist(),stable_coordinate=coordinates[:,1].tolist(),
        parent_period=parent['period'],parent_full_expansion=parent['full_multipliers'][0][0],
        parent_half_expansion=float(vals[0].real),local_expansion_checks=ratio_checks)


def main():
    forecast=json.loads((OUT/'predictions-frozen.json').read_text())
    search=json.loads((OUT/'search-results.json').read_text())
    audit=json.loads((OUT/'independent-checks.json').read_text())
    assert search['forecast_sha256']==hashlib.sha256((OUT/'predictions-frozen.json').read_bytes()).hexdigest()
    assert all(r['passed'] for r in audit['records'])
    registry=json.loads((ROOT/'results/chua-layer-atlas/family-registry.json').read_text())
    families=registry['families'];rows=[]
    ids=['S04','S05','S06','S07','P55','P65']
    for old,id in zip(old_layers(),ids):
        family=next(r for r in families if r['id']==id);core=family['core']
        rows.append(dict(n=old['n'],id=id,peaks=old['peaks'],fold_beta=old['fold']['beta'],
            width=old['width'],fold_period=old['fold']['period'],stage='known_before_forecast',
            quadratic_core_coordinate=(old['fold']['beta']-core['beta'])/old['width'],
            motion=section_word(core,family['regions'],family['targets'])))
    comparison=[];extended=list(families)
    for target,rec in zip(forecast['predictions'],search['records']):
        mp.mp.dps=55;bf=mp.mpf(rec['fold']['beta']);width=mp.mpf(rec['symmetric_core_width']);bc=mp.mpf(rec['core']['beta'])
        pred_fold=mp.mpf(str(target['fold_beta_corrected']));pred_width=mp.mpf(str(target['symmetric_core_width']))
        offset=(pred_fold-bf)/width;pred_left=offset-pred_width/width
        overlap=max(mp.mpf(0),min(mp.mpf(0),offset)-max(mp.mpf(-1),pred_left))
        actual_period=float(rec['fold']['period']);peaks=rec['peaks_requested'];id='P'+str(peaks)
        rows.append(dict(n=target['n'],id=id,peaks=peaks,fold_beta=float(bf),width=float(width),
            fold_period=actual_period,stage='new_target',quadratic_core_coordinate=float((bf-bc)/width),
            motion=section_word(rec['core'],rec['regions'],rec['targets'],True)))
        comparison.append(dict(peaks=peaks,predicted_fold_beta=target['fold_beta_corrected'],measured_fold_beta=str(bf),
            predicted_width=float(pred_width),measured_width=str(width),
            width_error_relative_to_prediction=float(abs(width/pred_width-1)),
            signed_fold_error_in_measured_widths=float(offset),
            predicted_interval_in_measured_width_units=[float(pred_left),float(offset)],
            interval_overlap_fraction_of_measured_core=float(overlap),
            predicted_midpoint_inside_measured_core=bool(-1<offset-pred_width/width/2<0),
            predicted_period=target['fold_period'],measured_period=actual_period,
            period_absolute_error=abs(target['fold_period']-actual_period),
            frozen_engineering_gates_passed=bool(abs(pred_fold-bf)<target['acceptance']['corrected_location_absolute_error'] and abs(width/pred_width-1)<target['acceptance']['width_relative_error'] and abs(target['fold_period']-actual_period)<target['acceptance']['period_absolute_error'])))
        f=FastFlights(Chua(),rec['regions'],rec['targets'],True)
        core=f.record(float(rec['core']['beta']),np.array(rec['core']['q'],float))
        extended.append(dict(id=id,symmetric=True,core=core,regions=rec['regions'],targets=rec['targets'],target=.5,
            precision_core=rec['core'],precision_evidence='search-results.json and independent-checks.json; core retains float64 compatibility only',
            peak_count=peaks,section_returns=target['section_returns'],
            sources=[dict(id='frozen_prediction_'+str(peaks),chain=[])],source_count=1,
            coverage='One independently checked reference root in the original box; no 3D sheet grid in this update.'))
    for i,row in enumerate(rows):
        if i:
            row['width_ratio']=row['width']/rows[i-1]['width']
            row['period_increment']=row['fold_period']-rows[i-1]['fold_period']
        if i>1:
            row['spacing_ratio']=(row['fold_beta']-rows[i-1]['fold_beta'])/(rows[i-1]['fold_beta']-rows[i-2]['fold_beta'])
    r=forecast['calibration']['ratio']
    result=dict(rows=rows,comparison=comparison,frozen_commit='0fc47ea',
        expected_spacing_ratio=r,expected_width_ratio=r*r,expected_period_increment=2*forecast['calibration']['parent']['period'],
        resolvable_members_per_decimal_digit=float(np.log(10)/(-2*np.log(r))),
        count_statement=dict(original_reference_catalogue=32,new_reference_roots=2,extended_reference_count=34,
            proved_infinite=False,global_number_unknown=True,
            scope='34 distinct roots under the existing normalization convention, not a count of globally disconnected components or simultaneous attractors.'),
        validation_scope='High precision is not interval certification; the global tangency and a uniform large-n remainder bound remain unproved.')
    (OUT/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'extended-family-registry.json').write_text(json.dumps(dict(family_count=34,families=extended,
        identity_rule=registry['identity_rule'],scope=result['count_statement']['scope']),indent=2)+'\n')
    print(json.dumps(dict(comparison=comparison,ratios=[{k:v for k,v in x.items() if k in ['n','width_ratio','spacing_ratio','period_increment','quadratic_core_coordinate']} for x in rows]),indent=2))


if __name__=='__main__':main()
