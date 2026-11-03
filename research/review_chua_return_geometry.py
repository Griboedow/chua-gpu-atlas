"""Audit preserved roots, forecast scores, scaling and independent controls."""
import ast
import json
from pathlib import Path
import numpy as np
from study_chua_return_geometry import ROOT, OUT, IDS, read, digest


def main():
    checks={}
    forecast=read(OUT/'forecasts-frozen.json')
    checks['forecast_sources_unchanged']=all(digest(ROOT/p)==h for p,h in forecast['source_sha256'].items())
    references=[]
    for path in sorted(OUT.glob('*-d70-s8.json')):
        row=read(path)
        assert digest(ROOT/row['source'])==row['source_sha256'],path
        assert row['passed'],path
        optimized=read(OUT/f"{row['id']}-optimized.json")
        assert optimized['success'] and optimized['maximum']['gain']>=row['maximum_gain']*(1-1e-8),path
        with np.load(path.with_suffix('.npz')) as z:
            assert abs(z['gain'].max()/row['maximum_gain']-1)<1e-12
            assert np.isfinite(z['local']).all()
        references.append(row)
    checks['all_36_roots_and_gain_records']=len(references)==36
    fine=[]
    for id in IDS+['P95','P105']:
        row=read(OUT/f'{id}-d85-s128.json')
        opt=read(OUT/f'{id}-optimized.json')['maximum']['gain']
        fine.append(dict(id=id,sampled_gain=row['maximum_gain'],refined_gain=opt,
                         relative_difference=abs(row['maximum_gain']/opt-1)))
    checks['fine_mesh_agrees_with_local_refinement']=max(r['relative_difference'] for r in fine)<.002
    independent=[]
    for p in sorted(OUT.glob('*-d*-s*-check.json')):
        r=read(p)
        assert r['passed'] and digest(OUT/r['source'])==r['source_sha256'],p
        independent.append(dict(file=p.name,maximum_relative_gain_error=max(c['relative_gain_error'] for c in r['checks'])))
    checks['independent_variation']=len(independent)>=4
    scores=[]
    for target in forecast['targets']:
        id=target['id'];edge=read(OUT/f'{id}-edges.json');saddle=read(OUT/f'{id}-companion-d70.json')
        observed=dict(width=float(edge['width']),maximum_gain=read(OUT/f'{id}-optimized.json')['maximum']['gain'],
                      companion_section_distance=float(saddle['section_distance']),
                      companion_half_multiplier=float(saddle['companion']['half_multiplier']))
        errors={k:observed[k]/float(target[k])-1 for k in target['relative_gates']}
        accepted=all(abs(errors[k])<g for k,g in target['relative_gates'].items()) and abs(observed['companion_half_multiplier']-1.5)<target['multiplier_absolute_gate']
        assert edge['forecast_sha256']==digest(OUT/'forecasts-frozen.json')
        scores.append(dict(id=id,predicted=target,observed=observed,signed_relative_errors=errors,passed=accepted))
    checks['held_out_forecasts']=all(r['passed'] for r in scores)
    orbit_checks=[]
    for id,peaks in [('P95',95),('P105',105)]:
        r=read(OUT/f'{id}-independent-orbit.json')
        whole=r.get('whole_period',r.get('whole_period_refinement'));section=r['section_variation']
        accepted=whole['success'] and whole['peaks']==peaks and whole['section_returns']==peaks//5 and section['closure']<1e-10 and section['multiplier_error']<2e-5
        orbit_checks.append(dict(id=id,whole_period=whole,section_checks={k:v for k,v in section.items() if k!='checks'},passed=accepted))
    checks['new_orbit_counts_and_section_variation']=all(r['passed'] for r in orbit_checks)
    precision=[read(OUT/f'{id}-precision-check.json') for id in ['P95','P105']]
    checks['new_orbit_eigen_flow_and_precision']=all(r['passed'] for r in precision)
    edges=[read(OUT/f'{id}-edge-precision-check.json') for id in ['P95','P105']]
    checks['new_window_endpoint_precision']=all(r['passed'] and digest(OUT/r['source'])==r['source_sha256'] for r in edges)
    sequence=read(OUT/'sequence-analysis.json')['rows']
    ratios=[];constants=[]
    for r in sequence:
        g=read(OUT/f"{r['id']}-optimized.json")['maximum']['gain']
        constants.append(r['width']*g*g)
        ratios.append(g)
    last_ratios=[ratios[-2]/ratios[-3],ratios[-1]/ratios[-2]]
    expansion=float(forecast['expansion'])
    checks['initial_gain_ratio_gate']=all(abs(r/expansion-1)<.2 for r in last_ratios)
    checks['initial_product_gate']=max(constants[-4:])/min(constants[-4:])<1.5
    low=read(OUT/'P85-companion-d70.json');high=read(OUT/'P85-companion-d90.json')
    # Compare decimal strings using high precision, not float64 subtraction.
    import mpmath as mp
    mp.mp.dps=95
    companions=[]
    for id in IDS+['P95','P105']:
        r=read(OUT/f'{id}-companion-d70.json');source=ROOT/r['source']
        assert r['passed'] and digest(source)==r['source_sha256']
        assert abs(mp.mpf(r['companion']['beta'])-mp.mpf(read(source)['beta']))<mp.mpf('1e-50')
        companions.append(id)
    checks['all_10_fixed_parameter_companions']=len(companions)==10
    diff=max(abs(mp.mpf(a)-mp.mpf(b)) for a,b in zip(low['companion']['q'],high['companion']['q']))
    checks['companion_precision_repeat']=diff<mp.mpf('1e-40')
    maps=[read(OUT/f'{id}-normal-form-d75.json') for id in ['P55','P85']]
    high_map=read(OUT/'P85-normal-form-d90.json')
    map_diff=max(abs(float(a['return_coordinate'])-float(b['return_coordinate'])) for a,b in zip(maps[-1]['rows'],high_map['rows']))
    checks['normal_form_precision_repeat']=map_diff<1e-10
    checks['normal_form_sample_paths']=all(r['event_checks']['maximum_sampled_region_violation']<1e-20 for r in high_map['rows'])
    result=dict(checks={k:bool(v) for k,v in checks.items()},passed=all(checks.values()),
        reference_count=len(references),forecast_scores=scores,mesh_comparisons=fine,
        independent_gain_checks=independent,new_orbit_checks=orbit_checks,
        new_orbit_precision_checks=precision,
        endpoint_precision_checks=[dict(id=r['id'],relative_width_change=r['relative_width_change']) for r in edges],
        sequence_width_gain_squared=constants,last_gain_ratios=last_ratios,
        companion_precision_difference=mp.nstr(diff,15),normal_form_repeat_difference=map_diff,
        normal_form_maximum_errors={r['id']:r['maximum_absolute_error'] for r in maps},
        scope='Numerical consistency, prediction and independent implementation audit; no interval certificate of gain maxima, basin boundaries or global atlas completeness')
    (OUT/'review.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(2)


if __name__=='__main__':main()
