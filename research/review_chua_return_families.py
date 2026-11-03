"""Replay record-level checks without treating forecast success as an audit gate."""
import subprocess
import mpmath as mp
import numpy as np
from test_chua_return_families import OUT,ROOT,read,save,digest


def main():
    forecast=read(OUT/'forecasts-frozen.json');cal=read(OUT/'calibration.json');scores=read(OUT/'scores.json')
    frozen=subprocess.check_output(['git','show','a2444a0:results/chua-return-families/forecasts-frozen.json'],cwd=ROOT)
    gates=dict(forecast_unchanged=frozen.replace(b'\r\n',b'\n')==(OUT/'forecasts-frozen.json').read_bytes().replace(b'\r\n',b'\n'),
        calibration_sources=all(digest(ROOT/r[k])==r[k+'_sha256'] for seq in cal['sequences'] for r in seq['calibration'] for k in ['source','gain_source']),
        protocol_sources=all(digest(ROOT/p)==d for p,d in forecast['source_sha256'].items()))
    observations=[]
    for t in forecast['targets']:
        id=t['id'];r=read(OUT/f'{id}-reference.json');a=read(OUT/f'{id}-audit.json');e=read(OUT/f'{id}-edges.json')
        g=read(OUT/f'{id}-gain-s8.json');fine=read(OUT/f'{id}-gain-s128.json');c=read(OUT/f'{id}-companion.json')
        s=next(v for v in scores['rows'] if v['id']==id)
        measured=dict(width=float(e['width']),gain=g['refined']['maximum']['gain'],period=float(r['period']))
        source=digest(OUT/f'{id}-reference.json');rel=abs(fine['maximum_gain']/measured['gain']-1)
        b=float(r['beta']);ends=sorted(float(e[k]['beta']) for k in ['fold','symmetry'])
        checks=dict(orbit_audit=a['passed'] and all(a['gates'].values()),
            references=a['source_sha256']==source==g['source_sha256']==fine['source_sha256'],
            forecast_digest=e['forecast_sha256']==digest(OUT/'forecasts-frozen.json')==scores['forecast_sha256'],
            interior_reference=ends[0]<b<ends[1],score=measured==s['measured'],
            gain_reconstruction=g['passed'] and fine['passed'],gain_refinement=g['refined']['success'] and rel<.002,
            companion=c['passed'] and c['half_multiplier_error']<.002,
            section_count=len(a['section_variation']['points'])==t['section_returns'])
        for key in ['width','gain']:
            err=measured[key]/t['corrected'][key]-1
            checks['score_'+key]=err==s['corrected_error'][key] and (abs(err)<=t['gates'][key+'_relative'])==s['passed'][key]
        gates[id]=all(checks.values())
        observations.append(dict(id=id,checks=checks,fine_gain_relative_difference=rel,
            relative_width_precision_change=a['relative_width_change'],
            eigen_closure=float(a['eigen_reconstruction']['closure']),
            full_field_maximum_closure=max(x['closure'] for x in a['whole_period']),
            full_field_maximum_symmetry_error=max(x['half_inversion_error'] for x in a['whole_period']),
            section_error=a['section_variation']['maximum_next_point_error'],
            independent_multiplier_error=a['section_variation']['multiplier_error'],
            companion_half_multiplier=float(c['companion']['signed_map_multiplier'])))
    normal=[]
    for id in ['S11-n5','E05-n5']:
        a,b=[read(OUT/f'{id}-normal-form-d{d}.json') for d in [75,90]];mp.mp.dps=60
        error=max(abs(mp.mpf(x['coordinate'])-mp.mpf(y['coordinate'])) for x,y in zip(a['samples'],b['samples']))
        control=read(OUT/f'{id}-gain-control.json')
        paths=all(x['checks']['minimum_duration']>0 and x['checks']['minimum_crossing_speed']>.01 and x['checks']['maximum_sampled_region_violation']<1e-30 for x in b['samples'])
        gates[id+'_nonlinear_repeat']=bool(error<mp.mpf('1e-20') and paths and b['maximum_error']<.001)
        gates[id+'_independent_gain']=control['passed']
        normal.append(dict(id=id,precision_difference=mp.nstr(error,15),maximum_return_error=b['maximum_error'],maximum_off_secant=b['maximum_off_secant'],maximum_gain_control_error=max(c['relative_gain_error'] for c in control['checks'])))
    for id in ['P55','S11-n5','E05-n5']:
        d=read(OUT/f'{id}-directions.json');gates[id+'_fixed_direction']=abs(d['floquet_final']-.25)<1e-8
    summary=dict(gates=gates,passed=all(gates.values()),observations=observations,normal_form=normal,
        forecast_passes=scores['passed_by_quantity'],scope='Numerical record audit. Failed forecast gates remain scientific results; no interval certification or global classification.')
    save(OUT/'review.json',summary)
    print('Audit gates:',gates,'FORECAST results:',scores['passed_by_quantity'],flush=True)
    for key in ['fine_gain_relative_difference','relative_width_precision_change','eigen_closure','full_field_maximum_closure','section_error','independent_multiplier_error']:
        print(key,max(r[key] for r in observations),flush=True)
    print('Nonlinear checks',normal,flush=True)
    if not summary['passed']:raise SystemExit('Record audit failed')


if __name__=='__main__':main()
