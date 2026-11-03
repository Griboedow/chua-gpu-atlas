"""Recompute nested-return record checks and retain every forecast miss."""
import subprocess
import mpmath as mp
import numpy as np
from test_chua_nested_returns import ROOT,OUT,read,save,digest
from chua_block_return import BlockReturn


def main():
    mp.mp.dps=70;frozen=read(OUT/'forecasts-frozen.json');scores=read(OUT/'scores.json');cal=read(OUT/'calibration.json')
    gates={}
    for name in ['forecasts-frozen.json','calibration.json']:
        original=subprocess.check_output(['git','show',f'eaa2eca:results/chua-nested-returns/{name}'],cwd=ROOT)
        gates[name+'_unchanged']=original.replace(b'\r\n',b'\n')==(OUT/name).read_bytes().replace(b'\r\n',b'\n')
    gates['protocol_sources']=all(digest(ROOT/p)==h for p,h in frozen['source_sha256'].items())
    gates['parent_sources']=all(digest(ROOT/c[k])==c[k+'_sha256'] for c in cal['parents'] for k in ['source','companion_source'])
    observations=[];counts=dict(location=0,width=0,period=0)
    for t in frozen['targets']:
        id=t['id'];r=read(OUT/f'{id}-reference.json');a=read(OUT/f'{id}-audit.json');e=read(OUT/f'{id}-edges.json');g=read(OUT/f'{id}-geometry.json');p=read(OUT/f'{id}-precision.json')
        s=next(v for v in scores['rows'] if v['id']==id);beta=mp.mpf(r['beta']);width=abs(mp.mpf(e['fold']['beta'])-mp.mpf(e['symmetry']['beta']))
        ends=sorted(mp.mpf(e[k]['beta']) for k in ['fold','symmetry'])
        err=dict(location_in_widths=float((beta-mp.mpf(t['beta']))/mp.mpf(e['width'])),
                 width_relative=float(mp.mpf(e['width'])/mp.mpf(t['width'])-1),period_absolute=float(mp.mpf(r['period'])-mp.mpf(t['period'])))
        checks=dict(audits=a['passed'] and all(a['gates'].values()),interior=ends[0]<beta<ends[1],
            # Original subtraction used 60-digit arithmetic on beta near 14.
            width=abs(width-mp.mpf(e['width']))<mp.mpf('1e-58'),
            hashes=a['source_sha256']==g['source_sha256']==p['source_sha256']==digest(OUT/f'{id}-reference.json'),
            forecasts=r['forecast_sha256']==e['forecast_sha256']==scores['forecast_sha256']==digest(OUT/'forecasts-frozen.json'),
            repeated_word=g['exact_parent_word_repetition'] and g['side_word_minimal_period']==7,
            same_parameter=mp.mpf(g['parent']['beta'])==beta==mp.mpf(g['continued_attracting_branch']['beta'])==mp.mpf(g['short_parent']['beta']),
            section_count=len(g['points'])==7*t['block_count'])
        for quantity,key,limit in [('location','location_in_widths',1),('width','width_relative',.3),('period','period_absolute',.5)]:
            passed=abs(err[key])<=limit;counts[quantity]+=passed
            checks['score_'+quantity]=err[key]==s['errors'][key] and passed==s['passed'][quantity]
        ev=BlockReturn(g['parent'],70).evaluate(beta,[mp.mpf(v) for v in g['parent']['q'][:2]])
        tr=ev['matrix'][0,0]+ev['matrix'][1,1];mu=(tr+mp.sqrt(tr*tr-4*ev['determinant']))/2
        checks['long_parent_is_saddle']=mu>1 and abs(ev['determinant']/mu)<1
        gates[id]=bool(all(checks.values()))
        observations.append(dict(id=id,checks={k:bool(v) for k,v in checks.items()},
            parent_period=float(g['parent']['period']),parent_half_multiplier=float(mu),parent_full_multiplier=float(mu**2),
            short_period=float(g['short_parent']['period']),minimum_block_distance=min(g['block_rms_distance_to_parent']),
            minimum_other_branch_distance=min(g['block_rms_distance_to_other_parent']),
            period_residual=g['period_minus_parent_multiple'],minimum_distinct_section_distance=a['minimum_distinct_section_distance'],
            maximum_beta_precision_change=max(float(v['beta_change']) for v in p['controls'].values()),
            relative_width_precision_change=p['relative_width_change'],
            independent_closure=float(a['exponential_reconstruction']['half_closure']),
            independent_multiplier_error=a['section_controls'][-1]['multiplier_error'],
            next_section_error=a['section_controls'][-1]['maximum_next_point_error'],
            whole_period_drift=max(w['closure'] for w in a['whole_period'])))
    controls=read(OUT/'negative-controls.json');c=controls['ordinary_105']
    gates['ordinary_105_control']=c['side_word_minimal_period']==21 and c['same_side_transitions']==1 and digest(ROOT/c['source'])==c['source_sha256']
    gates['repeated_parent_rejected']=all(c['rejected_as_nonprimitive'] and c['minimum_distinct_section_distance']==0 and abs(float(c['apparent_signed_half_multiplier'])-.5)<1e-15 and float(c['block_closure_residual'])<1e-30 for c in controls['repeated_parent'])
    increments=[]
    for parent in ['S06','S11-n3']:
        rows=[next(o for o in observations if o['id']==f'{parent}-k{k}') for k in [3,5,7]]
        d=[r['minimum_block_distance'] for r in rows];gates[parent+'_increasing_shadowing']=d[2]<d[1]<d[0]
        for k in [5,7]:
            a,b=[read(OUT/f'{parent}-k{j}-reference.json') for j in [k-2,k]]
            obs=next(o for o in rows if o['id']==b['id']);dt=float(b['period'])-float(a['period'])
            increments.append(dict(parent=parent,k=k,period_increment=dt,twice_long_parent_period=2*obs['parent_period'],
                twice_short_parent_period=2*obs['short_period'],long_period_error=dt-2*obs['parent_period']))
        gates[parent+'_long_clock']=all(abs(v['long_period_error'])<.002 for v in increments if v['parent']==parent)
    broad=read(OUT/'parent-block-map.json');gates['sampled_map_independent']=all(c['independent_error']<1e-45 and c['checks']['maximum_sampled_region_violation']<1e-30 and c['checks']['minimum_duration']>0 for row in broad['rows'] for c in row['independent_controls'])
    branches=read(OUT/'parent-branch-continuation.json')['rows'];gates['fold_connection']=all([float(r['target_multiplier']) for r in row['roots']]==[.5,.75,1,1.25,1.5,2,3,4] and max(float(r['residual']) for r in row['roots'])<1e-30 for row in branches)
    gates['score_totals']=counts==scores['passed_by_quantity']==dict(location=2,width=6,period=6)
    summary=dict(gates=gates,passed=all(gates.values()),observations=observations,period_increments=increments,forecast_passes=counts,
        scope='Record consistency and numerical controls, not interval existence or an arbitrary-depth theorem. Forecast misses are retained.')
    save(OUT/'review.json',summary);print('GATES',gates);print('FORECASTS',counts);print('PERIOD INCREMENTS',increments)
    if not summary['passed']:raise SystemExit('Review failed')


if __name__=='__main__':main()
