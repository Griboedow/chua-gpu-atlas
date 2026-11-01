"""Audit repeated-precision enclosures and all stored independent controls."""
import json
from flint import arb,ctx
from study_chua_hidden_boundary import OUT,dump


def read(name):
    return json.loads((OUT/name).read_text())


def ball(d):
    return arb(d.get('enclosure',d['ball']))


def main():
    ctx.dps=110
    overlap=0
    for base in ['saddle-certificate','global-sector']:
        low,high=[read(f'{base}-{d}.json') for d in [80,100]]
        assert low['passed'] and high['passed']
        if base=='saddle-certificate':
            pairs=[]
            for name in ['saddle','attractor']:
                a,b=low[name],high[name]
                pairs.append((a['period'],b['period']))
                pairs.extend(zip(a['root'],b['root']))
                pairs.extend(zip(a['multipliers'],b['multipliers']))
            assert low['shifted_root_rejected'] and high['shifted_root_rejected']
        else:
            pairs=[(low['constants'][k],high['constants'][k]) for k in low['constants']]
        for a,b in pairs:
            assert ball(a).overlaps(ball(b))
            overlap+=1
    controls=read('cycle-independent-controls.json')['controls']
    scaling=read('linear-contact-results.json')['records']
    independent=[r['independent'] for r in scaling if 'independent' in r]
    assert len(scaling)==62 and len(independent)==8
    for c in controls+independent:
        assert c['success'] and c['closure']<1e-9 and c['multiplier_error']<1e-8
    nearest=[r for r in scaling if r['delta']==1e-7]
    errors=[abs(r[k]-1) for r in nearest for k in ['amplitude_ratio','duration_ratio']]
    assert len(nearest)==2 and max(errors)<.000061
    basin=read('basin-separator.json')
    assert len(basin['curve'])==71 and len(basin['controls'])==16
    assert all(c['label']==c['expected'] for c in basin['controls'])
    maximum_difference=max(c['difference'] for c in basin['curve'])
    assert maximum_difference<1e-12
    algebra=read('global-sector-algebra.json')
    assert all(algebra['identities'].values()) and algebra['wrong_feedthrough_rejected']
    dump(OUT/'numerical-audit.json',dict(passed=True,overlapping_80_100_digit_enclosures=overlap,
        independent_periodic_controls=len(controls),independent_scaling_controls=len(independent),
        maximum_independent_closure=max(c['closure'] for c in controls+independent),
        maximum_independent_multiplier_error=max(c['multiplier_error'] for c in controls+independent),
        nearest_contact_maximum_relative_prediction_error=max(errors),
        basin_two_sided_controls=len(basin['controls']),maximum_four_six_iterate_curve_difference=maximum_difference,
        exact_symbolic_identity_count=len(algebra['identities']),
        scope='Checks recorded computations. Floating-point controls do not certify a global basin boundary or asymptotic remainder.'))
    print('Repeated precision, independent integration, local separator and scaling checks passed.')


if __name__=='__main__':main()
