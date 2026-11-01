"""Replay diagram identifications and audit resolution/precision agreement.

--replay-diagrams additionally repeats the reduction from every raw PD code.
Neither mode certifies the polygon-to-ODE isotopy or its geometric predicates.
"""
import argparse,hashlib,json,random
from collections import Counter,defaultdict
from pathlib import Path
import mpmath as mp
import numpy as np
from spherogram import Link
from spherogram.links import simplify
from chua_knot_diagrams import canonical_gauss,fox_coefficients,polynomial_coefficients,identify_summand
from reconstruct_chua_knot_orbits import rows

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/chua-knot-atlas'


def replay_identification(s):
    assert s['type']!='unresolved'
    k=Link(s['start_pd'])
    for descriptor in s['type_III_moves']:
        matches=[tri for tri in simplify.possible_type_III_moves(k)
                 if [[k.crossings.index(c),e] for c,e in tri]==descriptor]
        assert len(matches)==1
        simplify.reidemeister_III(k,matches[0]);k._rebuild()
    assert canonical_gauss(k)==canonical_gauss(Link(s['final_pd']))
    reference=Link(s['type'].replace(' mirror',''))
    if s['type'].endswith(' mirror'):reference=reference.mirror()
    assert canonical_gauss(k)==canonical_gauss(reference)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--replay-diagrams',action='store_true');args=parser.parse_args()
    groups=defaultdict(list);paths=sorted((OUT/'diagrams').glob('*.json'));replayed=0
    for path in paths:
        r=json.loads(path.read_text());groups[r['id']].append(r)
        assert hashlib.sha256((OUT/r['orbit_file']).read_bytes()).hexdigest()==r['orbit_sha256']
        link=Link(r['simplified_pd_code'])
        assert polynomial_coefficients(link)==fox_coefficients(link)==r['alexander_coefficients']==r['fox_alexander_coefficients']
        assert len(link.crossings)==r['simplified_crossings']
        pieces=link.deconnect_sum() if link.crossings else []
        assert Counter(canonical_gauss(k) for k in pieces)==Counter(canonical_gauss(Link(s['start_pd'])) for s in r['summands'])
        for s in r['summands']:replay_identification(s)
        if args.replay_diagrams:
            raw=Link(r['pd_code']);assert raw.is_planar() and len(raw.link_components)==1
            random.seed(20261101+r['seed']);raw.simplify('global')
            assert len(raw.crossings)==r['simplified_crossings']
            assert polynomial_coefficients(raw)==r['alexander_coefficients']
            names=sorted(identify_summand(k)['type'] for k in raw.deconnect_sum()) if raw.crossings else []
            assert (' # '.join(names) or 'unknot')==r['knot_type']
            replayed+=1
        assert r['minimum_crossing_height']>0 and r['minimum_vertex_fraction']>0
    assert len(groups)==34 and len(paths)==108
    summary=[];mp.mp.dps=90
    for source in rows():
        id=source['id'];rr=groups[id]
        assert len({r['knot_type'] for r in rr})==1
        required={f'{id}-d60-n96.npz',f'{id}-d70-n192.npz'}
        assert required.issubset({r['orbit_file'] for r in rr})
        assert {r['seed'] for r in rr if r['orbit_file']==f'{id}-d70-n192.npz'}=={0,2}
        if id in ('A04','A05','P85'):
            assert {r['seed'] for r in rr if r['orbit_file']==f'{id}-d80-n384.npz'}=={0,2}
            third=json.loads((OUT/f'{id}-d80-n384.json').read_text())
            assert float(third['residual'])<1e-69
        a=json.loads((OUT/f'{id}-d60-n96.json').read_text());b=json.loads((OUT/f'{id}-d70-n192.json').read_text())
        assert float(a['residual'])<1e-39 and float(b['residual'])<1e-49
        assert float(b['beta_change'])<3e-13
        assert b['regions']==source['regions'] and b['targets']==source['targets']
        assert min(mp.mpf(v) for v in b['q'][2:])>0
        assert b['sampling']['maximum_sampled_region_violation']<1e-39
        x=np.load(OUT/f'{id}-d60-n96.npz')['xyz'];y=np.load(OUT/f'{id}-d70-n192.npz')['xyz']
        error=float(np.max(np.abs(x-y[::2])));assert error<1e-13
        assert np.array_equal(x[0],x[-1]) and np.array_equal(y[0],y[-1])
        if b['symmetric']:
            half=(len(y)-1)//2;assert np.array_equal(y[half:-1],-y[:half])
        # The other normalized-map multiplier is det / target. Symmetric full
        # multipliers are their squares; the chosen target is already inside 1.
        other=abs(mp.mpf(b['determinant'])/mp.mpf(b['target_multiplier']));assert other<1
        summary.append(dict(id=id,knot_type=rr[0]['knot_type'],peaks=source.get('peak_count',int(id[1:]) if id.startswith('P') else None),
            beta=b['beta'],period=b['period'],symmetric=b['symmetric'],checks=len(rr),
            alexander_coefficients=rr[0]['alexander_coefficients'],determinant=rr[0]['determinant'],
            precision_grid_difference=error,minimum_crossing_height=min(r['minimum_crossing_height'] for r in rr),
            numerical_stability_other_map_multiplier=float(other)))
    thin=['S04','S05','S06','S07','P55','P65','P75','P85']
    assert all(next(r for r in summary if r['id']==id)['knot_type']=='unknot' for id in thin)
    controls=json.loads((OUT/'controls.json').read_text());assert len(controls['tests'])==22
    counter=json.loads((OUT/'sampling-counterexample.json').read_text())
    assert len(counter['uniform'][0]['alexander_coefficients'])==13
    assert counter['uniform'][1]['alexander_coefficients']==next(r for r in summary if r['id']=='A05')['alexander_coefficients']
    report=dict(families=sorted(summary,key=lambda r:r['id']),family_count=len(summary),diagram_count=len(paths),
        counts=dict(Counter(r['knot_type'] for r in summary)),thin_sequence=thin,
        analytic_curve_controls=len(controls['tests']),raw_diagrams_replayed=replayed,
        scope='Numerical orbit topology. Exact algebra and diagram matches do not enclose the ODE curve.')
    (OUT/'review.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='families'},indent=2))


if __name__=='__main__':main()
