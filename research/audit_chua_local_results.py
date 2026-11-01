"""Audit the finite certificates used by the two local-results arguments.

This is a consistency/provenance audit. The localization, graph transform and
Newhouse--Robinson theorem are mathematical arguments, not code assertions.
"""
import hashlib
import json
from collections import Counter
from flint import arb,ctx
from certify_chua_basin_separator import ROOT,OUT


def read(path):return json.loads((ROOT/path).read_text())
def ball(v):return arb(v.get('enclosure',v['ball']))
def leaves(flow):
    if flow is None:return []
    if flow['method']=='subdivision':return sum((leaves(x) for x in flow['children']),[])
    return [flow]
def digest(path):return hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    ctx.dps=120;basins=[];tangencies=[];inputs=[]
    for digits in [80,100]:
        path=f'results/chua-basin-proof/separator-{digits}.json';inputs.append(path);r=read(path)
        assert r['passed'] and not r['partial'] and all(r['local']['gates'].values())
        assert len(r['exits'])==128 and r['pieces_per_side']==64
        for side,label in [(-1,'origin'),(1,'cycle')]:
            rows=[e for e in r['exits'] if e['side']==side]
            assert [e['index'] for e in rows]==list(range(64))
            assert all(e['label']==label and e['error'] is None for e in rows)
            intervals=[side*ball(e['initial_u']) for e in rows]
            assert intervals[0].contains(arb(r['local']['delta']))
            assert intervals[-1].overlaps(ball(r['exit_u_max']))
            assert all(a.overlaps(b) for a,b in zip(intervals[:-1],intervals[1:]))
        flows=sum((leaves(e['continuous']) for e in r['exits']),[])
        assert all(x['passed'] and ball(x['energy_upper'])<ball(r['origin_capture']['level']) for x in flows)
        basins.append(dict(digits=digits,covered_rectangles=128,
            event_returns=sum(e['returns'] for e in r['exits']),
            event_time_cells=sum(e['path_cells'] for e in r['exits']),
            continuous_rectangles=sum(e['continuous'] is not None for e in r['exits']),
            continuous_leaves=len(flows),continuous_steps=sum(x['steps'] for x in flows),
            continuous_methods=dict(Counter(x['method'] for x in flows)),
            origin_returns=[min(e['returns'] for e in r['exits'] if e['side']<0),max(e['returns'] for e in r['exits'] if e['side']<0)],
            cycle_returns=[min(e['returns'] for e in r['exits'] if e['side']>0),max(e['returns'] for e in r['exits'] if e['side']>0)]))
        path=f'results/chua-layer-infinity/tangency-certificate-{digits}.json';inputs.append(path);t=read(path)
        assert t['status']=='passed' and all(t['gates'].values()) and all(t['local_gates'].values()) and all(t['graph_gates'].values())
        assert ball(t['full_return_determinant'])>0 and ball(t['full_return_determinant'])<1
        assert ball(t['normal_curvature'])>0 and ball(t['parameter_unfolding'])>0
        assert t['return_iterations']==27
        tangencies.append(dict(digits=digits,gates=len(t['gates']),
            flights=sum(len(x) for x in t['itinerary_cells']),time_cells=sum(map(sum,t['itinerary_cells'])),
            contact_for_ordinary_return='H^28 = P^14 after one local half return'))
    controls=read('results/chua-basin-proof/independent-checks.json');assert controls['passed'] and controls['shifted_root_rejected']
    inputs+=['research/certify_chua_basin_separator.py','research/check_chua_basin_proof.py',
             'research/certify_chua_homoclinic_tangency.py','research/chua_interval_map.py',
             'research/chua_interval_connection.py','research/certify_chua_layer_organizer.py',
             'docs/chua-basin-and-coexistence.md','results/chua-basin-proof/independent-checks.json']
    report=dict(passed=True,basin_certificates=basins,tangency_certificates=tangencies,
        theorem_A='Computer-assisted local basin separation; graph-transform argument in the report.',
        theorem_B='Theorem consequence: localized parameter Newhouse--Robinson, with the report\'s domain argument.',
        formal_scope='Neither complete theorem is verified in Lean. This program checks finite premises and coverage only.',
        source_sha256={p:digest(p) for p in inputs})
    (OUT/'audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=True,basins=basins,tangencies=tangencies),indent=2))


if __name__=='__main__':main()
