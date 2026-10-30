"""Compare the automatic return map with independently closed affine orbits."""
import argparse
import json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_event_map import EventMap
from bridge_chua_itineraries import local_transport


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    out = parser.parse_args().output
    families = {r['id']: r for r in json.loads((out/'family-registry.json').read_text())['families']}
    checks = []
    for name in ['S02', 'S08', 'P55', 'E09', 'E13']:
        item = families[name]; r = item['core']
        model = Chua(alpha=r['alpha'], m0=r['m0'], m1=r['m1'])
        closure, matrix, aux = EventMap(model, item['section_returns']).evaluate(np.r_[r['q'][:2], r['beta']])
        checks.append(dict(id=name, closure=float(max(abs(closure))), period_error=abs(aux['period']-r['period']),
            rho_error=float(abs(max(abs(np.linalg.eigvals(matrix)))-r['spectral_radius'])),
            half_inversion_error=aux['half_inversion_error']))
    passed = all(r['closure']<1e-8 and r['period_error']<1e-7 and r['rho_error']<2e-5 and r['half_inversion_error']<2e-7 for r in checks)
    report = dict(checks=checks, passed=passed, scope='Numerical agreement on selected reference orbits; not an interval event-detection proof.')
    (out/'event-map-checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    path = out/'off21-transport.json'
    if path.exists():
        steps = []
        for row in json.loads(path.read_text()):
            closure, matrix, aux = EventMap(Chua(alpha=row['alpha'],m0=row['m0']),3).evaluate(np.array(row['w']))
            steps.append(dict(alpha=row['alpha'],m0=row['m0'],beta=row['w'][-1],closure=float(max(abs(closure))),
                half_inversion_error=aux['half_inversion_error'],half_trace=aux['half_trace']))
        path_passed = all(r['closure']<2e-7 and r['half_inversion_error']<2e-7 and abs(r['half_trace']-.5)<2e-5 for r in steps)
        (out/'off21-path-symmetry-check.json').write_text(json.dumps(dict(checks=steps,passed=path_passed),indent=2)+'\n',encoding='utf-8',newline='\n')
        passed = passed and path_passed
    sheet=json.loads((out/'E09-sheet.json').read_text());nodes={tuple(n['index']):n for n in sheet['nodes']}
    template=dict(symmetric=True,target_multiplier=.5,regions=sheet['regions'],targets=sheet['targets']);local=[]
    for ij in [(4,3),(4,5),(3,4),(5,4),(4,8)]:
        target=nodes[ij];node,history=local_transport(template,nodes[4,4],target['alpha'],target['m0'])
        local.append(dict(index=list(ij),beta_error=abs(node['core']['beta']-target['core']['beta'])))
    (out/'local-bridge-checks.json').write_text(json.dumps(local,indent=2)+'\n',encoding='utf-8',newline='\n')
    passed=passed and max(r['beta_error'] for r in local)<1e-8
    print(json.dumps(report),flush=True)
    if not passed:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
