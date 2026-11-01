"""Continue the verified five- and seven-peak symmetric representatives in slope volumes.

Family labels are inherited by local continuation from closed reference cycles,
not assigned by the number of disconnected pieces in a sampled mask.
"""
import argparse
from collections import deque
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq
from chua_fast_flights import FastFlights
from chua_layer_bifurcations import from_orbit
from trace_chua_window import Chua

ROOT = Path(__file__).resolve().parents[1]


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8', newline='\n')


class CoordinateFlights(FastFlights):
    """Use any one parameter as the unknown in the existing boundary solver."""
    def __init__(self, parameters, regions, targets, free='m0', symmetric=True):
        self.parameters = {k:float(parameters[k]) for k in ['alpha', 'beta', 'm0', 'm1']}
        self.free = free
        super().__init__(Chua(**{k: parameters[k] for k in ['alpha', 'm0', 'm1']}),
                         regions, targets, symmetric)
        self._coordinate = None

    def evaluate(self, value, q, check=False):
        if value != self._coordinate:
            self.parameters[self.free] = float(value)
            self.model = Chua(**{k: self.parameters[k] for k in ['alpha', 'm0', 'm1']})
            self._beta = None
            self._coordinate = value
        return super().evaluate(self.parameters['beta'], q, check)

    def record(self, value, q, check=False):
        result = super().record(value, q, check)
        result.update(self.parameters)
        return result


def valid(record):
    if record['residual'] > 1e-8 or record['region_violation'] > 1e-7:
        raise RuntimeError('Closing residual or affine itinerary failed')
    if min(record['q'][2:]) <= 0 or record['minimum_breakpoint_speed'] < .02:
        raise RuntimeError('Nonpositive flight or nearly grazing breakpoint')


def extrema_check(f, record):
    """Check arc extrema, rather than only uniformly sampled region membership."""
    f.evaluate(record[f.free], np.array(record['q']))
    matrices = f.matrices(record['beta'])
    turns = 0
    violation = 0.
    for arc in record['arcs']:
        region, duration = arc['region'], arc['duration']
        A, c, eig, V, W = matrices[region]
        coefficients = V[0] * (W @ (np.array(arc['start'])-c))
        def x(t):
            return float(c[0]+np.real(np.sum(coefficients*np.exp(eig*t))))
        def dx(t):
            return float(np.real(np.sum(eig*coefficients*np.exp(eig*t))))
        grid = np.linspace(0, duration, max(3, int(np.ceil(duration/.025))+1))
        speeds = np.real(np.exp(grid[:, None]*eig) @ (eig*coefficients))
        roots = [brentq(dx, a, b, xtol=1e-13) for a, b, da, db in
                 zip(grid[:-1], grid[1:], speeds[:-1], speeds[1:]) if da*db < 0]
        turns += len(roots)
        values = np.array([x(t) for t in [0., *roots, duration]])
        violation = max(violation, float(np.max(abs(values)-1 if region == 0 else 1-region*values)))
    if violation > 1e-7:
        raise RuntimeError('An affine arc leaves its actual diode region')
    return dict(full_cycle_peaks=turns, extrema_region_violation=violation,
                derivative_bracket_max_step=.025)


def anchors(out):
    seeds = json.loads((ROOT/'results/chua-slope-plane/closed-cycles.json').read_text())['cycles']
    records = []
    for seed in seeds:
        r = seed['cycle']
        params = {k: r[k] for k in ['alpha', 'beta', 'm0', 'm1']}
        model = Chua(**{k: params[k] for k in ['alpha', 'm0', 'm1']})
        template, q = from_orbit(model, params['beta'], r['shooting_vector'], True)
        f = CoordinateFlights(params, template.regions, template.targets)
        item = dict(id='P'+str(seed['peak_count']), anchor=params,
                    regions=template.regions, targets=template.targets, reference_peak_count=seed['peak_count'])
        for name, target in [('fold', 1), ('symmetry', -1)]:
            m0, qq = f.boundary(params['m0'], q, target)
            edge = f.record(m0, qq, True)
            valid(edge)
            item[name] = edge
        a, b = item['fold']['m0'], item['symmetry']['m0']
        m0 = (a+b)/2
        qq = f.solve(m0, np.array(item['symmetry']['q']))
        item['representative'] = f.record(m0, qq, True)
        valid(item['representative'])
        item['representative'].update(extrema_check(f, item['representative']))
        assert item['representative']['spectral_radius'] < .9999
        item['width_m0'] = abs(b-a)
        records.append(item)
        print(item['id'], 'edges', a, b, 'width', item['width_m0'],
              'period', item['representative']['period'], flush=True)
    dump(out/'anchors.json', records)
    return records


def advance(template, source, parameters, depth=0):
    """Implicit tangent prediction, corrected boundaries, and an attracting representative."""
    try:
        f = CoordinateFlights(parameters, template['regions'], template['targets'])
        source_parameters = {k: source['representative'][k] for k in ['alpha', 'beta', 'm0', 'm1']}
        old = CoordinateFlights(source_parameters, template['regions'], template['targets'])
        row = {}
        for name, target in [('fold', 1), ('symmetry', -1)]:
            previous = source[name]
            w = np.r_[previous['q'], previous['m0']]
            forcing = np.zeros(len(w))
            for parameter in ['alpha', 'beta', 'm1']:
                delta = parameters[parameter]-source_parameters[parameter]
                if delta == 0:
                    continue
                h = 1e-6 if parameter != 'm1' else 1e-7
                values = []
                for sign in [-1, 1]:
                    perturbed = dict(source_parameters)
                    perturbed[parameter] += sign*h
                    trial = CoordinateFlights(perturbed, template['regions'], template['targets'])
                    values.append(trial.boundary_system(w, target))
                forcing += delta*(values[1]-values[0])/(2*h)
            prediction = w-np.linalg.solve(old.boundary_system(w, target, True), forcing)
            m0, q = f.boundary(prediction[-1], prediction[:-1], target)
            r = f.record(m0, q, True)
            valid(r)
            r.update(extrema_check(f, r))
            if abs(r['period']/previous['period']-1) > .08:
                raise RuntimeError('Period jump to another closing root')
            if not -1.5 < m0 < -.9:
                raise RuntimeError('Boundary left the original slope-volume m0 domain')
            row[name] = r
        width = row['fold']['m0']-row['symmetry']['m0']
        if width < 1e-10 or abs(width/source['width_m0']-1) > .7:
            raise RuntimeError('Stability width collapsed or jumped')
        m0 = (row['fold']['m0']+row['symmetry']['m0'])/2
        q = f.solve(m0, np.array(row['symmetry']['q']))
        row['representative'] = f.record(m0, q, True)
        valid(row['representative'])
        row['representative'].update(extrema_check(f, row['representative']))
        if row['representative']['spectral_radius'] >= .9999:
            raise RuntimeError('The continued representative is not attracting')
        if row['representative']['full_cycle_peaks'] != template['reference_peak_count']:
            raise RuntimeError('Peak count changed along the continued family')
        row['width_m0'] = width
        return row
    except (RuntimeError, ValueError, OverflowError, np.linalg.LinAlgError) as exc:
        if depth >= 4:
            raise RuntimeError(str(exc))
        middle = {k: (parameters[k]+source['representative'][k])/2 for k in ['alpha', 'beta', 'm0', 'm1']}
        halfway = advance(template, source, middle, depth+1)
        return advance(template, halfway, parameters, depth+1)


def sheets(out, nodes):
    templates = json.loads((out/'anchors.json').read_text())
    for template in templates:
        for kind, span in [('alpha', .4), ('beta', .8)]:
            path = out/(template['id']+'-'+kind+'.json')
            if path.exists() and json.loads(path.read_text()).get('complete'):
                continue
            started = time.monotonic()
            controls = template['anchor'][kind]+np.linspace(-span, span, nodes)
            slopes = template['anchor']['m1']+np.linspace(-.08, .08, nodes)
            center = nodes//2
            initial = advance(template, template, template['anchor'])
            initial['index'] = [center, center]
            found = {(center, center): initial}
            attempted = {(center, center)}
            queue = deque([(center, center)])
            failures = []
            while queue:
                i, j = queue.popleft()
                for ii, jj in [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]:
                    if not (0 <= ii < nodes and 0 <= jj < nodes) or (ii, jj) in attempted:
                        continue
                    attempted.add((ii, jj))
                    params = dict(template['anchor'])
                    params.update({kind: float(controls[ii]), 'm1': float(slopes[jj])})
                    try:
                        row = advance(template, found[i, j], params)
                        row['index'] = [ii, jj]
                        found[ii, jj] = row
                        queue.append((ii, jj))
                    except (RuntimeError, ValueError, OverflowError, np.linalg.LinAlgError) as exc:
                        failures.append(dict(index=[ii, jj], parameters=params, error=str(exc)))
                if len(attempted) % 40 < 4 or not queue:
                    report = dict(id=template['id'], kind=kind, reference_peak_count=template['reference_peak_count'],
                        regions=template['regions'], targets=template['targets'], controls=controls.tolist(), m1=slopes.tolist(),
                        nodes=list(found.values()), failures=failures, complete=len(found) == nodes*nodes,
                        elapsed_seconds=time.monotonic()-started,
                        scope='Locally continued attracting symmetric representatives and their two marginal-stability edges. '
                              'A failed continuation is not evidence of an orbit ending. Full stable asymmetric extensions are not included.')
                    dump(path, report)
                    print(template['id'], kind, len(found), '/', nodes*nodes, 'failures', len(failures), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=['anchors', 'sheets'], default='anchors')
    parser.add_argument('--nodes', type=int, default=21)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.nodes < 3 or args.nodes % 2 != 1:
        parser.error('--nodes must be odd and at least three')
    if args.stage == 'anchors':
        anchors(args.output)
    else:
        sheets(args.output, args.nodes)


if __name__ == '__main__':
    main()
