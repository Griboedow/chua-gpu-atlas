"""Resolve symmetric family boundaries with a separate itinerary for each root.

Fixed-itinerary roots are accepted only after checking all bracketed extrema.
When an affine arc is invalid, the actual piecewise field supplies a replacement
itinerary. A failed corrector remains a hole, never an inferred bifurcation.
"""
import argparse
from collections import deque
import json
from pathlib import Path
import time

import numpy as np
from continue_chua_slope_families import CoordinateFlights, extrema_check, dump
from trace_chua_window import Chua

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/chua-slope-refinement'
KEYS = ['alpha', 'beta', 'm0', 'm1']


def parameters(r):
    return {k: float(r[k]) for k in KEYS}


def attach(r, template):
    return {**r, 'regions': template['regions'], 'targets': template['targets']}


def flights(p, r):
    return CoordinateFlights(p, r['regions'], r['targets'])


def checked(f, value, q):
    r = f.record(value, q, True)
    if (r['residual'] > 2e-9 or r['region_violation'] > 2e-8
            or min(q[2:]) <= 0 or r['minimum_breakpoint_speed'] < 1e-7):
        raise RuntimeError('Invalid affine root or unresolved grazing crossing')
    r.update(extrema_check(f, r))
    if r['extrema_region_violation'] > 2e-8:
        raise RuntimeError('An affine extremum leaves its region')
    r.pop('arcs')
    return attach(r, dict(regions=f.regions, targets=f.targets))


def rebuild(p, r):
    """Recover event order from the actual field, not the extrapolated arcs."""
    model = Chua(alpha=p['alpha'], m0=p['m0'], m1=p['m1'], rtol=2e-12, step=.015)
    duration = r['period']/2
    sol = model.integrate(p['beta'], [0., *r['q'][:2]], duration)
    if not sol.success:
        raise RuntimeError('Actual-field itinerary integration failed')
    events = sorted((float(t), target) for times, target in zip(sol.t_events[:2], [-1, 1])
                    for t in times if 1e-8 < t < duration-1e-8)
    times = [0.]+[t for t, _ in events]+[duration]
    regions = []
    for left, right in zip(times[:-1], times[1:]):
        x = sol.sol((left+right)/2)[0]
        regions.append(0 if abs(x) < 1 else int(np.sign(x)))
    q = np.r_[r['q'][:2], np.diff(times)]
    # At a fold the old m0 may lie outside the new existence interval. Do not
    # insist on a fixed-m0 closing root before the augmented boundary solve.
    return CoordinateFlights(p, regions, [target for _, target in events]), q


def boundary(previous, p, target):
    oldp = parameters(previous)
    old = flights(oldp, previous)
    w = np.r_[previous['q'], previous['m0']]
    forcing = np.zeros(len(w))
    for parameter in ['alpha', 'beta', 'm1']:
        delta = p[parameter]-oldp[parameter]
        if not delta:
            continue
        h = 1e-7 if parameter == 'm1' else 1e-6
        values = []
        for sign in [-1, 1]:
            trialp = {**oldp, parameter: oldp[parameter]+sign*h}
            values.append(flights(trialp, previous).boundary_system(w, target))
        forcing += delta*(values[1]-values[0])/(2*h)
    prediction = w-np.linalg.solve(old.boundary_system(w, target, True), forcing)
    f = flights(p, previous)
    m0, q = f.boundary(prediction[-1], prediction[:-1], target)
    for attempt in range(5):
        try:
            r = checked(f, m0, q)
            if max(abs(f.boundary_system(np.r_[q, m0], target))) > 2e-9:
                raise RuntimeError('Boundary equation residual')
            if abs(r['period']/previous['period']-1) > .08 or not -1.5 < m0 < -.9:
                raise RuntimeError('Boundary jumped to another root or left domain')
            return r
        except RuntimeError:
            if attempt == 4:
                raise
            candidate = dict(q=q.tolist(), period=2*sum(q[2:]))
            f, q = rebuild({**p, 'm0': float(m0)}, candidate)
            m0, q = f.boundary(m0, q, target)


def representative(previous, p):
    f = flights(p, previous)
    q = f.solve(p['m0'], np.array(previous['q']))
    for attempt in range(5):
        try:
            r = checked(f, p['m0'], q)
            if r['spectral_radius'] >= .9999:
                raise RuntimeError('Representative is not attracting')
            return r
        except RuntimeError:
            if attempt == 4:
                raise
            f, q = rebuild(p, dict(q=q.tolist(), period=2*sum(q[2:])))
            q = f.solve(p['m0'], q)


def advance(source, p, peaks, depth=0):
    try:
        row = {name: boundary(source[name], p, target)
               for name, target in [('fold', 1), ('symmetry', -1)]}
        width = row['fold']['m0']-row['symmetry']['m0']
        if width < 1e-10 or abs(width/source['width_m0']-1) > .7:
            raise RuntimeError('Width collapsed or changed too rapidly')
        mid = (row['fold']['m0']+row['symmetry']['m0'])/2
        row['representative'] = representative(row['symmetry'], {**p, 'm0': mid})
        if row['representative']['full_cycle_peaks'] != peaks:
            raise RuntimeError('Peak count changed; requires a separate branch study')
        row['width_m0'] = width
        return row
    except (RuntimeError, ValueError, OverflowError, np.linalg.LinAlgError) as exc:
        if depth >= 5:
            raise RuntimeError(str(exc)) from exc
        middle = {k: (p[k]+source['representative'][k])/2 for k in KEYS}
        halfway = advance(source, middle, peaks, depth+1)
        return advance(halfway, p, peaks, depth+1)


def run(family, kind, count):
    started = time.monotonic()
    OUT.mkdir(exist_ok=True)
    templates = json.loads((ROOT/'results/chua-slope-families/anchors.json').read_text())
    t = next(x for x in templates if x['id'] == family)
    span = .4 if kind == 'alpha' else .8
    controls = t['anchor'][kind]+np.linspace(-span, span, count)
    slopes = t['anchor']['m1']+np.linspace(-.08, .08, count)
    old = json.loads((ROOT/f'results/chua-slope-families/{family}-{kind}.json').read_text())
    found = {}
    # Previously checked nodes anchor the branch on both coordinate directions.
    for row in old['nodes']:
        i, j = row['index']
        ii, jj = round(i*(count-1)/20), round(j*(count-1)/20)
        if abs(controls[ii]-old['controls'][i]) > 1e-10 or abs(slopes[jj]-old['m1'][j]) > 1e-10:
            continue
        fresh = {}
        for role in ['fold', 'symmetry', 'representative']:
            f = CoordinateFlights(parameters(row[role]), t['regions'], t['targets'])
            fresh[role] = checked(f, row[role]['m0'], np.array(row[role]['q']))
        fresh.update(width_m0=row['width_m0'], index=[ii, jj], origin='previous_mesh')
        found[ii, jj] = fresh
    queue = deque(found)
    attempts = {}
    errors = {}
    path = OUT/f'{family}-{kind}.json'
    last_saved = 0
    while queue:
        i, j = queue.popleft()
        for ii, jj in [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]:
            key = ii, jj
            if not (0 <= ii < count and 0 <= jj < count) or key in found:
                continue
            # Retry a failure from at most two different verified neighbours.
            if attempts.get(key, 0) >= 2:
                continue
            attempts[key] = attempts.get(key, 0)+1
            p = {**t['anchor'], kind: float(controls[ii]), 'm1': float(slopes[jj])}
            try:
                row = advance(found[i, j], p, t['reference_peak_count'])
                row.update(index=[ii, jj], origin='continued', parent_index=[i, j])
                found[key] = row
                errors.pop(key, None)
                queue.append(key)
            except (RuntimeError, ValueError, OverflowError, np.linalg.LinAlgError) as exc:
                errors[key] = dict(index=[ii, jj], parameters=p, error=str(exc))
        if time.monotonic()-last_saved > 30 or not queue:
            report = dict(id=family, kind=kind, controls=controls.tolist(), m1=slopes.tolist(),
                nodes=sorted(found.values(), key=lambda r: r['index']), failures=list(errors.values()),
                requested=count*count, complete=len(found) == count*count,
                elapsed_seconds=time.monotonic()-started,
                scope='Locally continued symmetric attracting representatives and their fold/symmetry edges; '
                      'each root has its own checked itinerary. Missing nodes are unresolved, not orbit deaths. '
                      'No exhaustive family census or interval proof is claimed.')
            dump(path, report)
            last_saved = time.monotonic()
            print(family, kind, len(found), '/', count*count, 'unresolved', len(errors), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family', choices=['P5', 'P7'], required=True)
    parser.add_argument('--kind', choices=['alpha', 'beta'], required=True)
    parser.add_argument('--nodes', type=int, default=41)
    args = parser.parse_args()
    if args.nodes < 21 or (args.nodes-1) % 20:
        parser.error('Use 21, 41, 61, ... nodes to include the original mesh')
    run(args.family, args.kind, args.nodes)
