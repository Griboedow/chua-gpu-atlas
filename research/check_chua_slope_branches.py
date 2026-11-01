"""Repeat branch switches at separated locations in both refined slope volumes."""
import json
from refine_chua_slope_atlases import ROOT, OUT, parameters, advance, dump, flights, checked
from continue_chua_slope_families import CoordinateFlights, extrema_check
from chua_layer_bifurcations import from_orbit
from trace_chua_window import Chua
from check_chua_slope_refinement import independent_segmented
from audit_chua_census import independent
import numpy as np


def half_error(f, r):
    mats = f.matrices(r['beta'])
    state = np.r_[0., r['q'][:2]]
    remaining = r['period']/2
    for duration, region in zip(r['q'][2:], f.regions):
        duration = min(duration, remaining)
        A, c, eig, V, W = mats[region]
        state = c+((V*np.exp(eig*duration))@W).real@(state-c)
        remaining -= duration
        if remaining < 1e-12:
            break
    return float(max(abs(state+np.r_[0., r['q'][:2]])))


def branch_control(t):
    result = dict(id=t['id'], parameters=t['anchor'])
    r = t['fold']; f = flights(parameters(r), r); q = np.array(r['q'])
    direction = np.linalg.svd(f.evaluate(r['m0'], q)[1])[2][-1]
    m0 = r['m0']-.02*t['width_m0']
    for amplitude in [.01, .02, .005, .05]:
        try:
            pair = [checked(f, m0, f.solve(m0, q+sign*amplitude*direction)) for sign in [-1, 1]]
            if min(a['spectral_radius'] for a in pair) < 1 < max(a['spectral_radius'] for a in pair):
                break
        except (RuntimeError, ValueError):
            continue
    else:
        raise RuntimeError('Stable/saddle fold pair not found')
    for a in pair:
        a['independent'] = independent_segmented(a)
    result['fold_partners'] = pair
    r = t['symmetry']
    model = Chua(alpha=r['alpha'], m0=r['m0'], m1=r['m1'])
    full, q = from_orbit(model, r['beta'], [*r['q'][:2], r['period']], False)
    f = CoordinateFlights(parameters(r), full.regions, full.targets, symmetric=False)
    direction = np.linalg.svd(f.evaluate(r['m0'], q)[1])[2][-1]
    attempts = []
    accepted = False
    for fraction in [.1, .03, .01, .003]:
        m0 = r['m0']-fraction*t['width_m0']
        for amplitude in [.05, .02, .1, .01, .005, .2]:
            try:
                pair = []
                for sign in [-1, 1]:
                    v = f.solve(m0, q+sign*amplitude*direction)
                    a = checked(f, m0, v)
                    a['full_cycle_peaks'] //= 2
                    a['half_period_symmetry_error'] = half_error(f, a)
                    pair.append(a)
                radii = [a['spectral_radius'] for a in pair]
                errors = [a['half_period_symmetry_error'] for a in pair]
                separation = float(np.linalg.norm(np.array(pair[0]['q'][:2])-pair[1]['q'][:2]))
                attempts.append(dict(fraction=fraction, amplitude=amplitude, radii=radii,
                                     symmetry_errors=errors, section_separation=separation))
                if max(radii) < 1 and min(errors) > 1e-3 and separation > 1e-4:
                    accepted = True
                    break
            except (RuntimeError, ValueError, OverflowError) as exc:
                attempts.append(dict(fraction=fraction, amplitude=amplitude, error=str(exc)))
        if accepted:
            break
    result['switch_attempts'] = attempts
    if not accepted:
        raise RuntimeError('Asymmetric attracting pair not resolved: '+str(attempts[-1]))
    for a in pair:
        a['independent'] = independent(a)
    # Event states at negative x=0 sections identify the inversion-related orbit.
    a, b = pair
    model = Chua(alpha=a['alpha'], m0=a['m0'], m1=a['m1'], rtol=2e-12, step=.02)
    from scipy.integrate import solve_ivp
    def section(time, state):
        return state[0]
    section.direction = -1
    sol = solve_ivp(lambda time, state: model.field(a['beta'], state), [0, a['period']],
        [0., *a['q'][:2]], events=section, rtol=2e-12, atol=2e-13, max_step=.02)
    result['inversion_pair_error'] = float(min(max(abs(-v-np.r_[0., b['q'][:2]])) for v in sol.y_events[0]))
    parentf = flights({**parameters(r), 'm0': m0}, r)
    parent = checked(parentf, m0, parentf.solve(m0, np.array(r['q'])))
    parent['independent'] = independent_segmented(parent)
    result.update(asymmetric_partners=pair, unstable_symmetric_parent=parent)
    checks = [a['independent'] for a in result['fold_partners']+pair+[parent]]
    if (parent['spectral_radius'] <= 1 or result['inversion_pair_error'] > 2e-6
            or not all(c['success'] and c['closure'] < 2e-7 and c['multiplier_error'] < 2e-5 for c in checks)):
        raise RuntimeError('Independent branch controls failed')
    return result


def main():
    anchors, heldout = [], []
    for family in ['P5', 'P7']:
        for kind in ['alpha', 'beta']:
            d = json.loads((OUT/f'{family}-{kind}.json').read_text())
            nodes = {tuple(r['index']): r for r in d['nodes']}
            # Offset the beta-volume high-m1 check so the common section is not
            # counted twice as two distinct branch-switch locations.
            for index in [(0, 20), (40, 20), (20 if kind == 'alpha' else 24, 40)]:
                row = nodes[index]
                if row['fold']['regions'] != row['symmetry']['regions']:
                    raise RuntimeError('Selected branch control straddles two itineraries')
                anchors.append(dict(id=f'{family}-{kind}-{index[0]}-{index[1]}',
                    anchor=parameters(row['core']), regions=row['core']['regions'],
                    targets=row['core']['targets'], reference_peak_count=int(family[1:]), **row))
            # Unused cell centres check interpolation on the new 41-point grid.
            for i, j in [(4, 8), (8, 24), (16, 30), (20, 33), (24, 35), (28, 12), (34, 26), (38, 38)]:
                corners = [(i, j), (i+1, j), (i, j+1), (i+1, j+1)]
                if not all(k in nodes for k in corners):
                    continue
                p = parameters(nodes[i, j]['core'])
                p[kind] = (d['controls'][i]+d['controls'][i+1])/2
                p['m1'] = (d['m1'][j]+d['m1'][j+1])/2
                exact = advance(nodes[i, j], p, int(family[1:]))
                new_error = max(abs(np.mean([nodes[k][role]['m0'] for k in corners])-exact[role]['m0'])
                                for role in ['fold', 'symmetry'])
                # The old 21-point interpolation uses the containing 2x2 cell.
                a, b = i//2*2, j//2*2
                old_corners = [(a, b), (a+2, b), (a, b+2), (a+2, b+2)]
                u, v = (i+.5-a)/2, (j+.5-b)/2
                weights = [(1-u)*(1-v), u*(1-v), (1-u)*v, u*v]
                old_error = max(abs(sum(w*nodes[k][role]['m0'] for w, k in zip(weights, old_corners))-exact[role]['m0'])
                                for role in ['fold', 'symmetry'])
                heldout.append(dict(family=family, kind=kind, cell=[i, j], parameters=p,
                    corrected=exact, old_edge_error=float(old_error), new_edge_error=float(new_error),
                    new_error_over_width=float(new_error/exact['width_m0'])))
    dump(OUT/'bifurcation-anchors.json', anchors)
    dump(OUT/'heldout-interpolation.json', dict(nodes=heldout,
        scope='32 newly corrected cell centres. Old-grid spacing is compared using corresponding '
              'coarse vertices from the refined mesh, including recovered vertices. No uniform error bound is claimed.'))
    records, failures = [], []
    for t in anchors:
        try:
            record = branch_control(t)
            records.append(record)
            print(t['id'], 'branch switches checked', flush=True)
        except (RuntimeError, ValueError) as exc:
            failures.append(dict(id=t['id'], error=str(exc)))
            print(t['id'], 'unresolved', str(exc), flush=True)
        dump(OUT/'bifurcation-controls.json', dict(records=records, failures=failures, complete=False))
    result = dict(records=records, failures=failures, complete=True,
                  status='passed' if len(records) == len(anchors) else 'incomplete')
    result['scope'] = ('Twelve separated slope-volume locations: numerical stable/saddle fold pairs, '
                       'inversion-related asymmetric attracting pairs, and unstable symmetric parents. '
                       'The controls are not a proof of an unchanging normal form everywhere in the volumes.')
    dump(OUT/'bifurcation-controls.json', result)


if __name__ == '__main__':
    main()
