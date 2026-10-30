"""Independent controls, bifurcation witnesses and overlap tests for slope families."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from audit_chua_census import independent
from chua_layer_bifurcations import from_orbit
from continue_chua_slope_families import CoordinateFlights, dump, extrema_check, valid
from trace_chua_window import Chua


def verify(r):
    check = independent(r)
    if not (check['success'] and check['closure'] < 2e-7 and check['multiplier_error'] < 2e-5):
        raise RuntimeError('Independent closure or multiplier failed: '+str(check))
    return check


def family(template, parameters, symmetric=True, regions=None, targets=None):
    return CoordinateFlights(parameters, regions or template['regions'], targets or template['targets'], symmetric=symmetric)


def boundaries(out):
    records = []
    for t in json.loads((out/'anchors.json').read_text()):
        item = dict(id=t['id'], fold_partners=[], asymmetric_partners=[], fold_coefficients=[])
        r = t['fold']; f = family(t, r); q = np.array(r['q'])
        U, S, V = np.linalg.svd(f.evaluate(r['m0'], q)[1]); left, right = U[:, -1], V[-1]
        item['fold_singular_values'] = S.tolist()
        for h in [1e-4, 3e-5, 1e-5]:
            residual = f.evaluate(r['m0'], q)[0]
            parameter = (f.evaluate(r['m0']+h, q)[0]-f.evaluate(r['m0']-h, q)[0])/(2*h)
            quadratic = (f.evaluate(r['m0'], q+h*right)[0]-2*residual+f.evaluate(r['m0'], q-h*right)[0])/h**2
            item['fold_coefficients'].append(dict(step=h, parameter=float(left@parameter), quadratic=float(left@quadratic)))
        for sign in [-1, 1]:
            m0 = r['m0']-.02*t['width_m0']
            v = f.solve(m0, q+sign*.01*right)
            a = f.record(m0, v, True); valid(a)
            a.update(extrema_check(f, a)); a['independent'] = verify(a)
            item['fold_partners'].append(a)
        radii = [a['spectral_radius'] for a in item['fold_partners']]
        assert min(radii) < 1 < max(radii)
        assert min(abs(c['parameter']) for c in item['fold_coefficients']) > 1e-4
        assert min(abs(c['quadratic']) for c in item['fold_coefficients']) > 1e-4
        r = t['symmetry']; model = Chua(alpha=r['alpha'], m0=r['m0'], m1=r['m1'])
        full, q = from_orbit(model, r['beta'], [*r['q'][:2], r['period']], False)
        f = family(t, r, False, full.regions, full.targets)
        direction = np.linalg.svd(f.evaluate(r['m0'], q)[1])[2][-1]
        m0 = r['m0']-.1*t['width_m0']
        for sign in [-1, 1]:
            v = f.solve(m0, q+sign*.05*direction)
            a = f.record(m0, v, True); valid(a); a['independent'] = verify(a)
            model = Chua(alpha=r['alpha'], m0=m0, m1=r['m1'], rtol=2e-12, step=.02)
            full_record = model.record(r['beta'], [*v[:2], a['period']])
            a['half_period_symmetry_error'] = full_record['half_period_symmetry_error']
            a['full_cycle_peaks'] = len(full_record['maxima'])
            assert a['spectral_radius'] < 1 and a['half_period_symmetry_error'] > 1e-3
            item['asymmetric_partners'].append(a)
        # Verify inversion pairing allowing a phase shift to the other positive section.
        a, b = item['asymmetric_partners']; model = Chua(alpha=a['alpha'], m0=m0, m1=a['m1'], rtol=2e-12, step=.02)
        sol = model.integrate(a['beta'], [0., *a['q'][:2]], a['period'])
        def negative_section(time, state): return state[0]
        negative_section.direction = -1
        cross = solve_ivp(lambda time, state: model.field(a['beta'], state), [0, a['period']],
                          [0., *a['q'][:2]], events=negative_section, rtol=2e-12, atol=2e-13, max_step=.02)
        candidates = cross.y_events[0]
        item['inversion_pair_error'] = float(min(np.max(abs(-v-[0., *b['q'][:2]])) for v in candidates))
        item['section_separation'] = float(np.linalg.norm(np.array(a['q'][:2])-b['q'][:2]))
        assert item['inversion_pair_error'] < 2e-6 and item['section_separation'] > 1e-4
        half = family(t, r); v = half.solve(m0, np.array(r['q']))
        parent = half.record(m0, v, True); valid(parent); parent['independent'] = verify(parent)
        assert parent['spectral_radius'] > 1
        item['unstable_symmetric_parent'] = parent
        records.append(item)
        print('bifurcations', t['id'], 'fold pair', radii, 'inversion error', item['inversion_pair_error'], flush=True)
    dump(out/'bifurcations.json', dict(records=records,
        scope='Numerical evidence for a generic saddle-node and symmetry-breaking pitchfork at the central section; not interval proofs.'))


def overlaps(out):
    templates = {t['id']:t for t in json.loads((out/'anchors.json').read_text())}
    records = []
    for kind in ['alpha', 'beta']:
        maps = [{tuple(r['index']):r for r in json.loads((out/(name+'-'+kind+'.json')).read_text())['nodes']} for name in ['P5','P7']]
        for index in sorted(maps[0].keys() & maps[1].keys()):
            low = max(d[index]['symmetry']['m0'] for d in maps)
            high = min(d[index]['fold']['m0'] for d in maps)
            if low >= high: continue
            parameters = {k:maps[0][index]['core'][k] for k in ['alpha','beta','m0','m1']}
            parameters['m0'] = (low+high)/2
            item = dict(kind=kind, index=list(index), interval_m0=[low, high], parameters=parameters, cycles=[])
            for name, rows in zip(['P5','P7'], maps):
                f = family(templates[name], parameters)
                q = f.solve(parameters['m0'], np.array(rows[index]['core']['q']))
                r = f.record(parameters['m0'], q, True); valid(r); r.update(extrema_check(f,r))
                assert r['spectral_radius'] < .9999 and r['full_cycle_peaks'] == templates[name]['reference_peak_count']
                r['id'] = name; r['independent'] = verify(r)
                item['cycles'].append(r)
            records.append(item)
            print('coexistence', kind, index, 'width', high-low, flush=True)
    assert records
    dump(out/'overlaps.json', dict(records=records,
        scope='Independent closure and attraction checked for both distinct cycles at identical parameters at every listed overlap witness. '
              'Overlap interval endpoints come from continuation; no whole-interval or whole-surface certification.'))
    witness = next(r for r in records if r['kind']=='alpha' and r['index']==[7,10])
    traces = {}
    for r in witness['cycles']:
        model = Chua(alpha=r['alpha'], m0=r['m0'], m1=r['m1'], rtol=2e-12, step=.02)
        time = np.linspace(0, r['period'], 2001)
        sol = model.integrate(r['beta'], [0., *r['q'][:2]], r['period'])
        traces[r['id']+'_time'] = time; traces[r['id']+'_state'] = sol.sol(time)
    np.savez_compressed(out/'coexisting-traces.npz', **traces)
    dump(out/'coexistence-witness.json', witness)


def sampled(out):
    checks = []
    for name in ['P5', 'P7']:
        for kind in ['alpha', 'beta']:
            data = json.loads((out/(name+'-'+kind+'.json')).read_text()); nodes=data['nodes']
            picks = {min(range(len(nodes)), key=lambda k:sum((i-10)**2 for i in nodes[k]['index']))}
            for coordinate in ['m1', kind]:
                picks.add(min(range(len(nodes)), key=lambda k:nodes[k]['core'][coordinate]))
                picks.add(max(range(len(nodes)), key=lambda k:nodes[k]['core'][coordinate]))
            for k in sorted(picks):
                for role in ['fold', 'symmetry', 'core']:
                    r = nodes[k][role]
                    checks.append(dict(id=name, kind=kind, index=nodes[k]['index'], role=role, **verify(r)))
            print('independent sheet samples', name, kind, len(picks), flush=True)
    dump(out/'independent-checks.json', dict(status='passed', checks=checks,
        scope='Selected full-field DOP853 state-and-variational controls; not all continuation nodes.'))


def itinerary_change(out):
    """Follow a real orbit beyond the fixed-itinerary sheet's missing edge."""
    template = next(t for t in json.loads((out/'anchors.json').read_text()) if t['id']=='P7')
    data = json.loads((out/'P7-alpha.json').read_text())
    source = next(r for r in data['nodes'] if r['index']==[10,16])
    previous = {key:source[key] for key in ['fold','symmetry']}
    guess = [*source['core']['q'][:2], source['core']['period']]
    records = []
    for m1 in np.linspace(source['core']['m1'], -5/7+.06, 13):
        params = dict(template['anchor']); params['m1'] = float(m1)
        f = family(template, params); algebraic = {}
        # Past grazing these affine roots supply only a parameter predictor.
        # They are explicitly not accepted as physical stability boundaries.
        for name, target in [('fold',1),('symmetry',-1)]:
            m0, q = f.boundary(previous[name]['m0'], np.array(previous[name]['q']), target)
            algebraic[name] = f.record(m0, q, True)
        previous = algebraic
        m0 = (algebraic['fold']['m0']+algebraic['symmetry']['m0'])/2
        model = Chua(alpha=9., m0=m0, m1=float(m1), rtol=2e-12, step=.02)
        guess = model.shoot(100/7, guess)
        r = model.record(100/7, guess)
        control = verify({**r, 'q':guess[:2].tolist(), 'full_multipliers':r['multipliers']})
        rho = max(abs(complex(*v)) for v in r['multipliers'])
        assert rho < .9999 and len(r['maxima']) == 7 and r['half_period_symmetry_error'] < 2e-6
        records.append(dict(cycle=r, spectral_radius=rho, independent=control,
            predictor_edge_violation=max(e['region_violation'] for e in algebraic.values())))
        print('actual orbit', m1, 'crossings', r['breakpoint_crossings'], 'rho', rho, flush=True)
    assert records[0]['cycle']['breakpoint_crossings']==8 and records[-1]['cycle']['breakpoint_crossings']==12
    dump(out/'itinerary-change.json', dict(records=records,
        scope='The full-field seven-peak attracting cycle persists past an invalid fixed-itinerary edge and changes from 8 to 12 breakpoint crossings. '
              'The transition is bracketed between sampled witnesses, not solved as an exact grazing bifurcation. '
              'Invalid affine boundary roots are used only to predict m0, never as verified edges.'))


def audit(out):
    summary = []; max_residual = 0.; max_violation = 0.
    for name in ['P5','P7']:
        for kind in ['alpha','beta']:
            d=json.loads((out/(name+'-'+kind+'.json')).read_text()); indices=[]
            for row in d['nodes']:
                i,j=row['index']; indices.append((i,j))
                assert row['width_m0'] > 0 and row['symmetry']['m0'] < row['core']['m0'] < row['fold']['m0']
                assert row['core']['full_cycle_peaks']==int(name[1:]) and row['core']['spectral_radius'] < .9999
                for key in ['fold','symmetry','core']:
                    r=row[key];valid(r)
                    assert r[kind]==d['controls'][i] and r['m1']==d['m1'][j]
                    assert r['extrema_region_violation'] < 1e-7
                    max_residual=max(max_residual,r['residual']);max_violation=max(max_violation,r['extrema_region_violation'])
                for key,target in [('fold',1),('symmetry',-1)]:
                    r=row[key]
                    assert abs(1-target*r['map_trace']+r['map_determinant_liouville']) < 1e-8
            assert len(set(indices))==len(indices)
            summary.append(dict(id=name,kind=kind,continued=len(indices),attempted_failures=len(d['failures']),
                unvisited=441-len(indices)-len(d['failures']),width_m0_range=[min(r['width_m0'] for r in d['nodes']),max(r['width_m0'] for r in d['nodes'])]))
    overlaps=json.loads((out/'overlaps.json').read_text())['records']
    for witness in overlaps:
        assert witness['interval_m0'][0] < witness['parameters']['m0'] < witness['interval_m0'][1]
        assert {r['full_cycle_peaks'] for r in witness['cycles']}=={5,7}
        for r in witness['cycles']:
            assert all(r[k]==v for k,v in witness['parameters'].items())
            assert r['spectral_radius']<1 and r['independent']['success']
    checks=json.loads((out/'independent-checks.json').read_text())['checks']
    itinerary=json.loads((out/'itinerary-change.json').read_text())['records']
    assert all(len(r['cycle']['maxima'])==7 and r['spectral_radius']<1 for r in itinerary)
    basin=np.load(out/'basin.npz');ref=np.load(out/'basin-refinement.npz');time=np.load(out/'basin-time-only.npz')
    assert np.array_equal(basin['initial'][ref['indices']],ref['initial'])
    assert np.array_equal(ref['indices'],time['indices'])
    old=basin['labels'][ref['indices']];resolved=old!=0
    report=dict(status='passed',sheets=summary,total_continued=sum(r['continued'] for r in summary),
        maximum_affine_residual=max_residual,maximum_extrema_region_violation=max_violation,
        independent_sheet_controls=len(checks),maximum_independent_sheet_closure=max(r['closure'] for r in checks),
        maximum_independent_sheet_multiplier_error=max(r['multiplier_error'] for r in checks),
        coexistence_witnesses=len(overlaps),overlap_width_range=[min(np.ptp(r['interval_m0']) for r in overlaps),max(np.ptp(r['interval_m0']) for r in overlaps)],
        basin_sampled_resolved=int(resolved.sum()),basin_resolved_changes_time_only=int(np.sum(old[resolved]!=time['labels'][resolved])),
        basin_resolved_changes_half_step=int(np.sum(old[resolved]!=ref['labels'][resolved])),
        scope='Numerical consistency and witness checks; no global layer count, complete basin certification or exact bifurcation theorem.')
    dump(out/'audit.json',report);print(json.dumps(report),flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--stage', choices=['boundaries','overlaps','sampled','itinerary','audit'], required=True)
    a=p.parse_args(); {'boundaries':boundaries,'overlaps':overlaps,'sampled':sampled,'itinerary':itinerary_change,'audit':audit}[a.stage](a.output)


if __name__ == '__main__': main()
