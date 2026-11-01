"""Locate an internal diode contact and compare a Taylor prediction with cycles.

The old outer flight develops a minimum at x=1. A signed half-return multiplier
of -0.43 selects the base contact. The successful comparison uses a fixed smooth
parameter line, after an unsuccessful constant-multiplier experiment. Neither
the revised comparison nor the floating-point contact is a rigorous proof.
"""
import json
import time
import numpy as np
from scipy.optimize import root, brentq
from refine_chua_slope_atlases import (ROOT, OUT, CoordinateFlights, parameters,
                                      attach, checked, boundary, representative, flights, dump)
from check_chua_slope_refinement import independent_segmented


def touch_system(w, p, template, target, details=False):
    q, m0, m1, dt = w[:-3], w[-3], w[-2], w[-1]
    f = CoordinateFlights({**p, 'm1': m1}, template['regions'], template['targets'])
    equations = f.boundary_system(np.r_[q, m0], target)
    mats = f.matrices(p['beta'])
    state = np.r_[0., q[:2]]
    # Contact is the second internal minimum of the first outer flight.
    for k, duration in [(0, q[2]), (1, dt)]:
        A, c, eig, V, W = mats[template['regions'][k]]
        state = c+((V*np.exp(eig*duration))@W).real@(state-c)
    velocity = A@(state-c)
    residual = np.r_[equations, state[0]-1, velocity[0]]
    if details:
        return f, dict(state=state.tolist(), acceleration=float((A@velocity)[0]),
                      contact_time_within_outer_flight=float(dt), residual=float(max(abs(residual))))
    return residual


def locate(p, template, target, seed=None):
    if seed is None:
        previous = json.loads((ROOT/'results/chua-slope-families/P7-alpha.json').read_text())['nodes']
        row = min(previous, key=lambda r: abs(r['representative']['alpha']-9)+abs(r['representative']['m1']+.655))
        guess = row['fold' if target == 1 else 'symmetry' if target == -1 else 'representative']
        f = CoordinateFlights({**p, 'm1': -.655}, template['regions'], template['targets'])
        m0, q = f.boundary(guess['m0'], np.array(guess['q']), target)
        seed = np.r_[q, m0, -.655, 4.03]
    sol = root(lambda w: touch_system(w, p, template, target), seed, tol=2e-11)
    f, info = touch_system(sol.x, p, template, target, True)
    if info['residual'] > 2e-9:
        raise RuntimeError('Contact equations did not close')
    q, m0, m1, dt = sol.x[:-3], sol.x[-3], sol.x[-2], sol.x[-1]
    if not 0 < dt < q[3] or info['acceleration'] <= 0:
        raise RuntimeError('Contact is not an interior quadratic minimum')
    r = checked(f, m0, q)
    return dict(target=target, orbit=r, **info, solution=sol.x.tolist())


def extended_minimum(p, template, contact, delta):
    r = contact['orbit']
    f = CoordinateFlights({**p, 'm1': r['m1']+delta}, template['regions'], template['targets'])
    m0, q = f.boundary(r['m0'], np.array(r['q']), contact['target'])
    record = f.record(m0, q, True)
    arc = record['arcs'][1]
    A, c, eig, V, W = f.matrices(p['beta'])[1]
    coeff = V[0]*(W@(np.array(arc['start'])-c))
    dx = lambda t: float(np.real(np.sum(eig*coeff*np.exp(eig*t))))
    tc = contact['contact_time_within_outer_flight']
    dt = brentq(dx, tc-.15, tc+.15, xtol=1e-14)
    minimum = float(c[0]+np.real(np.sum(coeff*np.exp(eig*dt))))
    return minimum


def main():
    OUT.mkdir(exist_ok=True)
    t = json.loads((ROOT/'results/chua-slope-families/anchors.json').read_text())[1]
    p = t['anchor']
    contacts = [locate(p, t, target) for target in [1., -.43, -1.]]
    contact = contacts[1]
    derivatives = []
    for h in [1e-4, 3e-5, 1e-5]:
        derivative = (extended_minimum(p, t, contact, h)-extended_minimum(p, t, contact, -h))/(2*h)
        derivatives.append(dict(step=h, minimum_derivative=derivative))
    opening = derivatives[-1]['minimum_derivative']
    coefficient = 2*np.sqrt(-2*opening/contact['acceleration'])
    forecast = dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        constraint='Signed half-return multiplier = -0.43; alpha=9, beta=100/7.',
        contact=contact, derivative_checks=derivatives,
        predicted_new_central_flight_coefficient=float(coefficient), exponent=.5,
        prediction='For positive delta=m1-m1_contact, a new central flight has duration '
                   'coefficient*sqrt(delta)+o(sqrt(delta)); full-period breakpoint crossings '
                   'increase from 8 to 12 while the same attracting symmetric P7 cycle persists.',
        scope='Numerical Taylor prediction. The base contact is solved before searching post-contact cycles.')
    forecast_path = OUT/'grazing-prediction.json'
    if forecast_path.exists():
        frozen = json.loads(forecast_path.read_text())
        if abs(frozen['predicted_new_central_flight_coefficient']-coefficient) > 1e-7:
            raise RuntimeError('Existing forecast differs; preserve it and investigate')
    else:
        dump(forecast_path, forecast)
    print('Original forecast retained', contact['orbit']['m1'], coefficient, flush=True)
    # A preliminary search on the constant-multiplier constraint did not
    # approach the predicted coefficient. Keep that unsuccessful control.
    original_controls = []
    for delta in [1e-3, 1e-4, 1e-5, 3e-6, 1e-6]:
        try:
            r = boundary(contact['orbit'], {**p, 'm1': contact['orbit']['m1']+delta}, -.43)
            predicted = coefficient*np.sqrt(delta)
            original_controls.append(dict(delta=delta, orbit=r,
                ratio=float(r['q'][4]/predicted), independent=independent_segmented(r)))
        except RuntimeError as exc:
            original_controls.append(dict(delta=delta, error=str(exc)))
    dump(OUT/'grazing-constraint-control.json', dict(
        status='prediction_not_validated_for_this_constraint', controls=original_controls,
        note='The initial protocol fixed a Floquet multiplier across contact. Its square-root '
             'coefficient was not validated. The revised experiment follows a fixed smooth '
             'parameter line; it is not a blind or preregistered forecast.'))
    r0 = contact['orbit']
    slopes = []
    for h in [1e-4, 3e-5, 1e-5]:
        values = []
        for sign in [-1, 1]:
            f = flights({**parameters(r0), 'm1': r0['m1']+sign*h}, r0)
            values.append(f.boundary(r0['m0'], np.array(r0['q']), -.43)[0])
        slopes.append(dict(step=h, dm0_dm1=(values[1]-values[0])/(2*h)))
    slope = slopes[-1]['dm0_dm1']
    dump(OUT/'grazing-linear-prediction.json', dict(contact=contact, derivative_checks=derivatives,
        slope_checks=slopes, line_slope=slope, coefficient=float(coefficient), exponent=.5,
        protocol='m1=m1_contact+delta; m0=m0_contact+line_slope*delta; alpha=9; beta=100/7.',
        status='Analytic prediction compared after preliminary numerical exploration; not a blind forecast.',
        original_protocol='Retained in grazing-prediction.json and grazing-constraint-control.json.'))
    witnesses = []
    previous = contact['orbit']
    for delta in [-1e-4, -1e-5, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6]:
        pp = {**p, 'm1': r0['m1']+delta, 'm0': r0['m0']+slope*delta}
        r = representative(previous, pp)
        previous = r
        check = independent_segmented(r)
        crossing_count = 2*len(r['targets'])
        new_duration = r['q'][4] if crossing_count == 12 else None
        predicted = float(coefficient*np.sqrt(delta)) if delta > 0 else None
        witnesses.append(dict(delta=delta, orbit=r, crossings=crossing_count,
            new_central_flight=new_duration, predicted_duration=predicted,
            relative_error=new_duration/predicted-1 if predicted and new_duration else None,
            independent=check))
        print('Search', delta, crossing_count, new_duration, check['closure'], flush=True)
    # A second integrator checks one orbit on each side of contact.
    radau = [dict(delta=w['delta'], **independent_segmented(w['orbit'], 'Radau'))
             for w in [witnesses[0], witnesses[4], witnesses[-1]]]
    curves = []
    for kind, span in [('alpha', .4), ('beta', .8)]:
        for target in [1., -.43, -1.]:
            center = next(c for c in contacts if c['target'] == target)
            rows = [dict(control=p[kind], contact=center)]
            for sign in [-1, 1]:
                seed = np.array(center['solution'])
                for offset in np.linspace(span/8, span, 8):
                    value = p[kind]+sign*offset
                    result = locate({**p, kind: value}, t, target, seed)
                    seed = np.array(result['solution'])
                    rows.append(dict(control=float(value), contact=result))
            curves.append(dict(kind=kind, target=target, nodes=sorted(rows, key=lambda r: r['control'])))
    gates = dict(contacts_close=max(c['residual'] for c in contacts) < 2e-9,
        nondegenerate_opening=opening < -.01 and contact['acceleration'] > .01,
        independent_closure=all(w['independent']['success'] and w['independent']['closure'] < 2e-7 for w in witnesses),
        independent_multipliers=all(w['independent']['multiplier_error'] < 2e-5 for w in witnesses),
        radau_controls=all(r['success'] and r['closure'] < 2e-7 and r['multiplier_error'] < 2e-5 for r in radau),
        predicted_crossings=all(w['crossings'] == (8 if w['delta'] < 0 else 12) for w in witnesses),
        same_seven_peak_attractor=all(w['orbit']['full_cycle_peaks'] == 7 and w['orbit']['spectral_radius'] < .2 for w in witnesses),
        square_root_limit=abs(witnesses[-1]['relative_error']) < .01)
    result = dict(status='passed' if all(gates.values()) else 'failed', gates=gates,
                  contacts=contacts, witnesses=witnesses, contact_curves=curves, radau_controls=radau,
                  scope='Floating-point contacts and independent ODE controls; not an interval proof or universal grazing classification.')
    dump(OUT/'grazing-results.json', result)
    print(result['status'], gates, flush=True)
    if not all(gates.values()):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
