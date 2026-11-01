"""Independent field integration and branch/interpolation checks of slope meshes."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from refine_chua_slope_atlases import ROOT, OUT, parameters, flights, advance, dump
from trace_chua_window import Chua


def independent_segmented(r, method='DOP853'):
    """Use the actual field on every segment, including very short new flights.

    Segment endpoints are time partitions, not imposed state/region resets.
    The state and variational matrix remain continuous throughout integration.
    """
    model = Chua(alpha=r['alpha'], m0=r['m0'], m1=r['m1'])
    initial = np.r_[0., r['q'][:2]]
    v = np.r_[initial, np.eye(3).ravel()]
    def rhs(t, state):
        return np.r_[model.field(r['beta'], state[:3]),
            (model.jacobian(r['beta'], abs(state[0]) < 1)@state[3:].reshape(3, 3)).ravel()]
    success = True
    for duration in r['q'][2:]*(2 if r['symmetric'] else 1):
        sol = solve_ivp(rhs, [0, duration], v, method=method, rtol=2e-12,
                        atol=2e-13, max_step=min(.02, duration/8))
        success &= sol.success
        v = sol.y[:, -1]
    velocity = model.field(r['beta'], v[:3])
    M = v[3:].reshape(3, 3)
    P = ((np.eye(3)-np.outer(velocity, [1., 0., 0.])/velocity[0])@M)[1:, 1:]
    eig = np.linalg.eigvals(P)
    dominant = eig[np.argmax(abs(eig))]
    expected = complex(*r['full_multipliers'][0])
    return dict(method=method, success=bool(success), closure=float(max(abs(v[:3]-initial))),
                multiplier_error=float(abs(dominant-expected)), dominant=float(dominant.real))


def fold_coefficients(r):
    f = flights(parameters(r), r)
    q = np.array(r['q'])
    m = r['m0']
    U, S, V = np.linalg.svd(f.evaluate(m, q)[1])
    left, direction = U[:, -1], V[-1]
    base = f.evaluate(m, q)[0]
    values = []
    for h in [1e-4, 3e-5, 1e-5]:
        parameter = left@((f.evaluate(m+h, q)[0]-f.evaluate(m-h, q)[0])/(2*h))
        quadratic = left@((f.evaluate(m, q+h*direction)[0]-2*base+f.evaluate(m, q-h*direction)[0])/h**2)
        values.append(dict(step=h, parameter=float(parameter), quadratic=float(quadratic)))
    return dict(singular_values=S.tolist(), coefficients=values)


def main():
    audits, controls, interpolation, folds = [], [], [], []
    for family in ['P5', 'P7']:
        for kind in ['alpha', 'beta']:
            data = json.loads((OUT/f'{family}-{kind}.json').read_text())
            nodes = {tuple(r['index']): r for r in data['nodes']}
            old = json.loads((ROOT/f'results/chua-slope-families/{family}-{kind}.json').read_text())
            old_indices = {tuple(r['index']) for r in old['nodes']}
            recovered = sum((2*i, 2*j) in nodes and (i, j) not in old_indices
                            for i in range(21) for j in range(21))
            allr = [r[role] for r in nodes.values() for role in ['fold', 'symmetry', 'representative']]
            audits.append(dict(family=family, kind=kind, accepted=len(nodes), requested=data['requested'],
                recovered_old_nodes=recovered,
                new_itinerary_representatives=sum(len(r['representative']['targets']) != 4 for r in nodes.values()),
                widths=[min(r['width_m0'] for r in nodes.values()), max(r['width_m0'] for r in nodes.values())],
                maximum_closure=max(r['residual'] for r in allr),
                maximum_region_violation=max(r['extrema_region_violation'] for r in allr)))
            # Include corners, centre, and both sides of the P7 itinerary change.
            picks = [(0, 20), (20, 0), (20, 20), (40, 20), (20, 30), (20, 35), (20, 40)]
            for index in picks:
                if index not in nodes:
                    continue
                row = nodes[index]
                for role in ['fold', 'symmetry', 'representative']:
                    result = independent_segmented(row[role])
                    controls.append(dict(family=family, kind=kind, index=index, role=role, **result))
                if index in [(20, 20), (20, 35), (20, 40)]:
                    folds.append(dict(family=family, kind=kind, index=index, **fold_coefficients(row['fold'])))
            # Old-cell centre: compare old bilinear edge estimates with the new root.
            for i in range(1, 40, 2):
                for j in range(1, 40, 2):
                    corners = [(i-1, j-1), (i-1, j+1), (i+1, j-1), (i+1, j+1)]
                    if (i, j) not in nodes or not all(k in nodes and (k[0]//2, k[1]//2) in old_indices for k in corners):
                        continue
                    row = nodes[i, j]
                    err = max(abs(np.mean([nodes[k][role]['m0'] for k in corners])-row[role]['m0'])
                              for role in ['fold', 'symmetry'])
                    interpolation.append(dict(family=family, kind=kind, index=[i, j],
                        edge_error=float(err), error_over_width=float(err/row['width_m0'])))
            print(family, kind, 'accepted', len(nodes), 'recovered', recovered, flush=True)
    gates = dict(accepted_roots_close=max(a['maximum_closure'] for a in audits) < 2e-9,
        physical_arcs=max(a['maximum_region_violation'] for a in audits) < 2e-8,
        independent_integrations=all(c['success'] and c['closure'] < 2e-7 and c['multiplier_error'] < 2e-5 for c in controls),
        sampled_fold_nondegeneracy=all(min(abs(c['parameter']) for c in f['coefficients']) > 1e-4
            and min(abs(c['quadratic']) for c in f['coefficients']) > 1e-4 for f in folds))
    result = dict(status='passed' if all(gates.values()) else 'failed', gates=gates, meshes=audits,
        independent_controls=controls, sampled_fold_coefficients=folds, old_cell_interpolation=interpolation,
        scope='Selected independent controls and finite-difference fold tests. Mesh interpolation is not certified, '
              'and symmetry labels away from tested branch switches refer to the signed half-map multiplier condition.')
    dump(OUT/'audit.json', result)
    print(result['status'], gates, flush=True)
    if not all(gates.values()):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
