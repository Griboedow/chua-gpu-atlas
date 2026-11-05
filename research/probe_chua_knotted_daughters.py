"""Probe symmetry-breaking daughters of knotted, inversion-symmetric cycles.

An amplitude constraint on the difference of the two residence-time blocks
prevents Newton iteration from returning to the symmetric parent.  Parameter
beta is solved together with the closed full orbit.  Results are numerical.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
import os

ROOT = Path(__file__).resolve().parents[1]
TASK_TEMP=ROOT.parent/'.agent-tmp';TASK_TEMP.mkdir(parents=True,exist_ok=True)
os.environ['TEMP']=os.environ['TMP']=str(TASK_TEMP)
os.environ['PIP_CACHE_DIR']='D:/_TMP/python-cache/pip'

import numpy as np
from scipy.optimize import root

from chua_fast_flights import FastFlights
from trace_chua_window import Chua

OUT = ROOT / 'results/chua-link-search/daughters'


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8', newline='\n')


def source_digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def probe(family, amplitudes=(1e-4, 3e-4, 1e-3, 3e-3)):
    source = ROOT / f'results/chua-layer-atlas/{family}-sheet.json'
    data = json.loads(source.read_text(encoding='utf-8'))
    node = next(n for n in data['nodes'] if n['index'] == [4, 4])
    edge = node['symmetry']
    model = Chua(alpha=edge['alpha'], m0=edge['m0'], m1=edge['m1'],
                 rtol=2e-12, step=.025)
    half = FastFlights(model, data['regions'], data['targets'], True)
    beta0, qhalf = half.boundary(edge['beta'], np.asarray(edge['q']), -1)
    count = len(half.regions)
    full = FastFlights(model, half.regions + [-r for r in half.regions],
                       half.targets + [0] + [-t for t in half.targets], False)
    q0 = np.r_[qhalf[:2], qhalf[2:], qhalf[2:]]
    _, jac, _, _ = full.evaluate(beta0, q0)
    _, singular, right = np.linalg.svd(jac)
    direction = right[-1]
    antisymmetric = direction[2:2+count] - direction[2+count:]
    antisymmetric /= np.linalg.norm(antisymmetric)
    gradient = np.r_[0., 0., antisymmetric, -antisymmetric]
    direction /= gradient @ direction
    rows = []
    print(family, 'boundary', beta0, 'singular values', singular[-2:], flush=True)
    for amplitude in amplitudes:
        started = time.perf_counter()
        row = dict(family=family, amplitude=float(amplitude), parent_beta=float(beta0),
                   parent_period=float(2*sum(qhalf[2:])), source=source.relative_to(ROOT).as_posix(),
                   source_sha256=source_digest(source),
                   node_index=[4, 4], smallest_parent_singular_values=singular[-2:].tolist())
        try:
            def equations(w):
                residual = full.evaluate(w[-1], w[:-1])[0]
                return np.r_[residual, gradient @ (w[:-1] - q0) - amplitude]

            def derivative(w):
                _, jacobian, _, _ = full.evaluate(w[-1], w[:-1])
                step = 1e-6
                beta_column = (full.evaluate(w[-1]+step, w[:-1])[0] -
                               full.evaluate(w[-1]-step, w[:-1])[0]) / (2*step)
                return np.block([[jacobian, beta_column[:, None]],
                                 [gradient[None, :], np.zeros((1, 1))]])

            initial = np.r_[q0 + amplitude*direction, beta0]
            answer = root(equations, initial, jac=derivative, tol=1e-10)
            w = answer.x
            error = float(np.max(np.abs(equations(w))))
            if error > 2e-9 or min(w[2:-1]) <= 0:
                raise RuntimeError(f'Augmented closure {error}: {answer.message}')
            beta, q = float(w[-1]), w[:-1]
            record = full.record(beta, q, True)
            if record['region_violation'] > 1e-7:
                raise RuntimeError(f"Wrong itinerary: {record['region_violation']}")
            period = record['period']
            state = np.r_[0., q[:2]]
            orbit = model.integrate(beta, state, period)
            independent = model.integrate(beta, state, period, method='Radau')
            times = np.linspace(0., period, 4001)
            xyz = orbit.sol(times).T
            opposite_phase = orbit.sol((times + period/2) % period).T
            symmetry_error = float(np.max(np.abs(xyz + opposite_phase)))
            closure = float(np.max(np.abs(independent.y[:, -1] - state)))
            row.update(status='closed', beta=beta, beta_offset=beta-beta0,
                       cycle=record, regions=full.regions, targets=full.targets,
                       amplitude_residual=error, solver_success=bool(answer.success),
                       half_period_symmetry_error=symmetry_error,
                       radau_closure=closure,
                       attracting=record['spectral_radius'] < 1,
                       asymmetric=symmetry_error > 1e-7,
                       scope='Numerical full-period closure and independent Radau check; no interval certificate')
            stem = f'{family}-a{amplitude:.0e}'
            np.savez_compressed(OUT / (stem+'.npz'), xyz=xyz, t=times, partner_xyz=-xyz)
            row['polygon_file'] = stem+'.npz'
            print(family, amplitude, 'beta offset', beta-beta0,
                  'rho', record['spectral_radius'], 'asymmetry', symmetry_error,
                  'Radau closure', closure, flush=True)
        except Exception as error:
            row.update(status='unresolved', error=f'{type(error).__name__}: {error}')
            print(family, amplitude, row['error'], flush=True)
        row['elapsed_seconds'] = time.perf_counter()-started
        rows.append(row)
        save(OUT / (family+'-probe.json'), dict(family=family, rows=rows))
    return rows


def refine_and_analyze(family, amplitude, digits=0, subdivisions=96, seeds=(0,)):
    """Refine at fixed decimal beta and compare the joint and split Jones data."""
    from search_chua_links import split_screen
    from compare_chua_joint_jones import compare
    import mpmath as mp
    from chua_precision_flights import PrecisionFlights
    from reconstruct_chua_knot_orbits import sample

    source = OUT / (family+'-probe.json')
    rows = json.loads(source.read_text(encoding='utf-8'))['rows']
    row = next(r for r in rows if r['amplitude'] == amplitude)
    if row['status'] != 'closed' or not row['attracting'] or not row['asymmetric']:
        print(family, amplitude, 'skipped: no attracting asymmetric root', flush=True)
        return
    label = f'{family}-a{amplitude:.0e}'
    stem = label if not digits else f'{label}-d{digits}-n{subdivisions}'
    if digits:
        target = OUT / (stem+'.json')
        if not target.exists():
            f = PrecisionFlights(row['regions'], row['targets'], False, digits)
            # The central atlas node has exactly alpha=9, m0=-8/7, m1=-5/7.
            beta = mp.mpf(str(row['beta']))
            q = mp.matrix([mp.mpf(str(v)) for v in row['cycle']['q']])
            history = []
            for _ in range(25):
                ev = f.evaluate(beta, list(q), False)
                error = max(abs(v) for v in ev['residual'])
                history.append(mp.nstr(error, 15))
                if error < mp.mpf(10)**(-(digits-20)):
                    break
                q += mp.lu_solve(ev['jacobian'], -ev['residual'])
                assert min(q[2:]) > 0
            else:
                raise RuntimeError('High-precision fixed-parameter closure failed')
            eigen = mp.polyroots([1, -ev['trace'], ev['determinant']], maxsteps=1000)
            radius = max(abs(v) for v in eigen)
            assert radius < 1
            record = dict(id=label, beta=mp.nstr(beta, digits), alpha='9',
                          m0='-8/7', m1='-5/7', gamma='0', symmetric=False,
                          q=[mp.nstr(v, digits-5) for v in q],
                          period=mp.nstr(sum(q[2:]), digits-5),
                          regions=f.regions, targets=f.targets,
                          trace=mp.nstr(ev['trace'], digits-5),
                          determinant=mp.nstr(ev['determinant'], digits-5),
                          spectral_radius=float(radius), residual=mp.nstr(error, 15),
                          working_digits=digits, newton_residuals=history,
                          source=source.name,
                          source_sha256=source_digest(source))
            xyz, times, diagnostic = sample(f, record, ev, subdivisions)
            record['sampling'] = diagnostic
            assert diagnostic['maximum_sampled_region_violation'] < 1e-25
            np.savez_compressed(OUT / (stem+'.npz'), xyz=xyz, t=times, partner_xyz=-xyz)
            save(target, record)
            print(stem, 'refined', 'rho', float(radius), 'residual', mp.nstr(error, 8), flush=True)
    xyz = np.load(OUT / (stem+'.npz'))['xyz']
    for seed in seeds:
        joint_target = OUT / f'{stem}-joint-s{seed}.json'
        jones_target = OUT / f'{stem}-jones-s{seed}.json'
        if joint_target.exists():
            joint = json.loads(joint_target.read_text(encoding='utf-8'))
        else:
            joint = split_screen([xyz, -xyz], seed)
            joint.update(ids=[label, label+'-inverted'], polygon_file=stem+'.npz',
                         parameters=dict(alpha=9, m0='-8/7', m1='-5/7', beta=str(row['beta'])),
                         digits=digits or 'double', subdivisions=subdivisions if digits else None)
            save(joint_target, joint)
        if jones_target.exists():
            comparison = json.loads(jones_target.read_text(encoding='utf-8'))
        else:
            comparison = compare(joint, [xyz, -xyz])
            comparison.update(ids=joint['ids'], joint_record=joint_target.name)
            save(jones_target, comparison)
        print(stem, seed, joint['crossings'], '->', joint['simplified_crossings'],
              joint['status'], comparison['conclusion'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--families', nargs='+', default=['S01', 'S02'])
    parser.add_argument('--amplitudes', nargs='+', type=float, default=[1e-4, 3e-4, 1e-3, 3e-3])
    parser.add_argument('--stage', choices=['probe', 'analyze'], default='probe')
    parser.add_argument('--digits', type=int, default=0)
    parser.add_argument('--subdivisions', type=int, default=96)
    parser.add_argument('--seeds', type=int, nargs='+', default=[0])
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for family in args.families:
        if args.stage == 'probe':
            probe(family, args.amplitudes)
        else:
            for amplitude in args.amplitudes:
                refine_and_analyze(family, amplitude, args.digits, args.subdivisions, args.seeds)
