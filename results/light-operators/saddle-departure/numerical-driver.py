"""Finite-time two-sided departure test at the previously closed optical saddle cycle.

The assumed logistic heating profile is unchanged. Small section perturbations
approximate the local unstable direction; they do not compute an exact global
manifold or establish that a stable manifold is a basin boundary.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from smooth_operators import rhs_and_jac

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/light-operators/cycle-arclength/unstable-w004.json'
STABLE = ROOT / 'results/light-operators/lifetimes/periodic-controls.json'


def run(output, wall=420.):
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'probe.json').exists():
        raise RuntimeError('Use a new output folder; existing evidence is not overwritten')
    begin = time.perf_counter()
    record = json.loads(SOURCE.read_text())
    p = np.array([74.1567393414624, 1.784180421282182,
                  .4139517339970917, .25683679452165964, .04])
    scale = np.array([1., np.sqrt(p[0]), 1., 1., np.sqrt(p[0]), 1.])
    base = np.array(record['state']) / scale
    period = record['period']
    stable_record = next(x for x in json.loads(STABLE.read_text()) if x['parameters'][-1] == .04)
    # The source controls use the same x_minus=delta_minus upward section.
    stable_keys = list(stable_record)
    def budget():
        if time.perf_counter() - begin > wall:
            raise TimeoutError('CPU wall budget reached; partial trajectories retained')
    def field(z):
        f, jac = rhs_and_jac(z * scale, p)
        return f / scale, jac / scale[:, None] * scale[None, :]
    def augmented(t, y):
        budget()
        f, jac = field(y[:6])
        tail = y[6:]
        return np.r_[f, (jac @ tail.reshape(6, -1)).ravel()]
    def section(t, y):
        return y[0] - p[2]
    section.direction = 1
    solutions = {}
    for method in ['DOP853', 'Radau']:
        sol = solve_ivp(augmented, (0, period), np.r_[base, np.eye(6).ravel()],
                        method=method, rtol=3e-12, atol=3e-14, max_step=.003)
        if not sol.success:
            raise RuntimeError(sol.message)
        end = sol.y[:6, -1]
        matrix = sol.y[6:, -1].reshape(6, 6)
        flow = field(end)[0]
        pm = ((np.eye(6) - np.outer(flow, np.eye(6)[0]) / flow[0]) @ matrix)[1:, 1:]
        values, vectors = np.linalg.eig(pm)
        index = int(np.argmax(abs(values)))
        vec = vectors[:, index].real
        vec /= np.linalg.norm(vec)
        if vec[np.argmax(abs(vec))] < 0:
            vec = -vec
        solutions[method] = dict(closure=float(max(abs(end-base))),
            multiplier=float(values[index].real), direction=vec.tolist(),
            multiplier_imag=float(values[index].imag))
    direction = np.r_[0., solutions['DOP853']['direction']]
    if abs(solutions['DOP853']['multiplier']-solutions['Radau']['multiplier']) > 1e-7:
        raise RuntimeError('Independent saddle multiplier check failed')
    # Obtain a comparable stable section state directly from its archived orbit.
    stable_initial = np.array(stable_record['section_state']) / scale
    stable_sol = solve_ivp(lambda t, y: field(y)[0], (0, 40), stable_initial,
                          method='DOP853', rtol=2e-12, atol=2e-14, max_step=.005, events=section)
    if not stable_sol.success or len(stable_sol.y_events[0]) < 3:
        raise RuntimeError('Stable section reference failed')
    stable_section = stable_sol.y_events[0][-1]
    rows = []
    metadata = dict(parameters=p.tolist(), coordinate_scale=scale.tolist(),
        state_columns=['x_minus','v_minus','temperature_minus','x_plus','v_plus','temperature_plus'],
        tangent_norm='Euclidean after scaling velocities by sqrt(k)',
        source_saddle='results/light-operators/cycle-arclength/unstable-w004.json',
        source_stable='results/light-operators/lifetimes/periodic-controls.json',
        source_sha256={q.relative_to(ROOT).as_posix(): hashlib.sha256(q.read_bytes()).hexdigest() for q in [SOURCE, STABLE]},
        solver_controls=solutions, stable_section_balanced=stable_section.tolist(),
        stable_source_keys=stable_keys, rows=rows,
        settings=dict(warm=100., observation=400., renormalization_interval=1.,
                      sample_interval=.2, method='DOP853', rtol=2e-10, atol=2e-12, max_step=.02, wall_budget=wall),
        scope='Four declared local section perturbations, two amplitudes and two signs. Not an exact unstable manifold, global basin survey, asymptotic proof, or novelty claim.')
    for amplitude in [1e-4, 1e-5]:
        for sign in [-1, 1]:
            start = base + sign*amplitude*direction
            tangent = np.array([1., .3, .4, -.2, .5, -.1]); tangent /= np.linalg.norm(tangent)
            z = np.r_[start, tangent]
            logs, times, samples, crossings, crossing_times = [], [], [], [], []
            row = dict(amplitude=amplitude, sign=sign, start_balanced=start.tolist(),
                       completed=False, completed_time=0.)
            try:
                for step in range(500):
                    budget()
                    sol = solve_ivp(augmented, (0, 1), z, method='DOP853',
                                    rtol=2e-10, atol=2e-12, max_step=.02,
                                    t_eval=np.linspace(0, 1, 6)[1:], events=section)
                    if not sol.success:
                        raise RuntimeError(sol.message)
                    z = sol.y[:, -1]
                    norm = np.linalg.norm(z[6:]); logs.append(float(np.log(norm)))
                    z[6:] /= norm
                    times.extend(step+sol.t); samples.extend(sol.y[:6].T*scale)
                    crossing_times.extend(step+sol.t_events[0]); crossings.extend(sol.y_events[0].reshape(-1, 12)[:, :6])
                    row['completed_time'] = float(step+1)
                row['completed'] = True
            except (TimeoutError, RuntimeError) as error:
                row['error'] = str(error)
            arr = np.array(samples)
            cross = np.array(crossings).reshape(-1, 6)
            ct = np.array(crossing_times)
            distance = np.linalg.norm(cross-stable_section, axis=1)
            row.update(crossings=len(cross),
                last_section_distances=distance[-20:].tolist(),
                min_observed_temperature=float(arr[:, [2,5]].min()) if len(arr) else None,
                max_observed_temperature=float(arr[:, [2,5]].max()) if len(arr) else None,
                max_balanced_abs_state=float(np.max(abs(arr/scale))) if len(arr) else None,
                end_balanced=z[:6].tolist(),
                last_return_intervals=np.diff(ct)[-20:].tolist())
            if row['completed']:
                post = np.array(logs[100:])
                row.update(leading_finite_time_mean=float(post.mean()),
                           leading_windows=np.mean(post.reshape(4, 100), axis=1).tolist(),
                           last_20_sections_near_known_cycle=bool(len(distance)>=20 and np.all(distance[-20:]<1e-6)))
            name = f"amplitude-{amplitude:g}-sign-{sign:+d}.npz"
            np.savez_compressed(output/name, times=times, states=arr, loggrowth=logs,
                section_times=ct, section_states_balanced=cross,
                section_distance_to_stable=distance, final_balanced_state_tangent=z)
            row['array_file'] = name
            rows.append(row)
            metadata['wall_seconds'] = time.perf_counter()-begin
            (output/'probe.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            print(json.dumps(row), flush=True)
            if not row['completed']:
                return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'results/light-operators/saddle-departure')
    parser.add_argument('--wall', type=float, default=420.)
    args = parser.parse_args()
    run(args.output, args.wall)
