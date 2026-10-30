"""Adaptive FP64 screening of the exact published constrained N-row model.

State and tangent both use the reduced (c,s,X[0:N-1]) coordinates. Physical
time is tau; computational s=min(1,rho)*tau, and exponents are converted back.
Source-like seeded starts obey the exact rate/density invariant disks. Broad
rho/r grids are mathematical exploration, not source apparatus calibration.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp

from motor_rows import dimension, field_jac, rhs, source_parameters, source_start, unpack

ROOT = Path(__file__).resolve().parents[1]


def default_tangent(n):
    vector = np.sin(np.arange(dimension(n))+1.)
    return vector/np.linalg.norm(vector)


class MotorRowsGPU:
    def __init__(self, n):
        self.n, self.d = int(n), dimension(n)
        devices = [d for platform in cl.get_platforms() for d in platform.get_devices()
                   if d.type & cl.device_type.GPU and d.double_fp_config]
        if not devices:
            raise RuntimeError('An FP64 GPU is required')
        self.device = devices[0]
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context)
        program = cl.Program(self.context, Path(__file__).with_suffix('.cl').read_text()).build(options=[f'-DN={self.n}'])
        self.kernel = cl.Kernel(program, 'advance')

    def scan(self, parameters, starts, transient=100., duration=300., rtol=1e-8,
             step_factor=1., max_steps=8000000, max_wall=200., tangents=None, time_mode='scaled'):
        params = np.ascontiguousarray(parameters, dtype=float)
        initial = np.ascontiguousarray(starts, dtype=float)
        lanes, d, n = len(params), self.d, self.n
        if params.shape != (lanes, 3) or initial.shape != (lanes, d) or lanes < 1:
            raise ValueError('Invalid parameter/state shape')
        if not np.isfinite(params).all() or np.any(params[:, :2] <= 0) or not np.isfinite(initial).all():
            raise ValueError('Nonfinite inputs or nonpositive rho/r')
        if transient < 0 or min(duration, rtol, step_factor, max_steps, max_wall) <= 0:
            raise ValueError('Invalid numerical settings')
        if time_mode not in ('scaled', 'original'):
            raise ValueError('Unknown time mode')
        p = np.zeros((lanes, 5)); p[:, :3] = params
        p[:, 3] = np.minimum(1., params[:, 0]) if time_mode == 'scaled' else 1.
        full_x = np.column_stack((initial[:, 2*n:], -initial[:, 2*n:].sum(axis=1)))
        p[:, 4] = np.maximum(abs(full_x).max(axis=1), 2*params[:, 1]*np.sqrt(1+params[:, 2]**2)/params[:, 0])
        initial_radius = (initial[:, :n]**2+initial[:, n:2*n]**2).max(axis=1)
        if np.any(initial_radius > 1+params[:, 2]**2+1e-10):
            raise ValueError('Initial state lies outside the declared invariant disk')
        state = np.zeros((lanes, 2*d+1)); state[:, :d] = initial
        if tangents is None:
            state[:, d:2*d] = default_tangent(n)
        else:
            tangent = np.asarray(tangents, dtype=float)
            if tangent.shape != (lanes, d) or not np.isfinite(tangent).all() or np.any(np.linalg.norm(tangent, axis=1) == 0):
                raise ValueError('Invalid carried tangent')
            state[:, d:2*d] = tangent/np.linalg.norm(tangent, axis=1)[:, None]
        clock = np.zeros((lanes, 4)); clock[:, 1] = np.minimum(.001, .1/(1+params[:, :2].sum(axis=1)))*p[:, 3]
        pairs = n*(n-1)//2
        logs = np.zeros((lanes, 4)); stats = np.zeros((lanes, 4*d+4+3*pairs))
        stats[:, :d] = np.inf; stats[:, d:2*d] = -np.inf
        stats[:, 2*d:3*d] = np.inf; stats[:, 3*d:4*d] = -np.inf
        stats[:, 4*d+4:4*d+4+pairs] = np.inf
        status = np.zeros(lanes, np.int32); counts = np.zeros((lanes, 2), np.int32)
        arrays = [state, p, clock, logs, stats, status, counts]
        buffers = [cl.Buffer(self.context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=a) for a in arrays]
        started = time.perf_counter(); last = started; calls = 0
        while True:
            self.kernel(self.queue, (lanes,), None, *buffers, np.int32(lanes), np.int32(32),
                        np.float64(transient), np.float64(duration), np.float64(rtol), np.float64(rtol*.001),
                        np.float64(step_factor), np.int32(max_steps)).wait()
            calls += 1
            if calls % 8 == 0:
                cl.enqueue_copy(self.queue, clock, buffers[2]).wait()
                cl.enqueue_copy(self.queue, status, buffers[5]).wait()
                active = (status == 0) & (clock[:, 0] < transient+duration)
                if time.perf_counter()-last > 25:
                    print(json.dumps(dict(N=n, mean_scaled_progress=float(np.mean(clock[:, 0])/(transient+duration)),
                                          complete=int(np.sum((status == 0) & ~active)), failed=int(np.sum(status != 0)),
                                          wall=round(time.perf_counter()-started, 1))), flush=True)
                    last = time.perf_counter()
                if not active.any():
                    break
                if time.perf_counter()-started > max_wall:
                    status[active] = 4
                    cl.enqueue_copy(self.queue, buffers[5], status).wait()
                    break
        for array, buffer in zip(arrays, buffers):
            cl.enqueue_copy(self.queue, array, buffer).wait()
        rates = p[:, 3]
        windows = logs/(duration/4.)*rates[:, None]
        windows[status != 0] = np.nan
        return dict(N=n, dimension=d, parameters=params, initial_states=initial,
                    final_states=state[:, :d].copy(), final_tangents=state[:, d:2*d].copy(),
                    lle=windows.mean(axis=1), lle_windows=windows, status=status,
                    elapsed_scaled=clock[:, 0].copy(), elapsed_original=clock[:, 0]/rates,
                    raw_state=state, raw_clock=clock, raw_log_sums=logs, raw_stats=stats, step_counts=counts,
                    observation_min=stats[:, :d].copy(), observation_max=stats[:, d:2*d].copy(),
                    late_min=stats[:, 2*d:3*d].copy(), late_max=stats[:, 3*d:4*d].copy(),
                    maximum_disk_radius_squared=stats[:, 4*d].copy(), maximum_abs_X=stats[:, 4*d+1].copy(),
                    mean_divergence_original=stats[:, 4*d+2]/duration*rates,
                    late_pair_min_distance=np.sqrt(stats[:, 4*d+4:4*d+4+pairs]),
                    late_pair_max_distance=np.sqrt(stats[:, 4*d+4+pairs:4*d+4+2*pairs]),
                    late_pair_rms_distance=np.sqrt(stats[:, 4*d+4+2*pairs:4*d+4+3*pairs]/(duration/4.)),
                    row_distance_definition='Euclidean distance between row states (c_j,s_j,X_j), including eliminated last X; unordered pairs in lexicographic zero-based order',
                    row_distance_quadrature='Last observation quarter only; RMS uses accepted-step right-endpoint quadrature. For incomplete lanes these statistics may cover a partial or absent last window.',
                    row_pairs=[[a, b] for a in range(n) for b in range(a+1, n)],
                    settings=dict(transient_scaled=transient, observation_scaled=duration, rtol=rtol, atol=rtol*.001,
                                  step_factor=step_factor, max_steps_per_lane=max_steps, max_wall_seconds=max_wall,
                                  time_mapping=('s=min(1,rho)*tau; lambda_tau=min(1,rho)*lambda_s' if time_mode == 'scaled' else 's=tau; all requested durations and exponents are in original tau'),
                                  time_mode=time_mode,
                                  state_norm='Euclidean in reduced (c[0:N],s[0:N],X[0:N-1]); last X eliminated exactly',
                                  cap='h_tau <= step_factor*min(.05,.5/conservative_local_Jacobian_row_bound)',
                                  invariants='No clipping; check motor disks and analytical displacement bound on every accepted step'),
                    status_codes={'0': 'completed', '1': 'invariant-domain violation', '2': 'numerical failure',
                                  '3': 'step budget incomplete', '4': 'wall budget incomplete'},
                    device=self.device.name, wall_seconds=time.perf_counter()-started,
                    caveat='Exact published model, exploratory parameters. Finite-time positive exponents are candidates only; incomplete lanes are unresolved.')


def save(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **{k: v for k, v in data.items() if isinstance(v, np.ndarray)})
    metadata = {k: v for k, v in data.items() if not isinstance(v, np.ndarray)}
    valid = data['status'] == 0
    candidate = valid & (data['lle'] > .005) & np.all(data['lle_windows'] > .001, axis=1)
    metadata.update(lanes=len(valid), completed=int(valid.sum()), preliminary_count=int(candidate.sum()),
                    preliminary_ids=np.flatnonzero(candidate).tolist(),
                    preliminary_rule='All four original-tau FTLE windows >.001 and mean >.005; screening rule, not a theorem',
                    largest_completed_lle_ids=np.flatnonzero(valid)[np.argsort(data['lle'][valid])[-10:][::-1]].tolist(),
                    weak_positive_count=int(np.sum(valid & (data['lle'] > 0))),
                    threshold_scope='The screening cutoff is arbitrary; all raw exponents are retained and subthreshold lanes do not prove absence of chaos.',
                    status_counts={str(i): int(np.sum(data['status'] == i)) for i in range(5)},
                    source_canonical_lf_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                                               for p in [ROOT/'research/motor_rows.py', Path(__file__), Path(__file__).with_suffix('.cl')]})
    path.with_suffix('.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k: metadata[k] for k in ['N', 'lanes', 'completed', 'preliminary_count', 'status_counts', 'wall_seconds']}), flush=True)


def validate():
    records = []; regular = []; instances = {}
    for n in (3, 4, 8):
        gpu = instances[n] = MotorRowsGPU(n)
        for rho, r, duration in ((3.5803912736843437, 4.755146047751339, 2.), (.03, 30.9, .1), (10., 330., .03)):
            p = np.array([rho, r, 1.]); initial, _ = source_start(n)
            tangent = default_tangent(n); d = dimension(n)
            def augmented(t, z):
                f, jac = field_jac(z[:d], p, n)
                return np.r_[f, jac@z[d:]]
            references = {}
            for method in ('DOP853', 'Radau'):
                sol = solve_ivp(augmented, (0., duration), np.r_[initial, tangent], method=method,
                                rtol=2e-12, atol=2e-14, max_step=min(.005, duration/100))
                if not sol.success:
                    raise RuntimeError(sol.message)
                end = sol.y[:, -1]; norm = np.linalg.norm(end[d:])
                references[method] = (end[:d], end[d:]/norm, np.log(norm)/duration)
            result = gpu.scan(p[None], initial[None], transient=0, duration=duration*min(1., rho),
                              rtol=1e-10, step_factor=.5, max_wall=30)
            a, b = references['DOP853'], references['Radau']
            row = dict(N=n, parameters=p.tolist(), duration_original=duration, status=int(result['status'][0]),
                       cpu_state_difference=float(max(abs(a[0]-b[0]))), cpu_lle_difference=float(abs(a[2]-b[2])),
                       gpu_state_difference=float(max(abs(a[0]-result['final_states'][0]))),
                       gpu_tangent_difference=float(max(abs(a[1]-result['final_tangents'][0]))),
                       gpu_lle_difference=float(abs(a[2]-result['lle'][0])), cpu_lle=float(a[2]))
            row['passed'] = bool(row['status'] == 0 and row['cpu_state_difference'] < 2e-8
                                 and row['gpu_state_difference'] < 2e-8 and row['gpu_tangent_difference'] < 2e-8
                                 and row['cpu_lle_difference'] < 1e-6*(1+abs(a[2]))
                                 and row['gpu_lle_difference'] < 2e-6*(1+abs(a[2])))
            records.append(row); print(json.dumps(row), flush=True)
    source = json.loads((ROOT/'results/motor-rows/source-audit.json').read_text(encoding='utf-8'))['controls']
    for case in source['cases']:
        sp = case['source_parameters']; p = np.array([sp['rho'], sp['r'], sp['b']])
        target = case['long']['Radau']; duration = target['end_time']
        initial = np.array(source['initial_state'])
        result = instances[8].scan(p[None], initial[None], transient=0., duration=duration*min(1., p[0]),
                                  rtol=1e-10, step_factor=.5, max_wall=60)
        save(ROOT/f"results/motor-rows/gpu-source-{'positive' if sp['delta'] > 0 else 'negative'}-control.npz", result)
        error = float(max(abs(result['final_states'][0]-np.array(target['final_state']))))
        record = dict(delta=sp['delta'], duration_original=duration, max_error_against_independent_Radau=error,
                      status=int(result['status'][0]), maximum_disk_radius_squared=float(result['maximum_disk_radius_squared'][0]),
                      passed=bool(result['status'][0] == 0 and error < 2e-5))
        regular.append(record); print(json.dumps(record), flush=True)
    report = dict(short_controls=records, source_controls=regular,
                  all_passed=all(row['passed'] for row in records+regular),
                  source_control_scope='Independent full-state agreement at tau=300 in published N=8 stable/oscillatory regimes; no chaos classification',
                  tangent_scope='Analytic constrained Jv compared with two independent CPU solvers for N=3,4,8, including stiff exploratory tuples')
    output = ROOT/'results/motor-rows/gpu-validation.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    if not report['all_passed']:
        raise RuntimeError('N-row GPU validation failed; do not launch a survey')


def survey(wall=600.):
    report = json.loads((ROOT/'results/motor-rows/gpu-validation.json').read_text())
    if not report['all_passed']:
        raise RuntimeError('Short and source controls must pass before a survey')
    rho_values, relative_values = np.geomspace(.03, 10., 16), np.geomspace(1.05, 30., 20)
    seeds = [20251030, 20251031, 20251032]
    for n in (3, 4, 8):
        tuples = np.array([[rho, rel*(rho+1), 1.] for rho in rho_values for rel in relative_values])
        parameters = np.repeat(tuples, 3, axis=0)
        starts = np.tile(np.array([source_start(n, seed)[0] for seed in seeds]), (len(tuples), 1))
        result = MotorRowsGPU(n).scan(parameters, starts, max_wall=wall/3.)
        result.update(tuple_id=np.repeat(np.arange(len(tuples)), 3), start_id=np.tile(np.arange(3), len(tuples)),
                      survey=dict(rho_range=[.03, 10.], relative_gain_range=[1.05, 30.], rho_points=16,
                                  relative_gain_points=20, relative_gain_definition='r/(rho+1)', b=1., seeds=seeds,
                                  initial_prescription='Source equilibrium Fourier densities, seeded absolute positions uniform[-.01,.01]nm, fixed endpoints; convert to constrained relative X',
                                  scope='Broader dimensionless family of the exact source equations. Source activity path does not independently vary rho and r.',
                                  geometry='Dimension 3N-1; any three-coordinate plot is a projection. Equal-size two-cluster states reduce to the previously known N=2 motor core.'))
        save(ROOT/f'results/motor-rows/gpu-scout-N{n}.npz', result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['validate', 'survey'])
    parser.add_argument('--wall', type=float, default=600.)
    args = parser.parse_args()
    validate() if args.command == 'validate' else survey(args.wall)
