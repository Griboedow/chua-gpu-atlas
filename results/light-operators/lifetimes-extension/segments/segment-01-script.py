"""Continue the saved optical cohort without resampling or changing its model.

The original data are read-only. Every completed 500-unit block is checkpointed
atomically in a new directory, with source hashes and the carried tangent. A wall
budget stops between blocks whenever possible; a partial block is retained but
never counted as a completed observation interval. No escape law is fitted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from smooth_operators import SmoothOperatorsGPU, rhs_and_jac

SOURCE = Path('results/light-operators/lifetimes')
TARGET = Path('results/light-operators/lifetimes-extension')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_npz(path, arrays):
    temporary = path.with_suffix('.temporary.npz')
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def atomic_json(path, value):
    temporary = path.with_suffix('.temporary.json')
    temporary.write_text(json.dumps(value, indent=2), encoding='utf-8')
    temporary.replace(path)


def run(target_time=100000., wall_budget=515., resume=False):
    TARGET.mkdir(parents=True, exist_ok=True)
    checkpoint = TARGET / 'ensemble.npz'
    metadata_path = TARGET / 'ensemble.json'
    if checkpoint.exists() and not resume:
        raise RuntimeError('Extension exists: use --resume; never overwrite a completed run silently')
    source = np.load(SOURCE / 'ensemble.npz')
    controls = json.loads((SOURCE / 'periodic-controls.json').read_text())
    paths = [SOURCE / name for name in ['ensemble.npz', 'ensemble.json', 'periodic-controls.json', 'capture-cpu-audit.json', 'chunk-validation.json']]
    paths += [Path('research') / name for name in ['smooth_operators.py', 'smooth_operators.cl', 'extend_optical_lifetimes.py']]
    hashes = {p.as_posix(): digest(p) for p in paths}
    if resume:
        metadata = json.loads(metadata_path.read_text())
        if metadata['source_hashes'][str(SOURCE / 'ensemble.npz').replace('\\', '/')] != hashes[(SOURCE / 'ensemble.npz').as_posix()]:
            raise RuntimeError('Original cohort hash changed')
        data = {key: array.copy() for key, array in np.load(checkpoint).items()}
        if np.any(data['status'] > 1):
            raise RuntimeError('Partial or failed lanes require an explicit separate recovery audit')
    else:
        data = {key: array.copy() for key, array in source.items()}
        n = len(data['parameters'])
        data['cohort_ids'] = np.arange(n)
        data['source_status'] = source['status'].copy()
        data['source_capture_time'] = source['capture_time'].copy()
        data['near_streak'] = np.zeros(n, np.int32)
        data['first_near_time'] = np.full(n, np.nan)
        data['endpoint_time'] = np.full(n, source['checkpoint_times'][-1])
        for index in np.flatnonzero(source['status'] == 1):
            first = np.flatnonzero(source['status_history'][:, index] == 1)[0]
            data['endpoint_time'][index] = source['checkpoint_times'][first]
        # Reconstruct pending confirmations rather than assuming they are absent.
        for index in np.flatnonzero(source['status'] == 0):
            for j in range(len(source['checkpoint_times'])-1, -1, -1):
                near = source['period_residuals'][j, index] < 1e-7 and np.max(abs(source['lle_windows'][j, index])) < .01
                if not near:
                    break
                data['near_streak'][index] += 1
                data['first_near_time'][index] = source['checkpoint_times'][j]
        data['extension_step_counts'] = np.empty((0, n, 2), np.int64)
        data['extension_elapsed'] = np.empty((0, n))
        data['extension_gpu_status'] = np.empty((0, n), np.int32)
        data['extension_stats'] = np.empty((0, n, 8))
        data['extension_probe_status'] = np.empty((0, n), np.int32)
        data['extension_wall_seconds'] = np.empty(0)
        metadata = dict(source_hashes=hashes, source_horizon=float(source['checkpoint_times'][-1]),
            source_checkpoint_count=len(source['checkpoint_times']), runs=[],
            model='Six-state coupled optical operators; assumed logistic illumination, not a measured beam profile',
            capture_rule='Known-period velocity-balanced max residual <1e-7 and all four chunk exponent magnitudes <.01 at two consecutive checkpoints. Capture time is the first near checkpoint; endpoint time is the confirmation checkpoint.',
            status_codes={'0': 'right-censored and still observed', '1': 'near independently checked stable periodic orbit', '11': 'thermal failure', '12': 'numerical failure', '13': 'step budget unresolved', '14': 'wall budget unresolved'},
            caveat='Finite persistence and positive finite-time exponents do not prove an attracting invariant chaotic set; exact long chaotic paths are sensitive to numerical perturbations.')
    if not all(c.get('stable_numerically') for c in controls):
        raise RuntimeError('Missing stable periodic controls')
    gpu = SmoothOperatorsGPU()
    wall = time.perf_counter()
    chunk = 500.
    rtol = 1e-9
    step_factor = .5
    block_durations = []
    record = dict(start_horizon=float(data['checkpoint_times'][-1]), target_horizon=target_time,
        wall_budget_seconds=wall_budget, chunk_duration=chunk, rtol=rtol, atol=1e-12,
        step_factor=step_factor, max_attempts_per_chunk=10000000,
        device=gpu.device.name, script_sha256=hashes['research/extend_optical_lifetimes.py'])
    metadata['runs'].append(record)
    reason = 'target reached'
    while data['checkpoint_times'][-1] < target_time:
        active = np.flatnonzero(data['status'] == 0)
        if not len(active):
            reason = 'all observed survivors captured'
            break
        remaining = wall_budget - (time.perf_counter()-wall)
        estimate = max(20., 1.5*max(block_durations[-3:], default=20.))
        if remaining < estimate:
            reason = 'wall budget: stopped before starting another block'
            break
        began = time.perf_counter()
        current = float(data['checkpoint_times'][-1])
        duration = min(chunk, target_time-current)
        result = gpu.scan(data['parameters'][active], data['final_states'][active],
            transient=0, duration=duration, rtol=rtol, step_factor=step_factor,
            max_attempts=10000000, max_wall_seconds=min(remaining-5., 45.),
            balanced_tangents=data['balanced_tangents'][active])
        n = len(data['parameters'])
        exponents = np.full(n, np.nan)
        windows = np.full((n, 4), np.nan)
        residuals = np.full(n, np.nan)
        counts = np.full((n, 2), -1, np.int64)
        elapsed = np.full(n, np.nan)
        stats = np.full((n, 8), np.nan)
        gpu_status = np.full(n, -1, np.int32)
        probe_status = np.full(n, -1, np.int32)
        data['final_states'][active] = result['final_states']
        data['balanced_tangents'][active] = result['balanced_tangents']
        data['endpoint_time'][active] += result['elapsed']
        failed = result['status'] != 0
        data['status'][active[failed]] = 10 + result['status'][failed]
        exponents[active] = result['lle']
        windows[active] = result['lle_windows']
        counts[active] = result['step_counts']
        elapsed[active] = result['elapsed']
        stats[active] = result['stats']
        gpu_status[active] = result['status']
        if failed.any():
            # Preserve every partial endpoint and its actual clock; do not append
            # a fictitious common checkpoint or classify it as a capture.
            partial = dict(cohort_ids=active, nominal_start=current, requested_duration=duration,
                **{k: v for k, v in result.items() if isinstance(v, np.ndarray)})
            atomic_npz(TARGET / 'partial-block.npz', partial)
            reason = 'partial block retained; common horizon not advanced'
            metadata['partial_block'] = dict(start=current, requested_duration=duration,
                failed_ids=active[failed].tolist(), gpu_settings=result['settings'])
            break
        for control in controls:
            width = control['parameters'][4]
            group = np.flatnonzero((data['status'] == 0) & (data['parameters'][:, 4] == width))
            if not len(group):
                continue
            probe = gpu.scan(data['parameters'][group], data['final_states'][group],
                transient=0, duration=control['period'], rtol=1e-10,
                step_factor=.5, max_wall_seconds=10)
            probe_status[group] = probe['status']
            scale = np.array([1, np.sqrt(control['parameters'][0]), 1, 1, np.sqrt(control['parameters'][0]), 1])
            residuals[group] = np.max(abs((probe['final_states']-data['final_states'][group])/scale), axis=1)
            near = (probe['status'] == 0) & (residuals[group] < 1e-7) & (np.max(abs(windows[group]), axis=1) < .01)
            fresh = group[near & (data['near_streak'][group] == 0)]
            data['first_near_time'][fresh] = current+duration
            data['near_streak'][group] = np.where(near, data['near_streak'][group]+1, 0)
            captured = group[data['near_streak'][group] >= 2]
            data['status'][captured] = 1
            data['capture_time'][captured] = data['first_near_time'][captured]
        data['checkpoint_times'] = np.r_[data['checkpoint_times'], current+duration]
        entries = dict(exponents=exponents, lle_windows=windows, period_residuals=residuals,
            status_history=data['status'].copy(), state_history=data['final_states'].copy(),
            extension_step_counts=counts, extension_elapsed=elapsed, extension_gpu_status=gpu_status,
            extension_stats=stats, extension_probe_status=probe_status)
        for key, value in entries.items():
            data[key] = np.concatenate([data[key], value[None]], axis=0)
        block_durations.append(time.perf_counter()-began)
        data['extension_wall_seconds'] = np.r_[data['extension_wall_seconds'], block_durations[-1]]
        record.update(completed_horizon=current+duration, wall_seconds=time.perf_counter()-wall)
        metadata['completed_extension_blocks'] = len(data['extension_elapsed'])
        atomic_npz(checkpoint, data)
        atomic_json(metadata_path, metadata)
        print(json.dumps(dict(time=current+duration, active=[int(np.sum((data['status']==0)&(data['parameters'][:,4]==c['parameters'][4]))) for c in controls],
            new_captures=int(np.sum((data['source_status']==0)&(data['status']==1))), wall=round(time.perf_counter()-wall, 1))), flush=True)
    record.update(stop_reason=reason, completed_horizon=float(data['checkpoint_times'][-1]),
        wall_seconds=time.perf_counter()-wall)
    metadata['completed_extension_blocks'] = len(data['extension_elapsed'])
    metadata['original_source_unchanged'] = all(digest(p) == hashes[p.as_posix()] for p in paths if p.parent == SOURCE)
    atomic_npz(checkpoint, data)
    atomic_json(metadata_path, metadata)
    print(json.dumps(record, indent=2), flush=True)


def audit():
    data = np.load(TARGET / 'ensemble.npz')
    controls = json.loads((SOURCE / 'periodic-controls.json').read_text())
    rows = []
    for index in np.flatnonzero((data['source_status'] == 0) & (data['status'] == 1)):
        p = data['parameters'][index]
        control = next(c for c in controls if c['parameters'][4] == p[4])
        scale = np.array([1, np.sqrt(p[0]), 1, 1, np.sqrt(p[0]), 1])
        def event(t, y):
            return y[0]-p[2]
        event.direction = 1
        period = control['period']
        sol = solve_ivp(lambda t, y: rhs_and_jac(y, p)[0], (0, 4*period), data['final_states'][index],
            method='DOP853', rtol=1e-12, atol=1e-14, max_step=.005, events=event, dense_output=True)
        row = dict(cohort_id=int(index), width=float(p[4]), capture_checkpoint=float(data['capture_time'][index]),
            endpoint_time=float(data['endpoint_time'][index]), solver_success=bool(sol.success))
        if sol.success and len(sol.t_events[0]) >= 3:
            errors = dict(section_error=float(np.max(abs((sol.y_events[0][-2:]-np.array(control['section_state']))/scale))),
                period_error=float(np.max(abs(np.diff(sol.t_events[0][-3:])-period))),
                scaled_flow_closure=float(np.max(abs((sol.sol(4*period)-sol.sol(3*period))/scale))))
            row.update(**errors, agrees_with_checked_cycle=max(errors.values()) < 1e-7)
        else:
            row.update(agrees_with_checked_cycle=False, message=sol.message)
        rows.append(row)
    result = dict(method='Independent DOP853, rtol=1e-12, atol=1e-14, max_step=.005, four periods per newly captured endpoint',
        scope='Only new captures; the 54 original captures retain their original CPU audit.',
        source_control_sha256=digest(SOURCE / 'periodic-controls.json'), captured_count=len(rows),
        all_passed=all(r['agrees_with_checked_cycle'] for r in rows), rows=rows)
    atomic_json(TARGET / 'capture-cpu-audit.json', result)
    print(json.dumps(result, indent=2), flush=True)


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data = np.load(TARGET / 'ensemble.npz')
    times = data['checkpoint_times']
    widths = np.unique(data['parameters'][:, 4])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), layout='constrained')
    groups = []
    for width, color in zip(widths, ['#2680b9', '#cf6945', '#238b75']):
        ids = np.flatnonzero(data['parameters'][:, 4] == width)
        active = ids[data['status'][ids] == 0]
        captured = ids[data['status'][ids] == 1]
        new = ids[(data['source_status'][ids] == 0) & (data['status'][ids] == 1)]
        counts = np.sum(data['status_history'][:, ids] == 0, axis=1)
        axes[0].step(np.r_[0, times], np.r_[len(ids), counts], where='post', color=color, label=f'w = {width:g}: {len(active)}/32 remain')
        if len(active):
            minimum = np.min(data['lle_windows'][-1, active], axis=1)
            axes[1].scatter(active, minimum, color=color, label=f'w = {width:g}', s=28)
        groups.append(dict(width=float(width), starts=len(ids), original_survivors=int(np.sum(data['source_status'][ids]==0)),
            newly_captured=len(new), total_captured=len(captured), still_observed=len(active),
            new_capture_ids=new.tolist(), surviving_ids=active.tolist(),
            late_min_window_range=None if not len(active) else [float(min(minimum)), float(max(minimum))],
            all_four_late_windows_above_point01=int(np.sum(np.all(data['lle_windows'][-1, active] > .01, axis=1))),
            last_capture=None if not len(new) else float(np.max(data['capture_time'][new]))))
    axes[0].axvline(20000, color='#666', linestyle='--', lw=1, label='Continuation begins')
    axes[0].set(xlabel='Time / thermal relaxation time', ylabel='Trajectories still under observation',
        ylim=(-.5, 33), title='The same 96 initial conditions')
    axes[0].legend(fontsize=8)
    axes[1].axhline(0, color='#666', lw=1)
    axes[1].set(xlabel='Original cohort ID', ylabel='Minimum of four late finite-time exponents',
        title=f'Last completed 500-unit block: t = {times[-1]:,.0f}')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Longer optical lifetimes: finite survival is not proof of a chaotic attractor\nSix-state model with assumed logistic illumination', fontsize=12)
    path = Path('docs/assets/optical_lifetime_extension.png')
    fig.savefig(path, dpi=170)
    plt.close(fig)
    summary = dict(observation_horizon=float(times[-1]), completed_extension_blocks=len(data['extension_elapsed']),
        groups=groups, total_new_captures=int(np.sum((data['source_status']==0)&(data['status']==1))),
        unresolved_ids=np.flatnonzero(data['status'] > 1).tolist(),
        caveat='Descriptive finite-cohort counts, not basin probabilities or asymptotic classification. Captures are checkpoint-resolved and surviving trajectories are right-censored.')
    atomic_json(TARGET / 'summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--audit', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--target', type=float, default=100000.)
    parser.add_argument('--budget', type=float, default=515.)
    args = parser.parse_args()
    if args.run:
        run(args.target, args.budget, args.resume)
    if args.audit:
        audit()
    if args.plot:
        plot()
