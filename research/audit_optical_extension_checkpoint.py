"""Incrementally verify newly captured optical endpoints with an independent CPU solver.

This can run beside the GPU continuation: captured endpoints are frozen. Cached
checks are reused only when the parameters, endpoint bytes, and periodic-control
content hash agree. Original pre-extension captures are never rerun here.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from smooth_operators import rhs_and_jac

ROOT = Path('results/light-operators/lifetimes-extension')
CONTROLS = Path('results/light-operators/lifetimes/periodic-controls.json')


def content_hash(p, state, control):
    payload = np.r_[p, state].astype('<f8').tobytes() + json.dumps(control, sort_keys=True, separators=(',', ':')).encode()
    return hashlib.sha256(payload).hexdigest()


def audit(width=None):
    # Freeze this complete checkpoint in memory while the GPU writes later ones.
    data = {k: v.copy() for k, v in np.load(io.BytesIO((ROOT / 'ensemble.npz').read_bytes())).items()}
    controls = json.loads(CONTROLS.read_text())
    path = ROOT / 'capture-cpu-audit.json'
    old = json.loads(path.read_text()) if path.exists() else {'rows': []}
    rows = {int(r['cohort_id']): r for r in old['rows']}
    selected = (data['source_status'] == 0) & (data['status'] == 1)
    if width is not None:
        selected &= data['parameters'][:, 4] == width
    performed = []
    for index in np.flatnonzero(selected):
        p = data['parameters'][index]
        state = data['final_states'][index]
        control = next(c for c in controls if c['parameters'][4] == p[4])
        checksum = content_hash(p, state, control)
        if int(index) in rows and rows[int(index)].get('input_content_sha256') == checksum:
            continue
        scale = np.array([1, np.sqrt(p[0]), 1, 1, np.sqrt(p[0]), 1])
        def event(t, y):
            return y[0]-p[2]
        event.direction = 1
        period = control['period']
        sol = solve_ivp(lambda t, y: rhs_and_jac(y, p)[0], (0, 4*period), state,
            method='DOP853', rtol=1e-12, atol=1e-14, max_step=.005, events=event, dense_output=True)
        row = dict(cohort_id=int(index), width=float(p[4]), parameters=p.tolist(), endpoint=state.tolist(),
            input_content_sha256=checksum, capture_checkpoint=float(data['capture_time'][index]),
            endpoint_time=float(data['endpoint_time'][index]), solver_success=bool(sol.success))
        if sol.success and len(sol.t_events[0]) >= 3:
            errors = dict(section_error=float(np.max(abs((sol.y_events[0][-2:]-np.array(control['section_state']))/scale))),
                period_error=float(np.max(abs(np.diff(sol.t_events[0][-3:])-period))),
                scaled_flow_closure=float(np.max(abs((sol.sol(4*period)-sol.sol(3*period))/scale))))
            row.update(**errors, agrees_with_checked_cycle=max(errors.values()) < 1e-7)
        else:
            row.update(agrees_with_checked_cycle=False, message=sol.message)
        rows[int(index)] = row
        performed.append(int(index))
    result = dict(method='DOP853, rtol=1e-12, atol=1e-14, max_step=.005, four periods per captured endpoint',
        scope='Independent checks only for new captures. Original 54 captures retain their existing audit. Cached checks require exact parameter/endpoint bytes and canonical periodic-control content.',
        checked_checkpoint_horizon=float(data['checkpoint_times'][-1]), selected_width=width,
        new_checks_this_call=performed, captured_count=len(rows),
        all_passed=all(r['agrees_with_checked_cycle'] for r in rows.values()),
        rows=[rows[i] for i in sorted(rows)])
    temporary = path.with_suffix('.temporary.json')
    temporary.write_text(json.dumps(result, indent=2), encoding='utf-8')
    temporary.replace(path)
    print(json.dumps(result, indent=2), flush=True)
    if not result['all_passed']:
        raise RuntimeError('A captured endpoint needs further investigation')


def integrity():
    source_path = Path('results/light-operators/lifetimes/ensemble.npz')
    source = np.load(source_path)
    data = np.load(io.BytesIO((ROOT / 'ensemble.npz').read_bytes()))
    metadata = json.loads((ROOT / 'ensemble.json').read_text())
    nold = len(source['checkpoint_times'])
    old_captured = source['status'] == 1
    history_keys = ['checkpoint_times', 'exponents', 'lle_windows', 'period_residuals', 'status_history', 'state_history']
    checks = dict(original_binary_hash_unchanged=hashlib.sha256(source_path.read_bytes()).hexdigest() == metadata['source_hashes'][source_path.as_posix()],
        original_history_exact=all(np.array_equal(data[k][:nold], source[k], equal_nan=True) for k in history_keys),
        original_parameters_exact=np.array_equal(data['parameters'], source['parameters']),
        original_initial_states_exact=np.array_equal(data['initial_states'], source['initial_states']),
        original_cohort_ids_exact=np.array_equal(data['cohort_ids'], np.arange(len(source['parameters']))),
        original_captured_states_untouched=np.array_equal(data['final_states'][old_captured], source['final_states'][old_captured]),
        original_captured_tangents_untouched=np.array_equal(data['balanced_tangents'][old_captured], source['balanced_tangents'][old_captured]),
        original_capture_times_untouched=np.array_equal(data['capture_time'][old_captured], source['capture_time'][old_captured]),
        all_checkpoint_states_finite=bool(np.isfinite(data['state_history']).all()),
        all_active_endpoint_times_match_horizon=bool(np.all(data['endpoint_time'][data['status']==0] == data['checkpoint_times'][-1])),
        all_extension_lanes_completed=bool(np.all(np.isin(data['extension_gpu_status'], [-1, 0]))),
        all_capture_probes_completed=bool(np.all(np.isin(data['extension_probe_status'], [-1, 0]))))
    archive = ROOT / 'segments/segment-01.npz'
    if archive.exists():
        first = np.load(archive)
        length = len(first['checkpoint_times'])
        checks['first_segment_history_exact'] = all(np.array_equal(data[k][:length], first[k], equal_nan=True) for k in history_keys)
        if len(metadata['runs']) > 1:
            checks['resume_input_matches_first_archive'] = hashlib.sha256(archive.read_bytes()).hexdigest() == metadata['runs'][1]['starting_checkpoint_npz_sha256']
    for filename, checksum in metadata.get('canonical_dependency_text_hashes', {}).items():
        value = hashlib.sha256(Path(filename).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        checks['canonical_dependency:' + filename] = value == checksum
    for index, run in enumerate(metadata['runs'], 1):
        script = ROOT / f'segments/segment-{index:02d}-script.py'
        if script.exists():
            content = script.read_bytes()
            expected = run.get('script_canonical_lf_sha256', run['script_sha256'])
            checks[f'archived_driver_{index}_matches_record'] = hashlib.sha256(content.replace(b'\r\n', b'\n')).hexdigest() == expected
    stats = data['extension_stats']
    counts = data['extension_step_counts']
    accepted = counts[:, :, 0] >= 0
    result = dict(horizon=float(data['checkpoint_times'][-1]), checks=checks, all_passed=all(checks.values()),
        completed_extension_blocks=len(data['extension_elapsed']),
        accepted_steps=int(counts[:, :, 0][accepted].sum()), rejected_steps=int(counts[:, :, 1][accepted].sum()),
        thermal_range_at_accepted_steps=[float(np.nanmin(stats[:, :, 4])), float(np.nanmax(stats[:, :, 5]))],
        displacement_ranges_at_accepted_steps=[[float(np.nanmin(stats[:, :, 0])), float(np.nanmax(stats[:, :, 1]))],
            [float(np.nanmin(stats[:, :, 2])), float(np.nanmax(stats[:, :, 3]))]],
        max_tangent_unit_norm_error=float(np.max(abs(np.linalg.norm(data['balanced_tangents'], axis=1)-1))),
        cumulative_gpu_segment_wall_seconds=sum(r['wall_seconds'] for r in metadata['runs']),
        scope='Provenance and continuation checks, not an independent proof of chaotic dynamics or exact reproduction of numerical escape times.')
    (ROOT / 'integrity-audit.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    artifacts = [ROOT / 'ensemble.npz', ROOT / 'ensemble.json', ROOT / 'capture-cpu-audit.json',
        Path('research/extend_optical_lifetimes.py'), Path('research/audit_optical_extension_checkpoint.py')]
    artifacts += sorted((ROOT / 'segments').glob('*')) if (ROOT / 'segments').exists() else []
    records = []
    for path in artifacts:
        content = path.read_bytes()
        records.append(dict(path=path.as_posix(), raw_sha256=hashlib.sha256(content).hexdigest(),
            canonical_lf_sha256=None if path.suffix == '.npz' else hashlib.sha256(content.replace(b'\r\n', b'\n')).hexdigest()))
    manifest = dict(policy='Binary NPZ hashes are exact bytes. Text has both historical raw-byte and CRLF-to-LF canonical hashes. Archived scripts identify the actual numerical runs; final drivers additionally add portable CLI and analysis features.',
        horizon=float(data['checkpoint_times'][-1]), artifacts=records)
    (ROOT / 'provenance-manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
    if not result['all_passed']:
        raise RuntimeError('Continuation integrity check failed')


def windows():
    data = np.load(io.BytesIO((ROOT / 'ensemble.npz').read_bytes()))
    groups = []
    for width in np.unique(data['parameters'][:, 4]):
        ids = np.flatnonzero((data['parameters'][:, 4] == width) & (data['status'] == 0))
        values = data['lle_windows'][:, ids]
        groups.append(dict(width=float(width), surviving_ids=ids.tolist(),
            windows_per_survivor=4*len(data['checkpoint_times']),
            all_recorded_values_finite=bool(np.isfinite(values).all()),
            minimum_across_full_history=None if not len(ids) else float(values.min()),
            maximum_across_full_history=None if not len(ids) else float(values.max()),
            survivors_with_every_window_above_point01=int(np.sum(np.all(values > .01, axis=(0, 2)))),
            survivor_full_history_mean_range=None if not len(ids) else [float(values.mean(axis=(0, 2)).min()), float(values.mean(axis=(0, 2)).max())],
            survivor_late_block_mean_range=None if not len(ids) else [float(data['exponents'][-1, ids].min()), float(data['exponents'][-1, ids].max())]))
    latest_captures = []
    for width in [.003, .005]:
        ids = np.flatnonzero((data['parameters'][:, 4] == width) & (data['status'] == 1))
        if not len(ids):
            continue
        index = ids[np.argmax(data['capture_time'][ids])]
        values = data['lle_windows'][:, index].ravel()
        first = np.flatnonzero(~(values > .01))
        count = int(first[0]) if len(first) else len(values)
        latest_captures.append(dict(width=width, cohort_id=int(index),
            selection='Chosen retrospectively as the latest confirmed capture in this width cohort',
            capture_checkpoint=float(data['capture_time'][index]),
            confirmation_checkpoint=float(data['endpoint_time'][index]),
            consecutive_initial_125_unit_windows_above_point01=count,
            duration_covered_by_those_windows=125.*count,
            smallest_estimate_in_that_prefix=None if not count else float(values[:count].min()),
            prefix_mean=None if not count else float(values[:count].mean())))
    result = dict(horizon=float(data['checkpoint_times'][-1]), window_duration=125.,
        scope='Surviving trajectories include the original initial observation. Selection by survival and retrospective selection of latest captures are explicit; these are finite-time tangent-growth values, not a full asymptotic spectrum.',
        groups=groups, latest_capture_examples=latest_captures)
    (ROOT / 'window-persistence-audit.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--width', type=float)
    parser.add_argument('--integrity', action='store_true')
    parser.add_argument('--windows', action='store_true')
    parser.add_argument('--root', type=Path, help='Alternate extension result directory')
    args = parser.parse_args()
    if args.root:
        ROOT = args.root
    if args.integrity:
        integrity()
    elif args.windows:
        windows()
    else:
        audit(args.width)
