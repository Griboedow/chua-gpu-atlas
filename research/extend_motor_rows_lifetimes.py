"""Fixed-cohort long continuation of the 24 retained N=3 motor-row candidates.

No reseeding or replacement of weak/failed trajectories. Every 500 original
tau block preserves state, tangent, all exponent windows, domain diagnostics,
and late-quarter row-pair ranges. The total device wall budget is bounded.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from motor_rows_gpu import MotorRowsGPU, save

ROOT = Path(__file__).resolve().parents[1]


def main(output, wall=900.):
    output = Path(output)
    if not output.is_absolute():
        output = ROOT/output
    if output.exists() and any(output.iterdir()):
        raise RuntimeError('Refusing to overwrite a lifetime study; select a new output directory')
    output.mkdir(parents=True, exist_ok=True)
    source_path = ROOT/'results/motor-rows/gpu-refined-N3.npz'
    source = np.load(source_path)
    mask = (source['status'] == 0) & (source['lle'] > .005) & np.all(source['lle_windows'] > .001, axis=1)
    selected = np.flatnonzero(mask)
    if len(selected) != 24:
        raise RuntimeError('The authorized fixed cohort changed; review the input rather than resampling')
    source_ids = source['source_scout_ids'][selected]
    parameters = source['parameters'][selected]
    state = source['final_states'][selected].copy()
    tangent = source['final_tangents'][selected].copy()
    scout = np.load(ROOT/'results/motor-rows/gpu-scout-N3.npz')
    prior_tau = scout['elapsed_original'][source_ids]+source['elapsed_original'][selected]
    manifest = dict(N=3, dimension=8, source_scout_ids=source_ids.tolist(),
                    parameters=parameters.tolist(), initial_states=state.tolist(), initial_tangents=tangent.tolist(),
                    prior_original_tau=prior_tau.tolist(),
                    source_file=source_path.relative_to(ROOT).as_posix(),
                    source_exact_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                    orchestrator=Path(__file__).relative_to(ROOT).as_posix(),
                    orchestrator_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
                    protocol=dict(additional_warmup=0., target_additional_original_tau=10000., block_original_tau=500.,
                                  windows_per_block=4, rtol=1e-10, atol=1e-13, step_cap_factor=.5,
                                  maximum_total_wall_seconds=wall, fixed_cohort_size=24,
                                  initial_tangents='Carried aligned vectors from the completed refinement',
                                  row_distances='Physical normalized (c,s,X), at accepted endpoints in each block\'s final125tau; these ranges do not continuously cover the entire10000tau'),
                    caveat='Long finite-time persistence can reduce a transient explanation; it cannot prove an asymptotic attractor or scientific novelty.')
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    gpu = MotorRowsGPU(3)
    start = time.perf_counter(); cumulative_logs = np.zeros(24); elapsed = np.zeros(24)
    blocks = []; stop_reason = 'target reached'; last_status = np.zeros(24, dtype=int)
    for number in range(1, 21):
        remaining = wall-(time.perf_counter()-start)
        if remaining <= 0:
            stop_reason = 'wall budget exhausted between blocks'
            break
        result = gpu.scan(parameters, state, transient=0., duration=500., time_mode='original',
                          rtol=1e-10, step_factor=.5, max_steps=8000000, max_wall=remaining, tangents=tangent)
        result.update(source_scout_ids=source_ids, lifetime_block=number,
                      lifetime_scope='Same fixed24 refined trajectories, no added transient or replacement starts',
                      additional_start_original_tau=elapsed.copy(),
                      cumulative_start_original_tau=prior_tau+elapsed)
        save(output/f'block-{number:03d}.npz', result)
        cumulative_logs += result['raw_log_sums'].sum(axis=1)
        elapsed += result['elapsed_original']
        state, tangent = result['final_states'].copy(), result['final_tangents'].copy()
        last_status = result['status'].copy()
        record = dict(number=number, file=f'block-{number:03d}.npz',
                      completed=int(np.sum(result['status'] == 0)),
                      original_tau_per_lane=result['elapsed_original'].tolist(),
                      window_lle=result['lle_windows'].tolist(), status=result['status'].tolist(),
                      cumulative_lle=(cumulative_logs/elapsed).tolist(),
                      wall_seconds=result['wall_seconds'])
        blocks.append(record)
        np.savez_compressed(output/'checkpoint.npz', parameters=parameters, source_scout_ids=source_ids,
                            states=state, tangents=tangent, cumulative_log_sums=cumulative_logs,
                            additional_elapsed_original=elapsed, prior_original_tau=prior_tau,
                            status=last_status, block_number=np.array(number))
        summary = dict(source_scout_ids=source_ids.tolist(), completed_blocks=len(blocks), blocks=blocks,
                       additional_elapsed_original=elapsed.tolist(), cumulative_lle=(cumulative_logs/elapsed).tolist(),
                       last_status=last_status.tolist(), all_target_reached=bool(np.all(elapsed >= 10000.)),
                       wall_seconds=time.perf_counter()-start,
                       scope='Fixed-cohort finite-time continuation, not a proof of an asymptotic attractor. Every block and incomplete status is retained.')
        (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=True)+'\n', encoding='utf-8')
        print(json.dumps(dict(lifetime_block=number, complete=record['completed'],
                              minimum_cumulative_lle=float(np.min(cumulative_logs/elapsed)),
                              maximum_cumulative_lle=float(np.max(cumulative_logs/elapsed)),
                              total_wall=round(time.perf_counter()-start, 1))), flush=True)
        if np.any(last_status != 0):
            stop_reason = 'a block was incomplete or failed; fixed cohort stopped without silently dropping lanes'
            break
    summary['stop_reason'] = stop_reason
    summary['wall_seconds'] = time.perf_counter()-start
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=True)+'\n', encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/motor-rows/lifetimes')
    parser.add_argument('--wall', type=float, default=900.)
    args = parser.parse_args()
    main(args.output, args.wall)
