"""Tighter continuation of all 40 N=3 positives from the declared scout.

This tests finite-time persistence, not scientific priority or asymptotic chaos.
The initial state and tangent are the saved scout endpoints. Warm-up 300 and
observation 1200 are original tau for every parameter tuple, not scaled time.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from motor_rows import field_jac, unpack
from motor_rows_gpu import MotorRowsGPU, save

ROOT = Path(__file__).resolve().parents[1]


def validate_extended_kernel(gpu, scout, selected):
    """Exercise original-time mode and new pair statistics at two time scales."""
    cases = [int(selected[np.argmax(scout['parameters'][selected, 0])]),
             int(selected[np.argmin(scout['parameters'][selected, 0])])]
    records = []
    for index in cases:
        p, initial, tangent = scout['parameters'][index], scout['final_states'][index], scout['final_tangents'][index]
        initial = np.r_[initial, tangent/np.linalg.norm(tangent)]
        def augmented(t, z):
            f, jac = field_jac(z[:8], p, 3)
            return np.r_[f, jac@z[8:]]
        references = {}
        for method in ('DOP853', 'Radau'):
            sol = solve_ivp(augmented, (0., 2.), initial, method=method, rtol=2e-12, atol=2e-14,
                            max_step=.005, dense_output=True)
            if not sol.success:
                raise RuntimeError(sol.message)
            endpoint = sol.y[:, -1]
            norm = np.linalg.norm(endpoint[8:])
            samples = sol.sol(np.linspace(1.5, 2., 5001))[:8]
            c, s, x = samples[:3], samples[3:6], np.vstack((samples[6:], -samples[6:].sum(axis=0)))
            distances = np.array([np.sqrt((c[a]-c[b])**2+(s[a]-s[b])**2+(x[a]-x[b])**2)
                                  for a, b in ((0, 1), (0, 2), (1, 2))])
            references[method] = (endpoint[:8], endpoint[8:]/norm, np.log(norm)/2., distances.max(axis=1))
        result = gpu.scan(p[None], initial[None, :8], transient=0., duration=2., rtol=1e-11,
                          step_factor=.5, max_wall=30., tangents=initial[None, 8:], time_mode='original')
        a, b = references['DOP853'], references['Radau']
        record = dict(source_scout_id=index, parameters=p.tolist(), status=int(result['status'][0]),
                      cpu_state_max_difference=float(max(abs(a[0]-b[0]))),
                      gpu_state_max_difference=float(max(abs(a[0]-result['final_states'][0]))),
                      gpu_tangent_max_difference=float(max(abs(a[1]-result['final_tangents'][0]))),
                      gpu_lle_difference=float(abs(a[2]-result['lle'][0])),
                      pair_max_distance_sample_difference=float(max(abs(a[3]-result['late_pair_max_distance'][0]))))
        record['passed'] = bool(record['status'] == 0 and record['cpu_state_max_difference'] < 1e-8
                                and record['gpu_state_max_difference'] < 1e-8
                                and record['gpu_tangent_max_difference'] < 1e-8
                                and record['gpu_lle_difference'] < 1e-7
                                and record['pair_max_distance_sample_difference'] < .02)
        records.append(record); print(json.dumps(record), flush=True)
    report = dict(all_passed=all(row['passed'] for row in records), rows=records,
                  note='Independent solver agreement and direct sampled row-distance checks after adding original-time mode and passive diagnostic accumulators; source field/Jv unchanged.')
    (ROOT/'results/motor-rows/gpu-refinement-validation.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    if not report['all_passed']:
        raise RuntimeError('Extended-kernel controls failed')


def main():
    source_path = ROOT/'results/motor-rows/gpu-scout-N3.npz'
    scout = np.load(source_path)
    selected = np.flatnonzero((scout['status'] == 0) & (scout['lle'] > .005)
                              & np.all(scout['lle_windows'] > .001, axis=1))
    if len(selected) != 40:
        raise RuntimeError('The declared scout or selection changed; review rather than silently expanding scope')
    gpu = MotorRowsGPU(3)
    validate_extended_kernel(gpu, scout, selected)
    result = gpu.scan(scout['parameters'][selected], scout['final_states'][selected],
                      transient=300., duration=1200., time_mode='original', rtol=1e-10,
                      step_factor=.5, max_steps=8000000, max_wall=600., tangents=scout['final_tangents'][selected])
    result.update(source_scout_ids=selected, source_tuple_ids=scout['tuple_id'][selected],
                  source_start_ids=scout['start_id'][selected], previous_lle=scout['lle'][selected],
                  previous_lle_windows=scout['lle_windows'][selected],
                  provenance=dict(source_npz=source_path.relative_to(ROOT).as_posix(),
                                  source_npz_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                                  orchestrator=Path(__file__).relative_to(ROOT).as_posix(),
                                  orchestrator_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
                                  selection='All 40 declared scout positives; no discarded failures or cherry-picked starts',
                                  continuation='State and tangent carried from the scout endpoint; an additional 300 original-tau warm-up precedes 1200 original-tau observation'))
    valid = result['status'] == 0
    late_collapsed = valid & np.any(result['late_pair_max_distance'] < 1e-6, axis=1)
    result['cluster_screen'] = dict(threshold=1e-6,
                                   description='A pair is descriptively collapsed only if its (c,s,X) distance stays below tolerance throughout the last observation quarter. This is a finite-window numerical classification.',
                                   completed_with_collapsed_pair_source_ids=selected[late_collapsed].tolist(),
                                   completed_with_all_pairs_resolved_source_ids=selected[valid & ~late_collapsed].tolist(),
                                   caveat='A small instantaneous minimum does not imply clustering. RMS uses right-endpoint quadrature; all pair ranges are retained.')
    save(ROOT/'results/motor-rows/gpu-refined-N3.npz', result)


if __name__ == '__main__':
    main()
