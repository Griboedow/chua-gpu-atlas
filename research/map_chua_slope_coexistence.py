"""Map initial-state capture by the two independently verified coexisting cycles."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from chua_gpu import GPU, DEFAULT_SEEDS
from detail_analysis import collect_peaks
from survey_chua_layers import recurrence
from trace_chua_window import Chua
from continue_chua_slope_families import dump


def signatures(witness):
    result = {}
    for r in witness['cycles']:
        model = Chua(alpha=r['alpha'], m0=r['m0'], m1=r['m1'], rtol=2e-12, step=.02)
        record = model.record(r['beta'], [*r['q'][:2], r['period']])
        result[r['full_cycle_peaks']] = np.array(record['maxima'])
    return result


def classify(peaks, status, classes, references):
    periods = np.array([recurrence(v) if s == 0 else 0 for v, s in zip(peaks, status)])
    labels = np.zeros(len(peaks), dtype=np.int8)
    errors = np.full(len(peaks), np.nan)
    for i, count in enumerate(periods):
        if count not in references or classes[i] != 2: continue
        v = peaks[i][np.isfinite(peaks[i])][-4*count:]
        reference = references[count]
        error = min(np.max(abs(v-np.tile(np.roll(reference, offset), 4))) for offset in range(count))
        errors[i] = error
        if error < 5e-4*max(1., max(abs(reference))): labels[i] = count
    return labels, periods, errors


def run(gpu, parameters, initial, refs, dt, transient, duration):
    params = np.tile(parameters, (len(initial), 1))
    result = gpu.scan(params, initial_states=initial[:, None, :], dt=dt, transient=transient, duration=duration)
    peaks, counts, status = collect_peaks(gpu, params, result['final_states'][:, 0], duration=800, dt=dt, capacity=384)
    labels, periods, errors = classify(peaks, status, result['classes'][:, 0], refs)
    arrays = {k:v for k,v in result.items() if isinstance(v, np.ndarray) and k != 'stats'}
    # Retain enough late maxima to inspect four repeats of the longer reference.
    tails = np.full((len(initial), 28), np.nan)
    for i, v in enumerate(peaks):
        finite = v[np.isfinite(v)][-28:]; tails[i, -len(finite):] = finite if len(finite) else np.nan
    return dict(initial=initial, labels=labels, peak_periods=periods, peak_counts=counts,
                peak_status=status, signature_error=errors, late_maxima=tails, **arrays), result['settings']


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--stage', choices=['basin','time-control'], default='basin')
    a = p.parse_args(); out = a.output
    w = json.loads((out/'coexistence-witness.json').read_text()); refs = signatures(w)
    parameters = [w['parameters'][k] for k in ['alpha','beta','m0','m1']]
    if a.stage == 'time-control':
        base = np.load(out/'basin.npz'); refined = np.load(out/'basin-refinement.npz')
        arrays, settings = run(GPU(), parameters, refined['initial'], refs, .0025, 2400, 1600)
        np.savez_compressed(out/'basin-time-only.npz', indices=refined['indices'], **arrays)
        old = base['labels'][refined['indices']]; new = arrays['labels']
        dump(out/'basin-time-only.json', dict(samples=len(old), settings=settings,
            all_label_agreement=float(np.mean(old==new)),
            resolved_label_agreement=float(np.mean(old[old!=0]==new[old!=0])),
            counts={str(x):{str(y):int(((old==x)&(new==y)).sum()) for y in [0,5,7]} for x in [0,5,7]},
            original_starts_baseline=base['labels'][-3:].tolist(),original_starts_time_only=new[-3:].tolist(),
            original_starts_step_and_time=refined['labels'][-3:].tolist(),
            scope='Same step as the baseline with doubled durations; separates transient-length sensitivity from the combined half-step check.'))
        return
    y = np.linspace(.19, .36, 121); z = np.linspace(-.60, -.30, 121)
    zz, yy = np.meshgrid(z, y, indexing='ij')
    initial = np.column_stack([np.zeros(yy.size), yy.ravel(), zz.ravel()])
    initial = np.vstack([initial, DEFAULT_SEEDS])
    gpu = GPU()
    arrays, settings = run(gpu, parameters, initial, refs, .0025, 1200, 800)
    np.savez_compressed(out/'basin.npz', y=y, z=z, parameters=parameters, **arrays)
    print('basin', {str(k):int((arrays['labels']==k).sum()) for k in [0,5,7]}, flush=True)
    rng = np.random.default_rng(20261030)
    selected = []
    for label in np.unique(arrays['labels']):
        ids = np.flatnonzero(arrays['labels']==label)
        selected.extend(rng.choice(ids, min(96, len(ids)), replace=False).tolist())
    indices = np.array(sorted(set(selected) | set(range(len(initial)-3, len(initial)))))
    refined, refinement_settings = run(gpu, parameters, initial[indices], refs, .00125, 2400, 1600)
    np.savez_compressed(out/'basin-refinement.npz', indices=indices, **refined)
    report = dict(parameters=w['parameters'], shape=[121,121], axes_order=['z','y'], fixed_initial_x=0,
        settings=settings, peak_observation=800, reference_maxima={str(k):v.tolist() for k,v in refs.items()},
        labels={'0':'Unresolved or another outcome','5':'Matched the verified five-peak cycle','7':'Matched the verified seven-peak cycle'},
        counts={str(k):int((arrays['labels'][:-3]==k).sum()) for k in [0,5,7]},
        original_three_starts=dict(states=np.asarray(DEFAULT_SEEDS).tolist(), labels=arrays['labels'][-3:].tolist()),
        refinement=dict(samples=len(indices), settings=refinement_settings,
                        label_agreement=float(np.mean(arrays['labels'][indices]==refined['labels'])),
                        changed_indices=indices[arrays['labels'][indices]!=refined['labels']].tolist()),
        device=gpu.info, scope='Finite-time initial-condition slice, not a full basin volume. Matching requires a repeated peak count, '
        'the ordered peak-height signature of the closed reference, and an earlier regular-candidate classification. '
        'Small basin boundaries and unresolved outcomes are not certified.')
    dump(out/'basin.json', report); print('refinement', report['refinement']['label_agreement'], flush=True)


if __name__ == '__main__': main()
