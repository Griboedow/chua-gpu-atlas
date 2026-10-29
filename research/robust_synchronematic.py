"""Controlled robustness experiments around a published-model chaos candidate.

Geometry perturbations retain the published equations. Coefficient perturbations
and harmonic ablations change the reduced model and are labeled as such.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from synchronematic import COEFFICIENTS, SynchronematicGPU, flow_and_jac, geometry, gradient_coefficients, potential
from search_synchronematic import save


def validate_four():
    data = np.load('results/synchronematic/n4-scout/scan.npz')
    p, start = data['positions'][43], data['starts'][43]
    pairs = geometry(p)
    rng = np.random.default_rng(41)
    coefficients = np.array([COEFFICIENTS, COEFFICIENTS * rng.uniform(.9, 1.1, (2, 4))])
    gpu = SynchronematicGPU(4)
    errors = []
    for dt in [.01, .005, .0025]:
        result = gpu.scan(np.repeat(p[None], 2, axis=0), np.repeat(start[None], 2, axis=0),
                          dt=dt, transient=0, duration=1, coefficients=coefficients)
        row = []
        for c, end in zip(coefficients, result['final_states']):
            reference = solve_ivp(lambda t, y: flow_and_jac(y, pairs, c)[0],
                                  (0, 1), start, method='DOP853', rtol=1e-12, atol=1e-14)
            if not reference.success:
                raise RuntimeError(reference.message)
            row.append(float(np.max(abs(np.angle(np.exp(1j * (end - reference.y[:, -1])))))))
        errors.append(row)
    errors = np.array(errors)
    result = dict(dt=[.01, .005, .0025], errors=errors.tolist(),
                  columns=['published coefficients', 'perturbed coefficients'],
                  passed=bool(np.all(errors[-1] < 1e-6) and np.all(errors[-1] < errors[0] / 100)))
    g=gradient_coefficients();div_errors=[];gradient_errors=[]
    for s in rng.uniform(-np.pi,np.pi,(20,8)):
        f,j=flow_and_jac(s,pairs);fg,jg=flow_and_jac(s,pairs,g)
        eps=1e-6
        grad=np.array([(potential(s+eps*e,pairs)-potential(s-eps*e,pairs))/(2*eps) for e in np.eye(8)])
        gradient_errors.append(float(np.max(abs(fg+grad))));div_errors.append(float(abs(np.trace(j-jg))))
    result.update(gradient_coefficients=g.tolist(),gradient_error=max(gradient_errors),solenoidal_divergence_error=max(div_errors))
    result['passed']=result['passed'] and max(gradient_errors)<1e-7 and max(div_errors)<1e-10
    Path('results/synchronematic/validation-n4.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)
    if not result['passed']:
        raise RuntimeError('N=4 validation failed')


def make_ensemble(mode, geometries, seeds):
    data = np.load('results/synchronematic/n4-scout/scan.npz')
    p = data['positions'][43]
    p = p - p.mean(axis=0)
    rng = np.random.default_rng(2026102904)
    positions, coefficients, labels, starts = [], [], [], []
    if mode == 'step':
        ids = np.r_[np.arange(40, 48), np.arange(376, 384), np.arange(1696, 1704)]
        return data['positions'][ids], data['starts'][ids], np.broadcast_to(COEFFICIENTS, (len(ids), 2, 4)).copy(), ids
    if mode == 'geometry':
        levels = [0., .01, .03, .1, .25]
    elif mode == 'coefficients':
        levels = [0., .03, .06, .15]
    else:
        levels = np.linspace(0., 1.5, 61)
        geometries = 1
    for level in levels:
        for sample in range(geometries):
            pp, cc = p.copy(), COEFFICIENTS.copy()
            if mode == 'geometry':
                # Absolute displacement in units of the original minimum spacing.
                pp += rng.uniform(-level, level, pp.shape)
                distances = np.linalg.norm(pp[:, None] - pp[None, :], axis=2)
                distances += np.eye(4) * 1e6
                pp /= distances.min()  # Factor out the uniform time rescaling.
            elif mode == 'coefficients':
                cc *= rng.uniform(1 - level, 1 + level, cc.shape)
            elif mode=='nonreciprocity':
                cc[0, [0, 2]] *= level  # Even-in-phase, nonreciprocal phase terms.
            else:
                g=gradient_coefficients();cc=g+level*(COEFFICIENTS-g)
            for seed in range(seeds):
                positions.append(pp)
                coefficients.append(cc)
                labels.append([level, sample, seed])
                starts.append(rng.uniform(-np.pi, np.pi, 8))
    return np.array(positions), np.array(starts), np.array(coefficients), np.array(labels)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['geometry', 'coefficients', 'nonreciprocity', 'solenoidal', 'step'], default='geometry')
    parser.add_argument('--geometries', type=int, default=32)
    parser.add_argument('--seeds', type=int, default=4)
    parser.add_argument('--dt', type=float, default=.005)
    parser.add_argument('--transient', type=float, default=1000)
    parser.add_argument('--duration', type=float, default=2000)
    parser.add_argument('--validate', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.validate:
        validate_four()
        return
    positions, starts, coefficients, labels = make_ensemble(args.mode, args.geometries, args.seeds)
    print(f'Starting {args.mode}: {len(starts)} trajectories', flush=True)
    result = SynchronematicGPU(4).scan(positions, starts, dt=args.dt, transient=args.transient,
                                      duration=args.duration, coefficients=coefficients,
                                      progress=lambda f, t: print(f'{f:.1%}, {t:.1f}s', flush=True))
    result['experiment_labels'] = labels
    result['experiment'] = dict(mode=args.mode, seed=2026102904, source_index=43,
                                changes_equations=args.mode in ['coefficients', 'nonreciprocity', 'solenoidal'],
                                note='Coefficient jitter is a sensitivity test, not a bound on omitted Fourier harmonics.')
    path = args.output or Path('results/synchronematic') / ('n4-' + args.mode)
    save(path, result)
    grouped = []
    levels = labels if labels.ndim == 1 else labels[:, 0]
    for level in np.unique(levels):
        mask = levels == level
        passed = (result['status'][mask] == 0) & (result['lle'][mask] > .002) & (result['lle_windows'][mask].min(axis=1) > .001)
        grouped.append(dict(level=float(level), count=int(mask.sum()), retained=int(passed.sum()),
                            median_lle=float(np.nanmedian(result['lle'][mask])),
                            min_lle=float(np.nanmin(result['lle'][mask])), max_lle=float(np.nanmax(result['lle'][mask]))))
    (path / 'groups.json').write_text(json.dumps(grouped, indent=2))
    print(json.dumps(grouped) if len(grouped)<10 else f'{len(grouped)} groups saved to {path / "groups.json"}', flush=True)


if __name__ == '__main__':
    main()
