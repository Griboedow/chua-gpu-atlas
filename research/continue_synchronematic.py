"""Resolve regular windows in the explicit gradient-to-full-field homotopy.

This changes the interaction law. It is a mechanism experiment, not another
parameter survey of the unmodified physical model. The same random starts are
used at each eta to make the finite ensemble comparison interpretable.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from search_synchronematic import save
from synchronematic import COEFFICIENTS, SynchronematicGPU, gradient_coefficients


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lower', type=float, default=.93)
    parser.add_argument('--upper', type=float, default=1.18)
    parser.add_argument('--levels', type=int, default=126)
    parser.add_argument('--seeds', type=int, default=8)
    parser.add_argument('--transient', type=float, default=2000)
    parser.add_argument('--duration', type=float, default=4000)
    parser.add_argument('--dt', type=float, default=.005)
    parser.add_argument('--output', type=Path, default=Path('results/synchronematic/n4-transition'))
    args = parser.parse_args()
    if args.levels < 2 or args.seeds < 1 or args.upper <= args.lower:
        parser.error('Require increasing limits, at least two levels and one seed')
    original = np.load('results/synchronematic/n4-scout/scan.npz')
    levels = np.linspace(args.lower, args.upper, args.levels)
    rng = np.random.default_rng(2026102919)
    initial = rng.uniform(-np.pi, np.pi, (args.seeds, 8))
    initial[0] = original['starts'][43]
    starts = np.tile(initial, (len(levels), 1))
    positions = np.broadcast_to(original['positions'][43], (len(starts), 4, 2)).copy()
    gradient = gradient_coefficients()
    coefficients = gradient + levels[:, None, None] * (COEFFICIENTS - gradient)
    coefficients = np.repeat(coefficients, args.seeds, axis=0)
    labels = np.column_stack((np.repeat(levels, args.seeds), np.tile(np.arange(args.seeds), len(levels))))
    print(f'Starting {len(starts)} transition trajectories', flush=True)
    result = SynchronematicGPU(4).scan(
        positions, starts, dt=args.dt, transient=args.transient, duration=args.duration,
        coefficients=coefficients,
        progress=lambda f, t: print(f'{f:.1%}, {t:.1f}s', flush=True))
    result['experiment_labels'] = labels
    result['experiment'] = dict(
        mode='gradient-to-full-field continuation', seed=2026102919,
        source_index=43, changes_equations=True, common_starts=True,
        equation='F_eta=G+eta*(F-G), G=-grad(V), F=literal printed Fourier model')
    save(args.output, result)
    groups = []
    for index, eta in enumerate(levels):
        ids = slice(index * args.seeds, (index + 1) * args.seeds)
        passed = ((result['status'][ids] == 0) & (result['lle'][ids] > .002)
                  & (result['lle_windows'][ids].min(axis=1) > .001)
                  & (result['stats'][ids, -1] > 1e-5))
        groups.append(dict(eta=float(eta), retained=int(passed.sum()), count=args.seeds,
                           median_lle=float(np.nanmedian(result['lle'][ids])),
                           minimum_lle=float(np.nanmin(result['lle'][ids])),
                           maximum_lle=float(np.nanmax(result['lle'][ids]))))
    (args.output / 'groups.json').write_text(json.dumps(groups, indent=2))
    print(f'Saved {len(groups)} parameter groups', flush=True)


if __name__ == '__main__':
    main()
