"""Compare literal article coefficients with the authors' NumPy-branch equations.

Only mathematical coefficients are transcribed; the implementation is our own.
The author's script is neither executed nor redistributed. This explicitly
does not reproduce its separate, asymmetric Numba branch or Euler discretization.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from synchronematic import COEFFICIENTS, SynchronematicGPU, REPOSITORY_NUMPY_COEFFICIENTS
from search_synchronematic import save

REPOSITORY_NUMPY = REPOSITORY_NUMPY_COEFFICIENTS
SOURCE_URL = 'https://github.com/slevinskygra/ReducedSynchronematic/blob/016614ad6debe11e0353f0a01b1fad937fc77461/simulation_numba.py'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['reference', 'survey', 'refine'], default='reference')
    parser.add_argument('--source', type=Path, default=Path('results/synchronematic/n4-repository-survey'))
    args = parser.parse_args()
    audit = dict(date='2026-10-29', article_doi='10.1038/s41467-026-68552-8', code_url=SOURCE_URL,
                 article_coefficients=COEFFICIENTS.tolist(), numpy_branch_coefficients=REPOSITORY_NUMPY.tolist(),
                 coefficient_order=['cos Delta', 'sin Delta', 'cos 3Delta', 'sin 3Delta'],
                 main_discrepancy='Orientation sin(3 Delta) coefficient: article -1.359; repository NumPy branch +1.359334.',
                 derivation='cos(a-b)+cos(a+b)=2 cos(a) cos(b); sin(a-b)+sin(a+b)=2 sin(a) cos(b).',
                 additional_observation='The separate Numba reverse-pair expressions omit the phase sin(3 Delta) term and orientation sin(Delta) term present in the forward-pair expressions. That branch is not reproduced here.',
                 interpretation='Sources disagree. This audit does not determine which version was physically intended or used for every published figure. No author contacted.')
    Path('results/synchronematic/source-audit.json').write_text(json.dumps(audit, indent=2))
    if args.mode == 'reference':
        data = np.load('results/synchronematic/n4-geometry/scan.npz')
        p, s = data['positions'][:128], data['starts'][:128]
        high_precision_negative = REPOSITORY_NUMPY.copy()
        high_precision_negative[1, 3] *= -1
        variants = [COEFFICIENTS, high_precision_negative, REPOSITORY_NUMPY]
        positions = np.tile(p, (3, 1, 1));starts = np.tile(s, (3, 1))
        coefficients = np.repeat(np.array(variants), 128, axis=0)
        labels = np.repeat(np.arange(3), 128)
        transient, duration = 2000, 4000
        output = Path('results/synchronematic/n4-source-variants')
    elif args.mode == 'survey':
        data = np.load('results/synchronematic/n4-scout/scan.npz')
        positions, starts = data['positions'], data['starts']
        coefficients = REPOSITORY_NUMPY;labels = np.arange(len(starts))
        transient, duration = 200, 400
        output = Path('results/synchronematic/n4-repository-survey')
    else:
        data = np.load(args.source / 'scan.npz')
        ids = json.loads((args.source / 'summary.json').read_text())['candidate_indices']
        if not ids:
            raise SystemExit('No candidates to refine')
        positions, starts = data['positions'][ids], data['starts'][ids]
        coefficients = data['coefficients'][ids];labels = np.array(ids)
        transient, duration = 2000, 4000
        output = Path('results/synchronematic/n4-repository-refined')
    print(f'Source audit {args.mode}: {len(starts)} trajectories', flush=True)
    result = SynchronematicGPU(4).scan(positions, starts, coefficients=coefficients, dt=.005,
                                      transient=transient, duration=duration,
                                      progress=lambda f, t: print(f'{f:.1%}, {t:.1f}s', flush=True))
    result['source_variant_labels'] = labels
    result['source_audit'] = dict(mode=args.mode, source_url=SOURCE_URL,
                                 variants=['article-rounded-negative', 'repository-precision-negative', 'repository-numpy-positive'])
    save(output, result)
    if args.mode == 'reference':
        groups=[]
        for label in range(3):
            mask=labels==label
            keep=(result['status'][mask]==0)&(result['lle'][mask]>.002)&(result['lle_windows'][mask].min(axis=1)>.001)
            groups.append(dict(variant=result['source_audit']['variants'][label],count=int(mask.sum()),retained=int(keep.sum()),
                               mean_lle=float(np.mean(result['lle'][mask])),min_lle=float(np.min(result['lle'][mask])),max_lle=float(np.max(result['lle'][mask]))))
        (output/'groups.json').write_text(json.dumps(groups,indent=2));print(json.dumps(groups),flush=True)


if __name__ == '__main__':
    main()
