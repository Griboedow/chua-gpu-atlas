"""Select off-center candidates, coalescing nearby fragments of the GPU screen.

This is a discovery heuristic, not a complete window enumeration. Retain the
largest fragment in each group; discard one-node groups and already mapped
symmetric representatives of the same candidate peak count.
"""
import argparse
import json
from pathlib import Path


def select(out):
    summary = json.loads((out / 'offcenter-summary.json').read_text())
    sheets = [json.loads(p.read_text()) for p in sorted(out.glob('S*-sheet.json'))]
    candidates = []
    for j, section in enumerate(summary['slices']):
        groups = []
        for window in section['windows']:
            if groups and window['start'] - groups[-1][-1]['end'] < .0002:
                groups[-1].append(window)
            else:
                groups.append([window])
        for group in groups:
            window = max(group, key=lambda w: w['nodes'])
            if window['nodes'] < 2:
                continue
            beta = (window['start'] + window['end']) / 2
            hint = int(min(window['periods']))
            known = False
            for sheet in sheets:
                node = min(sheet['nodes'], key=lambda n: abs(n['alpha'] - 9) + abs(n['m0'] - section['m0']))
                width = node['fold']['beta'] - node['symmetry']['beta']
                if hint == sheet['peak_count'] and node['symmetry']['beta'] - .4 * width <= beta <= node['fold']['beta'] + 1e-5:
                    known = True
            if not known:
                candidates.append(dict(slice=j, m0=section['m0'], beta=beta,
                    period_hint=hint, nodes=window['nodes'], start=window['start'], end=window['end']))
    return candidates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--check', action='store_true', help='Compare with the recorded candidate list without replacing it.')
    args = parser.parse_args()
    rows = select(args.output)
    path = args.output / 'offcenter-candidates.json'
    if args.check:
        old = json.loads(path.read_text())
        assert rows == old, 'Selection differs from the recorded candidates'
        print('Matched', len(rows), 'recorded candidates')
    else:
        path.write_text(json.dumps(rows, indent=2) + '\n', encoding='utf-8', newline='\n')
        print('Selected', len(rows), 'candidates')


if __name__ == '__main__':
    main()
