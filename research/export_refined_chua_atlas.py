"""Publish compact atlas evidence without the large GPU peak matrices."""
import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact(value):
    if isinstance(value, dict):
        return {k: compact(v) for k, v in value.items() if k != 'arcs'}
    if isinstance(value, list):
        return [compact(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert args.source.resolve() != args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    provenance = []
    for path in sorted(args.source.glob('*.json')):
        if path.name.endswith('-partial.json'):
            continue
        data = compact(json.loads(path.read_text()))
        target = args.output / path.name
        target.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8', newline='\n')
        provenance.append(dict(source_name=path.name, source_sha256=digest(path), published_name=target.name,
            transformation='Removed sampled arc geometry; retained flight times, itineraries and checks.'))
    for path in sorted(args.source.glob('*.npz')):
        arrays = dict(np.load(path))
        if path.name == 'offcenter-raw.npz':
            arrays.pop('peaks', None)
        elif path.name.startswith('overlap-basins'):
            arrays = {k: v for k, v in arrays.items() if k in ['y', 'z', 'indices', 'initial', 'labels', 'status', 'counts']}
        target = args.output / path.name
        np.savez_compressed(target, **arrays)
        provenance.append(dict(source_name=path.name, source_sha256=digest(path), published_name=target.name,
            retained_arrays=list(arrays), source_bytes=path.stat().st_size, published_bytes=target.stat().st_size))
    report = dict(files=provenance, python=platform.python_version(),
        libraries={name: importlib.metadata.version(name) for name in ['numpy', 'scipy', 'matplotlib', 'plotly', 'pyopencl']},
        scope='Compact evidence for a finite numerical catalog. Raw GPU scans can be regenerated with the research scripts.')
    (args.output / 'raw-provenance.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8', newline='\n')
    print('Exported', len(provenance), 'files;', sum(p.stat().st_size for p in args.output.iterdir() if p.is_file()), 'bytes')


if __name__ == '__main__':
    main()
