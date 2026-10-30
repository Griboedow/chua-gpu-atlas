"""Export compact census evidence while recording hashes of full raw scans."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,default=ROOT/'results/chua-layer-census')
    args=p.parse_args();src=args.source;out=args.output;out.mkdir(parents=True,exist_ok=True)
    if src.resolve()==out.resolve():raise ValueError('Export destination must differ from raw source')
    provenance=[]
    for name in ['line','tight','adding']:
        path=src/(name+'.npz');raw=np.load(path)
        arrays={k:raw[k] for k in raw.files if k!='peaks'}
        if name=='line':arrays['period']=np.load(src/'line-recurrence.npz')['period']
        np.savez_compressed(out/path.name,**arrays)
        provenance.append(dict(name=path.name,raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                               raw_bytes=path.stat().st_size,exported_arrays=list(arrays),omitted_arrays=['peaks'],
                               period_note='Recomputed at tolerance 5e-4 from saved peaks' if name=='line' else 'Original 5e-4 recurrence'))
    for path in src.glob('*.json'):
        # Preserve original reports; do not silently rewrite trial hashes.
        (out/path.name).write_text(path.read_text(encoding='utf-8'),encoding='utf-8',newline='\n')
    for name in ['cycles','cycles-extra','adding-cycles','coexistence','atlas-comparison','line-recurrence']:
        (out/(name+'.npz')).write_bytes((src/(name+'.npz')).read_bytes())
    (out/'raw-provenance.json').write_text(json.dumps(dict(scans=provenance,
        note='Full raw peak matrices are omitted from Git. Compact files retain all final states, Lyapunov windows, classifications, counts, and accepted peak periods. Full raw matrices can be regenerated with the documented survey commands.'),indent=2)+'\n',encoding='utf-8',newline='\n')
    print('exported bytes',sum(p.stat().st_size for p in out.iterdir() if p.is_file()))


if __name__=='__main__':main()
