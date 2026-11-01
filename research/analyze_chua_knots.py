"""Numerical knot classification of independently reconstructed atlas orbits."""
import argparse,json,time,hashlib,platform
from pathlib import Path
import numpy as np
import spherogram
from chua_knot_diagrams import analyze

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/chua-knot-atlas'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pattern',default='*-d70-n192.npz')
    p.add_argument('--seeds',type=int,nargs='+',default=[0,2])
    args=p.parse_args();dest=OUT/'diagrams';dest.mkdir(exist_ok=True)
    for path in sorted(OUT.glob(args.pattern)):
        xyz=np.load(path)['xyz']
        for seed in args.seeds:
            target=dest/f'{path.stem}-s{seed}.json'
            if target.exists():continue
            before=time.monotonic();record,link=analyze(xyz,seed)
            record.update(orbit_file=path.name,orbit_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                          id=path.stem.split('-')[0],sample_count=len(xyz),
                          spherogram_version=spherogram.__version__,python_version=platform.python_version(),
                          elapsed_seconds=time.monotonic()-before,
                          evidence='Numerical polygon; exact diagram invariants; no interval isotopy certificate.')
            target.write_text(json.dumps(record,separators=(',',':'))+'\n',encoding='utf-8')
            print(path.stem,'seed',seed,record['crossings'],'->',record['simplified_crossings'],
                  record['knot_type'],'seconds',round(record['elapsed_seconds'],1),flush=True)


if __name__=='__main__':main()
