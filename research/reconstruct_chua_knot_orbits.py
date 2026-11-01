"""Reconstruct atlas cycles with high-precision shooting and aligned flight samples.

All output is new evidence. Earlier reference records and display arrays remain
unchanged. The resulting polygon is a numerical approximation, not an enclosure.
"""
import argparse, json, time
from pathlib import Path
import mpmath as mp
import numpy as np
from chua_precision_flights import PrecisionFlights

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/chua-knot-atlas'


def rows():
    records=json.loads((ROOT/'results/chua-layer-atlas/family-registry.json').read_text())['families']
    later=json.loads((ROOT/'results/chua-layer-prediction/search-results.json').read_text())['records']
    return records+[dict(id='P'+str(r['peaks_requested']),symmetric=True,regions=r['regions'],
                        targets=r['targets'],representative=r['representative']) for r in later]


def reconstruct(row, digits):
    r=row['representative']; f=PrecisionFlights(row['regions'],row['targets'],row['symmetric'],digits)
    target='.5' if row['symmetric'] else '-.5'
    root=f.boundary(r['beta'],r['q'],target)
    beta=mp.mpf(root['beta']);q=[mp.mpf(x) for x in root['q']]
    ev=f.evaluate(beta,q,False)
    return f,root,ev


def sample(f,root,ev,subdivisions):
    beta=mp.mpf(root['beta']); arcs=[]
    for a in ev['arcs']:
        reg=a['region'];slope=f.m0 if reg==0 else f.m1
        A=mp.matrix([[-f.alpha*(1+slope),f.alpha,0],[1,-1,1],[0,-beta,0]])
        c=mp.matrix([reg*(f.m1-f.m0)/(1+f.m1),0,-reg*(f.m1-f.m0)/(1+f.m1)])
        start=a['start'];duration=a['duration']
        # Every full central crossing is split at x=0. Repeated visits then use
        # the same section-to-section phase grid, avoiding unrelated chord errors.
        cuts=[mp.mpf(0),duration]
        if reg==0 and start[0]*a['end'][0] < -mp.mpf('.5'):
            at=lambda t:(c+mp.expm(A*t)*(start-c))[0]
            cross=mp.findroot(at,(mp.mpf(0),duration),solver='bisect',maxsteps=f.dps*4)
            assert 0<cross<duration
            cuts.insert(1,cross)
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            arcs.append((reg,A,c,c+mp.expm(A*lo)*(start-c),hi-lo))
    xyz=[]; times=[]; total=mp.mpf(0); bound=mp.mpf(0)
    for reg,A,c,start,dt in arcs:
        E=mp.expm(A*(dt/subdivisions));u=start.copy()
        for k in range(subdivisions):
            xyz.append([float(v) for v in u]);times.append(float(total+k*dt/subdivisions))
            # Numerical validity diagnostic, not an all-time certificate.
            bound=max(bound,abs(u[0])-1 if reg==0 else 1-reg*u[0])
            u=c+E*(u-c)
        total+=dt
    if f.sign==-1:
        xyz += [[-v for v in p] for p in xyz.copy()]
        times += [float(total)+t for t in times.copy()]
        total*=2
    xyz.append(xyz[0]);times.append(float(total))
    return np.array(xyz),np.array(times),dict(half_or_full_arcs=len(arcs),subdivisions=subdivisions,
        maximum_sampled_region_violation=float(bound),closed_by='periodic shooting; inversion for symmetric second half')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--digits',type=int,default=60)
    p.add_argument('--subdivisions',type=int,default=96)
    p.add_argument('--ids',nargs='*')
    args=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    for row in rows():
        if args.ids and row['id'] not in args.ids:continue
        label=f'{row["id"]}-d{args.digits}-n{args.subdivisions}'
        target=OUT/(label+'.json')
        if target.exists() and (OUT/(label+'.npz')).exists():continue
        before=time.monotonic();f,r,ev=reconstruct(row,args.digits)
        xyz,t,diag=sample(f,r,ev,args.subdivisions)
        r.update(id=row['id'],symmetric=row['symmetric'],regions=row['regions'],targets=row['targets'],
                 alpha='9',m0='-8/7',m1='-5/7',gamma='0',sampling=diag,
                 source_beta=str(row['representative']['beta']),
                 beta_change=float(abs(mp.mpf(r['beta'])-mp.mpf(str(row['representative']['beta'])))),
                 elapsed_seconds=time.monotonic()-before)
        np.savez_compressed(OUT/(label+'.npz'),xyz=xyz,t=t)
        target.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
        print(row['id'],'samples',len(t),'residual',r['residual'],'seconds',round(r['elapsed_seconds'],2),flush=True)


if __name__=='__main__':main()
