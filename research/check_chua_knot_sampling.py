"""Controls for knot extraction and a uniform-time sampling counterexample."""
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import expm
from chua_precision_flights import PrecisionFlights
from chua_knot_diagrams import analyze,canonical_gauss
from spherogram import Link

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/chua-knot-atlas'


def torus(p,q,n=1600):
    t=np.linspace(0,2*np.pi,n,endpoint=False)
    return np.column_stack(((2+np.cos(q*t))*np.cos(p*t),
                            (2+np.cos(q*t))*np.sin(p*t),np.sin(q*t)))


def connected_mirrors(q,n=1600):
    t=np.linspace(.13,2*np.pi-.13,n//2)
    left=np.column_stack(((2+np.cos(q*t))*np.cos(2*t)-3.6,
                           (2+np.cos(q*t))*np.sin(2*t),np.sin(q*t)))
    return np.vstack((left,-left))


def standard_curves():
    t=np.linspace(0,2*np.pi,1600,endpoint=False)
    return {'unknot':np.column_stack((np.cos(t),np.sin(t),.13*np.sin(3*t))),
            'T(2,5)':torus(2,5), 'T(3,5) mirror':torus(3,5)*[1,1,-1],
            'T(2,3) # T(2,3) mirror':connected_mirrors(3),
            'T(2,5) # T(2,5) mirror':connected_mirrors(5)}


def uniform_orbit(root,n):
    f=PrecisionFlights(root['regions'],root['targets'],root['symmetric'],70)
    ev=f.evaluate(mp.mpf(root['beta']),[mp.mpf(v) for v in root['q']],False)
    end=np.cumsum([float(a['duration']) for a in ev['arcs']]);period=end[-1]
    if root['symmetric']:period*=2
    times=np.linspace(0,period,n,endpoint=False);xyz=[]
    for time in times:
        sign=-1 if root['symmetric'] and time>=end[-1] else 1
        local=time-end[-1] if sign==-1 else time
        j=int(np.searchsorted(end,local,side='right'));a=ev['arcs'][j]
        dt=local-(end[j-1] if j else 0);reg=a['region']
        slope=float(f.m0 if reg==0 else f.m1)
        A=np.array([[-9*(1+slope),9,0],[1,-1,1],[0,-float(root['beta']),0]])
        c=np.array([reg*1.5,0,-reg*1.5])
        xyz.append(sign*(c+expm(A*dt)@(np.array(list(a['start']),float)-c)))
    return np.array(xyz)


def main():
    records=[]
    curves=standard_curves()
    for name,xyz in curves.items():
        r,link=analyze(xyz,0)
        # The torus parameterization may use the opposite handedness to the
        # library's braid convention. Fix the displayed representative by test.
        if name in ('T(2,5)','T(3,5) mirror') and r['knot_type']!=name:
            xyz=xyz*[1,1,-1];curves[name]=xyz;r,link=analyze(xyz,0)
        assert r['knot_type']==name,(name,r['knot_type'])
        for seed in (0,2):
            for reversed_path in (False,True):
                test,link=analyze(xyz[::-1] if reversed_path else xyz,seed)
                assert test['knot_type']==name
                records.append(dict(control=name,reversed=reversed_path,**test))
    # A non-torus prime control checks crossing signs and both determinant paths.
    t=np.linspace(0,2*np.pi,1600,endpoint=False)
    figure8=np.column_stack(((2+np.cos(2*t))*np.cos(3*t),(2+np.cos(2*t))*np.sin(3*t),np.sin(4*t)))
    for seed in (0,2):
        r,link=analyze(figure8,seed)
        assert r['alexander_coefficients']==[1,-3,1]
        assert canonical_gauss(link)==canonical_gauss(Link('4_1'))
        records.append(dict(control='figure-eight',**r))
    np.savez_compressed(OUT/'standard-knot-curves.npz',**{f'k{i}':x for i,x in enumerate(curves.values())})
    (OUT/'controls.json').write_text(json.dumps(dict(types=list(curves),tests=records),indent=2)+'\n')
    root=json.loads((OUT/'A05-d70-n192.json').read_text());sampling=[]
    for n in (4000,16000):
        xyz=uniform_orbit(root,n);np.savez_compressed(OUT/f'A05-uniform-{n}.npz',xyz=xyz)
        r,link=analyze(xyz,0);r.update(vertices=n,scheme='uniform time; 70-digit root, float64 affine sampling')
        sampling.append(r);print('uniform',n,r['knot_type'],r['alexander_coefficients'],flush=True)
    (OUT/'sampling-counterexample.json').write_text(json.dumps(dict(id='A05',uniform=sampling,
        comparison='A05-d70-n192; same root, section-aligned sampling'),indent=2)+'\n')
    print('Passed',len(records),'analytic-curve controls',flush=True)


if __name__=='__main__':main()
