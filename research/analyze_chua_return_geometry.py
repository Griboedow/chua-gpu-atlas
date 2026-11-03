"""Compare the recorded predictions with return geometry measurements."""
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from chua_precision_flights import PrecisionFlights
from predict_chua_layers import parent_at
from study_chua_return_geometry import ROOT, OUT, IDS, read


def sections(row, root):
    f=PrecisionFlights(row['regions'],row['targets'],True,70)
    beta=mp.mpf(root['beta'])
    ev=f.evaluate(beta,[mp.mpf(v) for v in root['q']],False)
    points=[np.array([float(root['q'][0]),float(root['q'][1])])]
    times=[0.]
    elapsed=mp.mpf(0)
    for a in ev['arcs']:
        dt=a['duration']
        if a['region']==0 and a['start'][0]*a['end'][0]<-mp.mpf('.5'):
            A=mp.matrix([[-f.alpha*(1+f.m0),f.alpha,0],[1,-1,1],[0,-beta,0]])
            # Exact affine formula in high working precision. A double-precision
            # bracket seeds Newton; endpoints are not extra section returns.
            def state(t):return mp.expm(A*t)*a['start']
            def value(t):return state(t)[0]
            af=np.array(A.tolist(),float);uf=np.array(list(a['start']),float)
            initial=brentq(lambda t:(expm(af*t)@uf)[0],0,float(dt),xtol=1e-14)
            t=mp.findroot(value,mp.mpf(str(initial)),solver='newton',
                          df=lambda t:(A*state(t))[0],tol=mp.mpf('1e-60'))
            assert 0<t<dt
            u=state(t)
            sign=mp.sign((A*u)[0])
            points.append(np.array([float(sign*u[1]),float(sign*u[2])]))
            times.append(float(elapsed+t))
        elapsed+=dt
    return np.array(points),np.array(times)


def main():
    known=read(ROOT/'results/chua-layer-prediction/analysis.json')['rows']
    rows=[]
    for k,id in enumerate(IDS):
        reference=read(ROOT/f'results/chua-knot-atlas/{id}-d70-n192.json')
        companion_path=OUT/f'{id}-companion-d70.json'
        if not companion_path.exists():continue
        companion=read(companion_path)
        p,t=sections(reference,reference)
        s,st=sections(reference,companion['companion'])
        assert p.shape==s.shape
        parent=parent_at(float(reference['beta']))
        distance=np.linalg.norm(p-np.array(parent['q'][:2]),axis=1)
        separation=np.linalg.norm(p-s,axis=1)
        anchor=int(np.argmax(distance))
        order=np.roll(np.arange(len(p)),-anchor)
        r=dict(id=id,n=k+1,peaks=10*(k+1)+5,period=float(reference['period']),
               width=known[k]['width'],parent_period=parent['period'],
               parent_expansion=parent['full_multipliers'][0][0],
               companion_half_multiplier=float(companion['companion']['half_multiplier']),
               anchor_rule='Largest distance to the short saddle among normalized x=0 section crossings',
               original_anchor_index=anchor,section_points=p[order].tolist(),
               companion_section_points=s[order].tolist(),
               section_distance_to_parent=distance[order].tolist(),
               section_pair_separation=separation[order].tolist(),
               minimum_pair_separation=float(separation.min()),maximum_pair_separation=float(separation.max()),
               anchor_pair_separation=float(separation[anchor]))
        for subdivisions,digits in [(8,70),(32,85),(128,85)]:
            path=OUT/f'{id}-d{digits}-s{subdivisions}.json'
            if path.exists():
                g=read(path)
                r[f'gain_s{subdivisions}']=g['maximum_gain']
                r[f'phase_gain_ratio_s{subdivisions}']=g['phase_gain_ratio']
        rows.append(r)
        print(id,'pair min/max/anchor',r['minimum_pair_separation'],r['maximum_pair_separation'],r['anchor_pair_separation'],flush=True)
    (OUT/'sequence-analysis.json').write_text(json.dumps(dict(rows=rows),indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
