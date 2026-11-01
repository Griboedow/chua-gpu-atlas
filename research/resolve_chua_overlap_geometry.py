"""Resolve contacts of two independently continued Chua stability bands.

F and S denote the fold and symmetry-loss edges of the symmetric representative.
A contact is a simultaneous marginal event on TWO distinct periodic orbits;
it is not a collision of those orbits in state space.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import root, brentq

from continue_chua_slope_families import CoordinateFlights, valid, extrema_check, dump
from check_chua_slope_families import verify

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'results/chua-slope-families'
ROLES = {'F': ('fold', 1), 'S': ('symmetry', -1), 'D': ('doubling', -1)}
CONTACTS = ['FS', 'FF', 'SS', 'SF']
EXTENDED = ['FD', 'DF', 'SD', 'DS', 'DD']
PARAMETERS = ['alpha', 'beta', 'm0', 'm1']


def full_template(t):
    return {**t, 'regions':t['regions']+[-r for r in t['regions'][1:]],
            'targets':t['targets']+[-r for r in t['targets']], 'symmetric':False}


def checked_record(f, m0, q, expected):
    r=f.record(m0,q,True);valid(r);r.update(extrema_check(f,r))
    if not f.symmetric:
        assert r['full_cycle_peaks']%2==0
        r['full_cycle_peaks']//=2
    if r['full_cycle_peaks']!=expected:
        raise RuntimeError('A root changed its physical peak family')
    return r


def asymmetric_edge(t, row):
    """Switch at the symmetric branch point, then solve the first flip edge."""
    s=row['symmetry'];width=row['fold']['m0']-s['m0']
    full=full_template(t)
    f=CoordinateFlights(s,full['regions'],full['targets'],symmetric=False)
    half=np.array(s['q'][2:])
    q=np.r_[s['q'][:2],half[:-1],half[-1]+half[0],half[1:]]
    direction=np.linalg.svd(f.evaluate(s['m0'],q)[1])[2][-1]
    start=s['m0']-.3*width
    v=f.solve(start,q+.05*direction)
    trial=checked_record(f,start,v,t['reference_peak_count'])
    if trial['spectral_radius']>=1:
        raise RuntimeError('Branch switch did not find an attracting asymmetric orbit')
    m0,q=f.boundary(s['m0']-.52*width,v,-1)
    r=checked_record(f,m0,q,t['reference_peak_count'])
    if not .25 < (s['m0']-m0)/width < .9:
        raise RuntimeError('First doubling edge left the expected continued neighborhood')
    return r


class Contact:
    def __init__(self, kind, m1, code):
        self.kind, self.m1, self.code = kind, float(m1), code
        self.base = json.loads((SOURCE/'anchors.json').read_text())
        self.templates = [full_template(t) if role=='D' else t for t,role in zip(self.base,code)]
        self.n = [len(t['regions'])+2 for t in self.templates]

    def functions(self, w):
        p = dict(alpha=9., beta=100/7, m0=float(w[-2]), m1=self.m1)
        p[self.kind] = float(w[-1])
        return [CoordinateFlights(p, t['regions'], t['targets'], symmetric=t.get('symmetric',True)) for t in self.templates]

    def system(self, w, jacobian=False):
        fs = self.functions(w)
        qs = [w[:self.n[0]], w[self.n[0]:sum(self.n)]]
        blocks = []
        for f, q, role in zip(fs, qs, self.code):
            blocks.append(f.boundary_system(np.r_[q, w[-2]], ROLES[role][1], jacobian))
        if not jacobian:
            return np.concatenate(blocks)
        J = np.zeros((len(w), len(w)))
        offset = column = 0
        for n, block in zip(self.n, blocks):
            J[offset:offset+n+1, column:column+n] = block[:, :-1]
            J[offset:offset+n+1, -2] = block[:, -1]
            offset += n+1
            column += n
        h = 2e-6
        wp = w.copy(); wm = w.copy(); wp[-1] += h; wm[-1] -= h
        J[:, -1] = (self.system(wp)-self.system(wm))/(2*h)
        return J

    def seed(self):
        data = [json.loads((SOURCE/(t['id']+'-'+self.kind+'.json')).read_text()) for t in self.templates]
        j = min(range(len(data[0]['m1'])), key=lambda k: abs(data[0]['m1'][k]-self.m1))
        maps = [{tuple(r['index']): r for r in d['nodes']} for d in data]
        choices = []
        for index in maps[0].keys() & maps[1].keys():
            if index[1] != j:
                continue
            rows=[m[index] for m in maps]
            positions=[row[ROLES[role][0]]['m0'] if role!='D' else
                       row['symmetry']['m0']-.52*row['width_m0'] for row,role in zip(rows,self.code)]
            choices.append((abs(positions[0]-positions[1]), rows))
        _, rows = min(choices, key=lambda x:x[0])
        edges=[asymmetric_edge(t,row) if role=='D' else row[ROLES[role][0]]
               for t,row,role in zip(self.base,rows,self.code)]
        return np.r_[edges[0]['q'], edges[1]['q'],
                     np.mean([r['m0'] for r in edges]), edges[0][self.kind]]

    def solve(self, guess=None):
        guess = self.seed() if guess is None else np.array(guess)
        ans = root(self.system, guess, jac=lambda w:self.system(w, True), tol=1e-11)
        w = ans.x
        error = float(max(abs(self.system(w))))
        if error > 2e-9 or not np.isfinite(w).all():
            raise RuntimeError(f'Contact residual {error}: {ans.message}')
        qs = [w[:self.n[0]], w[self.n[0]:sum(self.n)]]
        records = []
        for t, f, q, role in zip(self.templates, self.functions(w), qs, self.code):
            r = checked_record(f,w[-2],q,t['reference_peak_count'])
            r.update(id=t['id'], role=ROLES[role][0])
            records.append(r)
        return dict(kind=self.kind, code=self.code, m1=self.m1,
                    control=float(w[-1]), m0=float(w[-2]), w=w.tolist(),
                    residual=error, cycles=records)

    def tangent(self, record, h=2e-6):
        w=np.array(record['w'])
        plus=Contact(self.kind,self.m1+h,self.code)
        minus=Contact(self.kind,self.m1-h,self.code)
        forcing=(plus.system(w)-minus.system(w))/(2*h)
        return -np.linalg.solve(self.system(w,True),forcing)


def edges_at(kind, control, m1, include_asymmetric=False):
    """Correct all four physical edges from nearby existing sheet records."""
    result = {}
    for t in json.loads((SOURCE/'anchors.json').read_text()):
        d=json.loads((SOURCE/(t['id']+'-'+kind+'.json')).read_text())
        source=min(d['nodes'],key=lambda r: ((r['representative'][kind]-control)/(.04 if kind=='alpha' else .08))**2
                   +((r['representative']['m1']-m1)/.008)**2)
        params={k:source['representative'][k] for k in PARAMETERS}
        params.update({kind:float(control),'m1':float(m1)})
        f=CoordinateFlights(params,t['regions'],t['targets'])
        row={}
        for role,target in [ROLES['F'],ROLES['S']]:
            old=source[role]
            m0,q=f.boundary(old['m0'],np.array(old['q']),target)
            r=f.record(m0,q,True);valid(r);r.update(extrema_check(f,r))
            if r['full_cycle_peaks']!=t['reference_peak_count']:
                raise RuntimeError('Edge changed peak family')
            row[role]=r
        if include_asymmetric:
            row['doubling']=asymmetric_edge(t,row)
        result[t['id']]=row
    return result


def trace(out, nodes, extended=False):
    center=-5/7
    slopes=np.unique(np.r_[np.linspace(-.79,-.674,nodes),center])
    for kind in ['alpha','beta']:
        for code in EXTENDED if extended else CONTACTS:
            path=out/(kind+'-'+code+'.json')
            if path.exists():
                print('exists',path.name,flush=True);continue
            started=time.monotonic()
            anchor=Contact(kind,center,code).solve()
            records=[anchor]; failures=[]
            for values in [slopes[slopes<center][::-1], slopes[slopes>center]]:
                previous=anchor
                for slope in values:
                    problem=Contact(kind,slope,code)
                    old=Contact(kind,previous['m1'],code)
                    guess=np.array(previous['w'])+(slope-previous['m1'])*old.tangent(previous)
                    try:
                        row=problem.solve(guess)
                        if max(abs(a['period']/b['period']-1) for a,b in zip(row['cycles'],previous['cycles']))>.04:
                            raise RuntimeError('Large period jump along contact curve')
                        row['tangent']=problem.tangent(row).tolist()
                        records.append(row);previous=row
                    except (RuntimeError,ValueError,np.linalg.LinAlgError) as e:
                        failures.append(dict(m1=float(slope),error=str(e)))
                        break
            anchor['tangent']=Contact(kind,center,code).tangent(anchor).tolist()
            records.sort(key=lambda r:r['m1'])
            dump(path,dict(kind=kind,code=code,records=records,failures=failures,
                          elapsed_seconds=time.monotonic()-started,
                          scope='Numerically solved simultaneous marginal events on distinct primary cycles; D refers to an asymmetric cycle. No orbit collision or exhaustive coexistence boundary is asserted.'))
            print(kind,code,len(records),'failures',len(failures),'seconds',round(time.monotonic()-started,1),flush=True)


def turning(out, extended=False):
    output=[]
    for kind in ['alpha','beta']:
        for code in EXTENDED if extended else CONTACTS:
            d=json.loads((out/(kind+'-'+code+'.json')).read_text())
            rows=d['records']
            i=min(range(len(rows)),key=lambda i: rows[i]['control'] if kind=='alpha' else -rows[i]['control'])
            if i==0 or i==len(rows)-1:
                raise RuntimeError('No interior turning point bracket')
            seed=rows[i]
            seed_tangent=Contact(kind,seed['m1'],code).tangent(seed)
            def evaluated(m1):
                problem=Contact(kind,m1,code)
                row=problem.solve(np.array(seed['w'])+(m1-seed['m1'])*seed_tangent)
                return row,problem.tangent(row)
            m1=brentq(lambda v:evaluated(v)[1][-1],rows[i-1]['m1'],rows[i+1]['m1'],xtol=5e-13)
            r,tangent=evaluated(m1)
            curvatures=[]
            for h in [2e-4,1e-4,5e-5]:
                a,ta=evaluated(m1-h);b,tb=evaluated(m1+h)
                curvatures.append(dict(step=h,second_derivative=float((tb[-1]-ta[-1])/(2*h)),
                    position_second_derivative=float((a['control']-2*r['control']+b['control'])/h**2)))
            r.update(tangent=tangent.tolist(),curvatures=curvatures)
            # Independent full-field checks of both marginal cycles at the tip.
            for cycle in r['cycles']:
                cycle['independent']=verify(cycle)
            output.append(r)
            print('turning',kind,code,'control',r['control'],'m1',m1,'curvature',curvatures[-1]['second_derivative'],flush=True)
    dump(out/('turning-extended.json' if extended else 'turning-points.json'),dict(records=output,
        scope='Extrema of the projection of contact curves onto alpha or beta; not new local orbit bifurcation types.'))


def predictions(out):
    """Use tip curvature only; independently correct held-out horizontal cuts."""
    tips=json.loads((out/'turning-points.json').read_text())['records']
    results=[]
    for tip in tips:
        kind,code=tip['kind'],tip['code']
        curv=tip['curvatures'][-1]['second_derivative']
        direction=np.sign(curv)
        for delta in [1e-5,4e-5,1.6e-4,6.4e-4]:
            target=tip['control']+direction*delta
            predicted=np.sqrt(2*delta/abs(curv))
            corrected=[]
            for sign in [-1,1]:
                def objective(slope):
                    return Contact(kind,slope,code).solve(tip['w'])['control']-target
                left=tip['m1'];right=tip['m1']+sign*2*predicted
                slope=brentq(objective,min(left,right),max(left,right),xtol=2e-13)
                row=Contact(kind,slope,code).solve(tip['w'])
                corrected.append(row)
            spacing=corrected[1]['m1']-corrected[0]['m1']
            results.append(dict(kind=kind,code=code,delta_control=delta,target_control=target,
                predicted_separation=2*predicted,corrected_separation=spacing,
                relative_separation_error=float(abs(spacing/(2*predicted)-1)),
                roots=corrected))
        print('predictions',kind,code,'maximum relative error',max(r['relative_separation_error'] for r in results if r['kind']==kind and r['code']==code),flush=True)
    dump(out/'predictions.json',dict(records=results,
        method='Local Taylor curvature at the independently solved tip predicts the separation 2*sqrt(2*delta/abs(curvature)). Horizontal-cut roots are corrected independently; they are not used to fit the prediction.',
        scope='A standard geometric square-root law checked for these particular contact curves; no new universal scaling law is claimed.'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/chua-overlap-geometry')
    parser.add_argument('--stage',choices=['contacts','turning','predictions'],required=True)
    parser.add_argument('--extended',action='store_true')
    parser.add_argument('--nodes',type=int,default=65)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if args.stage=='contacts':trace(args.output,args.nodes,args.extended)
    elif args.stage=='turning':turning(args.output,args.extended)
    elif args.stage=='predictions':predictions(args.output)


if __name__=='__main__':main()
