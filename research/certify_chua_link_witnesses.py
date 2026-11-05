"""Interval existence and attraction of link witnesses; no topological isotopy claim."""
from search_chua_links import ROOT,OUT,read,save
from fractions import Fraction
import argparse
import mpmath as mp
import numpy as np
from flint import arb,ctx
from flint import arb_mat
from scipy.linalg import expm
from scipy.optimize import brentq
from certify_chua_coexisting_cycles import ReturnMap
from chua_interval_map import col,eye,midpoint,bounds
from chua_precision_flights import PrecisionFlights


class KrawczykRejection(ValueError):
    """A proposed root box failed the required interval inclusion/contraction."""


def rational(text):
    value=Fraction(str(text));return arb(value.numerator)/arb(value.denominator)


class FullReturnMap(ReturnMap):
    """Split every central flight at all zero crossings, then validate each arc.

    Unlike the old four-cycle map, this return can contain intermediate upward
    section crossings. They are retained to exclude a multiple traversal.
    """
    def __init__(self,record,digits,parameters):
        ctx.dps=digits;self.alpha,self.beta,self.m0,self.m1=parameters
        self.record=record;self.schedule_rows=[];q=record['cycle']['q'];state=np.array([0.,*q[:2]])
        for reg,target,duration in zip(record['regions'],[*record['targets'],0],q[2:]):
            A,c,_=self.matrices(None,reg);An=np.array([[float(A[i,j].mid()) for j in range(3)] for i in range(3)]);cn=np.array([float(v.mid()) for v in c.entries()])
            at=lambda t:cn+expm(An*t)@(state-cn)
            cuts=[0.]
            if reg==0:
                grid=np.linspace(1e-8,duration-1e-8,max(3,int(duration/.01)+2));values=[at(t)[0] for t in grid]
                cuts += [brentq(lambda t:at(t)[0],a,b,xtol=1e-14) for a,b,x,y in zip(grid[:-1],grid[1:],values[:-1],values[1:]) if x*y<0]
            cuts.append(duration)
            for i,(a,b) in enumerate(zip(cuts[:-1],cuts[1:])):
                sign=reg if reg else int(np.sign(at((a+b)/2)[0]))
                self.schedule_rows.append((reg,target if i==len(cuts)-2 else 0,b-a,sign))
            state=at(duration);state[0]=target

    def evaluate(self,point,verify_path=False):
        state=col([0,*point]);D=arb_mat([[0,0],[1,0],[0,1]]);arcs=[]
        for reg,target,guess,sign in self.schedule_rows:
            A,c,_=self.matrices(None,reg);start=state;T=self.crossing(A,c,start,target,guess)
            end,E=self.flow(A,c,start,T);v=A*(end-c)
            assert not v[0,0].contains(0)
            F=E*D;D=arb_mat([[F[i,j]-v[i,0]*F[0,j]/v[0,0] for j in range(2)] for i in range(3)])
            end[0,0]=target;state=end
            arc=dict(region=reg,target=target,duration=T,start=start,end=end,speed=v[0,0],sign=sign)
            if verify_path:arc['path']=self.verify_arc(A,sign*c,dict(region=abs(reg),target=sign*target,duration=T,start=sign*start))
            arcs.append(arc)
        assert arcs[0]['sign']==1 and arcs[-1]['sign']==-1 and arcs[-1]['speed']>0
        return col([state[1,0],state[2,0]]),arb_mat([[D[i,j] for j in range(2)] for i in [1,2]]),arcs


def certify(record,digits,radius,parameters,shifted=False):
    f=FullReturnMap(record,digits,parameters);center=col([arb(str(v)) for v in record['cycle']['q'][:2]])
    for _ in range(8):
        value,J,_=f.evaluate(center.entries());center=midpoint(center-(J-eye(2)).inv()*(value-center))
    if shifted:center[0,0]+=arb('1e-12')
    box=col([arb(v,arb(radius)) for v in center.entries()]);value0,J0,_=f.evaluate(center.entries())
    _,J,arcs=f.evaluate(box.entries(),verify_path=True);C=midpoint((J0-eye(2)).inv());R=eye(2)-C*(J-eye(2))
    K=center-C*(value0-center)+R*(box-center);contraction=max(sum(abs(R[i,j]) for j in range(2)).abs_upper() for i in range(2))
    if not (all(box[i,0].contains_interior(K[i,0]) for i in range(2)) and contraction<1):
        raise KrawczykRejection('Krawczyk inclusion/contraction rejected')
    trace,det=J.trace(),J.det();schur=[1-trace+det,1+trace+det,1-det];assert all(s>0 for s in schur)
    upward=[arc for arc in arcs if arc['target']==0 and arc['speed']>0]
    for arc in upward[:-1]:assert any(not (arc['end'][i+1,0]-box[i,0]).contains(0) for i in range(2)),'Multiple traversal not excluded'
    return dict(id=record['id'],passed=True,root=[bounds(v) for v in box.entries()],krawczyk_image=[bounds(v) for v in K.entries()],
        contraction=bounds(contraction),period=bounds(sum(a['duration'] for a in arcs)),schur_margins=[bounds(v) for v in schur],
        return_trace=bounds(trace),return_determinant=bounds(det),primitive_period=True,upward_crossings=len(upward),
        upward_section_points=[[bounds(a['end'][i,0]) for i in [1,2]] for a in upward],
        section_itinerary=[dict(region=a['region'],target=a['target'],sign=a['sign']) for a in arcs],
        path_cells=[a['path']['accepted_cells'] for a in arcs])


def full_record(root,id=None):
    times=root['q'][2:];regions=root['regions'];targets=root['targets']
    if root['symmetric']:
        times=times+times;regions=regions+[-v for v in regions];targets=targets+[0]+[-v for v in targets]
    return dict(id=id or root['id'],cycle=dict(q=[float(v) for v in root['q'][:2]+times]),regions=regions,targets=targets)


def inverted_record(root):
    assert not root['symmetric'];mp.mp.dps=100
    f=PrecisionFlights(root['regions'],root['targets'],False,100)
    def number(s):
        value=Fraction(s);return mp.mpf(value.numerator)/value.denominator
    f.alpha=number(root['alpha']);f.m0=number(root['m0']);f.m1=number(root['m1'])
    q=list(map(mp.mpf,root['q']));ev=f.evaluate(number(root['beta']),q,False)
    cut=root['targets'].index(0)+1;state=ev['arcs'][cut-1]['end']
    regions=root['regions'][cut:]+root['regions'][:cut];targets=root['targets']+[0];targets=targets[cut:]+targets[:cut]
    times=q[2:];times=times[cut:]+times[:cut]
    return dict(id=root['id']+'-inverted',cycle=dict(q=[float(-state[1]),float(-state[2])]+[float(v) for v in times]),regions=[-r for r in regions],targets=[-v for v in targets[:-1]])


def main():
    p=argparse.ArgumentParser();p.add_argument('--digits',type=int,default=140);p.add_argument('--radius',default='1e-60');a=p.parse_args()
    pairs=[('S02-S08',[read(OUT/f'overlaps/S02-S08/{id}-d90-n384.json') for id in ['S02','S08']])]
    for name in ['S01','S02']:
        pairs.append((name+'-daughters',[read(OUT/f'daughters/{name}-a1e-03-d80-n192.json')]))
    results=[]
    for name,roots in pairs:
        ctx.dps=a.digits;params=[rational(roots[0][k]) for k in ['alpha','beta','m0','m1']]
        entries=[full_record(r) for r in roots]
        if len(roots)==1:entries.append(inverted_record(roots[0]))
        records=[]
        for entry in entries:
            result=certify(entry,a.digits,radius=a.radius,parameters=params);records.append(result)
            print(name,entry['id'],'interval sink passed',flush=True)
        # Arb's ball string includes outward rounding of both midpoint and
        # radius; parsing separately printed midpoint/radius need not do so.
        def enclosure(b):return arb(b['ball'])
        aroot=[enclosure(v) for v in records[0]['root']]
        separated=all(any(not (aroot[i]-enclosure(point[i])).contains(0) for i in range(2)) for point in records[1]['upward_section_points']);assert separated
        witness=dict(name=name,parameters={k:roots[0][k] for k in ['alpha','beta','m0','m1','gamma']},digits=a.digits,root_radius=a.radius,records=records,distinct_orbits_by_all_section_points=separated,passed=True)
        save(OUT/f'certificate-{name}-d{a.digits}.json',witness);results.append(witness)
    rejected=False;reason=None
    try:certify(entries[0],a.digits,radius=a.radius,parameters=params,shifted=True)
    except KrawczykRejection as error:rejected=True;reason=type(error).__name__+': '+str(error)
    assert rejected,'Shifted root control was not rejected'
    save(OUT/f'certificates-d{a.digits}.json',dict(witnesses=results,shifted_root_control=dict(rejected=rejected,reason=reason),passed=True,scope='Interval existence, distinctness, attraction and primitive periods; ambient isotopy and link topology are not certified.'))


if __name__=='__main__':main()
