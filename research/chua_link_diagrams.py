"""Joint polygonal link diagrams; all components share one spatial transform.

Floating-point crossing extraction is checked with projections, resolution and
an independent solid-angle Gauss integral. Diagram reductions are combinatorial;
they do not certify ambient isotopy to an exact ODE trajectory.
"""
import random
import itertools,math
from collections import Counter
import numpy as np
from numba import njit
from spherogram import Link
from chua_knot_diagrams import camera


def clean(curve):
    p=np.asarray(curve,float)
    return p[:-1] if np.linalg.norm(p[-1]-p[0])<1e-13 else p


@njit(cache=True)
def crossings(points,ends,components,local,sizes):
    rows=[]
    for i in range(len(points)):
        a=points[i];u=points[ends[i]]-a
        for j in range(i+1,len(points)):
            if components[i]==components[j] and (local[j]-local[i]==1 or local[j]-local[i]==sizes[components[i]]-1):continue
            c=points[j];v=points[ends[j]]-c
            if max(a[0],a[0]+u[0])<min(c[0],c[0]+v[0]) or max(c[0],c[0]+v[0])<min(a[0],a[0]+u[0]):continue
            if max(a[1],a[1]+u[1])<min(c[1],c[1]+v[1]) or max(c[1],c[1]+v[1])<min(a[1],a[1]+u[1]):continue
            det=u[0]*v[1]-u[1]*v[0]
            if abs(det)<1e-30:continue
            w=c-a;s=(w[0]*v[1]-w[1]*v[0])/det;t=(w[0]*u[1]-w[1]*u[0])/det
            if 0<s<1 and 0<t<1:
                height=a[2]+s*u[2]-c[2]-t*v[2]
                rows.append((i,j,s,t,height,det,components[i],components[j]))
    return rows


def diagram(curves,seed=0):
    curves=[clean(c) for c in curves];sizes=np.array([len(p) for p in curves]);starts=np.r_[0,np.cumsum(sizes)[:-1]]
    p=np.vstack(curves);scale=np.ptp(p,axis=0);scale[scale==0]=1
    p=(p-p.mean(axis=0))/scale@camera(seed)
    assert np.isfinite(p).all()
    components=np.repeat(np.arange(len(curves)),sizes);local=np.concatenate([np.arange(n) for n in sizes])
    ends=np.concatenate([s+(np.arange(n)+1)%n for s,n in zip(starts,sizes)])
    rows=crossings(p,ends,components,local,sizes);labels={};previous={};free=0;offset=0;spacings=[];matrix=np.zeros((len(curves),len(curves)),int)
    for component in range(len(curves)):
        events=[]
        for j,r in enumerate(rows):
            for side in [0,1]:
                if int(r[6+side])==component:events.append((r[side]-starts[component]+r[side+2],j,side))
        events.sort();n=len(events)
        if not n:free+=1;continue
        for i,(_,j,side) in enumerate(events):labels[j,side]=offset+i;previous[offset+i]=offset+(i-1)%n
        spacings.extend(events[i+1][0]-events[i][0] for i in range(n-1));offset+=n
    pd=[]
    for i,r in enumerate(rows):
        under,over=(1,0) if r[4]>0 else (0,1);u=labels[i,under];o=labels[i,over]
        det=r[5] if under==0 else -r[5]
        pd.append([previous[u],o,u,previous[o]] if det<0 else [previous[u],previous[o],u,o])
        a,b=int(r[6]),int(r[7])
        if a!=b:
            sign=int(np.sign(r[4]*r[5]));matrix[a,b]+=sign;matrix[b,a]+=sign
    assert np.all(matrix%2==0),'Nonintegral diagram linking number'
    matrix=matrix//2;link=Link(pd) if pd else Link([]);link.unlinked_unknot_components=free
    assert len(link.link_components)+free==len(curves)
    assert not pd or link.is_planar()
    return link,dict(seed=seed,components=len(curves),initial_free_components=free,pd_code=pd,
        crossings=len(rows),oriented_linking_matrix=matrix.tolist(),
        minimum_crossing_height=min((abs(r[4]) for r in rows),default=None),
        minimum_vertex_fraction=min((min(r[2],1-r[2],r[3],1-r[3]) for r in rows),default=None),
        minimum_event_spacing=min(spacings,default=None))


def classify(curves,seed=0):
    link,r=diagram(curves,seed);random.seed(9173+seed)
    r['library_initial_linking_matrix']=link.linking_matrix()
    stages=[]
    for mode in ['basic','global']:
        link.simplify(mode,type_III_limit=300)
        stages.append(dict(mode=mode,crossings=len(link.crossings),free_components=link.unlinked_unknot_components,
            remaining_components=len(link.link_components)))
        if not link.crossings:break
    assert len(link.link_components)+link.unlinked_unknot_components==len(curves),'A component was lost'
    zero=not link.crossings;nonzero=bool(np.any(np.array(r['oriented_linking_matrix'])))
    assert not (zero and nonzero),'Unlink reduction contradicts a nonzero linking number'
    polynomial=jones_A(link) if len(link.crossings)<=16 else None
    nontrivial_polynomial=polynomial is not None and polynomial!=unlink_jones(len(curves))
    r.update(reductions=stages,simplified_pd_code=link.PD_code(),simplified_crossings=len(link.crossings),
        final_free_components=link.unlinked_unknot_components,
        jones_A=polynomial,status='unlink' if zero else 'nontrivial_by_linking_number' if nonzero else 'nontrivial_by_jones' if nontrivial_polynomial else 'unresolved',
        scope='Numerical polygon diagram; a retained diagram is not evidence of nontrivial topology by itself')
    return r,link


def unlink_jones(components):
    k=components-1
    return {str(2*k-4*j):(-1)**k*math.comb(k,j) for j in range(k+1)}


def jones_A(link):
    """Exact Kauffman state sum, normalized by writhe; t=A^-4.

    Polynomial inequality proves nontriviality of the extracted diagram.
    Equality with a named link's polynomial does not identify that link.
    """
    pd=link.PD_code();n=len(pd);free=link.unlinked_unknot_components
    if not n:return unlink_jones(free)
    out=Counter()
    for state in itertools.product([0,1],repeat=n):
        parents=list(range(2*n))
        def find(a):
            while parents[a]!=a:a=parents[a]
            return a
        for (a,b,c,d),s in zip(pd,state):
            for u,v in ([(a,b),(c,d)] if s==0 else [(a,d),(b,c)]):parents[find(u)]=find(v)
        loops=len(set(find(a) for a in range(2*n)));base=n-2*sum(state);k=loops-1+free
        for j in range(k+1):out[base+2*k-4*j]+=(-1)**k*math.comb(k,j)
    w=link.writhe()
    return {str(e-3*w):v*(-1 if w%2 else 1) for e,v in sorted(out.items()) if v}


@njit(cache=True)
def triangle(a,b,c):
    a=a/np.sqrt(np.dot(a,a));b=b/np.sqrt(np.dot(b,b));c=c/np.sqrt(np.dot(c,c))
    det=a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0])
    return 2*np.arctan2(det,1+np.dot(a,b)+np.dot(b,c)+np.dot(c,a))


@njit(cache=True)
def polygon_gauss(a,b):
    total=0.
    for i in range(len(a)):
        p=a[i];q=a[(i+1)%len(a)]
        for j in range(len(b)):
            u=b[j];v=b[(j+1)%len(b)]
            total+=triangle(p-u,p-v,q-v)+triangle(p-u,q-v,q-u)
    return total/(4*np.pi)
