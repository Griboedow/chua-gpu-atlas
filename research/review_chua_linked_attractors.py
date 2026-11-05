"""Audit reconstructed roots, ODE trajectories and replayable link diagrams.

The numerical geometry checks do not enclose an ambient isotopy to the ODE.
"""
import itertools
import random
import subprocess
import importlib.metadata
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from numba import njit
from spherogram import Link
from study_chua_linked_attractors import ROOT, OUT, read, save, digest
from chua_link_diagrams import clean, jones_A, unlink_jones

IDS=['F1a','F1b','F2a','F2b']
LEVELS=[(60,96),(80,192),(90,384)]


@njit(cache=True)
def minimum_segment_distance(a,b):
    """Closest points of two closed polygons, in physical coordinates."""
    best=1e100;where=(-1,-1,0.,0.)
    for i in range(len(a)):
        p=a[i];q=a[(i+1)%len(a)];u=q-p;aa=np.dot(u,u)
        for j in range(len(b)):
            r=b[j];s=b[(j+1)%len(b)]
            lower=0.
            for k in range(3):
                gap=max(0.,min(p[k],q[k])-max(r[k],s[k]),min(r[k],s[k])-max(p[k],q[k]))
                lower+=gap*gap
            if lower>=best:continue
            v=s-r;w=p-r;bb=np.dot(u,v);cc=np.dot(v,v);dd=np.dot(u,w);ee=np.dot(v,w)
            denom=aa*cc-bb*bb
            along=0. if denom<=1e-30 else min(1.,max(0.,(bb*ee-cc*dd)/denom))
            across=(bb*along+ee)/cc
            if across<0:across=0.;along=min(1.,max(0.,-dd/aa))
            elif across>1:across=1.;along=min(1.,max(0.,(bb-dd)/aa))
            delta=w+along*u-across*v;distance=np.dot(delta,delta)
            if distance<best:best=distance;where=(i,j,along,across)
    return np.sqrt(best),where


def field(t,u):
    x,y,z=u;alpha=453/50;beta=71729/5000;m0=-8/7+3/500;m1=-5/7
    h=m1*x+.5*(m0-m1)*(abs(x+1)-abs(x-1))
    return [alpha*(y-x-h),x-y+z,-beta*y]


def main():
    mp.mp.dps=110
    certpath=ROOT/'results/chua-accessibility/four-cycle-certificate-100.json'
    cert=read(certpath);boxes={r['id']:r for r in cert['records']}
    source=ROOT/'results/chua-layer-atlas/four-cycles.json'
    protocol=subprocess.check_output(['git','show','39a84c8:results/chua-linked-attractors/protocol.md'],cwd=ROOT)
    assert protocol.replace(b'\r\n',b'\n')==(OUT/'protocol.md').read_bytes().replace(b'\r\n',b'\n')
    roots=[];curves={}
    for id in IDS:
        data=[read(OUT/f'{id}-d{d}-n{n}.json') for d,n in LEVELS]
        arrays=[np.load(OUT/f'{id}-d{d}-n{n}.npz') for d,n in LEVELS]
        final=data[-1];p=arrays[-1]['xyz'];t=arrays[-1]['t'];curves[id]=clean(p)
        assert final['source_sha256']==digest(source)
        for d in data:
            assert {k:d[k] for k in ['alpha','m0','m1','gamma']}=={k:cert['parameters'][k] for k in ['alpha','m0','m1','gamma']}
            assert mp.mpf(d['beta'])==mp.mpf(71729)/5000 and min(map(mp.mpf,d['q'][2:]))>0
            assert d['spectral_radius']<1
        ratios=[]
        for value,box in zip(final['q'][:2],boxes[id]['root']):
            radius=mp.mpf(box['radius'].strip('[]').split()[0]);ratios.append(float(abs(mp.mpf(value)-mp.mpf(box['midpoint']))/radius))
        assert max(ratios)<1
        refinements=[]
        for first,second,a,b in zip(data[:-1],data[1:],arrays[:-1],arrays[1:]):
            difference=max(abs(mp.mpf(x)-mp.mpf(y)) for x,y in zip(first['q'],second['q']))
            common=np.max(np.linalg.norm(a['xyz']-b['xyz'][::2],axis=1))
            midpoint=np.max(np.linalg.norm(b['xyz'][1::2]-(a['xyz'][:-1]+a['xyz'][1:])/2,axis=1))
            assert common<1e-12
            refinements.append(dict(from_digits=first['working_digits'],to_digits=second['working_digits'],root_max_difference=mp.nstr(difference,15),common_vertex_max_distance=float(common),maximum_midpoint_chord_deviation=float(midpoint)))
        integrations=[]
        for method in ['DOP853','Radau']:
            sol=solve_ivp(field,[0,t[-1]],p[0],method=method,rtol=2e-12,atol=2e-13,max_step=.01,dense_output=True)
            error=float(np.max(np.linalg.norm(sol.sol(t).T-p,axis=1)))
            closure=float(np.linalg.norm(sol.y[:,-1]-p[0]));peaks=len(find_peaks(sol.sol(np.linspace(0,t[-1],12001))[0])[0])
            assert sol.success and error<2e-7 and closure<2e-7 and peaks==5
            integrations.append(dict(method=method,maximum_trajectory_difference=error,closure_error=closure,x_maxima=peaks,nfev=sol.nfev))
        roots.append(dict(id=id,period=final['period'],spectral_radius=final['spectral_radius'],certificate_box_relative_offsets=ratios,refinements=refinements,independent_integrations=integrations))
        print(id,'root box',max(ratios),'ODE differences',[r['maximum_trajectory_difference'] for r in integrations],flush=True)
    allrecords=[];gauss=[];counts=[]
    for d,n in LEVELS:
        analysis=read(OUT/f'analysis-d{d}-n{n}.json');rows=analysis['records'];allrecords+=rows
        assert len(rows)==(45 if d==90 else 30)
        for r in rows:
            assert r['status']=='unlink' and r['simplified_crossings']==0 and r['final_free_components']==len(r['ids'])
            assert not np.any(r['oriented_linking_matrix'])
            for id,h in r['source_sha256'].items():assert digest(OUT/f'{id}-d{d}-n{n}.npz')==h
            # Replay from the saved integer diagram, with its crossing-free components.
            link=Link(r['pd_code']);link.unlinked_unknot_components=r['initial_free_components']
            random.seed(9173+r['seed'])
            for mode in ['basic','global']:
                link.simplify(mode,type_III_limit=300)
                if not link.crossings:break
            assert not link.crossings and link.unlinked_unknot_components==len(r['ids'])
            assert jones_A(link)==unlink_jones(len(r['ids']))
        for r in analysis['gauss_integrals']:
            assert abs(r['value'])<1e-8
            gauss.append(dict(digits=d,subdivisions=n,**r))
        counts.append(dict(digits=d,subdivisions=n,vertices_per_cycle=10*n,diagram_tests=len(rows),maximum_raw_crossings=max(r['crossings'] for r in rows)))
    controls=read(OUT/'controls.json');assert controls['passed'] and len(controls['rows'])==34
    control_replay=[]
    for r in controls['rows']:
        link=Link(r['simplified_pd_code']);link.unlinked_unknot_components=r['final_free_components']
        assert jones_A(link)==r['jones_A']
        if r['expected']=='borromean':assert r['jones_A']!=unlink_jones(3)
        control_replay.append(dict(control=r['control'],subset=r['subset'],seed=r['seed'],passed=True))
    gaps=[]
    for a,b in itertools.combinations(IDS,2):
        distance,location=minimum_segment_distance(curves[a],curves[b]);assert distance>0
        gaps.append(dict(ids=[a,b],minimum_polygon_segment_distance=float(distance),segment_pair_and_fractions=list(location)))
        print(a,b,'polygon gap',distance,flush=True)
    report=dict(protocol_commit='39a84c8',protocol_unchanged=True,source_sha256=digest(source),certificate_sha256=digest(certpath),
        reconstruction_checks=roots,levels=counts,diagram_replays=len(allrecords),control_polynomial_replays=len(control_replay),
        gauss_integrals=gauss,polygon_separations=gaps,
        versions={name:importlib.metadata.version(name) for name in ['spherogram','numpy','scipy','mpmath','numba']},
        conclusion='All 15 nonempty sublinks of the four numerical cycles reduce to unlinks at all tested resolutions and projections.',
        passed=True,scope='Numerical reconstruction and replay of extracted diagrams; no interval ambient-isotopy certificate, no global nonexistence theorem.')
    save(OUT/'review.json',report)
    print('Review passed:',len(allrecords),'Chua diagrams, 34 control polynomials, 8 independent integrations',flush=True)


if __name__=='__main__':main()
