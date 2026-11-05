"""Review the positive link witnesses, independent ODEs and interval records."""
from search_chua_links import ROOT,OUT,read,save,digest
from fractions import Fraction
from pathlib import Path
import json,subprocess
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from spherogram import Link
from chua_link_bracket import frontier_jones
from chua_link_invariants import split_union_jones,fox_nonsplit_test
from chua_link_diagrams import clean,polygon_gauss
from review_chua_linked_attractors import minimum_segment_distance
# The imported historical review sets a different D: temporary folder.
import os
os.environ['TEMP']=os.environ['TMP']=str(ROOT.parent/'.agent-tmp')


def value(text):return float(Fraction(str(text)))


def integration(root,xyz,t,reverse=False):
    alpha,beta,m0,m1=[value(root[k]) for k in ['alpha','beta','m0','m1']]
    if reverse:xyz=-xyz
    def field(time,u):
        x,y,z=u;h=m1*x+.5*(m0-m1)*(abs(x+1)-abs(x-1))
        return [alpha*(y-x-h),x-y+z,-beta*y]
    rows=[]
    for method in ['DOP853','Radau']:
        sol=solve_ivp(field,[0,t[-1]],xyz[0],method=method,rtol=5e-13,atol=5e-14,max_step=.0075,dense_output=True)
        err=float(np.max(np.linalg.norm(sol.sol(t).T-xyz,axis=1)));closure=float(np.linalg.norm(sol.y[:,-1]-xyz[0]));peaks=len(find_peaks(sol.sol(np.linspace(0,t[-1],18001))[0])[0])
        assert sol.success and err<1e-5 and closure<1e-5
        rows.append(dict(method=method,maximum_sampled_trajectory_error=err,closure_error=closure,x_maxima=peaks))
    return rows


def main():
    mp.mp.dps=180
    protocol=subprocess.check_output(['git','show','36ed4c5:results/chua-link-search/protocol.md'],cwd=ROOT)
    assert protocol.replace(b'\r\n',b'\n')==(OUT/'protocol.md').read_bytes().replace(b'\r\n',b'\n')
    certificates={d:read(OUT/f'certificates-d{d}.json') for d in [140,160]}
    for data in certificates.values():assert data['passed'] and data['shifted_root_control']['rejected']
    cases=[]
    for name in ['S02-S08','S01-daughters','S02-daughters']:
        if name=='S02-S08':
            folder=OUT/'overlaps'/name;rootfiles=[folder/f'{id}-d90-n384.json' for id in ['S02','S08']];polys=sorted(folder.glob('jones-d*.json'));pairs=[]
            for path in rootfiles:
                data=np.load(path.with_suffix('.npz'));pairs.append((read(path),data['xyz'],data['t'],False))
        else:
            id=name[:3];folder=OUT/'daughters';rootfiles=[folder/f'{id}-a1e-03-d90-n384.json'];data=np.load(rootfiles[0].with_suffix('.npz'));root=read(rootfiles[0]);pairs=[(root,data['xyz'],data['t'],False),(root,data['xyz'],data['t'],True)];polys=sorted(folder.glob(f'{id}-a1e-03-d*-jones-s*.json'))
        comparisons=[read(p) for p in polys];assert len(comparisons)==(6 if name=='S02-S08' else 7)
        actual=comparisons[0]['joint_polynomial'];expected=comparisons[0]['split_union_polynomial']
        assert actual!=expected
        for p,c in zip(polys,comparisons):
            assert c['joint_polynomial']==actual and c['split_union_polynomial']==expected
            assert split_union_jones(component['polynomial'] for component in c['components'])==expected
            joint=read(p.parent/c['joint_record']);k=Link(joint['simplified_pd_code']);k.unlinked_unknot_components=joint['final_free_components'];replayed,_=frontier_jones(k)
            assert replayed==actual and len(k.link_components)+k.unlinked_unknot_components==2
            assert not np.any(joint['oriented_linking_matrix']) and not any(fox_nonsplit_test(k)['evaluations'].values())
        curve=[-p[1] if p[3] else p[1] for p in pairs];gauss=float(polygon_gauss(clean(curve[0]),clean(curve[1])));assert abs(gauss)<1e-7
        distance,where=minimum_segment_distance(clean(curve[0]),clean(curve[1]));assert distance>0
        integrations=[]
        for root,xyz,t,reverse in pairs:
            if 'source_sha256' in root:assert digest(folder/root['source'])==root['source_sha256']
            assert root['spectral_radius']<1 and root['sampling']['maximum_sampled_region_violation']<1e-25
            integrations.append(dict(id=root['id']+('-inverted' if reverse else ''),checks=integration(root,xyz,t,reverse)))
        cert=next(w for w in certificates[160]['witnesses'] if w['name']==name)
        assert cert['distinct_orbits_by_all_section_points']
        for root in [read(p) for p in rootfiles]:
            cr=next(r for r in cert['records'] if r['id']==root['id'])
            for coordinate,box in zip(root['q'][:2],cr['root']):
                rad=mp.mpf(box['radius'].strip('[]').split()[0]);assert abs(mp.mpf(coordinate)-mp.mpf(box['midpoint']))<rad
            assert cr['primitive_period']
        for d in [140,160]:
            cc=next(w for w in certificates[d]['witnesses'] if w['name']==name)
            assert cc['parameters']==cert['parameters']
            for first,second in zip(cc['records'],cert['records']):
                assert first['upward_crossings']==second['upward_crossings']
                for a,b in zip(first['root'],second['root']):assert abs(mp.mpf(a['midpoint'])-mp.mpf(b['midpoint']))<mp.mpf('1e-60')
        case=dict(name=name,parameters=cert['parameters'],root_files=[p.relative_to(ROOT).as_posix() for p in rootfiles],
            root_sha256={p.name:digest(p) for p in rootfiles},joint_polynomial=actual,split_union_polynomial=expected,
            jones_checks=len(polys),gauss_integral=gauss,minimum_polygon_segment_distance=float(distance),
            integrations=integrations,certified_distinct_attractors=True,primitive_return_counts=[r['upward_crossings'] for r in cert['records']])
        cases.append(case);print(name,'passed','gap',distance,'gauss',gauss,flush=True)
    assert read(OUT/'overlaps/S02-S08/exhaustive-jones.json')['equals_frontier']
    controls=read(OUT/'invariant-controls.json');assert controls['passed']
    assert read(OUT/'bracket-controls.json')['passed']
    negatives=[read(p) for p in (OUT/'symmetry').glob('*.json')];assert len(negatives)==9 and all(r['status']=='split_union' for r in negatives)
    save(OUT/'review.json',dict(passed=True,witnesses=cases,protocol_commit='36ed4c5',protocol_unchanged=True,
        negative_symmetry_pairs=9,interval_precision_repeats=[140,160],jones_diagram_replays=sum(c['jones_checks'] for c in cases),independent_ODE_integrations=12,
        scope='Certified attracting periodic cycles and primitive periods; numerical non-split link types supported by exact diagram polynomials, density/precision/projection repeats. No certified ambient isotopy or Brunnian triple.'))


if __name__=='__main__':main()
