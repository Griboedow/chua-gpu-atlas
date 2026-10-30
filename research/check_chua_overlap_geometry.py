"""Physical cycles, projection cuts, and independent controls for the overlap atlas."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from resolve_chua_overlap_geometry import (ROOT,SOURCE,PARAMETERS,CONTACTS,EXTENDED,ROLES,
    Contact,edges_at,full_template,checked_record)
from continue_chua_slope_families import CoordinateFlights,dump
from check_chua_slope_families import verify
from trace_chua_window import Chua


def primary(t, row, m0, sign=1):
    s=row['symmetry'];f0=row['fold'];symmetric=m0>s['m0']
    template=t if symmetric else full_template(t)
    f=CoordinateFlights(s,template['regions'],template['targets'],symmetric=symmetric)
    if symmetric:
        q=f.solve(m0,np.array(s['q']))
    else:
        # The corrected flip edge already selects the asymmetric branch.
        # Continue inward from it rather than risk converging to the parent.
        q=np.array(row['doubling']['q']);start=row['doubling']['m0']
        for value in np.linspace(start,m0,9)[1:]:
            q=f.solve(value,q)
        if sign==-1:
            r=f.record(m0,q,True);k=len(t['regions'])-1
            arc=r['arcs'][k];A,c,eig,V,W=f.matrices(r['beta'])[0]
            def state(time):return c+np.real(V@(np.exp(eig*time)*(W@(np.array(arc['start'])-c))))
            tau=brentq(lambda time:state(time)[0],0,arc['duration'],xtol=1e-14)
            dt=q[2:]
            durations=np.r_[dt[k]-tau,dt[k+1:-1],dt[-1]+dt[0],dt[1:k],tau]
            q=f.solve(m0,np.r_[-state(tau)[1:],durations])
    r=checked_record(f,m0,q,t['reference_peak_count'])
    if r['spectral_radius']>=1:
        raise RuntimeError(f'Requested primary witness is not attracting: {t["id"]}, alpha={s["alpha"]}, m1={s["m1"]}, m0={m0}, symmetric={symmetric}, rho={r["spectral_radius"]}')
    r.update(id=t['id'],branch='symmetric' if symmetric else ('asymmetric A' if sign==1 else 'asymmetric B'))
    return r


def slice_roots(out,kind,control,code):
    rows=json.loads((out/(kind+'-'+code+'.json')).read_text())['records']
    roots=[]
    for a,b in zip(rows,rows[1:]):
        if (a['control']-control)*(b['control']-control)>=0:continue
        tangent=Contact(kind,a['m1'],code).tangent(a)
        def guess(m1):return np.array(a['w'])+(m1-a['m1'])*tangent
        def value(m1):return Contact(kind,m1,code).solve(guess(m1))['control']-control
        m1=brentq(value,a['m1'],b['m1'],xtol=2e-13)
        roots.append(Contact(kind,m1,code).solve(guess(m1)))
    return roots


def slices(out):
    templates=json.loads((SOURCE/'anchors.json').read_text())
    for alpha in [8.84,8.85,8.86,8.865]:
        path=out/f'slice-{alpha:g}.json'
        if path.exists():print('exists',path.name,flush=True);continue
        records=[]
        for slope in np.linspace(-.785,-.718,101):
            edges=edges_at('alpha',alpha,float(slope),True)
            checks=[]
            for t in templates:
                row=edges[t['id']]
                for lo,hi in [('symmetry','fold'),('doubling','symmetry')]:
                    m0=(row[lo]['m0']+row[hi]['m0'])/2
                    cycle=primary(t,row,m0)
                    checks.append(dict(id=t['id'],branch=cycle['branch'],m0=m0,
                        residual=cycle['residual'],spectral_radius=cycle['spectral_radius'],
                        period=cycle['period'],minimum_breakpoint_speed=cycle['minimum_breakpoint_speed']))
            records.append(dict(m1=float(slope),edges=edges,interior_checks=checks))
        roots={code:slice_roots(out,'alpha',alpha,code) for code in ['FS','SF','FD','DF','SD','DS']}
        dump(path,dict(alpha=alpha,beta=100/7,records=records,roots=roots,
            scope='Six corrected edges and four attracting interior cycles per sampled m1. Counts cover the two primary families through their first period doubling only.'))
        print('slice',alpha,'nodes',len(records),'contact roots',{k:len(v) for k,v in roots.items()},flush=True)


def witnesses(out):
    templates=json.loads((SOURCE/'anchors.json').read_text())
    records=[];traces={}
    for count,alpha,low_role,high_role in [(2,8.85,'symmetry','fold'),
                                          (3,8.86,'doubling','fold'),
                                          (4,8.854,'doubling','symmetry')]:
        m1=-.7521
        edges=edges_at('alpha',alpha,m1,True)
        lo=max(e[low_role]['m0'] for e in edges.values())
        hi=min(e[high_role]['m0'] for e in edges.values())
        assert lo<hi
        m0=(lo+hi)/2;cycles=[]
        for t in templates:
            row=edges[t['id']]
            signs=[1] if m0>row['symmetry']['m0'] else [1,-1]
            partners=[]
            for sign in signs:
                r=primary(t,row,m0,sign);r['independent']=verify(r)
                model=Chua(alpha=alpha,m0=m0,m1=m1,rtol=2e-12,step=.01)
                sol=model.integrate(100/7,[0.,*r['q'][:2]],r['period'])
                times=np.linspace(0,r['period'],2401)
                states=sol.sol(times)
                r['half_period_inversion_error']=float(np.max(abs(sol.sol(r['period']/2)+[0.,*r['q'][:2]])))
                if not r['symmetric']:
                    assert r['half_period_inversion_error']>1e-3
                key=f'n{count}_{t["id"]}_{sign}'
                r['trace_key']=key
                traces[key+'_time']=times;traces[key+'_state']=states
                partners.append(r);cycles.append(r)
            if len(partners)==2:
                a,b=partners
                def down(t,state):return state[0]
                down.direction=-1
                model=Chua(alpha=alpha,m0=m0,m1=m1)
                sol=solve_ivp(lambda t,u:model.field(100/7,u),[0,a['period']],[0.,*a['q'][:2]],
                    events=down,rtol=2e-12,atol=2e-13,max_step=.01)
                error=float(min(np.max(abs(-u-[0.,*b['q'][:2]])) for u in sol.y_events[0]))
                assert error<2e-7
                for r in partners:r['inversion_partner_error']=error
        assert len(cycles)==count
        records.append(dict(count=count,parameters=dict(alpha=alpha,beta=100/7,m0=m0,m1=m1),
            interval_m0=[lo,hi],cycles=cycles,edges=edges))
        print('witness',count,'alpha',alpha,'m0',m0,'radii',[r['spectral_radius'] for r in cycles],flush=True)
    dump(out/'witnesses.json',dict(records=records,
        scope='At least 2, 3 or 4 attracting cycles at the respective parameter points. Other families and doubled descendants have not been excluded.'))
    np.savez_compressed(out/'witness-traces.npz',**traces)


def controls(out):
    checks=[];transversality=[]
    for kind in ['alpha','beta']:
        for code in CONTACTS+EXTENDED:
            rows=json.loads((out/(kind+'-'+code+'.json')).read_text())['records']
            for index in [0,len(rows)//4,len(rows)//2,3*len(rows)//4,len(rows)-1]:
                r=rows[index];problem=Contact(kind,r['m1'],code)
                slopes=[]
                for f,cycle,role in zip(problem.functions(np.array(r['w'])),r['cycles'],code):
                    checks.append(dict(kind=kind,code=code,m1=r['m1'],id=cycle['id'],**verify(cycle)))
                    w=np.r_[cycle['q'],cycle['m0']]
                    h=2e-6;values=[]
                    for sign in [-1,1]:
                        params={k:cycle[k] for k in PARAMETERS};params[kind]+=sign*h
                        shifted=CoordinateFlights(params,f.regions,f.targets,symmetric=f.symmetric)
                        values.append(shifted.boundary_system(w,ROLES[role][1]))
                    derivative=-np.linalg.solve(f.boundary_system(w,ROLES[role][1],True),(values[1]-values[0])/(2*h))
                    slopes.append(float(derivative[-1]))
                transversality.append(dict(kind=kind,code=code,m1=r['m1'],edge_slopes=slopes,
                    difference=slopes[0]-slopes[1],smallest_joint_jacobian_singular_value=float(np.linalg.svd(problem.system(np.array(r['w']),True),compute_uv=False)[-1])))
            print('independent contacts',kind,code,flush=True)
    dump(out/'independent-checks.json',dict(checks=checks,transversality=transversality,
        scope='Selected full-field state-and-variational integrations and local surface-transversality controls.'))


def branches(out):
    """Resolve stable segments at selected controls instead of only their midpoint."""
    templates=json.loads((SOURCE/'anchors.json').read_text());parameters=[]
    for kind in ['alpha','beta']:
        curves=[json.loads((out/(kind+'-'+code+'.json')).read_text())['records'] for code in ['FS','SF']]
        for requested in [-.782,-.754,-.69]:
            rs=[min(rows,key=lambda r:abs(r['m1']-requested)) for rows in curves]
            parameters.append((kind,float(np.mean([r['control'] for r in rs])),rs[0]['m1']))
    parameters += [('alpha',alpha,-.7521) for alpha in [8.85,8.86,8.854]]
    records=[]
    for kind,control,m1 in parameters:
        edges=edges_at(kind,control,m1,True)
        for t in templates:
            row=edges[t['id']];s=row['symmetry']
            for symmetric in [True,False]:
                template=t if symmetric else full_template(t)
                f=CoordinateFlights(s,template['regions'],template['targets'],symmetric=symmetric)
                low=s['m0'] if symmetric else row['doubling']['m0']
                high=row['fold']['m0'] if symmetric else s['m0']
                q=np.array(s['q'] if symmetric else row['doubling']['q'])
                samples=[]
                for fraction in np.linspace(.01,.99,65):
                    m0=low+fraction*(high-low);q=f.solve(m0,q)
                    r=checked_record(f,m0,q,t['reference_peak_count'])
                    assert r['spectral_radius']<1
                    item=dict(fraction=float(fraction),m0=m0,period=r['period'],
                              dominant=r['map_multipliers'][0][0],spectral_radius=r['spectral_radius'])
                    if len(samples)==32:item['independent']=verify(r)
                    samples.append(item)
                eigenvalues=np.array([r['dominant'] for r in samples])
                assert np.min(np.diff(eigenvalues))>-1e-7
                records.append(dict(kind=kind,control=control,m1=m1,id=t['id'],symmetric=symmetric,
                    interval_m0=[low,high],samples=samples))
        print('branch segments',kind,control,m1,flush=True)
    dump(out/'branch-controls.json',dict(records=records,
        scope='65 interior points on each of 36 selected stable branch segments, with signed multiplier monotonicity checks and independent midpoint controls. Sampling is not a proof of stability at every intervening parameter.'))


def audit(out):
    contacts=[]
    for kind in ['alpha','beta']:
        for code in CONTACTS+EXTENDED:
            d=json.loads((out/(kind+'-'+code+'.json')).read_text())
            assert not d['failures']
            for r in d['records']:
                assert r['residual']<2e-9
                assert abs(r['cycles'][0]['period']-r['cycles'][1]['period'])>5
                for c,role in zip(r['cycles'],code):
                    assert c['full_cycle_peaks']==int(c['id'][1:])
                    assert abs(c['map_multipliers'][0][0]-ROLES[role][1])<1e-7
                    assert c['extrema_region_violation']<1e-7 and c['minimum_breakpoint_speed']>.02
                contacts.append(r)
    checks=json.loads((out/'independent-checks.json').read_text())
    assert all(x['closure']<2e-7 and x['multiplier_error']<2e-5 for x in checks['checks'])
    assert min(abs(x['difference']) for x in checks['transversality'])>1e-3
    assert min(x['smallest_joint_jacobian_singular_value'] for x in checks['transversality'])>1e-5
    tips=sum([json.loads((out/f).read_text())['records'] for f in ['turning-points.json','turning-extended.json']],[])
    for tip in tips:
        assert abs(tip['tangent'][-1])<1e-6
        curv=[c['second_derivative'] for c in tip['curvatures']]
        assert np.ptp(curv)/abs(curv[-1])<1e-3
    predictions=json.loads((out/'predictions.json').read_text())['records']
    assert max(r['relative_separation_error'] for r in predictions)<1e-3
    total=0;roots={};components={}
    for alpha in [8.84,8.85,8.86,8.865]:
        d=json.loads((out/f'slice-{alpha:g}.json').read_text());total+=len(d['records'])
        for row in d['records']:
            for e in row['edges'].values():assert e['doubling']['m0']<e['symmetry']['m0']<e['fold']['m0']
            assert all(c['spectral_radius']<1 for c in row['interior_checks'])
        roots[str(alpha)]={code:len(values) for code,values in d['roots'].items()}
        counts=[]
        for low in ['symmetry','doubling']:
            mask=np.array([min(e['fold']['m0'] for e in r['edges'].values())>max(e[low]['m0'] for e in r['edges'].values()) for r in d['records']])
            counts.append(int(np.sum(np.diff(np.r_[False,mask].astype(int))==1)))
        components[str(alpha)]=counts
    assert list(components.values())==[[0,0],[1,1],[2,1],[2,2]]
    branches=json.loads((out/'branch-controls.json').read_text())['records']
    assert len(branches)==36 and all(len(r['samples'])==65 for r in branches)
    assert all(s['spectral_radius']<1 for r in branches for s in r['samples'])
    w=json.loads((out/'witnesses.json').read_text())['records']
    for row in w:
        assert len(row['cycles'])==row['count']
        assert all(c['spectral_radius']<1 and c['independent']['closure']<2e-7 for c in row['cycles'])
    report=dict(status='passed',contact_points=len(contacts),contact_curves=18,
        independent_contact_cycle_checks=len(checks['checks']),
        maximum_independent_closure=max(r['closure'] for r in checks['checks']),
        maximum_independent_multiplier_error=max(r['multiplier_error'] for r in checks['checks']),
        minimum_contact_transversality=min(abs(r['difference']) for r in checks['transversality']),
        minimum_joint_jacobian_singular_value=min(r['smallest_joint_jacobian_singular_value'] for r in checks['transversality']),
        projection_turning_points=len(tips),held_out_cut_predictions=len(predictions),
        maximum_prediction_relative_error=max(r['relative_separation_error'] for r in predictions),
        slice_cross_sections=total,stable_slice_interior_checks=total*4,slice_root_counts=roots,
        slice_symmetric_and_primary_component_counts=components,
        refined_stable_branch_samples=sum(len(r['samples']) for r in branches),
        independent_stable_branch_midpoints=len(branches),
        coexistence_witness_counts=[r['count'] for r in w],
        maximum_witness_closure=max(c['independent']['closure'] for r in w for c in r['cycles']),
        scope='Numerical validation of the two selected primary families. No completeness of the full attractor set, interval proof, or world-first claim.')
    dump(out/'audit.json',report);print(json.dumps(report,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/chua-overlap-geometry')
    p.add_argument('--stage',choices=['slices','witnesses','controls','branches','audit'],required=True)
    args=p.parse_args();globals()[args.stage](args.output)


if __name__=='__main__':main()
