"""Independent full-field integrations and precision checks for the new cycles."""
import argparse
import json
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

from chua_precision_flights import PrecisionFlights
from chua_layer_bifurcations import Flights
from trace_chua_window import Chua
from predict_chua_layers import OUT


def full_field(rec,method,step,variational=False):
    representative=rec['representative'];beta=float(representative['beta']);q=np.array(representative['q'],float)
    u=np.r_[0.,q[:2]];period=float(representative['period'])
    model=Chua(rtol=2.5e-14,step=step,method=method)
    if not variational:
        sol=model.integrate(beta,u,period)
        ts=np.linspace(0,period/2,1001)
        crossings=sol.y_events[2][sol.t_events[2]<period-1e-5]
        distances=np.linalg.norm(crossings[:,None,:]-crossings[None,:,:],axis=2)
        distances+=np.eye(len(crossings))*1e10
        return dict(method=method,max_step=step,rtol=2.5e-14,success=sol.success,
            closure=float(max(abs(sol.y[:,-1]-u))),
            peaks=int(len(sol.t_events[3])),section_returns=int(len(crossings)),
            minimum_distinct_section_distance=float(distances.min()),
            half_inversion_error=float(np.max(abs(sol.sol(ts)+sol.sol(ts+period/2)))))
    def rhs(t,v):
        return np.r_[model.field(beta,v[:3]),(model.jacobian(beta,abs(v[0])<1)@v[3:].reshape(3,3)).ravel()]
    sol=solve_ivp(rhs,[0,period],np.r_[u,np.eye(3).ravel()],method=method,
        rtol=2.5e-14,atol=2.5e-15,max_step=step)
    end=sol.y[:3,-1];M=sol.y[3:,-1].reshape(3,3);v=model.field(beta,end)
    P=((np.eye(3)-np.outer(v,[1.,0.,0.])/v[0])@M)[1:,1:]
    eig=np.linalg.eigvals(P);large=eig[np.argmax(abs(eig))]
    return dict(method=method,max_step=step,rtol=2.5e-14,variational=True,
        success=sol.success,closure=float(max(abs(end-u))),
        dominant_multiplier=float(large.real),expected=.25,multiplier_error=float(abs(large-.25)))


def local_return_variation(rec,step):
    """Integrate the actual field between successive section crossings.

    Project away the neutral flow direction after every crossing before
    composing 2-D derivatives. This avoids carrying a poorly conditioned
    full-period 3-D monodromy product. Reference points are inputs, while
    each next crossing and its derivative are independently integrated.
    """
    from analyze_chua_prediction_evidence import section_word
    representative=rec['representative'];beta=float(representative['beta'])
    motion=section_word(representative,rec['regions'],rec['targets'],high=True)
    points=np.array(motion['normalized_section_points']);model=Chua()
    mp.mp.dps=60;product=mp.eye(2);determinant=mp.mpf(1);checks=[];duration=0.
    for i,point in enumerate(points):
        initial=np.r_[0.,point]
        def rhs(t,v):
            return np.r_[model.field(beta,v[:3]),(model.jacobian(beta,abs(v[0])<1)@v[3:].reshape(3,3)).ravel()]
        def section(t,v):return v[0]
        section.direction=-1;section.terminal=True
        sol=solve_ivp(rhs,[0,20.],np.r_[initial,np.eye(3).ravel()],method='DOP853',
            rtol=2.5e-14,atol=2.5e-15,max_step=step,events=section)
        assert sol.success and len(sol.t_events[0])==1
        end=sol.y[:3,-1];M=sol.y[3:,-1].reshape(3,3);v=model.field(beta,end)
        H=-((np.eye(3)-np.outer(v,[1.,0.,0.])/v[0])@M)[1:,1:]
        local=mp.matrix([[mp.mpf(str(x)) for x in row] for row in H])
        product=local*product;determinant*=mp.det(local)
        closure=float(max(abs(-end[1:]-points[(i+1)%len(points)])))
        duration+=sol.t[-1];checks.append(dict(index=i,duration=sol.t[-1],next_point_error=closure,half_return_matrix=H.tolist()))
    # Multiply local determinants: subtracting entries of the ill-conditioned
    # long product would erase its extremely small contracting multiplier.
    trace=product[0,0]+product[1,1]
    eigenvalues=[(trace+mp.sqrt(trace*trace-4*determinant))/2,(trace-mp.sqrt(trace*trace-4*determinant))/2]
    half=max(eigenvalues,key=abs);dominant=float(mp.re(half*half))
    return dict(method='DOP853 section-by-section full-field variational integration',
        max_step=step,rtol=2.5e-14,variational=True,success=True,
        closure=max(r['next_point_error'] for r in checks),closure_meaning='Maximum independent next-section error, not full-period drift',
        dominant_multiplier=dominant,expected=.25,multiplier_error=abs(dominant-.25),
        signed_half_multiplier=float(mp.re(half)),half_determinant=mp.nstr(determinant,15),
        period_difference=abs(2*duration-float(representative['period'])),
        checks=checks)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUT)
    p.add_argument('--refine-failed',action='store_true')
    p.add_argument('--section-controls',action='store_true')
    args=p.parse_args();out=args.output;search=json.loads((out/'search-results.json').read_text())
    if args.refine_failed or args.section_controls:
        path=out/'independent-checks.json';report=json.loads(path.read_text())
        for row in report['records']:
            if row['passed']:continue
            rec=next(r for r in search['records'] if r['peaks_requested']==row['peaks'])
            for step in [.0025,.00125]:
                check=local_return_variation(rec,step) if args.section_controls else full_field(rec,'DOP853',step,True)
                row['integrations'].append(check)
                row['gates']['independent_attraction']=check['closure']<2e-7 and check['multiplier_error']<2e-5 and abs(check['dominant_multiplier'])<1
                row['passed']=all(row['gates'].values())
                path.write_text(json.dumps(report,indent=2)+'\n')
                print(row['peaks'],'refinement',{k:v for k,v in check.items() if k!='checks'},'passed',row['passed'],flush=True)
                if row['passed'] and not args.section_controls:break
        if not all(r['passed'] for r in report['records']):raise SystemExit(2)
        return
    rows=[];plot={};started=time.monotonic()
    def save():
        (out/'independent-checks.json').write_text(json.dumps(dict(records=rows,elapsed_seconds=time.monotonic()-started),indent=2)+'\n')
    for rec in search['records']:
        if rec['status']!='high_precision_candidate':continue
        row=dict(peaks=rec['peaks_requested'],integrations=[]);rows.append(row);save()
        high=PrecisionFlights(rec['regions'],rec['targets'],True,75)
        controls={}
        for key,target in [('fold',1),('symmetry',-1),('representative',.5)]:
            old=rec[key];new=high.boundary(old['beta'],old['q'],target)
            controls[key]=dict(root=new,beta_change=mp.nstr(abs(mp.mpf(new['beta'])-mp.mpf(old['beta'])),12),
                period_change=mp.nstr(abs(mp.mpf(new['period'])-mp.mpf(old['period'])),12))
            print(row['peaks'],key,'75-digit comparison',controls[key]['beta_change'],flush=True)
        row['precision_checks']=controls
        f=Flights(Chua(),rec['regions'],rec['targets'],True);c=rec['representative']
        rr=f.record(float(c['beta']),np.array(c['q'],float),True)
        row['independent_matrix_exponential']=dict(closure=rr['residual'],
            full_multiplier=rr['full_multipliers'][0][0],
            minimum_breakpoint_speed=rr['minimum_breakpoint_speed'],sampled_region_violation=rr['region_violation'])
        save()
        for method,step,var in [('DOP853',.005,False),('Radau',.005,False),('DOP853',.005,True)]:
            check=full_field(rec,method,step,var);row['integrations'].append(check);save();print(row['peaks'],check,flush=True)
        state=[r for r in row['integrations'] if not r.get('variational')]
        var=row['integrations'][-1]
        row['gates']=dict(precision_stable=max(float(v['beta_change']) for v in controls.values())<1e-28,
            independent_affine=rr['residual']<2e-8 and abs(rr['full_multipliers'][0][0]-.25)<2e-5,
            transverse_itinerary=rr['minimum_breakpoint_speed']>.02 and rr['region_violation']<1e-7,
            full_field_closure=all(r['closure']<2e-7 for r in state),
            independent_counts=all(r['peaks']==rec['peaks_requested'] and r['section_returns']==(rec['peaks_requested']//5) for r in state),
            primitive_orbit=all(r['minimum_distinct_section_distance']>1e-8 for r in state),
            full_field_symmetry=all(r['half_inversion_error']<2e-7 for r in state),
            independent_attraction=var['closure']<2e-7 and var['multiplier_error']<2e-5 and abs(var['dominant_multiplier'])<1)
        row['passed']=all(row['gates'].values());save()
        # Sample affine arcs for illustration; the full-field solvers above
        # count the extrema independently of this display sampling.
        times=[];xyz=[];t0=0.
        # Use the independent matrix-exponential implementation for display.
        for arc in rr['arcs']:
            from scipy.linalg import expm
            reg=arc['region'];A=Chua().jacobian(float(c['beta']),reg==0)
            B=np.zeros((4,4));B[:3,:3]=A;B[0,3]=-9*(-8/7+5/7)*reg
            for dt in np.linspace(0,arc['duration'],max(8,int(arc['duration']/.012)),endpoint=False):
                E=expm(B*dt);xyz.append(E[:3,:3]@arc['start']+E[:3,3]);times.append(t0+dt)
            t0+=arc['duration']
        half=np.array(xyz);tt=np.array(times)
        plot[str(row['peaks'])+'_xyz']=np.vstack([half,-half,half[:1]])
        plot[str(row['peaks'])+'_t']=np.r_[tt,tt+t0,2*t0]
        print(row['peaks'],'passed',row['passed'],row['gates'],flush=True)
    np.savez_compressed(out/'predicted-orbits.npz',**plot)
    if not all(r['passed'] for r in rows):raise SystemExit(2)


if __name__=='__main__':main()
