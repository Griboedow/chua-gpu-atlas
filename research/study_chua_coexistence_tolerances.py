"""Numerical slope section and exact component image of a certified box."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from flint import arb,ctx
from study_chua_circuit_constraints import LossyChua,LossyFlights,parameters
from chua_interval_map import bounds
from chua_research_followup import ROOT,OUT,dump
from certify_chua_coexistence_box import nominal,enclosure

DIR=OUT/'coexistence-box'


def solve(record,p):
    f=LossyFlights(LossyChua(**{k:v for k,v in p.items() if k!='beta'}),record['regions'],record['targets'],False)
    q=f.solve(p['beta'],np.array(record['cycle']['q']))
    rec=f.record(p['beta'],q,True)
    if rec['region_violation']>1e-7 or min(q[2:])<=0:raise ValueError('Inadmissible itinerary')
    if np.linalg.norm(q-np.array(record['cycle']['q']))>.5:raise ValueError('Continuation too far from the initial branch')
    return rec


def main():
    ctx.dps=80
    rows=json.loads((ROOT/'results/chua-layer-atlas/four-cycles.json').read_text())['records']
    center={k:float(v) for k,v in zip(['alpha','beta','m0','m1','gamma'],nominal())}
    offsets=np.linspace(-.0004,.0004,31)
    status=np.zeros((31,31,4),np.int8);rhos=np.full((31,31,4),np.nan)
    points=np.full((31,31,4,2),np.nan);distinct=np.zeros((31,31),np.int8)
    failures=[]
    for i,dm0 in enumerate(offsets):
        for j,dm1 in enumerate(offsets):
            p=dict(center,m0=center['m0']+dm0,m1=center['m1']+dm1)
            for k,row in enumerate(rows):
                try:
                    rec=solve(row,p);rhos[i,j,k]=rec['spectral_radius']
                    points[i,j,k]=rec['q'][:2]
                    status[i,j,k]=2 if rec['spectral_radius']<1 else 1
                except (ValueError,RuntimeError,np.linalg.LinAlgError) as err:
                    failures.append(dict(index=[i,j,k],message=str(err)))
            representatives=[]
            for k in np.flatnonzero(status[i,j]==2):
                point=points[i,j,k]
                if all(np.linalg.norm(point-old)>1e-6 for old in representatives):representatives.append(point)
            distinct[i,j]=len(representatives)
        print('slope row',i,flush=True)
    np.savez_compressed(DIR/'slope-section.npz',offsets=offsets,status=status,spectral_radii=rhos,
                        section_points=points,distinct_attracting_counts=distinct)
    checks=[]
    for row in rows:
        for sign in [-1,1]:
            X=np.r_[0.,row['cycle']['q'][:2]]+np.array([0,sign*1e-5,-sign*1e-5])
            model=LossyChua(alpha=center['alpha'],m0=center['m0'],m1=center['m1'])
            def event(t,x):return x[0]
            event.direction=1
            sol=solve_ivp(lambda t,x:model.field(center['beta'],x),[0,100],X,events=event,
                method='DOP853',rtol=2e-11,atol=2e-12,max_step=.02)
            distance=float(np.linalg.norm(sol.y_events[0][-1,1:]-np.array(row['cycle']['q'][:2])))
            assert sol.success and distance<2e-7,(row['id'],distance)
            checks.append(dict(id=row['id'],initial_state=X.tolist(),last_section_distance=distance,passed=True))
    dump(DIR/'slope-section.json',dict(parameters=center,grid_shape=list(status.shape),
        status={'0':'unresolved','1':'continued but not attracting','2':'continued attracting'},
        failed_attempts=failures,selection_checks=checks,
        four_attracting_nodes=int(np.sum(distinct==4)),section_clustering_tolerance=1e-6,
        nodes_with_merged_attracting_references=int(np.sum(distinct<np.sum(status==2,axis=2))),
        scope='Finite numerical continuation screen. A failed solve is unresolved, not evidence of cycle absence.'))
    components()


def components():
    ctx.dps=100
    # All component formulas are evaluated as rational/ball arithmetic. This
    # is an ideal realization, not an assertion of attainable hardware tolerances.
    alpha,beta,m0,m1,gamma=nominal();R=arb(1000);K1=arb(40);K2=arb(2);C2=arb(1)/10**7
    delta=m1-m0
    c=dict(R=R,R1=R*K1/delta,R2=R*(K2-1)/(delta/K1-m1),R3=arb(39000),R4=arb(1000),
           R5=arb(1000),R6=arb(1000),C1=C2/alpha,C2=C2,L=C2*R*R/beta,RL=arb(0))
    relative=arb('1e-10');cp={k:v*arb(1,relative) for k,v in c.items() if k!='RL'}
    cp['RL']=arb(0).union(arb('3e-7'))
    image,saturation=parameters(cp)
    cert=json.loads((DIR/'certificate.json').read_text())
    containment={k:enclosure(cert['parameter_box'][k],True).contains(v) for k,v in image.items() if k!='gamma'}
    # Arb balls slightly over-enclose the zero endpoint of RL. Physical RL is
    # restricted to [0,3e-7], so gamma>=0 follows from the component signs.
    assert all(cp[k]>0 for k in ['C2','R','L'])
    containment['gamma']=bool(enclosure(cert['parameter_box']['gamma'],True).contains(0)
        and image['gamma'].upper()<arb(cert['parameter_box']['gamma']['upper']).lower())
    assert all(containment.values()),containment
    assert all(arb(row['x_abs_upper']['ball'])<saturation for row in cert['records'])
    dump(DIR/'component-image.json',dict(nominal_components_SI={k:bounds(v) for k,v in c.items()},
        relative_tolerance=str(relative),inductor_resistance_ohm='[0, 3e-7]',
        component_box={k:bounds(v) for k,v in cp.items()},parameter_image={k:bounds(v) for k,v in image.items()},
        component_constraints=['All components other than RL are strictly positive.','RL is restricted to [0,3e-7] ohm.'],
        gamma_physical_lower='0',
        contained=containment,outer_breakpoint=bounds(saturation),passed=True,
        scope='Conservative interval image inside the mathematical certificate. These extremely small tolerances are not a practical hardware specification. The five-segment ideal diode remains unsaturated; op-amp bandwidth and noise are not modeled.',
        source='https://doi.org/10.1007/s11071-022-08078-y'))


if __name__=='__main__':main()
