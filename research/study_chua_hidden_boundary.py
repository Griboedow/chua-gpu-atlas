"""Follow a saddle periodic orbit, its local basin separator, and linear contacts.

The saddle curve is a numerical finite-iterate stable-manifold approximation.
Separate interval certificates establish the two periodic roots, not a global
classification of the basin boundary or the absence of other attractors.
"""
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import root, brentq
from scipy.integrate import solve_ivp
from continue_chua_slope_families import CoordinateFlights, extrema_check, dump
from chua_layer_bifurcations import from_orbit
from chua_event_map import EventMap
from trace_chua_window import Chua
from check_chua_slope_refinement import independent_segmented

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/chua-hidden-boundary'
REGIONS, TARGETS = [0,1,0], [1,1]


def record(m0, m1, q, alpha=8.4, beta=12):
    f = CoordinateFlights(dict(alpha=alpha,beta=beta,m0=float(m0),m1=float(m1)),REGIONS,TARGETS)
    q = f.solve(m0,np.array(q)); r = f.record(m0,q,True)
    r.update(extrema_check(f,r))
    assert r['residual']<2e-9 and r['extrema_region_violation']<2e-8
    arc = r['arcs'][1]; A,c,e,V,W = f.matrices(beta)[1]
    coeff = V[0]*(W@(np.array(arc['start'])-c))
    dx = lambda t: float(np.real(np.sum(e*coeff*np.exp(e*t))))
    ts = np.linspace(0,arc['duration'],101)
    turns = [brentq(dx,l,h,xtol=1e-14) for l,h in zip(ts[:-1],ts[1:]) if dx(l)*dx(h)<0]
    r['amplitude'] = max(float(c[0]+np.real(np.sum(coeff*np.exp(e*t)))) for t in turns)
    r.update(regions=REGIONS,targets=TARGETS)
    return r


def full_record(r):
    q = r['q']; durations = [q[2],q[3],q[4]+q[2],q[3],q[4]]
    f = CoordinateFlights(r,[0,1,0,-1,0],[1,1,-1,-1],symmetric=False)
    qq = f.solve(r['m0'],np.r_[q[:2],durations])
    result = f.record(r['m0'],qq,True)
    extrema_check(f,result)
    return dict(id='Saddle' if r['spectral_radius']>1 else 'HiddenCycle',
                regions=f.regions,targets=f.targets,cycle=result)


def seeds():
    a,b = 8.4,12.
    critical = np.sort(np.roots([a*a,a+a*a,b]))
    omega = np.sqrt(b+a*critical[0])
    base = np.array([omega/a,b*(1+critical[0])/omega])
    inner = record(-.95,-.05,np.r_[1.1*base,.65,.3,.65])
    assert inner['spectral_radius']>1
    source = min(json.loads((ROOT/'results/chua-accessibility/hidden-inner-continuation.json').read_text()),
                 key=lambda r:abs(r['m0']+.95))
    f,q = from_orbit(Chua(alpha=a,m0=-.95,m1=-.05),b,[*source['q'][:2],source['period']],True)
    assert f.regions==REGIONS
    outer = record(-.95,-.05,q)
    dump(OUT/'cycle-seeds.json',dict(saddle=inner,attractor=outer,full_records=[full_record(inner),full_record(outer)]))
    return critical,inner,outer


def scaling(critical,inner,outer):
    forecasts=[]
    for name,mstar in zip(['saddle','attractor'],critical):
        omega=np.sqrt(12+8.4*mstar)
        amp=(3*np.pi/(8*np.sqrt(2)*(-.05-mstar)))**(2/3)
        forecasts.append(dict(branch=name,mstar=float(mstar),omega=float(omega),
            amplitude_excess_coefficient=float(amp),outer_duration_coefficient=float(2*np.sqrt(2*amp)/omega),
            amplitude_exponent=2/3,outer_duration_exponent=1/3))
    prediction=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),forecasts=forecasts,
        derivation='First-harmonic balance: 1-N(1+eta) ~ 8*sqrt(2)/(3*pi)*eta^(3/2); '
                   'mstar-m0 ~ (m1-mstar)*(1-N). This is a local asymptotic prediction, not a rigorous remainder bound.',
        protocol='Analytic coefficients fixed after coarse exploratory continuation and before the fine-distance sweep; not an externally published blind forecast.')
    path=OUT/'linear-contact-prediction.json'
    if not path.exists():dump(path,prediction)
    else:
        old=json.loads(path.read_text());assert old['forecasts']==forecasts
    rows=[]
    for guess,pred in zip([inner,outer],forecasts):
        # Reach the first distance with a conventional m0 continuation.
        target=pred['mstar']-.01
        for m0 in np.linspace(guess['m0'],target,120):guess=record(m0,-.05,guess['q'])
        for delta in np.geomspace(.01,1e-7,31):
            guess=record(pred['mstar']-delta,-.05,guess['q'])
            excess=guess['amplitude']-1
            row=dict(branch=pred['branch'],delta=float(delta),cycle=guess,
                predicted_excess=float(pred['amplitude_excess_coefficient']*delta**(2/3)),
                predicted_outer_duration=float(pred['outer_duration_coefficient']*delta**(1/3)))
            row['amplitude_ratio']=excess/row['predicted_excess']
            row['duration_ratio']=guess['q'][3]/row['predicted_outer_duration']
            if len(rows)%10==0 or delta==1e-7:row['independent']=independent_segmented(guess)
            rows.append(row)
        print(pred['branch'],'last ratios',rows[-1]['amplitude_ratio'],rows[-1]['duration_ratio'],flush=True)
    dump(OUT/'linear-contact-results.json',dict(records=rows,
        scope='Actual piecewise-linear periodic roots. Limits approach amplitude one, not zero; '
              'a generic smooth Hopf amplitude formula cannot be substituted for this contact calculation.'))


def branches(inner,outer):
    result=[]
    for name,seed,stop in [('saddle',inner,-.938),('attractor',outer,-.182)]:
        q=seed['q']
        for m0 in np.linspace(-.95,-.999,60):q=record(m0,-.05,q)['q']
        for m0 in np.linspace(-.999,stop,180):
            r=record(m0,-.05,q);q=r['q'];result.append(dict(branch=name,cycle=r))
    dump(OUT/'branches.json',dict(records=result))


def basin(inner,outer):
    f=CoordinateFlights(inner,REGIONS,TARGETS)
    _,_,H,_=f.evaluate(inner['m0'],np.array(inner['q']))
    values,V=np.linalg.eig(H@H);order=np.argsort(-abs(values));V=V[:,order].real;values=values[order].real
    V/=np.linalg.norm(V,axis=0)
    if V[:,0]@np.array(inner['q'][:2])<0:V[:,0]*=-1
    if V[0,1]<0:V[:,1]*=-1
    W=np.linalg.inv(V);centre=np.array(inner['q'][:2])
    emap=EventMap(Chua(alpha=8.4,m0=-.95,m1=-.05),1,.04)
    curve=[]
    for stable in np.linspace(-.035,.035,71):
        approximations=[]
        for count in [4,6]:
            def residual(unstable):
                point=centre+stable*V[:,1]+unstable[0]*V[:,0];D=np.eye(2)
                for _ in range(count):
                    shift,J,_=emap.evaluate([*point,12]);point+=shift;D=J@D
                return np.array([W[0]@(point-centre)]),np.array([[W[0]@D@V[:,0]]])
            sol=root(lambda v:residual(v)[0],[0.],jac=lambda v:residual(v)[1],tol=1e-10)
            error=abs(residual(sol.x)[0][0]);assert error<1e-8
            point=centre+stable*V[:,1]+sol.x[0]*V[:,0]
            approximations.append(dict(iterations=count,point=point.tolist(),unstable_coordinate=float(sol.x[0]),residual=float(error)))
        curve.append(dict(stable_coordinate=float(stable),approximations=approximations,
            difference=float(np.linalg.norm(np.array(approximations[0]['point'])-approximations[1]['point']))))
    controls=[]
    for index in [0,10,20,30,40,50,60,70]:
        point=np.array(curve[index]['approximations'][-1]['point'])
        for sign in [-1,1]:
            initial=np.r_[0.,point+sign*1e-4*V[:,0]]
            model=Chua(alpha=8.4,m0=-.95,m1=-.05,rtol=2e-11,step=.08)
            sol=model.integrate(12,initial,700,dense=False)
            norm=float(np.linalg.norm(sol.y[:,-1]))
            distance=float(np.linalg.norm(sol.y_events[2][-1,1:]-outer['q'][:2]))
            label='origin' if norm<1e-7 else 'outer_cycle' if distance<1e-7 else 'unresolved'
            expected='origin' if sign<0 else 'outer_cycle'
            controls.append(dict(curve_index=index,side=sign,initial=initial.tolist(),label=label,expected=expected,
                final_norm=norm,last_section_distance=distance,duration=700,method='DOP853'))
            print('separator',index,sign,label,flush=True)
    assert all(c['label']==c['expected'] for c in controls)
    dump(OUT/'basin-separator.json',dict(cycle=inner,section_eigenvalues=values.tolist(),basis=V.tolist(),
        curve=curve,controls=controls,
        scope='Numerical local stable-manifold approximation and sixteen two-sided trajectory outcomes. '
              'This does not prove the complete global basin boundary or exclude other attractors.'))


if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    critical,inner,outer=seeds()
    scaling(critical,inner,outer)
    branches(inner,outer)
    basin(inner,outer)
