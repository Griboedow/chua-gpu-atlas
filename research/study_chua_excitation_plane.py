"""Initial-state basin section, stored-energy thresholds and rectangular forcing.

Numerical classification is deliberately separate from the earlier interval
proof on the positive y ray. A finite grid cannot exclude unresolved islands.
"""
import argparse
import json
import time
import numpy as np
from numba import njit, prange, set_num_threads
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from study_chua_trigger import setup, PARS
from chua_research_followup import ROOT, OUT, dump

DIR=OUT/'excitation-plane'


@njit(cache=True)
def field(x,y,z,u):
    h=-.05*x-.9*min(1.,max(-1.,x))
    return 8.4*(y-x-h),x-y+z+u,-12*y


@njit(cache=True)
def step(x,y,z,h,u):
    a,b,c=field(x,y,z,u)
    d,e,f=field(x+h*a/2,y+h*b/2,z+h*c/2,u)
    g,i,j=field(x+h*d/2,y+h*e/2,z+h*f/2,u)
    k,l,m=field(x+h*g,y+h*i,z+h*j,u)
    return x+h*(a+2*d+2*g+k)/6,y+h*(b+2*e+2*i+l)/6,z+h*(c+2*f+2*j+m)/6


@njit(cache=True)
def classify(x,y,z,h,pc,P,level,max_time=600.):
    hits=0
    for k in range(int(max_time/h)):
        if k%20==0:
            V=x*(P[0,0]*x+P[0,1]*y+P[0,2]*z)+y*(P[1,0]*x+P[1,1]*y+P[1,2]*z)+z*(P[2,0]*x+P[2,1]*y+P[2,2]*z)
            if V < .8*level: return 1,k*h
        xx,yy,zz=step(x,y,z,h,0.)
        if x<0 and xx>=0 and yy>0:
            s=-x/(xx-x)
            sy=y+s*(yy-y);sz=z+s*(zz-z)
            if (sy-pc[0])**2+(sz-pc[1])**2<4e-6:
                hits+=1
                if hits>=2:return 2,k*h
            else:hits=0
        x,y,z=xx,yy,zz
    return 0,max_time


@njit(cache=True,parallel=True)
def batch(points,h,pc,P,level):
    labels=np.zeros(len(points),np.int8);times=np.zeros(len(points))
    for i in prange(len(points)):
        labels[i],times[i]=classify(points[i,0],points[i,1],points[i,2],h,pc,P,level)
    return labels,times


@njit(cache=True)
def threshold(theta,lo,hi,h,pc,P,level):
    for _ in range(34):
        mid=(lo+hi)/2
        label,_=classify(0.,mid*np.cos(theta),np.sqrt(12.)*mid*np.sin(theta),h,pc,P,level)
        if label==0:return np.nan
        if label==1:lo=mid
        else:hi=mid
    return (lo+hi)/2


@njit(cache=True,parallel=True)
def ray_thresholds(angles,radii,labels,h,pc,P,level):
    roots=np.full(len(angles),np.nan)
    for i in prange(len(angles)):
        for j in range(len(radii)-1):
            if labels[i,j]==1 and labels[i,j+1]==2:
                roots[i]=threshold(angles[i],radii[j],radii[j+1],h,pc,P,level)
                break
    return roots


@njit(cache=True)
def pulse_end(amplitude,duration,h):
    x=y=z=0.
    n=max(1,int(np.ceil(duration/h)));hh=duration/n
    for _ in range(n): x,y,z=step(x,y,z,hh,amplitude)
    return x,y,z


@njit(cache=True,parallel=True)
def pulse_batch(durations,impulses,h,pc,P,level):
    labels=np.zeros((len(durations),len(impulses)),np.int8)
    for i in prange(len(durations)):
        for j in range(len(impulses)):
            point=pulse_end(impulses[j]/durations[i],durations[i],h)
            labels[i,j]=classify(*point,h,pc,P,level)[0]
    return labels


@njit(cache=True)
def pulse_threshold(tau,pc,P,level):
    lo=.2;hi=.9
    for _ in range(29):
        mid=(lo+hi)/2
        label=classify(*pulse_end(mid/tau,tau,.0005),.001,pc,P,level)[0]
        if label==0:return np.nan
        if label==1:lo=mid
        else:hi=mid
    return (lo+hi)/2


def independent(point,pc,P,level,forcing=None):
    def rhs(t,X,u=0.):return field.py_func(*X,u)
    rows=[]
    for method in ['DOP853','Radau']:
        start=np.array(point,float)
        if forcing:
            a,tau=forcing
            sol=solve_ivp(lambda t,X:rhs(t,X,a),[0,tau],start,method=method,
                          rtol=2e-11,atol=2e-12,max_step=.02)
            assert sol.success
            start=sol.y[:,-1]
        def crossing(t,X):return X[0]
        crossing.direction=1
        sol=solve_ivp(rhs,[0,600],start,events=crossing,method=method,
                      rtol=2e-11,atol=2e-12,max_step=.03)
        assert sol.success
        X=sol.y[:,-1]
        label=1 if X@P@X<level else 0
        crossings=sol.y_events[0]
        distance=float(np.linalg.norm(crossings[-1,1:]-pc)) if len(crossings) else None
        if distance is not None and distance<1e-6:label=2
        rows.append(dict(method=method,label=label,last_section_distance=distance))
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['plane','pulses','check']);args=ap.parse_args()
    set_num_threads(8)
    p,pc,P,level,*_=setup();start=time.monotonic()
    DIR.mkdir(parents=True,exist_ok=True)
    if args.phase=='plane':
        angles=np.linspace(-np.pi/2,np.pi/2,181)
        radii=np.linspace(0,1.5,241)
        aa,rr=np.meshgrid(angles,radii,indexing='ij')
        points=np.c_[np.zeros(aa.size),(rr*np.cos(aa)).ravel(),(np.sqrt(12)*rr*np.sin(aa)).ravel()]
        labels,times=batch(points,.004,pc,P,level)
        labels2,times2=batch(points,.002,pc,P,level)
        disagreements=labels!=labels2
        safe=labels2.copy();safe[disagreements]=0
        safe=safe.reshape(aa.shape)
        roots=ray_thresholds(angles,radii,safe,.002,pc,P,level)
        minima=[]
        # The sampled angle interval is periodic modulo pi by odd symmetry.
        # Its endpoints are neighbors, so a pure-current minimum must not be
        # discarded merely because it sits at the edge of the plot.
        if roots[0]<roots[1] and roots[0]<roots[-2]:
            fun=lambda a:threshold(a,max(0.,roots[0]-.04),roots[0]+.04,.001,pc,P,level)
            res=minimize_scalar(fun,bounds=(angles[0]-(angles[1]-angles[0]),angles[1]),
                                method='bounded',options={'xatol':2e-7})
            r=float(res.fun);a=float(res.x)
            minima.append(dict(angle=a,radius=r,normalized_energy=r*r/2,
                initial_state=[0.,r*np.cos(a),np.sqrt(12)*r*np.sin(a)],success=bool(res.success)))
        for i in range(1,len(roots)-1):
            if np.isfinite(roots[i]) and roots[i]<roots[i-1] and roots[i]<roots[i+1]:
                fun=lambda a:threshold(a,max(0.,roots[i]-.04),roots[i]+.04,.001,pc,P,level)
                res=minimize_scalar(fun,bounds=(angles[i-1],angles[i+1]),method='bounded',options={'xatol':2e-7})
                r=float(res.fun);a=float(res.x)
                minima.append(dict(angle=a,radius=r,normalized_energy=r*r/2,
                    initial_state=[0.,r*np.cos(a),np.sqrt(12)*r*np.sin(a)],success=bool(res.success)))
        transitions=[]
        for i in range(len(angles)):
            changes=[dict(radius_bracket=[float(radii[j]),float(radii[j+1])],
                          labels=[int(safe[i,j]),int(safe[i,j+1])])
                     for j in range(len(radii)-1) if safe[i,j]!=safe[i,j+1]]
            transitions.append(dict(angle=float(angles[i]),changes=changes))
        np.savez_compressed(DIR/'plane.npz',angles=angles,radii=radii,labels=safe,
                            labels_coarse=labels.reshape(aa.shape),labels_fine=labels2.reshape(aa.shape),
                            times=times2.reshape(aa.shape),thresholds=roots)
        dump(DIR/'plane.json',dict(parameters=PARS,grid_shape=list(aa.shape),
            steps=[.004,.002],max_time=600,cycle_section_tolerance=.002,required_consecutive_returns=2,
            disagreement_count=int(sum(disagreements)),
            labels={'0':'unresolved','1':'origin','2':'known attracting cycle'},
            energy_definition='E/(C2*Bp^2) = (x^2/alpha + y^2 + z^2/beta)/2; here x=0.',
            symmetry='The negative-y half-plane follows from the odd vector field.',
            minima=minima,transitions=transitions,elapsed_seconds=time.monotonic()-start,
            scope='Finite numerical basin section and local minima; no exhaustive proof of global energy optimality or star-shapedness.'))
        print('plane',minima,'disagreements',int(sum(disagreements)),flush=True)
    elif args.phase=='pulses':
        durations=np.geomspace(.02,10,31);impulses=np.linspace(0,3,151)
        lab=pulse_batch(durations,impulses,.002,pc,P,level)
        rows=[]
        for i,tau in enumerate(durations):
            intervals=[]
            for j in range(len(impulses)-1):
                if lab[i,j]==lab[i,j+1] or 0 in [lab[i,j],lab[i,j+1]]:continue
                lo,hi=impulses[j:j+2];left=lab[i,j]
                for _ in range(29):
                    mid=(lo+hi)/2
                    label=classify(*pulse_end(mid/tau,tau,.001),.001,pc,P,level)[0]
                    if label==0:raise RuntimeError('Pulse threshold remained unresolved within the integration horizon')
                    if label==left:lo=mid
                    else:hi=mid
                val=(lo+hi)/2
                intervals.append(dict(impulse=float(val),amplitude=float(val/tau),
                    squared_input_cost=float(val*val/tau),left=int(left),right=int(lab[i,j+1]),bracket=[float(lo),float(hi)]))
            rows.append(dict(duration=float(tau),transitions=intervals))
            dump(DIR/'pulses.json',dict(forcing='Add constant u=amplitude to dy/dt for 0 <= t <= duration, then u=0; start at the origin.',
                cost='Integral u(t)^2 dt = amplitude^2*duration; this is an input cost, not stored circuit energy.',
                impulse_grid=impulses.tolist(),rows=rows,
                scope='Numerical thresholds among rectangular pulses in one input channel; no optimum over arbitrary waveforms.'))
        np.savez_compressed(DIR/'pulses.npz',durations=durations,impulses=impulses,labels=lab)
        optimal=minimize_scalar(lambda tau:pulse_threshold(tau,pc,P,level)**2/tau,
                                bounds=(.8,1.6),method='bounded',options={'xatol':2e-5})
        impulse=pulse_threshold(optimal.x,pc,P,level)
        dump(DIR/'pulse-optimum.json',dict(duration=float(optimal.x),impulse=float(impulse),
            amplitude=float(impulse/optimal.x),squared_input_cost=float(optimal.fun),
            scope='Local numerical minimum over first-entry thresholds for rectangular y-channel pulses, duration restricted to [0.8,1.6].'))
        print('pulses complete',time.monotonic()-start,flush=True)
    else:
        data=json.loads((DIR/'plane.json').read_text());out=DIR/'independent-checks.json'
        rows=json.loads(out.read_text())['records'] if out.exists() else []
        for rec in data['minima']:
            for factor in [1-1e-5,1+1e-5]:
                point=factor*np.array(rec['initial_state'])
                if any(r.get('initial_state')==point.tolist() for r in rows):continue
                checks=independent(point,pc,P,level)
                expected=1 if factor<1 else 2
                assert all(c['label']==expected for c in checks),(rec,checks)
                rows.append(dict(initial_state=point.tolist(),expected=expected,checks=checks))
                dump(out,dict(records=rows,passed=True))
        pulses=json.loads((DIR/'pulses.json').read_text())
        for i in [0,10,20,23,30]:
            row=pulses['rows'][i]
            for tr in row['transitions'][:2]:
                for factor in [1-1e-5,1+1e-5]:
                    force=(tr['amplitude']*factor,row['duration'])
                    if any(r.get('forcing')==list(force) for r in rows):continue
                    checks=independent([0,0,0],pc,P,level,force)
                    expected=tr['left'] if factor<1 else tr['right']
                    assert all(c['label']==expected for c in checks),(force,checks)
                    rows.append(dict(forcing=force,expected=expected,checks=checks))
                    dump(out,dict(records=rows,passed=True))
        if (DIR/'pulse-optimum.json').exists():
            opt=json.loads((DIR/'pulse-optimum.json').read_text())
            for factor in [1-1e-5,1+1e-5]:
                force=(opt['amplitude']*factor,opt['duration'])
                if any(r.get('forcing')==list(force) for r in rows):continue
                checks=independent([0,0,0],pc,P,level,force);expected=1 if factor<1 else 2
                assert all(c['label']==expected for c in checks),(force,checks)
                rows.append(dict(forcing=force,expected=expected,checks=checks));dump(out,dict(records=rows,passed=True))
        dump(DIR/'independent-checks.json',dict(records=rows,passed=True))
        print('Independent initial-state/pulse checks',len(rows),flush=True)


if __name__=='__main__':main()
