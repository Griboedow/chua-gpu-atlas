"""Numerical pulse thresholds and saddle residence times at the hidden-cycle point."""
import argparse
import json
import time
import numpy as np
from scipy.optimize import brentq
from chua_research_followup import ROOT, OUT, dump, half_map

PARS=dict(alpha=8.4,beta=12.,m0=-.95,m1=-.05)


def setup():
    seeds=json.loads((ROOT/'results/chua-hidden-boundary/cycle-seeds.json').read_text())
    cert=json.loads((ROOT/'results/chua-basin-proof/separator-100.json').read_text())
    p=np.array(seeds['saddle']['q'][:2]);pc=np.array(seeds['attractor']['q'][:2])
    P=np.array(cert['origin_capture']['lyapunov']['P'],float)
    level=float(cert['origin_capture']['level']['midpoint'])
    J=np.array(seeds['saddle']['section_matrix'])
    e,V=np.linalg.eig(J);order=np.argsort(-abs(e));e,V=e[order],V[:,order]
    if V[:,0]@p<0:V[:,0]*=-1
    return p,pc,P,level,e,V,np.linalg.inv(V),seeds


def outcome(point,f,pc,P,level,max_returns=450,trace=False):
    state=np.array(point,float);total=0.;path=[]
    for k in range(max_returns):
        X=np.r_[0,state]
        if X@P@X<.9*level:return dict(label='origin',returns=k,time=total,trace=path)
        if np.linalg.norm(state-pc)<1e-5:return dict(label='cycle',returns=k,time=total,trace=path)
        state,J,aux=f.apply(state,12.);total+=aux['period']
        if trace:path.append([total,*state])
    return dict(label='unresolved',returns=max_returns,time=total,trace=path)


def boundary_on_ray(direction,bracket,iterations,f,p,W):
    def function(amplitude):
        state=amplitude*direction
        for _ in range(iterations):state,_,_=f.apply(state,12.)
        return W[0]@(state-p)
    value=brentq(function,*bracket,xtol=3e-15,rtol=1e-14)
    return value,float(function(value))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['screen','thresholds','delay'],default='screen')
    args=ap.parse_args();p,pc,P,level,e,V,W,seeds=setup();f=half_map(PARS)
    out=OUT/'trigger';out.mkdir(parents=True,exist_ok=True)
    if args.phase=='screen':
        rows=[]
        for angle in [-1.2,-.8,-.4,0.,.4,.8,1.2]:
            direction=np.array([np.cos(angle),np.sin(angle)])
            row=dict(angle=angle,direction=direction.tolist(),samples=[])
            for amplitude in np.linspace(0,4,81):
                r=outcome(amplitude*direction,f,pc,P,level)
                row['samples'].append(dict(amplitude=float(amplitude),**r))
            rows.append(row)
            dump(out/'screen.json',dict(parameters=PARS,records=rows,scope='Finite ray samples; not exhaustive basin coverage.'))
            labels=[r['label'] for r in row['samples']]
            print('angle',angle,{k:labels.count(k) for k in set(labels)},flush=True)
    elif args.phase=='thresholds':
        screen=json.loads((out/'screen.json').read_text());rows=[]
        for line in screen['records']:
            direction=np.array(line['direction']);a=line['samples']
            changes=[(u,v) for u,v in zip(a[:-1],a[1:]) if u['label']!=v['label']]
            for left,right in changes:
                lo,hi=left['amplitude'],right['amplitude']
                # Find a narrow basin bracket before solving a finite-manifold equation.
                for _ in range(20):
                    mid=(lo+hi)/2;label=outcome(mid*direction,f,pc,P,level)['label']
                    if label==left['label']:lo=mid
                    elif label==right['label']:hi=mid
                    else:raise ValueError('Unresolved bracket')
                approximants=[]
                for n in [4,6,8]:
                    try:
                        root,residual=boundary_on_ray(direction,(lo,hi),n,f,p,W)
                        approximants.append(dict(returns=n,amplitude=root,residual=residual))
                    except ValueError as exc:approximants.append(dict(returns=n,error=str(exc)))
                rows.append(dict(angle=line['angle'],direction=direction.tolist(),bracket=[lo,hi],
                    left=left['label'],right=right['label'],approximants=approximants))
                dump(out/'thresholds.json',dict(records=rows,scope='Numerical basin brackets and finite-iterate stable-manifold approximants.'))
                print(rows[-1],flush=True)
    else:
        rows=json.loads((out/'thresholds.json').read_text())['records']
        # Freeze the saddle-rate prediction before measuring any delay curve.
        half_period=seeds['saddle']['period']/2
        forecast=dict(unstable_half_multiplier=float(e[0]),saddle_half_period=half_period,
            predicted_time_per_log_distance=half_period/np.log(abs(e[0])),
            predicted_time_per_decimal_digit=half_period*np.log(10)/np.log(abs(e[0])),
            neighborhood_unstable_radius=.002,neighborhood_stable_radius=.002,
            formula='t_local = -T_half/log(abs(lambda_half))*log(abs(A-Astar)) + O(1)',
            scope='Asymptotic local residence time, not the full capture time.')
        path=out/'delay-forecast.json'
        if not path.exists():dump(path,forecast)
        result=[]
        for row in rows:
            if row['angle'] not in [0.,.8]:continue
            direction=np.array(row['direction'])
            astar=row['approximants'][-1].get('amplitude')
            if astar is None:continue
            for sign in [-1,1]:
                for distance in np.geomspace(1e-4,1e-10,25):
                    state=(astar+sign*distance)*direction;t=0.;inside=[];times=[]
                    for k in range(100):
                        coord=W@(state-p)
                        inside.append(bool(max(abs(coord))<.002));times.append(t)
                        if k>15 and not inside[-1] and any(inside):break
                        state,J,aux=f.apply(state,12.);t+=aux['period']
                    inds=np.flatnonzero(inside)
                    residence=times[inds[-1]+1]-times[inds[0]] if len(inds) and inds[-1]+1<len(times) else None
                    result.append(dict(angle=row['angle'],side=sign,distance=float(distance),amplitude=float(astar+sign*distance),
                        residence=residence,entry_time=times[inds[0]] if len(inds) else None,local_returns=len(inds)))
        fits=[]
        for angle in [0.,.8]:
            for sign in [-1,1]:
                data=[r for r in result if r['angle']==angle and r['side']==sign and r['residence'] is not None]
                if len(data)>3:
                    x=-np.log10([r['distance'] for r in data]);y=np.array([r['residence'] for r in data])
                    slope,intercept=np.polyfit(x,y,1)
                    fits.append(dict(angle=angle,side=sign,slope=float(slope),intercept=float(intercept),
                        relative_error=float(abs(slope/forecast['predicted_time_per_decimal_digit']-1))))
        dump(out/'delay.json',dict(prediction=forecast,records=result,fits=fits))
        print('delay fits',fits,flush=True)


if __name__=='__main__':main()
