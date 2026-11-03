"""Locally refine sampled start/end phases; retain lower estimates of maxima."""
import argparse
import json
import time
import mpmath as mp
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
from chua_precision_flights import PrecisionFlights
from study_chua_return_geometry import ROOT, OUT, IDS, read, norm2


def frame(v):
    v=v/np.linalg.norm(v)
    axis=np.eye(3)[np.argmin(abs(v))]
    a=np.cross(v,axis);a/=np.linalg.norm(a)
    return np.column_stack([a,np.cross(v,a)])


def arcs_for(row):
    f=PrecisionFlights(row['regions'],row['targets'],row['symmetric'],75)
    beta=mp.mpf(row['beta'])
    ev=f.evaluate(beta,[mp.mpf(v) for v in row['q']],False)
    arcs=[];elapsed=0.
    for a in ev['arcs']:
        r=a['region'];s=f.m0 if r==0 else f.m1
        A=np.array([[-float(f.alpha*(1+s)),float(f.alpha),0],[1,-1,1],[0,-float(beta),0]])
        c=np.array([r*float((f.m1-f.m0)/(1+f.m1)),0,-r*float((f.m1-f.m0)/(1+f.m1))])
        dt=float(a['duration']);u=np.array(list(a['start']),float)
        arcs.append(dict(start=elapsed,end=elapsed+dt,A=A,c=c,u=u,E=expm(A*dt)))
        elapsed+=dt
    if row['symmetric']:
        arcs += [dict(a,start=a['start']+elapsed,end=a['end']+elapsed,c=-a['c'],u=-a['u']) for a in arcs.copy()]
        elapsed*=2
    return arcs,elapsed


def transport(arcs,period,start,end,details=False):
    start=float(start);end=float(end)
    def locate(t):
        phase=t%period
        i=int(np.searchsorted([a['end'] for a in arcs],phase,side='right'))
        i=min(i,len(arcs)-1)
        a=arcs[i]
        u=a['c']+expm(a['A']*(phase-a['start']))@(a['u']-a['c'])
        return i,u,frame(a['A']@(u-a['c']))
    i,u,B=locate(start)
    initial_B=B.copy();initial_u=u.copy()
    product=np.eye(2);t=start
    while t<end-1e-12:
        a=arcs[i%len(arcs)]
        cycle=np.floor((t+1e-10)/period)
        endpoint=min(end,a['end']+cycle*period)
        dt=endpoint-t
        if dt < -1e-9:
            raise ValueError('Invalid interval')
        E=expm(a['A']*dt)
        # Evaluate the endpoint on the reference arc, avoiding accumulated drift.
        local_end=endpoint-cycle*period-a['start']
        v=a['c']+expm(a['A']*local_end)@(a['u']-a['c'])
        C=frame(a['A']@(v-a['c']))
        product=(C.T@E@B)@product
        B=C;t=endpoint;i+=1
    gain=float(norm2(product))
    if not details:return gain
    _,_,vt=np.linalg.svd(product)
    return dict(gain=gain,start=start,end=end,duration=end-start,phase=start/period,
                initial_state=initial_u.tolist(),final_state=v.tolist(),
                initial_basis=initial_B.tolist(),final_basis=B.tolist(),
                input_direction=(initial_B@vt[0]).tolist(),matrix=product.tolist())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ids',nargs='*')
    args=p.parse_args()
    for path in sorted(OUT.glob('*-d70-s8.json')):
        rough=read(path);id=rough['id']
        if args.ids and id not in args.ids:continue
        target=OUT/f'{id}-optimized.json'
        if target.exists():
            print('Retained',target.name,flush=True);continue
        started=time.monotonic()
        row=read(ROOT/rough['source'])
        arcs,period=arcs_for(row)
        with np.load(OUT/f'{id}-d70-s8.npz') as z:
            t=z['t'];i=rough['maximizing_start_index'];k=rough['maximizing_steps'];n=len(t)
            tt=np.r_[t,t+period,2*period]
            start=tt[i];end=tt[i+k]
            ds=2*max(tt[i+1]-tt[i],tt[i]-tt[i-1] if i else tt[n]-tt[n-1])
            de=2*max(tt[i+k+1]-tt[i+k],tt[i+k]-tt[i+k-1]) if i+k+1<len(tt) else ds
        bounds=[(max(0,start-ds),min(period,start+ds)),(end-de,end+de)]
        def objective(x):
            if not 0<x[1]-x[0]<=period:return 1e3
            return -np.log(transport(arcs,period,*x))
        sol=minimize(objective,[start,end],method='Powell',bounds=bounds,
                     options={'xtol':1e-10,'ftol':1e-12,'maxiter':80})
        result=dict(id=id,beta=row['beta'],period=period,source=path.name,
                    sampled_gain=rough['maximum_gain'],optimizer='Bounded Powell in neighboring sampled phase cells',
                    success=bool(sol.success),message=str(sol.message),bounds=bounds,
                    maximum=transport(arcs,period,*sol.x,True),
                    elapsed_seconds=time.monotonic()-started,
                    scope='Local refinement of a sampled maximum, not a certified global supremum')
        target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(id,result['maximum']['gain'],'seconds',round(result['elapsed_seconds'],1),flush=True)


if __name__=='__main__':main()
