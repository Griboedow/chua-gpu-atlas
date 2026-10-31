"""Finite-iterate approximations to the hypothesized organizing tangency.

The initial line is tangent to the saddle's unstable direction. Several local
iterations align it with the unstable manifold; a global excursion folds it
back towards the saddle. Extra terminal iterations reduce stable-coordinate
error before comparison with the tangent stable line. This numerical procedure
does NOT certify an exact manifold tangency or an infinite set of windows.
"""
import json
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq,root
from chua_event_map import EventMap
from trace_chua_window import Chua
from predict_chua_layers import OUT,parent_at
from analyze_chua_prediction_evidence import section_word


class HalfMap(EventMap):
    def __init__(self):super().__init__(Chua(),1,.08)

    def next_event(self,state,region,mat):
        A,c,eig,V,W=mat;coeff=V[0]*(W@(state-c))
        def x(t):return float(c[0]+np.real(np.sum(coeff*np.exp(eig*t))))
        def dx(t):return float(np.real(np.sum(eig*coeff*np.exp(eig*t))))
        eps=1e-8
        for horizon in [8.,24.]:
            grid=np.arange(eps,horizon+self.bracket_step/2,self.bracket_step)
            speeds=np.real(np.exp(grid[:,None]*eig[None,:])@(eig*coeff))
            turning=[brentq(dx,a,b,xtol=1e-13,rtol=1e-14) for a,b,da,db in zip(grid[:-1],grid[1:],speeds[:-1],speeds[1:]) if da*db<0]
            points=[eps]+turning+[horizon]
            for left,right in zip(points[:-1],points[1:]):
                xl,xr=x(left),x(right);events=[]
                for target in ([-1.,0.,1.] if region==0 else [float(region)]):
                    if (xl-target)*(xr-target)>=0:continue
                    dt=brentq(lambda t:x(t)-target,left,right,xtol=1e-13,rtol=1e-14)
                    speed=dx(dt)
                    if target==0 and speed>=0:continue
                    if target and ((region==0 and target*speed<=0) or (region and region*speed>=0)):continue
                    events.append((dt,target))
                if events:return min(events)
        raise RuntimeError('No half-return event')

    def apply(self,point,beta):
        mats=self.matrices(beta);state=np.r_[0.,point];M=np.eye(3);region=0;time=0.;itinerary=[];minspeed=1e30
        for _ in range(60):
            A,c,eig,V,W=mats[region];dt,target=self.next_event(state,region,mats[region])
            E=((V*np.exp(eig*dt))@W).real;state=c+E@(state-c);state[0]=target;M=E@M;time+=dt
            itinerary.append(region);v=self.model.field(beta,state);minspeed=min(minspeed,abs(v[0]))
            if target==0:
                H=-((np.eye(3)-np.outer(v,[1.,0.,0.])/v[0])@M)[1:,1:]
                return -state[1:],H,dict(period=time,itinerary=itinerary,minimum_speed=minspeed)
            region=int(target) if region==0 else 0
        raise RuntimeError('Half-map itinerary exceeded limit')


@lru_cache(maxsize=128)
def saddle(beta):
    p=parent_at(beta);values,V=np.linalg.eig(np.array(p['section_matrix']))
    order=np.argsort(-abs(values));V=np.real_if_close(V[:,order]);values=np.real_if_close(values[order])
    if np.iscomplexobj(V) or np.iscomplexobj(values):raise ValueError('Expected real saddle eigendirections')
    if V[0,0]<0:V[:,0]*=-1
    return np.array(p['q'][:2]),values,V,np.linalg.inv(V),p


def main():
    frozen=json.loads((OUT/'predictions-frozen.json').read_text())
    records=json.loads((OUT/'search-results.json').read_text())['records'];r=records[-1]
    motion=section_word(r['core'],r['regions'],r['targets'],True)
    b0=frozen['calibration']['limit'];p,vals,V,W,parent=saddle(b0)
    a0=float(W[0]@(np.array(motion['normalized_section_points'][11])-p))
    hmap=HalfMap();rows=[]
    for extra in [0,2,4]:
        for post in [1,2]:
            iterations=8+extra+post-1
            def quantities(a,b):
                p,ev,V,W,pr=saddle(float(b));state=p+a*V[:,0];tangent=V[:,0].copy();codes=[];minspeed=1e30
                for _ in range(iterations):
                    state,H,aux=hmap.apply(state,b);tangent=H@tangent;codes.append(aux['itinerary']);minspeed=min(minspeed,aux['minimum_speed'])
                return float(W[0]@(state-p)),float(W[0]@tangent),dict(final_state=state.tolist(),
                    stable_displacement=float(W[1]@(state-p)),itineraries=codes,minimum_speed=minspeed)
            def fun(w):
                g,d,_=quantities(w[0],w[1]);return [g,d*1e-4]
            sol=root(fun,[a0/float(vals[0]**extra),b0],tol=1e-10)
            a,b=sol.x;g,d,aux=quantities(a,b);controls=[]
            for step in [1e-6,3e-7,1e-7]:
                left=quantities(a-step,b);right=quantities(a+step,b)
                curvature=(right[1]-left[1])/(2*step)
                param=(quantities(a,b+step)[0]-quantities(a,b-step)[0])/(2*step)
                controls.append(dict(step=step,curvature=curvature,parameter_derivative=param))
            row=dict(extra_local_iterations=extra,terminal_iterations=post,total_half_returns=iterations,
                starting_coordinate=float(a),beta=float(b),normal_displacement=g,tangency_derivative=d,
                solver_success=bool(sol.success),controls=controls,**aux)
            rows.append(row);print({k:v for k,v in row.items() if k not in ['itineraries','controls']},flush=True)
            (OUT/'organizer-approximants.json').write_text(json.dumps(dict(records=rows,
                scope='Finite-iterate unstable-line and stable-line approximations; not an exact or interval-certified homoclinic tangency.'),indent=2)+'\n')
    # Show the normal-distance fold at the first approximation, with fixed axes.
    first=rows[0];a=first['starting_coordinate'];b=first['beta'];curves=[]
    for offset in [-1e-7,0.,1e-7]:
        xx=np.linspace(a-2e-5,a+2e-5,101);yy=[]
        for aa in xx:
            pp,ev,VV,WW,pr=saddle(b+offset);state=pp+aa*VV[:,0]
            for _ in range(first['total_half_returns']):state,_,_=hmap.apply(state,b+offset)
            yy.append(float(WW[0]@(state-pp)))
        curves.append(dict(beta_offset=offset,initial_coordinate=xx.tolist(),normal_distance=yy))
    rootspread=float(np.ptp([r['beta'] for r in rows]));report=dict(records=rows,curves=curves,
        beta_spread=rootspread,approximate_beta=float(np.mean([r['beta'] for r in rows])),
        frozen_fit_limit=frozen['calibration']['limit'],
        scope='Finite-iterate approximations with nonzero sampled curvature and unfolding derivative; not a certified exact manifold connection or a proof of infinitely many windows.')
    (OUT/'organizer-approximants.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Beta spread',rootspread,flush=True)


if __name__=='__main__':main()
