"""Numerical first return to x=1, dx/dt>0, within a positive Chua lobe."""
import numpy as np
from scipy.optimize import brentq
from trace_chua_window import Chua
from chua_event_map import EventMap


class LobeMap(EventMap):
    def __init__(self,parameters):
        super().__init__(Chua(parameters['alpha'],parameters['m0'],parameters['m1']),1,.04)

    def event(self,state,region,mat):
        A,c,e,V,W=mat;coeff=V[0]*(W@(state-c))
        x=lambda t:float(c[0]+np.real(np.sum(coeff*np.exp(e*t))))
        dx=lambda t:float(np.real(np.sum(e*coeff*np.exp(e*t))))
        grid=np.arange(1e-8,12.001,.04);speeds=np.real(np.exp(grid[:,None]*e)@(e*coeff))
        turns=[brentq(dx,a,b,xtol=1e-13) for a,b,da,db in zip(grid[:-1],grid[1:],speeds[:-1],speeds[1:]) if da*db<0]
        for a,b in zip([1e-8]+turns,turns+[12.]):
            hits=[]
            for target in ([1] if region else [-1,1]):
                if (x(a)-target)*(x(b)-target)>=0:continue
                t=brentq(lambda t:x(t)-target,a,b,xtol=1e-13)
                if region and dx(t)>=0:continue
                if not region and target*dx(t)<=0:continue
                hits.append((t,target))
            if hits:return min(hits)
        raise ValueError('No lobe return before horizon')

    def apply(self,point,beta):
        if not np.isfinite(point).all() or max(abs(point))>100:raise ValueError('Invalid section state')
        state=np.r_[1.,point];M=np.eye(3);mats=self.matrices(beta);arcs=[]
        if self.model.field(beta,state)[0]<=0:raise ValueError('Wrong initial section direction')
        for reg in [1,0]:
            A,c,e,V,W=mats[reg];t,target=self.event(state,reg,mats[reg])
            if target!=1:raise ValueError('Trajectory left the positive-lobe itinerary')
            E=((V*np.exp(e*t))@W).real;old=state.copy();state=c+E@(state-c);state[0]=1.;M=E@M
            arcs.append(dict(region=reg,duration=t,start=old.tolist(),end=state.tolist()))
        velocity=self.model.field(beta,state)
        J=((np.eye(3)-np.outer(velocity,[1,0,0])/velocity[0])@M)[1:,1:]
        return state[1:],J,dict(period=sum(a['duration'] for a in arcs),arcs=arcs)

    def fixed(self,point,beta):
        p=np.array(point,float)
        for _ in range(30):
            y,J,aux=self.apply(p,beta);r=y-p;err=max(abs(r))
            if err<1e-11:return p,J,aux
            change=np.linalg.solve(J-np.eye(2),-r)
            for factor in [1,.5,.25,.125,.0625,.03125,.015625]:
                trial=p+factor*change
                try:improved=max(abs(self.apply(trial,beta)[0]-trial))<err
                except (ValueError,RuntimeError,OverflowError):improved=False
                if improved:p=trial;break
            else:raise ValueError('Lobe fixed-point corrector stalled')
        raise ValueError('Lobe fixed-point iteration exhausted')
