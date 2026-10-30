"""Affine Poincare map with automatic breakpoint itinerary changes.

Stationary points bracket boundary crossings, including excursions narrower
than the sampling interval. The coarse derivative grid only locates brackets;
flight times are refined by scalar root finding. This is numerical, not an
interval-certified event enumerator.
"""
import numpy as np
from scipy.optimize import brentq,root


class EventMap:
    def __init__(self,model,returns,bracket_step=.1):
        self.model=model;self.returns=returns;self.bracket_step=bracket_step

    def matrices(self,beta):
        if getattr(self,'beta',None)!=beta:
            self.beta=beta;self.mats={}
            for region in [-1,0,1]:
                A=self.model.jacobian(beta,region==0);c=np.zeros(3)
                if region:c[0]=region*(self.model.m1-self.model.m0)/(1+self.model.m1);c[2]=-c[0]
                eig,V=np.linalg.eig(A);W=np.linalg.inv(V);self.mats[region]=(A,c,eig,V,W)
        return self.mats

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
                targets=[-1.,0.,1.] if region==0 else [float(region)]
                for target in targets:
                    if (xl-target)*(xr-target)>=0:continue
                    time=brentq(lambda t:x(t)-target,left,right,xtol=1e-13,rtol=1e-14)
                    speed=dx(time)
                    if target==0 and speed<=0:continue
                    if target and ((region==0 and target*speed<=0) or (region and region*speed>=0)):continue
                    events.append((time,target))
                if events:return min(events)
        raise RuntimeError('No next event within 24 time units')

    def evaluate(self,w):
        y,z,beta=w
        if not np.isfinite(w).all() or y<=0 or abs(z)>5 or not 10<beta<20:raise ValueError('Invalid map trial')
        mats=self.matrices(beta);u=np.array([0.,y,z]);state=u.copy();region=0;M=np.eye(3);time=0.;logdet=0.;count=0;events=[];segments=[]
        for _ in range(400):
            A,c,eig,V,W=mats[region];dt,target=self.next_event(state,region,mats[region]);E=((V*np.exp(eig*dt))@W).real
            segments.append((time,dt,region,state.copy(),M.copy()))
            state=c+E@(state-c);state[0]=target;M=E@M;time+=dt;logdet+=np.trace(A)*dt
            events.append(dict(time=time,target=target,region=region))
            if target==0:
                count+=1
                if count==self.returns:
                    f=self.model.field(beta,state);P=((np.eye(3)-np.outer(f,[1.,0.,0.])/f[0])@M)[1:,1:]
                    t0,dt,reg,st,mono=next(s for s in segments if s[0]+s[1]>=time/2)
                    A,c,eig,V,W=mats[reg];E=((V*np.exp(eig*(time/2-t0)))@W).real
                    hs=c+E@(st-c);hm=E@mono;hf=self.model.field(beta,hs)
                    half=(-((np.eye(3)-np.outer(hf,[1.,0.,0.])/hf[0])@hm))[1:,1:]
                    return state[1:]-u[1:],P,dict(period=time,logdet=logdet,events=events,final=state.tolist(),
                        half_inversion_error=float(max(abs(hs+u))),half_trace=float(np.trace(half)))
            else:region=int(target) if region==0 else 0
            if time>250 or max(abs(state))>100:raise RuntimeError('Map escaped or exceeded time limit')
        raise RuntimeError('Too many map events')

    def residual(self,w,target):
        r,P,aux=self.evaluate(w);return np.r_[r,target*target-target*np.trace(P)+np.exp(aux['logdet'])]

    def jacobian(self,w,target):
        r,P,aux=self.evaluate(w);J=np.empty((3,3));J[:2,:2]=P-np.eye(2)
        for k in range(3):
            h=1e-7 if k<2 else 1e-6;e=np.eye(3)[k]*h
            d=(self.residual(w+e,target)-self.residual(w-e,target))/(2*h)
            J[-1,k]=d[-1]
            if k==2:J[:,k]=d
        return J

    def solve(self,guess,target):
        ans=root(lambda w:self.residual(w,target),guess,jac=lambda w:self.jacobian(w,target),tol=2e-10)
        error=float(max(abs(self.residual(ans.x,target))))
        if not np.isfinite(error) or error>2e-7:raise RuntimeError(f'Event-map residual {error}')
        r,P,aux=self.evaluate(ans.x);return ans.x,aux
