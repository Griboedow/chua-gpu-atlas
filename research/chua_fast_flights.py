"""Exact affine flights with analytic duration derivatives for continuation.

The flow is evaluated by 3x3 eigendecomposition and the affine equilibrium.
The closing equations and half/full section conventions match Flights. The
boundary Jacobian differentiates the section projection as well as monodromy.
"""
import numpy as np
from scipy.optimize import root
from chua_layer_bifurcations import Flights


class FastFlights(Flights):
    def matrices(self,beta):
        if getattr(self,'_beta',None)!=beta:
            self._beta=beta;self._matrices={}
            for region in set(self.regions):
                A=self.model.jacobian(beta,region==0)
                c=np.array([0.,0.,0.])
                if region:
                    c[0]=region*(self.model.m1-self.model.m0)/(1+self.model.m1);c[2]=-c[0]
                eig,V=np.linalg.eig(A);W=np.linalg.inv(V)
                self._matrices[region]=(A,c,eig,V,W)
        return self._matrices

    def evaluate(self,beta,q,check=False):
        if not np.isfinite(q).all() or not np.isfinite(beta) or max(abs(q[2:]))>50:
            raise ValueError('Invalid trial flight duration or parameter')
        mats=self.matrices(beta);n=len(q);initial=np.r_[0.,q[:2]];state=initial.copy()
        D=np.zeros((3,n));D[1,0]=D[2,1]=1.;D0=D.copy();M=np.eye(3)
        residual=[];jac=[];arcs=[];violation=0.;speeds=[];logdet=0.;flows=[];prefix=[]
        for k,(dt,region) in enumerate(zip(q[2:],self.regions)):
            A,c,eig,V,W=mats[region];E=((V*np.exp(eig*dt))@W).real;old=state.copy()
            state=c+E@(old-c);D=E@D;D[:,k+2]+=A@(state-c)
            prefix.append(M);flows.append(E);M=E@M;logdet+=np.trace(A)*dt
            if check:
                ts=np.linspace(0,dt,101)
                x=c[0]+np.real((V[0,:]*(W@(old-c)))[None,:]*np.exp(ts[:,None]*eig[None,:])).sum(axis=1)
                error=abs(x)-1 if region==0 else 1-region*x
                violation=max(violation,float(np.max(error)))
                arcs.append(dict(region=region,duration=float(dt),start=old.tolist(),end=state.tolist()))
            if k<len(self.targets):
                residual.append(state[0]-self.targets[k]);jac.append(D[0].copy());speeds.append(abs((A@(state-c))[0]))
        residual.extend(state-self.sign*initial);jac.extend(D-self.sign*D0)
        f=self.model.field(beta,state);Jend=self.model.jacobian(beta,abs(state[0])<1)
        projection=np.eye(3)-np.outer(f,[1.,0.,0.])/f[0]
        P=(self.sign*projection@M)[1:,1:]
        # dM/d(duration_k), evaluated using prefix/suffix products.
        derivatives=[np.zeros((3,3)) for _ in range(n)];suffix=np.eye(3)
        for k in range(len(flows)-1,-1,-1):
            derivatives[k+2]=suffix@mats[self.regions[k]][0]@flows[k]@prefix[k]
            suffix=suffix@flows[k]
        grad=np.empty(n)
        for k in range(n):
            df=Jend@D[:,k]
            dproj=-np.outer(df,[1.,0.,0.])/f[0]+np.outer(f,[1.,0.,0.])*df[0]/f[0]**2
            dP=self.sign*(dproj@M+projection@derivatives[k])
            grad[k]=np.trace(dP[1:,1:])
        loggrad=np.r_[0.,0.,[np.trace(mats[r][0]) for r in self.regions]]
        return np.array(residual),np.array(jac),P,dict(arcs=arcs,violation=violation,minimum_speed=min(speeds),
                    logdet=logdet,mono=M,trace_gradient=grad,logdet_gradient=loggrad)

    def boundary_system(self,w,target,jacobian=False):
        q=w[:-1];beta=w[-1];r,J,P,aux=self.evaluate(beta,q)
        determinant=self.sign*np.exp(aux['logdet'])
        residual=np.r_[r,target*target-target*np.trace(P)+determinant]
        if not jacobian:return residual
        matrix=np.empty((len(w),len(w)));matrix[:-1,:-1]=J
        matrix[-1,:-1]=-target*aux['trace_gradient']+determinant*aux['logdet_gradient']
        h=1e-6
        wp=w.copy();wm=w.copy();wp[-1]+=h;wm[-1]-=h
        matrix[:,-1]=(self.boundary_system(wp,target)-self.boundary_system(wm,target))/(2*h)
        return matrix

    def boundary(self,beta,guess,target):
        w0=np.r_[guess,beta]
        sol=root(lambda w:self.boundary_system(w,target),w0,
                 jac=lambda w:self.boundary_system(w,target,True),tol=1e-12)
        w=sol.x
        for _ in range(10):
            residual=self.boundary_system(w,target);error=max(abs(residual))
            if error<3e-10:break
            change=np.linalg.solve(self.boundary_system(w,target,True),-residual)
            for scale in [1.,.5,.25,.125]:
                trial=w+scale*change
                if max(abs(self.boundary_system(trial,target)))<error:w=trial;break
            else:break
        error=max(abs(self.boundary_system(w,target)))
        if not np.isfinite(error) or not np.isfinite(w).all() or error>1e-8 or min(w[2:-1])<=0:
            raise RuntimeError(f'Fast boundary residual {error}; {sol.message}')
        return w[-1],w[:-1]
