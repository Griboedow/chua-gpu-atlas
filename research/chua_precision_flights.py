"""Arbitrary-precision affine shooting for a prescribed transverse itinerary.

Uses mpmath matrix exponentials, independently of FastFlights' double-precision
eigendecompositions. Duration derivatives are analytic; the parameter derivative
uses a precision-scaled symmetric difference. This is not interval arithmetic.
"""
import mpmath as mp


class PrecisionFlights:
    def __init__(self, regions, targets, symmetric=True, dps=55):
        mp.mp.dps=dps
        self.dps=dps; self.regions=list(regions); self.targets=list(targets)
        self.sign=-1 if symmetric else 1
        self.alpha=mp.mpf(9);self.m0=-mp.mpf(8)/7;self.m1=-mp.mpf(5)/7

    def exponential(self, A, duration):
        return mp.expm(A*duration)

    def evaluate(self, beta, q, gradients=True):
        n=len(q);initial=mp.matrix([0,q[0],q[1]]);state=initial.copy()
        D=mp.zeros(3,n);D[1,0]=D[2,1]=1;D0=D.copy();M=mp.eye(3)
        residual=[];jac=[];flows=[];prefix=[];logdet=mp.mpf(0);mats={}
        for reg in set(self.regions):
            slope=self.m0 if reg==0 else self.m1
            A=mp.matrix([[-self.alpha*(1+slope),self.alpha,0],[1,-1,1],[0,-beta,0]])
            c=mp.matrix([reg*(self.m1-self.m0)/(1+self.m1),0,-reg*(self.m1-self.m0)/(1+self.m1)])
            mats[reg]=(A,c)
        arcs=[]
        for k,(dt,reg) in enumerate(zip(q[2:],self.regions)):
            A,c=mats[reg];E=self.exponential(A,dt);old=state.copy();state=c+E*(state-c)
            D=E*D
            v=A*(state-c)
            for i in range(3):D[i,k+2]+=v[i]
            prefix.append(M);flows.append(E);M=E*M
            logdet+=sum(A[i,i] for i in range(3))*dt
            arcs.append(dict(region=reg,duration=dt,start=old,end=state.copy()))
            if k<len(self.targets):
                residual.append(state[0]-self.targets[k]);jac.append([D[0,j] for j in range(n)])
        for i in range(3):
            residual.append(state[i]-self.sign*initial[i]);jac.append([D[i,j]-self.sign*D0[i,j] for j in range(n)])
        A,c=mats[self.regions[-1]];v=A*(state-c)
        project=mp.eye(3)
        for i in range(3):project[i,0]-=v[i]/v[0]
        P=self.sign*project*M;tr=P[1,1]+P[2,2];det=self.sign*mp.exp(logdet)
        result=dict(residual=mp.matrix(residual),jacobian=mp.matrix(jac),trace=tr,
            determinant=det,logdet=logdet,arcs=arcs)
        if not gradients:return result
        dM=[mp.zeros(3) for _ in range(n)];suffix=mp.eye(3)
        for k in range(len(flows)-1,-1,-1):
            dM[k+2]=suffix*mats[self.regions[k]][0]*flows[k]*prefix[k]
            suffix=suffix*flows[k]
        grad=[]
        for k in range(n):
            dv=A*D[:,k];dproj=mp.zeros(3)
            for i in range(3):dproj[i,0]=-dv[i]/v[0]+v[i]*dv[0]/v[0]**2
            dP=self.sign*(dproj*M+project*dM[k])
            grad.append(dP[1,1]+dP[2,2])
        result['trace_gradient']=grad
        result['logdet_gradient']=[mp.mpf(0),mp.mpf(0)]+[sum(mats[r][0][i,i] for i in range(3)) for r in self.regions]
        return result

    def system(self,w,target,jacobian=False):
        ev=self.evaluate(w[len(w)-1],list(w[:-1]),gradients=jacobian)
        f=mp.matrix(list(ev['residual'])+[target*target-target*ev['trace']+ev['determinant']])
        if not jacobian:return f
        size=len(w);J=mp.zeros(size)
        for i in range(size-1):
            for j in range(size-1):J[i,j]=ev['jacobian'][i,j]
        for j in range(size-1):
            J[size-1,j]=-target*ev['trace_gradient'][j]+ev['determinant']*ev['logdet_gradient'][j]
        h=mp.mpf(10)**(-self.dps//3)
        plus=w.copy();minus=w.copy();plus[size-1]+=h;minus[size-1]-=h
        deriv=(self.system(plus,target)-self.system(minus,target))/(2*h)
        for i in range(size):J[i,size-1]=deriv[i]
        return f,J

    def boundary(self,beta,q,target):
        w=mp.matrix([mp.mpf(str(v)) for v in q]+[mp.mpf(str(beta))]);target=mp.mpf(str(target))
        history=[];gate=mp.mpf(10)**(-(self.dps-20))
        for iteration in range(18):
            f,J=self.system(w,target,True);err=max(abs(v) for v in f)
            history.append(mp.nstr(err,12))
            if err<gate:break
            change=mp.lu_solve(J,-f)
            for scale in [1,mp.mpf('.5'),mp.mpf('.25'),mp.mpf('.125'),mp.mpf('.0625')]:
                trial=w+scale*change
                if min(trial[2:-1])<=0:continue
                if max(abs(v) for v in self.system(trial,target))<err:
                    w=trial;break
            else:raise RuntimeError('High precision Newton line search stalled: '+str(history))
        else:raise RuntimeError('High precision Newton exhausted: '+str(history))
        ev=self.evaluate(w[len(w)-1],list(w[:-1]),False)
        return dict(beta=mp.nstr(w[len(w)-1],self.dps-5),q=[mp.nstr(v,self.dps-5) for v in w[:-1]],
            period=mp.nstr((2 if self.sign==-1 else 1)*sum(w[2:-1]),self.dps-5),
            target_multiplier=str(target),residual=mp.nstr(max(abs(v) for v in self.system(w,target)),12),
            determinant=mp.nstr(ev['determinant'],12),trace=mp.nstr(ev['trace'],20),
            decimal_digits=self.dps,newton_residuals=history)
