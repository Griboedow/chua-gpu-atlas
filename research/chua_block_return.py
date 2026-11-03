"""High-precision event returns and multiple shooting across parent blocks.

The prescribed half-itinerary is followed with event times solved anew. Local
section derivatives remove the neutral flow direction at every diode crossing.
Eigen-flow evaluation is independent of PrecisionFlights' matrix exponentials.
These are numerical maps, not interval enclosures.
"""
import mpmath as mp


class BlockReturn:
    def __init__(self,row,digits=60):
        mp.mp.dps=digits;self.digits=digits;self.row=row;self.cache={}
        self.regions=row['regions'];self.targets=row['targets']+[0]
        self.guesses=[mp.mpf(v) for v in row['q'][2:]]

    def matrices(self,beta):
        key=mp.nstr(beta,self.digits)
        if key not in self.cache:
            mats={}
            for r in set(self.regions):
                slope=-mp.mpf(8)/7 if r==0 else -mp.mpf(5)/7
                A=mp.matrix([[-9*(1+slope),9,0],[1,-1,1],[0,-beta,0]])
                c=mp.matrix([r*mp.mpf('1.5'),0,-r*mp.mpf('1.5')]);eig,V=mp.eig(A)
                mats[r]=(A,c,eig,V,V**-1)
            if len(self.cache)>8:self.cache.clear()
            self.cache[key]=mats
        return self.cache[key]

    def evaluate(self,beta,point):
        mats=self.matrices(beta);state=mp.matrix([0,*point]);initial_speed=9*state[1]
        product=mp.eye(2);logdiv=mp.mpf(0);durations=[];speeds=[]
        tolerance=mp.mpf(10)**(-(self.digits-12))
        for reg,target,guess in zip(self.regions,self.targets,self.guesses):
            A,c,eig,V,W=mats[reg];old=state.copy();coeff=W*(old-c)
            dt=guess
            for j in range(20):
                ex=[mp.exp(l*dt) for l in eig]
                state=c+(V*mp.matrix([ex[i]*coeff[i] for i in range(3)])).apply(mp.re)
                velocity=A*(state-c);error=state[0]-target
                if abs(error)<tolerance:break
                delta=error/velocity[0];dt-=delta
                if dt<=0 or dt>max(10,5*guess):raise ValueError('Invalid event duration in prescribed block')
            else:raise RuntimeError('Event Newton exhausted')
            E=(V*mp.diag(ex)*W).apply(mp.re)
            local=mp.matrix([[E[i,j]-velocity[i]*E[0,j]/velocity[0] for j in [1,2]] for i in [1,2]])
            product=local*product;logdiv+=mp.trace(A)*dt if hasattr(mp,'trace') else sum(A[i,i] for i in range(3))*dt
            state[0]=target;durations.append(dt);speeds.append(abs(velocity[0]))
        determinant=mp.exp(logdiv)*initial_speed/velocity[0]
        return dict(point=-state[1:,:],matrix=-product,determinant=determinant,
            durations=durations,time=sum(durations),minimum_speed=min(speeds))


class BlockShooting:
    def __init__(self,calibration,k,digits=60):
        mp.mp.dps=digits;self.digits=digits;self.k=k;self.calibration=calibration
        self.map=BlockReturn(calibration['representative'],digits)
        self.beta_fold=mp.mpf(calibration['fold']['beta']);self.width=mp.mpf(calibration['width'])
        self.origin=mp.matrix([mp.mpf(v) for v in calibration['representative']['q'][:2]])
        d=mp.matrix([mp.mpf(v) for v in calibration['companion']['q'][:2]])-self.origin
        self.S=mp.matrix([[d[0],-d[1]],[d[1],d[0]]]);self.W=self.S**-1

    def beta(self,M):return self.beta_fold-self.width*(M+mp.mpf('.25'))

    def block(self,w,M):
        ev=self.map.evaluate(self.beta(M),self.origin+self.S*mp.matrix(w))
        ev['normalized']=self.W*(ev['point']-self.origin)
        ev['derivative']=self.W*ev['matrix']*self.S
        return ev

    def components(self,w):
        return [self.block(w[2*i:2*i+2],w[2*self.k]) for i in range(self.k)]

    def residual(self,w,target,parts=None):
        parts=self.components(w) if parts is None else parts
        product=mp.eye(2);det=mp.mpf(1);res=[]
        for i,e in enumerate(parts):
            j=(i+1)%self.k;res.extend(e['normalized']-mp.matrix(w[2*j:2*j+2]))
            product=e['derivative']*product;det*=e['determinant']
        trace=product[0,0]+product[1,1]
        return mp.matrix(res+[target*target-target*trace+det]),product,det

    def jacobian(self,w,target,parts):
        size=2*self.k+1;J=mp.zeros(size);prefix=[mp.eye(2)]
        for e in parts:prefix.append(e['derivative']*prefix[-1])
        suffix=[mp.eye(2) for _ in parts];post=mp.eye(2)
        for i in range(self.k-1,-1,-1):suffix[i]=post;post=post*parts[i]['derivative']
        _,_,det=self.residual(w,target,parts);h=mp.mpf(10)**(-self.digits//3)
        for i,e in enumerate(parts):
            for a in range(2):
                for b in range(2):J[2*i+a,2*i+b]+=e['derivative'][a,b]
                J[2*i+a,2*((i+1)%self.k)+a]-=1
            for j in range(2):
                p=list(w[2*i:2*i+2]);m=p.copy();p[j]+=h;m[j]-=h
                plus=self.block(p,w[2*self.k]);minus=self.block(m,w[2*self.k])
                d=suffix[i]*((plus['derivative']-minus['derivative'])/(2*h))*prefix[i]
                dd=det/e['determinant']*(plus['determinant']-minus['determinant'])/(2*h)
                J[size-1,2*i+j]=-target*(d[0,0]+d[1,1])+dd
        plus=w.copy();minus=w.copy();plus[size-1]+=h;minus[size-1]-=h
        deriv=(self.residual(plus,target)[0]-self.residual(minus,target)[0])/(2*h)
        for i in range(size):J[i,size-1]=deriv[i]
        return J

    def solve(self,guess,target='.5',max_iterations=24):
        w=mp.matrix([mp.mpf(str(v)) for v in guess]);target=mp.mpf(str(target));history=[]
        gate=mp.mpf(10)**(-(self.digits-20))
        for iteration in range(max_iterations):
            parts=self.components(w);r,P,det=self.residual(w,target,parts);error=max(abs(x) for x in r)
            history.append(mp.nstr(error,16))
            if error<gate:break
            delta=mp.lu_solve(self.jacobian(w,target,parts),-r)
            for fraction in [1,mp.mpf('.5'),mp.mpf('.25'),mp.mpf('.125'),mp.mpf('.0625')]:
                trial=w+fraction*delta
                try:err=max(abs(x) for x in self.residual(trial,target)[0])
                except (ValueError,RuntimeError,ZeroDivisionError):continue
                if err<error:w=trial;break
            else:raise RuntimeError('Block Newton line search stalled: '+str(history))
        else:raise RuntimeError('Block Newton exhausted: '+str(history))
        points=[self.origin+self.S*mp.matrix(w[2*i:2*i+2]) for i in range(self.k)]
        regions=[];durations=[]
        for i,e in enumerate(parts):
            rr=[(-1)**i*r for r in self.map.regions];dd=e['durations']
            if not regions:regions=rr.copy();durations=dd.copy()
            else:
                assert regions[-1]==rr[0]==0
                durations[-1]+=dd[0];durations.extend(dd[1:]);regions.extend(rr[1:])
        targets=[regions[j] or regions[j+1] for j in range(len(regions)-1)]
        q=[*points[0],*durations];fmt=lambda v:mp.nstr(v,self.digits-5)
        return dict(beta=fmt(self.beta(w[2*self.k])),q=[fmt(x) for x in q],period=fmt(2*sum(durations)),
            regions=regions,targets=targets,symmetric=True,target_multiplier=fmt(target),
            normalized_coordinate=fmt(w[2*self.k]),block_coordinates=[[fmt(w[2*i]),fmt(w[2*i+1])] for i in range(self.k)],
            block_points=[[fmt(v) for v in p] for p in points],
            block_durations=[[fmt(v) for v in e['durations']] for e in parts],
            block_times=[fmt(e['time']) for e in parts],block_count=self.k,
            residual=fmt(error),trace=fmt(P[0,0]+P[1,1]),determinant=fmt(det),
            working_digits=self.digits,newton_residuals=history)
