"""Validate an entire local saddle exit strip into two invariant capture sets.

All-time event checks and mean-value propagation enclose continua of starts.
The analytic graph-transform/exit argument is given in the companion report.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from flint import arb,arb_mat,ctx
from chua_interval_map import col,eye,midpoint,upper
from derive_chua_global_sector import bounds
from certify_chua_coexisting_cycles import ReturnMap
from certify_chua_hidden_cycle import positive_definite
from trace_chua_window import Chua

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/chua-basin-proof'
PARS=(8.4,12.,-.95,-.05)


def matrix(rows):
    return arb_mat([[arb(str(v)) for v in row] for row in rows])


def npmat(M):
    return np.array([[float(M[i,j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])


class TightReturn(ReturnMap):
    """Sign-reduced half return; symmetry preserves both destination basins."""
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        end=next(i for i,r in enumerate(self.schedule_rows) if r[1]==0)
        self.schedule_rows=self.schedule_rows[:end+1]

    def evaluate(self,point,verify_path=False,centered=True):
        point=list(point)
        centers=None
        delta=col([v-v.mid() for v in point])
        if centered and any(v.rad()>arb('1e-50') for v in point):
            centers=self.evaluate([v.mid() for v in point],False,False)[2]
        state=col([0,*point]);D=arb_mat([[0,0],[1,0],[0,1]]);arcs=[]
        DD=[[col([0,0,0]) for _ in range(2)] for _ in range(2)]
        for k,(reg,target,guess,sign) in enumerate(self.schedule_rows):
            A,c,_=self.matrices(None,reg);start=state
            T=self.crossing(A,c,start,target,guess)
            end,E=self.flow(A,c,start,T);v=A*(end-c)
            assert not v[0,0].contains(0)
            F=E*D
            dt=[-F[0,j]/v[0,0] for j in range(2)]
            Ft=A*F;Ftt=A*v;newDD=[]
            for i in range(2):
                row=[]
                for j in range(2):
                    mixed=E*DD[i][j]+col([Ft[k,i] for k in range(3)])*dt[j]+col([Ft[k,j] for k in range(3)])*dt[i]+Ftt*dt[i]*dt[j]
                    dd=mixed-v*(mixed[0,0]/v[0,0]);dd[0,0]=0;row.append(dd)
                newDD.append(row)
            DD=newDD
            D=arb_mat([[F[i,j]+v[i,0]*dt[j] for j in range(2)] for i in range(3)])
            D[0,0]=D[0,1]=0
            end[0,0]=target
            if centers is not None:
                Dc=centers[k]['D']
                for i in range(3):
                    for j in range(2):
                        enclosure=Dc[i,j]+sum(DD[j][z][i,0]*delta[z,0] for z in range(2))
                        D[i,j]=D[i,j].intersection(enclosure)
                tighter=centers[k]['end']+Dc*delta
                for i in range(2):
                    for j in range(2):tighter+=DD[i][j]*delta[i,0]*delta[j,0]/2
                end=col([end[i,0].intersection(tighter[i,0]) for i in range(3)])
            state=end
            arc=dict(region=reg,target=target,duration=T,start=start,end=end,speed=v[0,0],sign=sign,D=D,DD=DD)
            if verify_path:
                arc['path']=self.verify_arc(A,sign*c,dict(region=abs(reg),target=sign*target,duration=T,start=sign*start))
            arcs.append(arc)
        assert arcs[0]['sign']==1 and arcs[-1]['sign']==1
        zeros=[a for a in arcs if a['target']==0]
        assert len(zeros)==1 and zeros[0]['speed']<0
        return -col([state[1,0],state[2,0]]),-arb_mat([[D[i,j] for j in range(2)] for i in [1,2]]),arcs


class DynamicReturn(TightReturn):
    """Numerical guesses select flights; interval path checks certify them."""
    def __init__(self,point,digits):
        ctx.dps=digits
        self.alpha,self.beta,self.m0,self.m1=arb(42)/5,arb(12),-arb(19)/20,-arb(1)/20
        self.schedule_rows=[]
        self.mats={}
        for reg in [-1,0,1]:
            A,c,_=self.matrices(None,reg)
            An,cn=npmat(A),npmat(c).ravel()
            eig,V=np.linalg.eig(An)
            self.mats[reg]=(An,cn,eig,V,np.linalg.inv(V))
        state=np.r_[0.,[float(v.mid()) for v in point]];reg=0;zeros=0
        for _ in range(30):
            A,c,e,V,W=self.mats[reg];coeff=V[0]*(W@(state-c))
            x=lambda t:float(c[0]+np.real(np.sum(coeff*np.exp(e*t))))
            dx=lambda t:float(np.real(np.sum(e*coeff*np.exp(e*t))))
            grid=np.arange(1e-8,12.001,.035)
            speeds=np.real(np.exp(grid[:,None]*e[None,:])@(e*coeff))
            turns=[brentq(dx,l,r,xtol=1e-13) for l,r,a,b in zip(grid[:-1],grid[1:],speeds[:-1],speeds[1:]) if a*b<0]
            found=[]
            for l,r in zip([1e-8]+turns,turns+[12.]):
                for target in ([-1,0,1] if reg==0 else [reg]):
                    if (x(l)-target)*(x(r)-target)>=0:continue
                    t=brentq(lambda t:x(t)-target,l,r,xtol=1e-13)
                    speed=dx(t)
                    if target and ((reg==0 and target*speed<0) or (reg and reg*speed>0)):continue
                    found.append((t,target))
                if found:break
            if not found:raise ValueError('No guessed return event')
            time,target=min(found)
            sign=reg if reg else int(np.sign(x(time/2)))
            self.schedule_rows.append((reg,target,time,sign))
            state=c+((V*np.exp(e*time))@W).real@(state-c);state[0]=target
            if target==0:
                zeros+=1
                if zeros==1:return
            else:reg=int(target) if reg==0 else 0
        raise ValueError('No complete first return')


def refine_root(record,digits,shift='0'):
    f=TightReturn(record,digits,(arb(42)/5,arb(12),-arb(19)/20,-arb(1)/20))
    p=col([arb(str(x)) for x in record['cycle']['q'][:2]])
    for _ in range(9):
        y,J,_=f.evaluate(p.entries());p=midpoint(p-(J-eye(2)).inv()*(y-p))
    p[0,0]+=arb(shift)
    radius=arb('1e-40');box=col([arb(v,radius) for v in p.entries()])
    y,J0,_=f.evaluate(p.entries());_,J,arcs=f.evaluate(box.entries(),True)
    C=midpoint((J0-eye(2)).inv());R=eye(2)-C*(J-eye(2));K=p-C*(y-p)+R*(box-p)
    assert all(box[i,0].contains_interior(K[i,0]) for i in range(2))
    assert max(sum(abs(R[i,j]) for j in range(2)).abs_upper() for i in range(2))<1
    return f,p,box,J0


def eigenbasis(J):
    eig,V=np.linalg.eig(npmat(J));V=V[:,np.argsort(-abs(eig))].real
    if V[:,0]@np.array([1.,1.])<0:V[:,0]*=-1
    if V[0,1]<0:V[:,1]*=-1
    return matrix(V)


def local_block(f,p,pbox,J0,delta='1e-7'):
    V=eigenbasis(J0);W=V.inv();d=arb(delta)
    box=pbox+V*col([arb(0,d),arb(0,d)])
    _,J,arcs=f.evaluate(box.entries(),True);A=W*J*V
    a=A[0,0].lower();b=abs(A[0,1]).upper();c=abs(A[1,0]).upper();e=abs(A[1,1]).upper();L=arb('.02')
    gates=dict(expands=a-b*L>1,exit_faces=a-b>1,stable_slope=(b+L*e)/(a-L*c)<L,
               stable_contracts=1/(a-L*c)<1,stable_image=c*L+e<1,
               vertical_exit=c+e<1,transverse_growth=a-L*c>1)
    assert all(gates.values()),(gates,A)
    return V,W,d,dict(passed=True,delta=delta,slope_bound=str(L),
        derivative=[[bounds(A[i,j]) for j in range(2)] for i in range(2)],
        a=bounds(a),b=bounds(b),c=bounds(c),d=bounds(e),gates=gates,
        path_cells=[r['path']['accepted_cells'] for r in arcs]),(abs(A[0,0]).upper()+b)*d,(c+e)*d


def capture_cycle(f,p,J0):
    V=eigenbasis(J0);W=V.inv();radius=arb('.01')
    domain=p+V*col([arb(0,radius),arb(0,radius)])
    y,J,arcs=f.evaluate(domain.entries(),True);yc,_,_=f.evaluate(p.entries())
    M=W*J*V
    norm=max(sum(abs(M[i,j]) for j in range(2)).abs_upper() for i in range(2))
    residual=W*(yc-p)
    assert norm<1 and all(abs(v)+norm*radius<radius for v in residual.entries())
    return (p,V,W,radius),dict(radius='.01',contraction=bounds(norm),
        residual=[bounds(x) for x in residual.entries()],path_cells=[r['path']['accepted_cells'] for r in arcs])


def origin_ellipsoid(A):
    eig,V=np.linalg.eig(npmat(A))
    real=int(np.argmin(abs(eig.imag)));complex_index=int(np.argmax(eig.imag))
    S=np.column_stack([V[:,real].real,V[:,complex_index].real,V[:,complex_index].imag])
    W=np.linalg.inv(S);Pn=W.T@W;Pn=(Pn+Pn.T)/2
    entries=[[format(Pn[i,j],'.15g') for j in range(3)] for i in range(3)]
    P=arb_mat(entries);Q=-(A.transpose()*P+P*A)
    certificate=dict(P=entries,P_minors=positive_definite(P),Q_minors=positive_definite(Q))
    return P,certificate


def captured(centre,V,radii,target,P0,level):
    p,Vs,W,r=target
    offsets=W*(centre-p);M=W*V
    if all(abs(offsets[i,0])+sum(abs(M[i,j])*radii[j] for j in range(2))<r for i in range(2)):
        return 'cycle'
    box=centre+V*col([arb(0,v) for v in radii]);u=col([0,*box.entries()])
    energy=(u.transpose()*P0*u)[0,0]
    if energy<level:return 'origin'
    return None


def propagate(centre,V,radii,digits):
    f=DynamicReturn(centre.entries(),digits)
    physical=centre+V*col([arb(0,v) for v in radii])
    _,J,arcs=f.evaluate(physical.entries(),True)
    value,Jc,_=f.evaluate(centre.entries())
    U,_,_=np.linalg.svd(npmat(Jc*V)@np.diag([float(r.mid()) for r in radii]))
    Vnext=matrix(U);Wnext=Vnext.inv();M=Wnext*J*V
    newcentre=midpoint(value);err=Wnext*(value-newcentre)
    newr=[(abs(err[i,0])+sum(abs(M[i,j])*radii[j] for j in range(2))).upper() for i in range(2)]
    return newcentre,Vnext,newr,sum(a['path']['accepted_cells'] for a in arcs),float(sum(a['duration'] for a in arcs).mid())


def affine_origin_capture(centre,V,radii,P,level,step='.005'):
    """Picard tubes and affine variation of constants through diode contacts.

    The transported basis retains correlations during fixed-time integration.
    Mixed-region steps use an interval forcing and require no event itinerary.
    """
    f=DynamicReturn(centre.entries(),ctx.dps)
    z=col([0,*centre.entries()]);B=arb_mat([[1,0,0],[0,V[0,0],V[0,1]],[0,V[1,0],V[1,1]]])
    r=[arb(0),*radii];h=arb(step);elapsed=arb(0);mixed_steps=0;trace=[]
    A0,c0,_=f.matrices(None,0);gap=f.alpha*(f.m1-f.m0)
    E0=(A0*h).exp();kernel=(A0*arb(h/2,h/2)).exp()*col([1,0,0])*h
    Ao,co,_=f.matrices(None,1);Eo=(Ao*h).exp()
    def forcing(x):
        positive=lambda v: arb(max(arb(0),v.lower())).union(max(arb(0),v.upper()))
        return -gap*(positive(x-1)-positive(-x-1))
    def field(X):
        v=A0*X;v[0,0]+=forcing(X[0,0]);return v
    for count in range(4000):
        Q=B.transpose()*P*B
        error_norm=sum(Q[i,j].abs_upper()*r[i]*r[j] for i in range(3) for j in range(3)).upper().sqrt().upper()
        norm=(z.transpose()*P*z)[0,0].sqrt()+error_norm
        if norm<level.sqrt():
            return dict(passed=True,method='transported_basis',steps=count,step=step,time=bounds(elapsed),
                        mixed_steps=mixed_steps,centre=[bounds(x) for x in z.entries()],
                        norm_radius=bounds(error_norm),energy_upper=bounds(norm**2),trace=trace)
        X=z+B*col([arb(0,v) for v in r]);Y=X+field(X)*arb(h/2,h/2)
        for attempt in range(30):
            Y=col([arb(v.mid(),v.rad()*arb('1.1')+arb('1e-30')) for v in Y.entries()])
            Z=X+field(Y)*arb(h/2,h/2)
            if all(Y[i,0].contains_interior(Z[i,0]) for i in range(3)):break
            Y=col([Y[i,0].union(Z[i,0]) for i in range(3)])
        else:raise ValueError('Continuous-flow Picard enclosure failed')
        if Y[0,0]>1 or Y[0,0]<-1:
            c=co if Y[0,0]>1 else -co;E=Eo
            value=c+E*(z-c)
        else:
            q=forcing(Y[0,0]);mixed_steps+=int(not q.is_zero());E=E0
            value=E*z+kernel*q
        Bnext=midpoint(E*B);W=Bnext.inv();zn=midpoint(value);err=W*(value-zn);M=W*E*B
        r=[(err[i,0].abs_upper()+sum(M[i,j].abs_upper()*r[j] for j in range(3))).upper() for i in range(3)]
        z=zn;B=Bnext;elapsed+=h
        if count%20==0:trace.append(dict(time=float(elapsed),point=[float(v.mid()) for v in z.entries()],radius_before_step=float(error_norm)))
    raise ValueError('Continuous flow did not enter origin ellipsoid')


def norm_origin_capture(centre,V,radii,P,level,step='.005'):
    """Alternative P-norm enclosure; useful when a transported basis stretches."""
    f=DynamicReturn(centre.entries(),ctx.dps)
    z=col([0,*centre.entries()]);h=arb(step);elapsed=arb(0);mixed_steps=0;trace=[]
    A,_,_=f.matrices(None,0);gap=f.alpha*(f.m1-f.m0)
    kappa=arb('.02');Q=-(A.transpose()*P+P*A)
    positive_definite(Q-2*kappa*P)
    inv=P.inv();coordinate=[inv[i,i].sqrt().upper() for i in range(3)]
    def norm_box(X):
        return sum(P[i,j].abs_upper()*X[i,0].abs_upper()*X[j,0].abs_upper() for i in range(3) for j in range(3)).upper().sqrt().upper()
    err=col([0,*(V*col([arb(0,v) for v in radii])).entries()]);rho=norm_box(err)
    E=(A*h).exp();kernel=(A*arb(h/2,h/2)).exp()*col([1,0,0])*h;decay=(-kappa*h).exp()
    Ao,co,_=f.matrices(None,1);S=Ao.transpose()*P+P*Ao
    mu=arb(int(np.ceil(max(np.linalg.eigvals(npmat(P.inv()*S)).real)/2))+1)
    positive_definite(2*mu*P-S);Eo=(Ao*h).exp();growth=(mu*h).exp()
    def forcing(x):
        positive=lambda v: arb(max(arb(0),v.lower())).union(max(arb(0),v.upper()))
        return -gap*(positive(x-1)-positive(-x-1))
    def field(X):
        v=A*X;v[0,0]+=forcing(X[0,0]);return v
    for count in range(4000):
        energy=(z.transpose()*P*z)[0,0]
        if energy.sqrt()+rho<level.sqrt():
            return dict(passed=True,method='quadratic_norm',steps=count,step=step,time=bounds(elapsed),
                        mixed_steps=mixed_steps,centre=[bounds(x) for x in z.entries()],
                        norm_radius=bounds(rho),energy_upper=bounds((energy.sqrt()+rho)**2),
                        kappa=bounds(kappa),outer_log_norm_bound=bounds(mu),trace=trace)
        if rho>level.sqrt():raise ValueError('Norm enclosure exceeded capture scale')
        X=col([z[i,0]+arb(0,rho*coordinate[i]) for i in range(3)])
        Y=X+field(X)*arb(h/2,h/2)
        for attempt in range(30):
            Y=col([arb(v.mid(),v.rad()*arb('1.1')+arb('1e-30')) for v in Y.entries()])
            Z=X+field(Y)*arb(h/2,h/2)
            if all(Y[i,0].contains_interior(Z[i,0]) for i in range(3)):break
            Y=col([Y[i,0].union(Z[i,0]) for i in range(3)])
        else:raise ValueError('Continuous-flow Picard enclosure failed')
        if Y[0,0]>1 or Y[0,0]<-1:
            c=co if Y[0,0]>1 else -co;value=c+Eo*(z-c);factor=growth
        else:
            q=forcing(Y[0,0]);mixed_steps+=int(not q.is_zero());value=E*z+kernel*q;factor=decay
        zn=midpoint(value);rho=(factor*rho+norm_box(value-zn)).upper();z=zn;elapsed+=h
        if count%20==0:trace.append(dict(time=float(elapsed),point=[float(v.mid()) for v in z.entries()],radius=float(rho)))
    raise ValueError('Continuous flow did not enter origin ellipsoid')


def continuous_origin_capture(centre,V,radii,P,level,depth=0):
    failures=[]
    for method,step in [(norm_origin_capture,'.005'),(affine_origin_capture,'.005'),(affine_origin_capture,'.001')]:
        try:
            result=method(centre,V,radii,P,level,step)
            result['rejected_enclosures']=failures
            return result
        except (ValueError,AssertionError) as exc:
            failures.append(dict(method=method.__name__,step=step,reason=str(exc)))
    if depth<4:
        half=radii[0]/2;children=[]
        for side in [-1,1]:
            child=centre+V*col([side*half,0])
            children.append(continuous_origin_capture(child,V,[half,radii[1]],P,level,depth+1))
        return dict(passed=True,method='subdivision',split_coordinate=0,children=children,
                    rejected_enclosures=failures)
    raise ValueError('All continuous enclosures failed: '+str(failures))


def run(digits=80,pieces=64,only=None):
    ctx.dps=digits;OUT.mkdir(exist_ok=True)
    seeds=json.loads((ROOT/'results/chua-hidden-boundary/cycle-seeds.json').read_text())['full_records']
    fs,ps,pbox,Js=refine_root(seeds[0],digits)
    fc,pc,cbox,Jc=refine_root(seeds[1],digits)
    V,W,d,local,umax,smax=local_block(fs,ps,pbox,Js)
    root_error=W*(pbox-ps)
    assert all(v.abs_upper()<arb('1e-38') for v in root_error.entries())
    print('local bounds',local['a']['ball'],local['b']['ball'],local['c']['ball'],local['d']['ball'],flush=True)
    target,targetcert=capture_cycle(fc,pc,Jc)
    A0,_,_=fs.matrices(None,0);P0,eqcert=origin_ellipsoid(A0);level=arb('.99')/P0.inv()[0,0]
    results=[]
    for side in [-1,1]:
        for index in range(pieces):
            if only is not None and index!=only:continue
            lower=d+(umax-d)*index/pieces;upperu=d+(umax-d)*(index+1)/pieces
            ua=lower.union(upperu)*side
            centre=ps+V*col([ua.mid(),arb(0)])
            basis=V;radii=[ua.rad()+arb('1e-38'),smax+arb('1e-38')]
            trace=[];totalcells=0;duration=0.;error=None;continuous=None
            try:
                for step in range(400):
                    label=captured(centre,basis,radii,target,P0,level)
                    if label:break
                    try:
                        centre,basis,radii,cells,dt=propagate(centre,basis,radii,digits)
                    except (ValueError,AssertionError):
                        if side>0:raise
                        continuous=continuous_origin_capture(centre,basis,radii,P0,level)
                        label='origin';break
                    totalcells+=cells;duration+=dt
                    trace.append(dict(point=[float(v.mid()) for v in centre.entries()],
                                      radii=[float(v.mid()) for v in radii],time=duration))
                else:raise ValueError('No certified capture within 400 returns')
                expected='origin' if side<0 else 'cycle'
                assert label==expected
            except (ValueError,AssertionError) as exc:
                label='unresolved';error=str(exc)
            row=dict(side=side,index=index,initial_u=bounds(ua),initial_s_radius=bounds(smax),
                     label=label,error=error,returns=len(trace),path_cells=totalcells,trace=trace,continuous=continuous)
            results.append(row)
            print('exit',side,index,label,'returns',len(trace),error,flush=True)
    report=dict(passed=all(r['label']!='unresolved' for r in results),working_digits=digits,
        parameters=dict(alpha='42/5',beta='12',m0='-19/20',m1='-1/20',gamma='0'),
        pieces_per_side=pieces,partial=only is not None,local=local,
        saddle_root=[bounds(v) for v in pbox.entries()],basis=[[bounds(V[i,j]) for j in range(2)] for i in range(2)],
        exit_u_max=bounds(umax),exit_s_radius=bounds(smax),cycle_capture=targetcert,
        origin_capture=dict(lyapunov=eqcert,level=bounds(level)),exits=results,
        scope='Uniform exit-strip capture plus local graph-transform inequalities. Read the report for the topological basin-boundary conclusion.')
    path=OUT/f'separator-{digits}{"-pilot" if only is not None else ""}.json'
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    return report


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--digits',type=int,default=80)
    ap.add_argument('--pieces',type=int,default=64);ap.add_argument('--only',type=int)
    a=ap.parse_args();r=run(a.digits,a.pieces,a.only)
    if not r['passed']:raise SystemExit(2)
