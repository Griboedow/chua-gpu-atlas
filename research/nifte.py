"""Literal smooth NIFTE NTP (5D) / four-pair NTP-TL (9D) models.

Wang et al. 2019, DOI 10.1186/s42480-019-0019-3, Eqs. 7-8.
Source pressure/flow scalings are retained; integration time is physical seconds.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]/'results/nifte'
COLUMNS=['F','K','Rf','Rl','Rth','Ld','Lf','Ll','Lp','Cad','Cd','Cp','Rtl0',
         'Rtl1','Rtl2','Rtl3','Rtl4','Ctl1','Ctl2','Ctl3','Ctl4','P0','U0']


def parameters(benchmark=False):
    return np.array([1.37,330,2.13e6,4.08e6,5.02e8,1.80e5,4.74e5 if benchmark else 4.74e6,
        1.27e7,3.77e5,1.76e-9,7.35e-8,7.43e-8,3.15e7,8.82e6,1.72e7,3.83e7,2.81e8,
        8.61e-10,5.95e-9,1.79e-8,2.66e-8,1.013e5,8e-4])


def matrix(p,loss=False):
    _,_,rf,rl,rth,ld,lf,ll,lp,cad,cd,cp=p[:12]
    p0,u0=p[-2:];a=ld+lf;den=ll*a+lp*(a+ll);c=u0/p0;b=p0/u0
    n=9 if loss else 5;A=np.zeros((n,n))
    A[0,0]=-1/(rth*cad);A[0,3:5]=c/cad
    A[1,3]=c/cd;A[2,4]=c/cp
    A[3,:5]=[-b*lp,-b*(lp+ll),b*ll,-ll*rf-lp*(rf+rl),-lp*rl];A[3,:5]/=den
    A[4,:5]=[-b*a,b*ll,-b*(a+ll),ll*rf-a*rl,-a*rl];A[4,:5]/=den
    if loss:
        r0=p[12];rs=p[13:17];cs=p[17:21]
        A[0,0]-=1/(r0*cad);A[0,5:]=1/(r0*cad)
        for i in range(4):
            A[5+i,0]=1/(r0*cs[i]);A[5+i,5:]=-1/(r0*cs[i]);A[5+i,5+i]-=1/(rs[i]*cs[i])
    return A,1/(rth*cad)


def rhs_jac(y,p,loss=False,linear=None):
    A,g=matrix(p,loss) if linear is None else linear
    th=np.tanh(p[1]*y[1]);f=A@y;f[0]+=g*p[0]*th
    j=A.copy();j[0,1]+=g*p[0]*p[1]*(1-th*th)
    return f,j


def mass_rhs(y,p,loss=False):
    """Independent dimensional pressure/flow balance, then rescale state only."""
    rf,rl,rth,ld,lf,ll,lp,cad,cd,cp=p[2:12];p0,u0=p[-2:]
    scale=np.r_[p0,p0,p0,u0,u0,np.repeat(p0,4)] if loss else np.array([p0,p0,p0,u0,u0])
    z=y*scale;P,D,Q,U,V=z[:5];a=ld+lf
    f=np.zeros(len(y),dtype=np.result_type(y,float))
    heat=p0*p[0]*np.tanh(p[1]*D/p0)
    f[0]=(heat-P)/(rth*cad)+(U+V)/cad;f[1]=U/cd;f[2]=V/cp
    mass=np.array([[a+ll,ll],[ll,lp+ll]])
    force=-np.array([P+D+(rf+rl)*U+rl*V,P+Q+rl*(U+V)])
    f[3:5]=np.linalg.solve(mass,force)
    if loss:
        r0=p[12];rs=p[13:17];cs=p[17:21];current=(P-z[5:].sum())/r0
        f[0]-=current/cad;f[5:]=(current-z[5:]/rs)/cs
    return f/scale


def save(name,data):
    name=Path(name);name.parent.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(str(name)+'.npz',**arrays)
    meta={k:v for k,v in data.items() if k not in arrays};meta['array_shapes']={k:list(v.shape) for k,v in arrays.items()}
    Path(str(name)+'.json').write_text(json.dumps(meta,indent=2,allow_nan=False),encoding='utf-8')


def validate():
    rng=np.random.default_rng(2026103007);err={'mass_equivalence':0.,'jacobian':0.,'odd_symmetry':0.,'constant_divergence':0.}
    for loss in [False,True]:
        n=9 if loss else 5
        for _ in range(64):
            p=parameters(True);p[:21]*=np.exp(rng.uniform(-.5,.5,21));y=rng.normal(size=n)*.03
            f,j=rhs_jac(y,p,loss);ref=mass_rhs(y,p,loss)
            err['mass_equivalence']=max(err['mass_equivalence'],float(np.max(abs(f-ref)/(1+abs(f)))))
            delta=1e-20;fd=np.column_stack([mass_rhs(y.astype(complex)+1j*np.eye(n)[i]*delta,p,loss).imag/delta for i in range(n)])
            err['jacobian']=max(err['jacobian'],float(np.max(abs(j-fd)/(1+abs(j)))))
            err['odd_symmetry']=max(err['odd_symmetry'],float(np.max(abs(rhs_jac(-y,p,loss)[0]+f))))
            err['constant_divergence']=max(err['constant_divergence'],float(abs(np.trace(j)-np.trace(matrix(p,loss)[0]))))
    p=parameters(True);y=np.array([.01,.005,-.003,.02,-.01]);lin=matrix(p)
    a=solve_ivp(lambda t,y:rhs_jac(y,p,linear=lin)[0],(0,10),y,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01)
    b=solve_ivp(lambda t,y:mass_rhs(y,p),(0,10),y,method='Radau',rtol=1e-11,atol=1e-13,max_step=.01)
    err['dop853_radau_state']=float(np.max(abs(a.y[:,-1]-b.y[:,-1])))
    passed=bool(max(err.values())<1e-7 and a.success and b.success)
    save(ROOT/'validation',dict(errors=err,passed=passed,random_states_per_model=64,models=[5,9],time_unit='physical s',
        scope='Independent dimensional mass balances, complex-step Jacobian and two CPU methods'))
    print(dict(errors=err,passed=passed),flush=True)
    if not passed:raise RuntimeError('NIFTE validation failed')


def control(loss=False,duration=400.):
    p=parameters(True);n=9 if loss else 5;y=np.zeros(n);y[:5]=[.01,.005,-.003,.02,-.01];lin=matrix(p,loss)
    def event(t,y):return y[1]
    event.direction=1
    wall=time.perf_counter();sol=solve_ivp(lambda t,y:rhs_jac(y,p,loss,lin)[0],(0,duration),y,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.02,events=event,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    times=np.arange(duration/2,duration,.001);states=sol.sol(times).T
    periods=np.diff(sol.t_events[0]);last=periods[-10:]
    out=dict(parameters=p,initial_state=y,times=times,states=states,section_times=sol.t_events[0],section_states=sol.y_events[0],
        last_periods=last,model='NTP-TL 9D' if loss else 'NTP 5D',source_override='Lf=4.74e5, all other parameters nominal',
        method='DOP853',rtol=1e-11,atol=1e-13,max_step=.02,duration=duration,time_unit='physical s',
        state_scaling='Source pressure/P0 and flow/U0; time is not divided by tau',
        initial_condition_scope='Chosen by us; source does not provide this complete initial state',wall_seconds=time.perf_counter()-wall)
    save(ROOT/('source-9d' if loss else 'source-5d'),out)
    print({'model':out['model'],'periods':last.tolist(),'ranges':np.ptp(states,axis=0).tolist(),'pressure_range':[float(states[:,0].min()),float(states[:,0].max())],'wall_seconds':out['wall_seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--control',action='store_true');parser.add_argument('--loss',action='store_true');a=parser.parse_args()
    if a.control:control(a.loss)
    else:validate()
