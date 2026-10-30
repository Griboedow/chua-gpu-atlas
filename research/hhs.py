"""Exact autonomous unsaturated single-cell S2, Mondal et al. arXiv2504.19655v3.

State M,E,D,y,x; parameters mu,tauD,tauL, with tauR=1,a=c=1,b=3,D0=.1.
No continuum cubic saturation or spatial forcing is inserted. Time is source model time.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

ROOT=Path(__file__).resolve().parents[1]/'results/hhs'
MREF=.6137372763943942
EREF=3.8626272360560576
DEFAULT=np.array([0.,1.1,1.2])


def field_jac(z,p):
    m,e,d,y,x=z;mu,td,tl=p
    f=np.array([1-4*m+m*m*e,3*m-m*m*e-d*e,(.1-d+mu*d*x)/td,(e-EREF-y)/tl,y-x])
    j=np.array([[-4+2*m*e,m*m,0,0,0],[3-2*m*e,-m*m-d,-e,0,0],
        [0,0,(-1+mu*x)/td,0,mu*d/td],[0,1/tl,0,-1/tl,0],[0,0,0,1,-1]])
    return f,j


def original_field(z,p,alpha):
    m,e,d,l0,l=z;mu,td,tl=p;beta=mu/alpha
    return np.array([1-4*m+m*m*e,3*m-m*m*e-d*e,(.1-d-beta*d*(l-1))/td,
        (1-l0-alpha*(e-EREF))/tl,l0-l])


def save(name,data):
    name=Path(name);name.parent.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(str(name)+'.npz',**arrays)
    meta={k:v for k,v in data.items() if k not in arrays};meta['array_shapes']={k:list(v.shape) for k,v in arrays.items()}
    Path(str(name)+'.json').write_text(json.dumps(meta,indent=2,allow_nan=False),encoding='utf-8')


def validate():
    rng=np.random.default_rng(2026103015);errors=dict(jacobian=0.,original_coordinates=0.,mass_identity=0.,reference_equilibrium=0.)
    for _ in range(128):
        p=np.exp(rng.uniform(-3,2,3));z=np.r_[np.exp(rng.uniform(-3,2,3)),rng.normal(size=2)]
        f,j=field_jac(z,p);ref=np.column_stack([field_jac(z.astype(complex)+1e-30j*np.eye(5)[i],p)[0].imag/1e-30 for i in range(5)])
        errors['jacobian']=max(errors['jacobian'],float(np.max(abs(j-ref)/(1+abs(j)))))
        alpha=.2;q=np.r_[z[:3],1-alpha*z[3:]];actual=original_field(q,p,alpha);actual[3:]/=-alpha
        errors['original_coordinates']=max(errors['original_coordinates'],float(np.max(abs(actual-f)/(1+abs(f)))))
        expected=1-z[0]-z[2]*z[1]
        errors['mass_identity']=max(errors['mass_identity'],float(abs(f[0]+f[1]-expected)/(1+abs(expected))))
        errors['reference_equilibrium']=max(errors['reference_equilibrium'],float(np.max(abs(field_jac([MREF,EREF,.1,0,0],p)[0]))))
    p=np.array([.7,1.1,1.2]);z=np.array([.4,4.,.1,0.,0.]);alpha=.2
    a=solve_ivp(lambda t,z:field_jac(z,p)[0],(0,30),z,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.02)
    b=solve_ivp(lambda t,q:original_field(q,p,alpha),(0,30),np.r_[z[:3],1-alpha*z[3:]],method='Radau',rtol=1e-11,atol=1e-13,max_step=.02)
    converted=np.r_[b.y[:3,-1],(1-b.y[3:,-1])/alpha]
    errors['two_solvers_two_coordinate_systems']=float(np.max(abs(a.y[:,-1]-converted)))
    passed=bool(a.success and b.success and max(errors.values())<1e-8)
    save(ROOT/'validation',dict(passed=passed,errors=errors,random_states=128,time_unit='source model time, not calibrated seconds'))
    print(dict(passed=passed,errors=errors),flush=True)
    if not passed:raise RuntimeError('CPU validation failed')


def source_cycle():
    p=DEFAULT.copy();initial=np.array([.4,4.,.1,0.,0.])
    def section(t,z):return z[0]-MREF
    section.direction=1
    warm=solve_ivp(lambda t,z:field_jac(z,p)[0],(0,1000),initial,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.05,events=section)
    if not warm.success or len(warm.t_events[0])<3:raise RuntimeError('Uncoupled periodic control not found')
    y0=warm.y_events[0][-1];period=warm.t_events[0][-1]-warm.t_events[0][-2]
    def shooting(q,augmented=False):
        start=np.r_[MREF,q[:4]];T=q[4]
        if not (1<T<100):return np.ones(5)*1e3
        if not augmented:
            sol=solve_ivp(lambda t,z:field_jac(z,p)[0],(0,T),start,method='DOP853',rtol=2e-12,atol=2e-14,max_step=.02)
            if not sol.success:raise RuntimeError(sol.message)
            return sol.y[:,-1]-start
        def rhs(t,z):
            f,j=field_jac(z[:5],p)
            return np.r_[f,(j@z[5:30].reshape(5,5)).ravel(),np.trace(j)]
        return solve_ivp(rhs,(0,T),np.r_[start,np.eye(5).ravel(),0.],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.02,dense_output=True)
    solved=root(shooting,np.r_[y0[1:],period],tol=1e-10)
    residual=float(np.max(abs(shooting(solved.x))))
    if residual>1e-9:raise RuntimeError('Cycle shooting failed')
    y0=np.r_[MREF,solved.x[:4]];period=solved.x[-1];sol=shooting(solved.x,True);monodromy=sol.y[5:30,-1].reshape(5,5)
    mult=np.linalg.eigvals(monodromy);exponents=np.sort(np.log(abs(mult))/period)[::-1]
    times=np.linspace(0,period,4001);states=sol.sol(times)[:5].T
    chemical_trace=-4+2*states[:,0]*states[:,1]-states[:,0]**2-.1
    radial=float(np.trapezoid(chemical_trace,times)/period)
    expected=np.sort([0.,radial,-1/p[1],-1/p[2],-1.])[::-1]
    alpha=.2;orig=np.r_[y0[:3],1-alpha*y0[3:]]
    other=solve_ivp(lambda t,q:original_field(q,p,alpha),(0,period),orig,method='Radau',rtol=2e-12,atol=2e-14,max_step=.02)
    closure=float(np.max(abs(other.y[:,-1]-orig)))
    errors=dict(shooting=residual,independent_original_radau_closure=closure,
        floquet_vs_triangular=float(np.max(abs(exponents-expected))),
        flow_identity=float(np.max(abs(monodromy@field_jac(y0,p)[0]-field_jac(y0,p)[0]))),
        spectrum_sum_vs_divergence=float(abs(exponents.sum()-sol.y[-1,-1]/period)))
    passed=bool(other.success and closure<1e-9 and errors['floquet_vs_triangular']<1e-4)
    save(ROOT/'uncoupled-cycle',dict(passed=passed,parameters=p,initial_state=y0,period=float(period),times=times,states=states,
        multipliers=mult,exponents=exponents,expected_triangular_exponents=expected,monodromy=monodromy,errors=errors,
        scope='Exact mu=0 five-state control: planar chemical cycle and three stable driven-filter modes; initial phase chosen here.',
        time_unit='source model time, not calibrated seconds',alpha_for_original_coordinate_control=alpha))
    print(dict(passed=passed,period=float(period),exponents=exponents.tolist(),expected=expected.tolist(),errors=errors),flush=True)
    if not passed:raise RuntimeError('Uncoupled cycle validation failed')


if __name__=='__main__':
    validate();source_cycle()
