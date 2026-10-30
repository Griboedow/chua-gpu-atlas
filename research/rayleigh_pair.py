"""Exact canonical 4D quadratic-drag thermal spring pair, Li et al. 2023 Eq.14.

Canonical time sigma and translated/scaled displacements; not physical seconds.
The source has already eliminated the thermal field. No extra feedback is added.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]/'results/rayleigh-pair'
SOURCE=np.array([.00880947943,.557378006])


def field_jac(y,p):
    mu,k=p; x,v,z,w=y
    f=np.array([v,-(1+k)*x-k*z+mu*(1-abs(v))*v,w,-k*x-(1+k)*z+mu*(1-abs(w))*w])
    j=np.array([[0,1,0,0],[-1-k,mu*(1-2*abs(v)),-k,0],[0,0,0,1],[-k,0,-1-k,mu*(1-2*abs(w))]])
    return f,j


def modal_energy(y,p):
    x,v,z,w=np.asarray(y).T;s=(x+z)/2;sv=(v+w)/2;a=(x-z)/2;av=(v-w)/2
    return np.stack([sv*sv+(1+2*p[1])*s*s,av*av+a*a],axis=-1)


def save(name,data):
    name=Path(name);name.parent.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(str(name)+'.npz',**arrays)
    meta={k:v for k,v in data.items() if k not in arrays};meta['array_shapes']={k:list(v.shape) for k,v in arrays.items()}
    Path(str(name)+'.json').write_text(json.dumps(meta,indent=2,allow_nan=False),encoding='utf-8')


def validate():
    rng=np.random.default_rng(2026103014);errors=dict(jacobian=0.,energy=0.,modal_restriction=0.,vdp_conjugacy=0.)
    for _ in range(256):
        p=np.exp(rng.uniform(-3,3,2));mu,k=p;y=rng.normal(size=4)
        f,j=field_jac(y,p);h=1e-6
        fd=np.column_stack([(field_jac(y+h*np.eye(4)[i],p)[0]-field_jac(y-h*np.eye(4)[i],p)[0])/(2*h) for i in range(4)])
        errors['jacobian']=max(errors['jacobian'],float(np.max(abs(j-fd)/(1+abs(j)))))
        x,v,z,w=y;grad=np.array([(1+k)*x+k*z,v,k*x+(1+k)*z,w]);ref=mu*(v*v+w*w-abs(v)**3-abs(w)**3)
        errors['energy']=max(errors['energy'],float(abs(grad@f-ref)/(1+abs(ref))))
        for sign in [1,-1]:
            q=np.array([x,v,sign*x,sign*v]);fq=field_jac(q,p)[0]
            errors['modal_restriction']=max(errors['modal_restriction'],float(np.max(abs(fq[2:]-sign*fq[:2]))))
        # Independent reflected nonsmooth van der Pol coordinates u=2*z'.
        pos=np.array([x,-z]);vel=np.array([v,-w]);M=np.array([[1+k,-k],[-k,1+k]])
        u=2*vel;b=2*(-M@pos+mu*(1-abs(vel))*vel)
        vdp_acc=-M@u+mu*(1-abs(u))*b
        ref_acc=2*np.array([j[1]@f,-j[3]@f])
        errors['vdp_conjugacy']=max(errors['vdp_conjugacy'],float(np.max(abs(vdp_acc-ref_acc)/(1+abs(ref_acc)))))
    y=np.array([.1,.7,-.3,.2]);p=np.array([3.,1.5])
    a=solve_ivp(lambda t,y:field_jac(y,p)[0],(0,10),y,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.01)
    b=solve_ivp(lambda t,y:field_jac(y,p)[0],(0,10),y,method='Radau',jac=lambda t,y:field_jac(y,p)[1],rtol=1e-11,atol=1e-13,max_step=.01)
    errors['two_cpu_methods']=float(np.max(abs(a.y[:,-1]-b.y[:,-1])))
    passed=bool(a.success and b.success and max(errors.values())<1e-7)
    save(ROOT/'validation',dict(passed=passed,errors=errors,random_states=256,time_unit='canonical sigma',
        note='C1 drag with absolute value; no saltation. Ordinary finite differences away from zero; exact analytic derivative at zero.'))
    print(dict(passed=passed,errors=errors),flush=True)
    if not passed:raise RuntimeError('Canonical model validation failed')


if __name__=='__main__':validate()
