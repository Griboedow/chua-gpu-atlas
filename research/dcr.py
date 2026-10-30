"""Independent literal FP64 residual for the published DCR author-code variant.

No contact capture/sliding convention is added: both predicates remain strict.
This permits short source-solver controls, not an event-resolved chaos claim.
Source: Warmington et al. 2026, DOI10.1098/rsif.2025.0438;
pinned author Julia commit796c5e39468b06fd6fc2c378437d618593948e98, CC0.
"""
from dataclasses import dataclass,asdict
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]/'results'/'dcr'
TWOPI=2*np.pi


@dataclass(frozen=True)
class Parameters:
    alpha:float=1.
    beta:float=1.
    c:float=.14
    lr:float=1.204
    lt:float=1.092
    d:float=.02
    n:int=100


def force(w):
    return 1-w/(4*np.pi)*(1+np.tanh(40*w))


def force_derivative(w):
    z=40*w;t=np.tanh(z)
    return -(1+t+z*(1-t*t))/(4*np.pi)


def geometry(y,p=Parameters(),literal=False):
    x=y[0];theta=y[1:];cc=np.abs(np.cos(theta));ss=np.sin(theta)
    wrapped=np.mod(theta,TWOPI);half=np.arcsin(p.c)
    gate=(wrapped>np.pi-half)&(wrapped<np.pi+half)
    rotor=np.arange(p.n)*p.lr-ss;nt=int(np.ceil(p.n*p.lr/p.lt))
    if literal:
        teeth=np.mod(np.arange(nt)*p.lt+x,nt*p.lt)
        distances=np.min(abs(rotor[:,None]-teeth[None,:]),axis=1)
    else:
        delta=np.mod(x,p.lt)
        index=np.clip(np.rint((rotor-delta)/p.lt),0,nt-1)
        distances=abs(rotor-(delta+index*p.lt))
    return gate&(distances<p.d),gate,cc,distances


def residual(y,yp,p=Parameters(),literal=False):
    contact,gate,cc,_=geometry(y,p,literal=literal)
    res=np.empty(p.n+1)
    res[0]=np.dot(force(yp[1:])*cc,gate)-p.alpha*yp[0]-p.beta*y[0]
    res[1:]=TWOPI-yp[1:]
    res[1:][contact]=yp[0]/cc[contact]-yp[1:][contact]
    return res


def residual_jacobian(y,yp,cj,p=Parameters()):
    contact,gate,cc,_=geometry(y,p)
    dc=-np.sin(y[1:])*np.sign(np.cos(y[1:]))
    j=np.zeros((p.n+1,p.n+1));j[0,0]=-p.beta-cj*p.alpha
    j[0,1:]=gate*(dc*force(yp[1:])+cj*cc*force_derivative(yp[1:]))
    j[np.arange(1,p.n+1),np.arange(1,p.n+1)]=-cj
    ci=np.flatnonzero(contact)+1
    j[ci,ci]-=yp[0]*dc[contact]/cc[contact]**2
    j[ci,0]=cj/cc[contact]
    return j


def field(y,p=Parameters()):
    contact,gate,cc,_=geometry(y,p)
    cc=cc[contact]
    if not len(cc):v=-p.beta*y[0]/p.alpha
    else:
        # Fs(2*pi) evaluates exactly to zero in the source's FP64 arithmetic.
        def balance(v):return p.alpha*v+p.beta*y[0]-np.dot(cc,force(v/cc))
        bound=(p.beta*abs(y[0])+len(cc)+1)/p.alpha+TWOPI
        v=brentq(balance,-bound,bound,xtol=5e-15,rtol=1e-14)
    out=np.full(p.n+1,TWOPI);out[0]=v;out[1:][contact]=v/cc
    return out


def save(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    meta={k:v for k,v in data.items() if k not in arrays}
    if arrays:np.savez_compressed(path.with_suffix('.npz'),**arrays)
    meta['arrays']={k:dict(shape=list(v.shape),dtype=str(v.dtype)) for k,v in arrays.items()}
    path.with_suffix('.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def algebra_controls():
    rng=np.random.default_rng(2026103020);worst=0.;distance_error=0.;root_error=0.;jac_error=0.
    records=[]
    for i in range(128):
        p=Parameters(alpha=1. if i%2 else 4.,beta=1. if i%3 else 4.)
        y=np.r_[rng.uniform(0,10),rng.uniform(-10,20,p.n)];yp=rng.normal(size=p.n+1)*3
        a=residual(y,yp,p);b=residual(y,yp,p,literal=True)
        worst=max(worst,float(np.max(abs(a-b))))
        c,g,cc,d=geometry(y,p);c2,g2,cc2,d2=geometry(y,p,literal=True)
        if not np.array_equal(c,c2):raise AssertionError('Contact identity failed away from interfaces')
        distance_error=max(distance_error,float(np.max(abs(d-d2))))
        root_error=max(root_error,float(np.max(abs(residual(y,field(y,p),p)))))
        # Fixed-branch central directional derivative, avoiding gate/distance interfaces.
        direction=rng.normal(size=p.n+1);cj=1.7;h=1e-6
        analytic=residual_jacobian(y,yp,cj,p)@direction
        plus=residual(y+h*direction,yp+h*cj*direction,p)
        minus=residual(y-h*direction,yp-h*cj*direction,p)
        if np.array_equal(geometry(y+h*direction,p)[0],c) and np.array_equal(geometry(y-h*direction,p)[0],c):
            jac_error=max(jac_error,float(np.max(abs((plus-minus)/(2*h)-analytic)/(1+abs(analytic)))))
    passed=worst<1e-12 and distance_error<1e-12 and root_error<1e-10 and jac_error<1e-7
    if not passed:raise AssertionError((worst,distance_error,root_error,jac_error))
    data=dict(passed=passed,trials=128,seed=2026103020,optimized_vs_literal_residual_max_error=worst,
              nearest_distance_max_error=distance_error,consistent_initial_derivative_residual=root_error,
              fixed_branch_jacobian_directional_relative_error=jac_error,
              scope='Random off-interface residual/geometry/root/Jacobian identities. Strict-boundary predicates and event continuation require a separate audit.',
              nominal_parameters=asdict(Parameters()),force_free_rate_fp64=float(force(TWOPI)))
    save(ROOT/'algebra-controls',data);print(data,flush=True)


if __name__=='__main__':algebra_controls()
