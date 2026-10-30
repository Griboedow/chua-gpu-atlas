"""Affine-flight shooting and symmetry analysis for a Chua atlas layer.

Every affine arc is checked against its region, then selected solutions are
independently integrated with the full vector field. This prevents an algebraic
closing solution with a wrong itinerary from being accepted as a circuit orbit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import root

from trace_chua_window import Chua


class Flights:
    def __init__(self, model, regions, targets, symmetric):
        self.model=model
        self.regions=list(regions)
        self.targets=list(targets)
        self.symmetric=symmetric
        self.sign=-1 if symmetric else 1

    def evaluate(self, beta, q, check=False):
        n=len(q)
        initial=np.array([0.,q[0],q[1]])
        state=initial.copy()
        derivative=np.zeros((3,n)); derivative[1,0]=derivative[2,1]=1.
        initial_derivative=derivative.copy()
        mono=np.eye(3)
        residual=[]; jac=[]; arcs=[]; violation=0.; speeds=[]; logdet=0.
        mats={}
        for region in set(self.regions):
            A=self.model.jacobian(beta,region==0)
            b=np.array([-self.model.alpha*(self.model.m0-self.model.m1)*region,0.,0.])
            B=np.zeros((4,4));B[:3,:3]=A;B[:3,3]=b
            mats[region]=(A,b,B)
        for k,(duration,region) in enumerate(zip(q[2:],self.regions)):
            A,b,B=mats[region]
            E=expm(B*duration)
            old=state.copy()
            state=E[:3,:3]@state+E[:3,3]
            derivative=E[:3,:3]@derivative
            derivative[:,k+2]+=A@state+b
            mono=E[:3,:3]@mono
            logdet+=np.trace(A)*duration
            if check:
                for dt in np.linspace(0.,duration,101):
                    e=expm(B*dt); x=(e[:3,:3]@old+e[:3,3])[0]
                    error=abs(x)-1 if region==0 else 1-region*x
                    violation=max(violation,error)
                arcs.append(dict(region=region,duration=float(duration),start=old.tolist(),end=state.tolist()))
            if k<len(self.targets):
                residual.append(state[0]-self.targets[k]); jac.append(derivative[0].copy())
                speeds.append(abs((A@state+b)[0]))
        residual.extend(state-self.sign*initial)
        jac.extend(derivative-self.sign*initial_derivative)
        f=self.model.field(beta,state)
        projection=np.eye(3)-np.outer(f,[1.,0.,0.])/f[0]
        P=(self.sign*projection@mono)[1:,1:]
        return np.array(residual),np.array(jac),P,dict(arcs=arcs,violation=float(violation),
                 minimum_speed=float(min(speeds)),logdet=float(logdet),mono=mono)

    def solve(self,beta,guess):
        ans=root(lambda q:self.evaluate(beta,q)[0],guess,
                 jac=lambda q:self.evaluate(beta,q)[1],tol=1e-10)
        error=np.max(np.abs(self.evaluate(beta,ans.x)[0]))
        if error>2e-10 or np.min(ans.x[2:])<=0:
            raise RuntimeError(f'flight closure failed beta={beta}: {error} {ans.message}')
        return ans.x

    def boundary(self,beta,guess,target):
        def fun(w):
            r,j,P,aux=self.evaluate(w[-1],w[:-1])
            determinant=self.sign*np.exp(aux['logdet'])
            return np.r_[r,target*target-target*np.trace(P)+determinant]
        ans=root(fun,np.r_[guess,beta],tol=1e-9)
        error=np.max(np.abs(fun(ans.x)))
        if error>1e-8: raise RuntimeError(f'boundary failed {error} {ans.message}')
        return ans.x[-1],ans.x[:-1]

    def record(self,beta,q,check=False):
        r,j,P,aux=self.evaluate(beta,q,check)
        # Direct determinants lose the strongly contracted multiplier after
        # several turns. Liouville's identity retains it without cancellation.
        determinant=self.sign*np.exp(aux['logdet'])
        trace=np.trace(P)
        discriminant=trace*trace-4*determinant
        if discriminant>=0:
            large=(trace+np.copysign(np.sqrt(discriminant),trace))/2
            vals=[large,determinant/large]
        else:
            vals=[(trace+1j*np.sqrt(-discriminant))/2,(trace-1j*np.sqrt(-discriminant))/2]
        vals=sorted(vals,key=abs,reverse=True)
        full=[v*v if self.symmetric else v for v in vals]
        factor=2 if self.symmetric else 1
        return dict(beta=float(beta),alpha=self.model.alpha,m0=self.model.m0,m1=self.model.m1,
                    symmetric=self.symmetric,q=q.tolist(),period=float(factor*sum(q[2:])),
                    residual=float(np.max(np.abs(r))),
                    map_multipliers=[[float(v.real),float(v.imag)] for v in vals],
                    full_multipliers=[[float(v.real),float(v.imag)] for v in full],
                    spectral_radius=float(max(abs(v) for v in full)),
                    section_matrix=P.tolist(),minimum_breakpoint_speed=aux['minimum_speed'],
                    map_determinant_liouville=float(determinant),map_trace=float(trace),
                    full_logdet=factor*aux['logdet'],region_violation=aux['violation'],arcs=aux['arcs'])


def from_orbit(model,beta,v,symmetric):
    duration=v[2]/(2 if symmetric else 1)
    sol=model.integrate(beta,[0.,*v[:2]],duration)
    events=sorted((float(t),target) for times,target in zip(sol.t_events[:2],[-1,1]) for t in times if 1e-8<t<duration-1e-8)
    times=[0.]+[t for t,target in events]+[duration]
    regions=[]
    for l,r in zip(times[:-1],times[1:]):
        x=sol.sol((l+r)/2)[0]
        regions.append(0 if abs(x)<1 else int(np.sign(x)))
    flights=Flights(model,regions,[target for t,target in events],symmetric)
    q=np.r_[v[:2],np.diff(times)]
    return flights,flights.solve(beta,q)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results/chua-layers')
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    model=Chua(rtol=2e-12,step=.025)
    v=json.loads((out/'initial.json').read_text())['shooting_vector']
    half,q0=from_orbit(model,14.32,v,True)
    result=dict(regions=half.regions,targets=half.targets,central=[],boundaries={})
    for target,beta in [(-1,14.3173),(1,14.3305)]:
        b,q=half.boundary(beta,q0,target)
        r=half.record(b,q,True)
        result['boundaries'][str(target)]=r
        print('BOUNDARY',target,json.dumps(r),flush=True)
    for direction,betas in [(1,np.arange(14.32,14.33051,.00025)),(-1,np.arange(14.31975,14.2949,-.00025))]:
        q=q0.copy()
        for beta in betas:
            try:q=half.solve(float(beta),q)
            except RuntimeError as e: print(e,flush=True);break
            r=half.record(float(beta),q)
            result['central'].append(r)
    (out/'affine-central.json').write_text(json.dumps(result,indent=2)+'\n')
    print('saved central',len(result['central']),flush=True)


if __name__=='__main__': main()
