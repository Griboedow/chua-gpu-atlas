"""Shooting continuation of the regular orbit in the explicit F_eta homotopy.

This changes the interaction law; eta is not a fitted device parameter.
Phase differences remove common phase. Fixed integer windings keep the same
torus orbit branch. CPU adaptive integration is appropriate for this small
dependent sequence; bulk basin searches use the existing OpenCL solver.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from synchronematic import COEFFICIENTS, flow_and_jac, geometry, gradient_coefficients


class Orbit:
    def __init__(self, positions, rtol=2e-11, max_step=.05, method='DOP853'):
        self.pairs=geometry(positions); self.rtol=rtol; self.max_step=max_step; self.method=method
        self.E=np.zeros((8,7)); self.E[1:4,:3]=np.eye(3); self.E[4:,3:]=np.eye(4)
        self.D=self.E.T.copy(); self.D[:3,0]=-1
        self.gradient=gradient_coefficients(); self.residual=COEFFICIENTS-self.gradient

    def field(self,q,eta):
        f,j=flow_and_jac(self.E@q,self.pairs,self.gradient+eta*self.residual)
        return self.D@f,self.D@j@self.E

    def integrate(self,q,period,eta,sensitivity=False):
        if not 0 < period < 100: raise ValueError('Period outside continuation guard')
        def rhs(t,z):
            f,j=self.field(z[:7],eta)
            result=np.r_[f,(j@z[7:56].reshape(7,7)).ravel()]
            if sensitivity:
                feta=self.D@flow_and_jac(self.E@z[:7],self.pairs,self.residual)[0]
                result=np.r_[result,j@z[56:63]+feta]
            return result
        start=np.r_[q,np.eye(7).ravel(),np.zeros(7) if sensitivity else []]
        sol=solve_ivp(rhs,(0,period),start,method=self.method,rtol=self.rtol,
                      atol=self.rtol*.01,max_step=self.max_step)
        if not sol.success: raise RuntimeError(sol.message)
        z=sol.y[:,-1]
        return z[:7],z[7:56].reshape(7,7),z[56:63]

    def shoot(self,guess,eta,reference,winding,tol=2e-9):
        u=guess.copy(); normal=self.field(reference,eta)[0]; normal/=np.linalg.norm(normal)
        for iteration in range(18):
            end,mono,_=self.integrate(u[:7],u[7],eta)
            r=np.r_[end-u[:7]-2*np.pi*winding,normal@(u[:7]-reference)]
            error=float(np.linalg.norm(r,np.inf))
            if error<tol:return u,mono,error,iteration+1
            jac=np.zeros((8,8));jac[:7,:7]=mono-np.eye(7)
            jac[:7,7]=self.field(end,eta)[0];jac[7,:7]=normal
            update=np.linalg.solve(jac,-r)
            # Bound an individual correction; failures are saved, not relabeled.
            factor=min(1.,.5/max(np.linalg.norm(update[:7]),1e-15),.5/max(abs(update[7]),1e-15))
            u+=factor*update
        raise RuntimeError(f'Shooting did not converge: residual={error:g}')


def record(model,u,eta,mono,error,iterations,winding):
    eig=np.linalg.eigvals(mono); neutral=int(np.argmin(abs(eig-1)))
    transverse=np.delete(eig,neutral); transverse=transverse[np.argsort(-abs(transverse))]
    return dict(eta=float(eta),period=float(u[7]),state=u[:7].tolist(),
                winding=winding.tolist(),closure_residual=error,iterations=iterations,
                neutral_multiplier=dict(real=float(eig[neutral].real),imag=float(eig[neutral].imag)),
                transverse_multipliers=[dict(real=float(z.real),imag=float(z.imag),absolute=float(abs(z))) for z in transverse],
                spectral_radius=float(max(abs(transverse))))


def main():
    p=argparse.ArgumentParser();p.add_argument('--step',type=float,default=.00025)
    p.add_argument('--lower',type=float,default=1.068);p.add_argument('--upper',type=float,default=1.084)
    p.add_argument('--output',type=Path,default=Path('results/synchronematic/n4-periodic-branch'))
    args=p.parse_args();p0=Path('results/synchronematic/n4-regular-window-cpu')
    data=np.load(p0/'trace.npz');rec=json.loads((p0/'recurrence.json').read_text())
    model=Orbit(data['positions']);q=model.D@data['state'][-1]
    q=np.angle(np.exp(1j*q));T=rec['returns'][0]['time'];eta0=1.074
    end,_,_=model.integrate(q,T,eta0);winding=np.rint((end-q)/(2*np.pi)).astype(int)
    seed,mono,error,it=model.shoot(np.r_[q,T],eta0,q,winding)
    records=[record(model,seed,eta0,mono,error,it,winding)];failures=[];wall=time.perf_counter()
    args.output.mkdir(parents=True,exist_ok=True)
    def save():
        ordered=sorted(records,key=lambda x:x['eta'])
        jumps=[dict(lower_eta=a['eta'],upper_eta=b['eta'],period_jump=abs(a['period']-b['period']))
               for a,b in zip(ordered,ordered[1:]) if abs(a['period']-b['period'])>.1]
        result=dict(model='F_eta=G+eta*(F-G); explicit modification of the printed Fourier model',
                    changes_equations=True,positions=data['positions'].tolist(),
                    coordinates='theta[1:]-theta[0], chi[:]; common phase removed',
                    method='Fixed-parameter Newton shooting with analytic monodromy and phase condition; DOP853',
                    rtol=model.rtol,atol=model.rtol*.01,max_step=model.max_step,
                    caveat='Numerical orbit solutions, not validated existence or global basin proof. Fixed-parameter Newton can jump branches at folds; flagged period jumps must not be connected. Tiny multipliers may be unresolved.',
                    large_adjacent_period_jumps=jumps,
                    records=ordered,failures=failures,wall_seconds=time.perf_counter()-wall)
        (args.output/'branch.json').write_text(json.dumps(result,indent=2))
    for sign,limit in [(-1,args.lower),(1,args.upper)]:
        u=seed.copy();eta=eta0;step=args.step
        while sign*(limit-eta)>1e-10:
            trial=eta+sign*min(step,abs(limit-eta))
            try:
                new,mono,error,it=model.shoot(u,trial,u[:7],winding)
            except (RuntimeError,ValueError,np.linalg.LinAlgError) as e:
                failures.append(dict(eta=trial,step=step,message=str(e)));save()
                step*=.5
                if step<1e-7:break
                continue
            u=new;eta=trial
            row=record(model,u,eta,mono,error,it,winding);records.append(row);save()
            print(f'eta={eta:.9f} T={u[7]:.8f} radius={row["spectral_radius"]:.8f} closure={error:.2g}',flush=True)
            step=min(args.step,step*1.5)
    save();print(f'Saved {len(records)} orbit points in {time.perf_counter()-wall:.1f}s',flush=True)


if __name__=='__main__':main()
