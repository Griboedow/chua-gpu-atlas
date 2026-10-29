"""Numerical orbit-closure check for the independently verified regular window."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from scipy.linalg import block_diag,null_space
from synchronematic import geometry,flow_and_jac


def main():
    path=Path('results/synchronematic/n4-regular-window-cpu')
    data=np.load(path/'trace.npz');meta=json.loads((path/'spectrum.json').read_text())
    state=data['state'][-1];pairs=geometry(data['positions']);c=np.array(meta['coefficients'])
    sol=solve_ivp(lambda t,y:flow_and_jac(y,pairs,c)[0],(0,40),state,
                  method='DOP853',rtol=1e-12,atol=1e-14,max_step=.025,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    def residual(t):
        difference=sol.sol(float(np.asarray(t).ravel()[0]))-state
        difference[:4]-=difference[0]
        return np.angle(np.exp(1j*difference))
    def derivative(t):
        f=flow_and_jac(sol.sol(float(np.asarray(t).ravel()[0])),pairs,c)[0]
        f[:4]-=f[0]
        return f[:,None]
    times=np.linspace(.5,40,1581);distances=np.array([np.linalg.norm(residual(t)) for t in times])
    minima=find_peaks(-distances)[0];returns=[]
    for i in minima:
        result=least_squares(residual,[times[i]],jac=derivative,bounds=([times[i-1]],[times[i+1]]),
                             xtol=1e-13,ftol=1e-13,gtol=1e-13)
        error=float(np.linalg.norm(result.fun))
        if error<1e-5:returns.append(dict(time=float(result.x[0]),residual=error))
    output=dict(eta=1.074,source_index=576,method='DOP853 dense trajectory and least-squares circular orbit closure; common phase removed',
                returns=returns,interpretation='Numerical periodic-orbit evidence, not a validated existence or stability proof')
    if returns:
        basis=block_diag(null_space(np.ones((1,4))),np.eye(4))
        def variational(t,y):
            f,j=flow_and_jac(y[:8],pairs,c)
            return np.r_[f,(j@y[8:].reshape(8,7)).ravel()]
        mono=solve_ivp(variational,(0,returns[0]['time']),np.r_[state,basis.ravel()],
                       method='DOP853',rtol=1e-12,atol=1e-14,max_step=.025)
        if not mono.success:raise RuntimeError(mono.message)
        eigenvalues=np.linalg.eigvals(basis.T@mono.y[8:,-1].reshape(8,7))
        eigenvalues=eigenvalues[np.argsort(-abs(eigenvalues))]
        output['floquet_multipliers']=[{'real':float(z.real),'imag':float(z.imag),'absolute':float(abs(z))} for z in eigenvalues]
        output['floquet_caveat']='Direct monodromy integration resolves the dominant multipliers; extremely small contracting multipliers are limited by floating-point precision.'
    (path/'recurrence.json').write_text(json.dumps(output,indent=2));print(json.dumps(output,indent=2))


if __name__=='__main__':main()
