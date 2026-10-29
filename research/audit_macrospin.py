"""Independent axial reduction of the deterministic SIT macrospin model.

Source: Arun et al., arXiv:2601.07300v1, Eqs. 1-3.
Checks algebraic identities, not a parameter fit or noise simulation.
"""
import json
from pathlib import Path
import numpy as np


def main():
    rng=np.random.default_rng(2026102944);errors=[];norm_errors=[];minimum=np.inf
    for _ in range(128):
        m=rng.normal(size=3);m/=np.linalg.norm(m);z=m[2];pole=np.array([0.,0.,1.])
        alpha=rng.uniform(.01,.2);kappa=rng.uniform(0.,2.);beta=rng.uniform(-1.,1.)
        gamma=rng.uniform(.2,2.);hz=rng.uniform(-2.,2.)*z;hs=rng.uniform(.1,2.)/(1+.3*z)
        cross=np.array([[0,-m[2],m[1]],[m[2],0,-m[0]],[-m[1],m[0],0]])
        w=np.cross(m,pole);q=np.cross(m,w)
        # kappa=gamma*G. R=p dot (m cross m_dot)=-w dot m_dot.
        matrix=np.eye(3)-alpha*cross+kappa*np.outer(q,w)
        direct=np.linalg.solve(matrix,gamma*((beta*hs-hz)*w+hs*q))
        c=alpha-kappa*(1-z*z);denominator=1+alpha*c
        a=gamma*((beta-alpha)*hs-hz)/denominator
        b=gamma*((1+beta*c)*hs-c*hz)/denominator
        reduced=a*w+b*q
        errors.append(float(np.max(abs(direct-reduced))));norm_errors.append(float(abs(m@direct)))
        minimum=min(minimum,denominator)
    result=dict(source='https://arxiv.org/html/2601.07300v1',samples=128,
                maximum_cartesian_reduction_error=max(errors),maximum_norm_derivative_half=max(norm_errors),
                minimum_sampled_denominator=float(minimum),
                conclusion='On a regular deterministic axial branch, z=m_z is a scalar autonomous state; the azimuth is slaved. No autonomous strange-attractor search is warranted for this model.',
                scope='Constant current, no thermal noise, single macrospin, axial field and pinned layer. No claim for coupled spins, delayed circuits, forcing, or singular branches.')
    if max(errors)>1e-12 or max(norm_errors)>1e-12:raise RuntimeError('Axial reduction check failed')
    folder=Path('results/macrospin-exclusion');folder.mkdir(parents=True,exist_ok=True)
    (folder/'algebra.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
