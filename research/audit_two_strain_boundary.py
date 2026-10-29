"""Check the single-strain equilibrium factorization used in the boundary note.

The no-homoclinic argument itself is analytic; these finite checks only audit
the algebra and parameter substitutions, not a theorem by sampling.
"""
import json
import numpy as np
from scipy.optimize import linear_sum_assignment
from two_strain import ROOT,parameters,rhs_jac


def boundary_spectrum(p):
    mu,b1,b2,g1,g2,g12,sigma,eta=p
    a=b1/(g1+mu);b=b2/(g2+mu)
    point=np.array([1/b,0,mu*(b-1)/b2,0,0])
    expected=np.r_[(g1+mu)*(a/b-1),-mu*(1+sigma*(b-1)),-mu-g12,
        np.roots([1,mu*b,mu*(g2+mu)*(b-1)])]
    actual=np.linalg.eigvals(rhs_jac(point,p)[1])
    left,right=linear_sum_assignment(abs(expected[:,None]-actual[None,:]))
    return point,expected,float(np.max(abs(expected[left]-actual[right])))


if __name__=='__main__':
    rng=np.random.default_rng(2026102972);errors=[]
    for _ in range(128):
        p=parameters(rng.uniform(1.05,2.4),rng.uniform(.7,2),r1=2.5,mu=10**rng.uniform(-5,-2))
        point,eig,error=boundary_spectrum(p);errors.append(error)
    p=parameters(1.87125);point,eig,error=boundary_spectrum(p)
    result=dict(random_checks=128,max_eigenvalue_matching_error=max(errors),reference_parameters=p.tolist(),
        reference_equilibrium=point.tolist(),reference_eigenvalues=[[float(z.real),float(z.imag)] for z in eig],
        analytic_statement='For R1>R2>1, any physical trajectory with I1(t0)>0 cannot converge to the strain-2-only equilibrium as t tends to positive infinity: d(log I1)/dt tends to (gamma1+mu)*(R1/R2-1)>0. This excludes a physical interior homoclinic return to that boundary equilibrium.',
        scope='Does not exclude long near-boundary excursions, an asymptotic singular-limit mechanism, interior homoclinic orbits, or chaos. The source says Shilnikov-like, not a proved physical homoclinic orbit.',
        passed=bool(max(errors)<1e-10))
    ROOT.mkdir(parents=True,exist_ok=True);(ROOT/'boundary-audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
