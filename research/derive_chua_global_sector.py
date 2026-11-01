"""Exact symbolic storage identities and numerical/interval checks of a sector.

The accompanying report supplies the sign, properness and invariance argument.
SymPy checks algebraic identities; it is not a formal proof of LaSalle's theorem.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from flint import arb, arb_mat, ctx
from chua_interval_map import bounds as display_bounds, col
from certify_chua_hidden_cycle import positive_definite
from trace_chua_window import Chua

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/chua-hidden-boundary'


def bounds(x):
    result=display_bounds(x)
    # The separately rounded midpoint/radius strings are display fields, not a
    # lossless serialization. Arb's joint ball string includes decimal rounding.
    result['enclosure']=x.str(ctx.dps+10)
    return result


def controllability(A, B):
    return sp.Matrix.hstack(B, A*B, A*A*B)


def symbolic():
    alpha, r, s = sp.symbols('alpha r s', positive=True)
    beta = r*(alpha+1-r)
    gap = 2*r-alpha-1
    omega2 = (r-1)*(alpha+1-r)
    mstar = (r-alpha-1)/alpha
    q = (2*r-1)/((r-1)*gap)
    p = alpha**2*r/((r-1)*gap)
    D = alpha*q
    K = alpha*(1+q*(alpha+1-r))/r
    A = sp.Matrix([[-alpha*(1+mstar), alpha, 0], [1, -1, 1], [0, -beta, 0]])
    B = sp.Matrix([alpha, 0, 0]); C = sp.Matrix([[1, 0, 0]])
    Ac = sp.Matrix([[-r, 0, 0], [0, 0, 1], [0, -omega2, 0]])
    Bc = sp.Matrix([1, 0, 1]); Cc = sp.Matrix([[-p, 0, K]])
    T = (controllability(Ac, Bc)*controllability(A, B).inv()).applyfunc(sp.factor)
    zero = lambda M: all(sp.cancel(v) == 0 for v in M)
    transfer = alpha*(1+q*s)*(s*s+s+beta)/((s+r)*(s*s+omega2))
    identities = dict(characteristic=sp.factor((s*sp.eye(3)-A).det()-(s+r)*(s*s+omega2)) == 0,
        similarity=zero(T*A-Ac*T), input=zero(T*B-Bc), output=zero(C+q*C*A-Cc*T),
        feedthrough=sp.cancel((q*C*B)[0]-D) == 0,
        positive_real_decomposition=sp.cancel(transfer-(D-p/(s+r)+K*s/(s*s+omega2))) == 0,
        strict_dissipation_margin=sp.cancel(D*r-p-alpha*r/(r-1)) == 0,
        physical_controllability=sp.cancel(controllability(A, B).det()+alpha**3*beta) == 0)
    z, v, phi, x, xdot = sp.symbols('z v phi x xdot', real=True)
    edot = p*z*(-r*z-phi)-K*v*phi
    # x + q*xdot = -p*z + K*v - D*phi.
    lhs = (edot+q*phi*xdot).subs(xdot, (-p*z+K*v-D*phi-x)/q)
    rhs = -x*phi-D*phi**2-2*p*z*phi-r*p*z*z
    identities['storage_derivative'] = sp.cancel(lhs-rhs) == 0
    assert all(identities.values()), identities
    assert sp.cancel(transfer-(-D-p/(s+r)+K*s/(s*s+omega2))) != 0
    return dict(passed=True, identities=identities, sympy_version=sp.__version__,
        wrong_feedthrough_rejected=True,
        parameterization='beta=r*(alpha+1-r), r>1, alpha+1-r>0, 2*r-alpha-1>0, alpha>0',
        coordinate_matrix=[[str(T[i,j]) for j in range(3)] for i in range(3)],
        scope='Exact rational-function identities. Positivity assumptions and the global invariance argument are stated in the report.')


def numerical(digits=80):
    ctx.dps = digits
    a, b = arb(42)/5, arb(12)
    gap = ((a+1)**2-4*b).sqrt(); r = (a+1+gap)/2
    mstar = (r-a-1)/a; omega2 = b+a*mstar
    q = (2*r-1)/((r-1)*gap); p = a*a*r/((r-1)*gap)
    D = a*q; K = a*(1+q*(b-omega2))/r
    A = arb_mat([[-a*(1+mstar), a, 0], [1,-1,1], [0,-b,0]])
    B = col([a,0,0]); Ac = arb_mat([[-r,0,0],[0,0,1],[0,-omega2,0]])
    Bc = col([1,0,1]); C = arb_mat([[1,0,0]]); Cc = arb_mat([[-p,0,K]])
    def ctr(M, v):
        columns = [v, M*v, M*M*v]
        return arb_mat([[columns[j][i,0] for j in range(3)] for i in range(3)])
    T = ctr(Ac,Bc)*ctr(A,B).inv()
    P = T.transpose()*arb_mat([[p,0,0],[0,K*omega2,0],[0,0,K]])*T
    minors = positive_definite(P)
    identities = [T*A-Ac*T,T*B-Bc,C+q*C*A-Cc*T]
    assert all(v.contains(0) for M in identities for v in M.entries())
    assert all(v>0 for v in [q,p,D,K,r,omega2,D*r-p])
    number = lambda x: float(x.mid())
    Tn = np.array([[number(T[i,j]) for j in range(3)] for i in range(3)])
    Pn = np.array([[number(P[i,j]) for j in range(3)] for i in range(3)])
    rng = np.random.default_rng(20261101)
    errors, derivatives = [], []
    for m0,m1 in [(-.18,-.18),(-.18,-.05),(-.05,1.),(1.,-.18),(2.,2.)]:
        model = Chua(alpha=8.4,m0=m0,m1=m1)
        for state in rng.normal(size=(200,3))*10**rng.uniform(-4,3,size=(200,1)):
            x = state[0]; phi = (m1-number(mstar))*x+(m0-m1)*np.clip(x,-1,1)
            field = model.field(12,state); z = (Tn@state)[0]
            direct = (Pn@state)@field+number(q)*phi*field[0]
            identity = -x*phi-number(D)*phi*phi-2*number(p)*z*phi-number(r*p)*z*z
            errors.append(abs(direct-identity)/max(1.,abs(identity)))
            derivatives.append(direct)
    assert max(errors)<1e-9 and max(derivatives)<0
    return dict(passed=True,working_digits=digits,alpha='42/5',beta='12',gamma='0',
        constants={name:bounds(value) for name,value in dict(mstar=mstar,r=r,omega2=omega2,q=q,p=p,D=D,K=K).items()},
        coordinate_matrix=[[bounds(T[i,j]) for j in range(3)] for i in range(3)],
        storage_matrix=[[bounds(P[i,j]) for j in range(3)] for i in range(3)],
        positive_storage_minors=minors,strict_dissipation_margin=bounds(D*r-p),
        numerical_identity_controls=dict(samples=len(errors),maximum_relative_error=max(errors),maximum_derivative=max(derivatives)),
        conclusion='If m0>mstar and m1>=mstar, every trajectory converges to the unique origin. '
                   'At m0=mstar and m1>=mstar the central linear region contains nonzero periodic orbits, so global asymptotic stability fails.',
        scope='The theorem is derived analytically in the report. Interval numbers bound coefficients at the stated alpha,beta, not a finite sampling substitute for the global argument.')


if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    for name,data in [('global-sector-algebra.json',symbolic()),
                      ('global-sector-80.json',numerical(80)),('global-sector-100.json',numerical(100))]:
        (OUT/name).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(name,'passed',flush=True)
