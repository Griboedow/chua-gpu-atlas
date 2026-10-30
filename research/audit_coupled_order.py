"""Algebraic order audit for the explicitly smooth two-operator model.

No trajectory integration, GPU use, or parameter search.  The Poincare-Bendixson
conclusion is imported with its equilibrium-free hypothesis, not inferred from
the finite matrix checks below.
"""
from __future__ import annotations

from itertools import combinations
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.sparse.csgraph import connected_components
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/coupled-operators/order-audit.json"


def rates(k, b):
    discriminant = b*b-4*k
    if discriminant < -1e-12*max(b*b, 4*k):
        raise ValueError("The positive real-chain transformation needs b^2 >= 4k")
    r2 = (b+np.sqrt(max(0., discriminant)))/2
    r1 = k/r2
    return r1, r2


def transform(k, b):
    r1, _ = rates(k, b)
    # Physical order: xm,vm,ym,xp,vp,yp.
    # Cyclic order: yp,zp,xp,ym,zm,xm; z=v+r1*x.
    T = np.zeros((6, 6))
    T[0, 5] = 1
    T[1, 3] = r1; T[1, 4] = 1
    T[2, 3] = 1
    T[3, 2] = 1
    T[4, 0] = r1; T[4, 1] = 1
    T[5, 0] = 1
    return T


def physical_field_jac(state, p):
    k, b, dm, dp, width = p
    xm, vm, ym, xp, vp, yp = state
    am = (xp-dp)/width
    ap = (dm-xm)/width
    field = np.array([vm, k*(ym-xm)-b*vm, expit(am)-ym,
                      vp, k*(yp-xp)-b*vp, expit(ap)-yp])
    jac = np.zeros((6, 6))
    for j in [0, 3]:
        jac[j,j+1] = 1
        jac[j+1,j] = -k; jac[j+1,j+1] = -b; jac[j+1,j+2] = k
        jac[j+2,j+2] = -1
    # Positive in exact arithmetic; expit(-a) avoids 1-expit(a) cancellation.
    jac[2,3] = expit(am)*expit(-am)/width
    jac[5,0] = -expit(ap)*expit(-ap)/width
    return field, jac


def chain_field_jac(state, p):
    k, b, dm, dp, width = p
    r1, r2 = rates(k, b)
    yp, zp, xp, ym, zm, xm = state
    am = (xp-dp)/width; ap = (dm-xm)/width
    field = np.array([expit(ap)-yp, k*yp-r2*zp, zp-r1*xp,
                      expit(am)-ym, k*ym-r2*zm, zm-r1*xm])
    jac = np.diag([-1.,-r2,-r1,-1.,-r2,-r1])
    jac[0,5] = -expit(ap)*expit(-ap)/width
    jac[1,0] = k; jac[2,1] = 1
    jac[3,2] = expit(am)*expit(-am)/width
    jac[4,3] = k; jac[5,4] = 1
    return field, jac


def second_compound(A):
    pairs = list(combinations(range(len(A)), 2))
    lookup = {pair:i for i, pair in enumerate(pairs)}
    C = np.zeros((len(pairs), len(pairs)))
    for column, (j, k) in enumerate(pairs):
        for i in range(len(A)):
            if i != k:
                C[lookup[tuple(sorted((i,k)))],column] += (1 if i < k else -1)*A[i,j]
            if i != j:
                C[lookup[tuple(sorted((j,i)))],column] += (1 if j < i else -1)*A[i,k]
    return C, pairs


def equilibrium(p):
    k, b, dm, dp, width = p
    xm = brentq(lambda x:expit((expit((dm-x)/width)-dp)/width)-x,
                0, 1, xtol=5e-15)
    xp = expit((dm-xm)/width)
    y = np.array([xm, 0., xm, xp, 0., xp])
    f, J = physical_field_jac(y, p)
    gain = -J[2,3]*J[5,0]
    determinant = np.linalg.det(-J)
    expected = k*k*(1+gain)
    return {"parameters_k_b_dm_dp_width":p.tolist(), "state":y.tolist(),
            "field_residual":float(np.max(np.abs(f))), "gain_product":float(gain),
            "det_minus_J":float(determinant), "expected_k_squared_times_one_plus_gain":float(expected),
            "relative_determinant_error":float(abs(determinant-expected)/expected),
            "eigenvalues":[[float(z.real),float(z.imag)] for z in np.linalg.eigvals(J)]}


def main():
    started = time.perf_counter()
    rng = np.random.default_rng(2026103001)
    errors = {"factorization":0.,"inverse":0.,"field_conjugacy":0.,"jacobian_similarity":0.,
              "second_compound_exponential":0.,"box_outward_violation":0.}
    minimum_offdiag = np.inf
    max_components = 0
    critical_count = 0
    equilibrium_records = []
    for index in range(128):
        k = 10**rng.uniform(-1, 2)
        b = 2*np.sqrt(k)*(1 if index % 4 == 0 else rng.uniform(1.05,4.))
        critical_count += int(index % 4 == 0)
        p = np.array([k,b,*rng.uniform(.1,.9,2),rng.uniform(.08,.4)])
        r1,r2 = rates(k,b)
        errors["factorization"] = max(errors["factorization"],abs(r1+r2-b)/b,abs(r1*r2-k)/k)
        T = transform(k,b); inverse = np.linalg.inv(T)
        errors["inverse"] = max(errors["inverse"],float(np.max(np.abs(T@inverse-np.eye(6)))))
        y = rng.uniform(-.2,1.2,6)
        f,J = physical_field_jac(y,p); g,K = chain_field_jac(T@y,p)
        errors["field_conjugacy"] = max(errors["field_conjugacy"],float(np.max(np.abs(g-T@f)/(1+np.abs(g)))))
        errors["jacobian_similarity"] = max(errors["jacobian_similarity"],float(np.max(np.abs(K-T@J@inverse)/(1+np.abs(K)))))
        compound, pairs = second_compound(K)
        offdiag = compound.copy(); np.fill_diagonal(offdiag,0)
        minimum_offdiag = min(minimum_offdiag,float(offdiag.min()))
        max_components = max(max_components,connected_components(offdiag>0,directed=True,connection="strong",return_labels=False))
        if index < 4:
            dt = .01
            transition = expm(K*dt)
            minors = np.array([[np.linalg.det(transition[np.ix_(i,j)]) for j in pairs] for i in pairs])
            errors["second_compound_exponential"] = max(errors["second_compound_exponential"],float(np.max(np.abs(minors-expm(compound*dt)))))
            equilibrium_records.append(equilibrium(p))
        # Exact invariant box 0<=yp,ym,xp,xm<=1; 0<=zp,zm<=r1.
        upper = np.array([1.,r1,1.,1.,r1,1.])
        q = rng.uniform(0,1,6)*upper
        for coordinate in range(6):
            low = q.copy(); low[coordinate]=0
            high = q.copy(); high[coordinate]=upper[coordinate]
            errors["box_outward_violation"] = max(errors["box_outward_violation"],
                -chain_field_jac(low,p)[0][coordinate],chain_field_jac(high,p)[0][coordinate])
    references=[]
    source=ROOT/"results/light-operators/coupled-scout.npz"
    if source.exists():
        data=np.load(source)
        for row in [1026,3054]:
            k,b,dm,dp=map(float,data['parameters'][row])
            references.append({"source_row":row,"parameters":[k,b,dm,dp],
                               "b_squared_over_4k":b*b/(4*k),"critical_b":2*np.sqrt(k),
                               "order_audit_applies":bool(b*b>=4*k)})
    assert max(errors.values()) < 1e-11 and minimum_offdiag >= 0 and max_components == 1
    assert all(q['relative_determinant_error'] < 1e-11 for q in equilibrium_records)
    report={
        "scope":"Explicit logistic optical profile of positive width; continuous heat balance; identical positive k,b. No theorem applied to the ideal relay.",
        "assumption":"b^2>=4k, including equality",
        "transformation":"r1+r2=b, r1*r2=k; z=v+r1*x; cyclic state=(yp,zp,xp,ym,zm,xm)",
        "cyclic_link_signs":[-1,1,1,1,1,1],"sign_product":-1,
        "invariant_box":"0<=yp,ym,xp,xm<=1; 0<=zp,zm<=r1",
        "verified_theorem_conclusion":"Every compact omega-limit set containing no equilibrium is a single periodic orbit (Weiss-Margaliot Theorem 12).",
        "equilibrium":"Unique for every positive logistic width; det(-J)=k^2*(1+a_minus*a_plus)>0",
        "stronger_extension_status":"Hori-Hara-Kim Proposition 1 states the original cyclic theorem's positive-determinant exclusion of homoclinic omega-limits for its gene-chain model. The original 1990 full hypotheses were not recovered; this audit does not promote that restatement alone to an unconditional global no-chaos theorem for the present six-state family.",
        "checks":{"seed":2026103001,"random_algebraic_states":128,"critical_damping_states":critical_count,
                  "errors":errors,"minimum_second_compound_offdiagonal":minimum_offdiag,
                  "maximum_strongly_connected_components_of_compound_graph":int(max_components),"passed":True},
        "equilibrium_controls":equilibrium_records,"existing_underdamped_references":references,
        "literature":[
            {"title":"The Poincare-Bendixson theorem for monotone cyclic feedback systems","authors":"Mallet-Paret and Smith","year":1990,"doi":"10.1007/BF01054041","url":"https://link.springer.com/article/10.1007/BF01054041","access":"Publisher abstract/metadata; full PDF redirected to subscription preview; no purchase."},
            {"title":"A Generalization of Linear Positive Systems with Applications to Nonlinear Systems: Invariant Sets and the Poincare-Bendixson Property","authors":"Weiss and Margaliot","year":2021,"doi":"10.1016/j.automatica.2020.109358","url":"https://arxiv.org/abs/1902.01630v4","access":"Full primary author manuscript; Theorem 5/6, Corollary 4, Theorem 12 inspected."},
            {"title":"Existence Criteria of Periodic Oscillations in Cyclic Gene Regulatory Networks","authors":"Hori, Hara and Kim","year":2010,"url":"https://www.keisu.t.u-tokyo.ac.jp/data/2010/METR10-11.pdf","access":"Full author/institutional research report; model equation 1, Proposition 1, Lemmas 2/3 and Appendix B inspected. This is an application/restatement, not the original 1990 proof."}
        ],"wall_seconds":time.perf_counter()-started,
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"checks":report['checks'],"references":references,"wall_seconds":report['wall_seconds']},indent=2))


if __name__ == '__main__':
    main()
