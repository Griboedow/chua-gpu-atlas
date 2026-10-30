"""Independent algebraic audit of the literal Case et al. four-state closure.

No trajectories, fitted-coefficient changes, or GPU calls. Equilibrium roots
are bracketed numerically, not certified to be an exhaustive algebraic census.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
L1, L2, L3, L4, L5 = 24., 144., 24., 240., 120.
ALPHA, BETA, ZETA = 31., .2565, 1.08
LB, GAMMA, RB, RS, QS = 42., 7.35, 4.1, 5.5, 25.
EPS, K3, K5, P0, Q0 = 10.9, 15.75, 12., .12, 8e-6
A0, B0 = 1/L1+1/L5, 1/L2+1/L4
MAT = np.array([[A0+ZETA/L3, -1/L3], [-ZETA/L3, B0+1/L3]])


def field(s, pressures):
    x, y, z, r = s
    eta = P0+Q0*x
    xi = RS*RS*eta-QS/RS
    resistance = LB+GAMMA*(r-RB)**2
    mix = z*(x-z)/x
    a3 = -resistance*z/L3+K3*mix
    a4 = -(ALPHA+BETA*(y+z))*(y+z)/12
    a5 = -(x-z)+K5*mix
    p3, p4 = np.linalg.solve(MAT, [pressures[0]/L1-x-a3-a5,
                                 pressures[1]/L2-y-a4+a3])
    return np.array([12*((pressures[0]-p3)/L1-x),
                     12*((pressures[1]-p4)/L2-y),
                     12*((ZETA*p3-p4)/L3+a3),
                     EPS*(z-QS-eta*(r-RS)**3+xi*(r-RS))])


def jacobian(s, pressures):
    return np.column_stack([field(np.asarray(s, dtype=complex)+1e-25j*np.eye(4)[i], pressures).imag/1e-25
                            for i in range(4)])


def energy_gradient(s):
    x, y, z, r = s
    return np.array([ZETA*(L1*x+L5*(x-z)), L2*y+L4*(y+z),
                     L3*z+L4*(y+z)-ZETA*L5*(x-z), 0.])/12


def energy_power(s, pressures):
    x, y, z, r = s
    w, u = x-z, y+z
    mix = z*w/x
    resistance = LB+GAMMA*(r-RB)**2
    return (ZETA*pressures[0]*x+pressures[1]*y-ZETA*L1*x*x-L2*y*y-resistance*z*z
            -L4*(ALPHA+BETA*u)*u*u/12-ZETA*L5*w*w
            +L3*K3*z*mix+ZETA*L5*K5*w*mix)


def conditional_bound_constants(pressures):
    # Coordinates of this quadratic form are (w=Q5,y=Q2,z=Q3).
    gram = np.array([[ZETA*(L1+L5), 0., ZETA*L1],
                     [0., L2+L4, L4], [ZETA*L1, L4, ZETA*L1+L3+L4]])
    a = ZETA*L5
    cross = ZETA*L5*K5
    C = L3*K3+cross*cross/(2*a)
    d = L4*BETA/12
    forcing = ZETA*pressures[0]
    zp = (2*(C+1)+np.sqrt(4*(C+1)**2+12*d*forcing))/(6*d)
    Kz = forcing*zp+(C+1)*zp*zp-d*zp**3
    K = forcing**2/a+pressures[1]**2/(2*L2)+Kz
    mineig, maxeig = np.linalg.eigvalsh(gram)[[0, -1]]
    decay = 24*min(a/4, L2/2, 1.)/maxeig
    eventual_energy = K/decay
    eventual_component_bound = np.sqrt(24*eventual_energy/mineig)
    return {"energy_gram_w_y_z": gram.tolist(), "gram_eigenvalue_min": float(mineig),
            "gram_eigenvalue_max": float(maxeig), "quadratic_z_coefficient_C": C,
            "cubic_damping_d": d, "forcing_constant_K": float(K), "energy_decay_k": float(decay),
            "limsup_energy_bound_if_forward_complete": float(eventual_energy),
            "limsup_each_of_w_y_z_bound_if_forward_complete": float(eventual_component_bound),
            "limsup_r_bound_if_forward_complete": float(max(4*RS, (8*eventual_component_bound/(3*P0))**(1/3))),
            "qualification": "Conditional on remaining in the regular forward-flow cone for all forward time; not an invariant-set certificate."}


def z_max(P2):
    a, b = L4*BETA/12, L4*ALPHA/12
    return 2*P2/(b+np.sqrt(b*b+4*a*P2))


def from_z(z, branch_sign, pressures):
    P1, P2 = pressures
    b = P1+L5*(1+K5)*z
    disc = b*b-4*(L1+L5)*L5*K5*z*z
    if disc < 0:
        return None
    x = (b+np.sqrt(disc))/(2*(L1+L5))
    c4, b4 = L4*BETA/12, L2+L4*ALPHA/12
    u = 2*(P2+L2*z)/(b4+np.sqrt(b4*b4+4*c4*(P2+L2*z)))
    y = u-z
    mix = z*(x-z)/x
    required_R = (ZETA*(P1-L1*x)-(P2-L2*y)+L3*K3*mix)/z
    if required_R < LB:
        return None
    r = RB+branch_sign*np.sqrt((required_R-LB)/GAMMA)
    if x <= 0 or y < -1e-10 or z < 0 or z > x or r < 0:
        return None
    return np.array([x, y, z, r])


def equilibria(pressures, samples=1801):
    P1, P2 = pressures
    bound = z_max(P2)
    if bound >= P1/L1:
        raise ValueError('This root scan omits the second x quadratic root when it can enter the forward cone')
    found = []
    for sign in [-1, 1]:
        previous = None
        for z in np.linspace(1e-7, bound, samples):
            state = from_z(z, sign, pressures)
            if state is None:
                previous = None
                continue
            value = field(state, pressures)[3]
            if previous is not None and previous[1]*value < 0:
                root_z = brentq(lambda zz: field(from_z(zz, sign, pressures), pressures)[3],
                                 previous[0], z, xtol=1e-12)
                state = from_z(root_z, sign, pressures)
                if all(np.linalg.norm(state-v) > 1e-7 for v in found):
                    found.append(state)
            previous = z, value
    return sorted(found, key=lambda s: s[3])


def equilibrium_record(s, pressures):
    eig = np.linalg.eigvals(jacobian(s, pressures))
    eta = P0+Q0*s[0]
    return {"state": s.tolist(), "field_residual_inf": float(max(abs(field(s, pressures)))),
            "eigenvalues_real_imag": [[float(v.real), float(v.imag)] for v in eig],
            "maximum_eigenvalue_real_part": float(max(eig.real)),
            "eta": float(eta), "xi": float(RS*RS*eta-QS/RS),
            "forward_flow_domain": bool(s[0] > 0 and s[1] >= 0 and 0 <= s[2] <= s[0] and s[3] >= 0)}


def main():
    pressures = np.array([62500., 62500.])
    boundary_y = brentq(lambda y: pressures[1]/L2-y+(ALPHA+BETA*y)*y/12
                       -ZETA*B0*pressures[0]/(L1*A0), 0., 1000.)
    boundary = np.array([1000., 100., 0., 4.1])
    rng = np.random.default_rng(2026103017)
    errors = []
    for _ in range(128):
        x = rng.uniform(1, 3000)
        state = np.array([x, rng.uniform(0, 300), rng.uniform(0, 1)*x, rng.uniform(0, 12)])
        p = rng.uniform(1000, 150000, 2)
        direct = energy_gradient(state)@field(state, p)
        expected = energy_power(state, p)
        errors.append(abs(direct-expected)/(1+abs(expected)))
    sampled_pressures = np.unique(np.r_[np.linspace(40000, 100000, 25), 62500])
    linear_controls = [{"equal_inlet_dimensionless_pressure": float(p),
                        "equal_inlet_pressure_Pa": .004*float(p),
                        "bracketed_equilibria": [equilibrium_record(s, [p, p]) for s in equilibria([p, p])]}
                       for p in sampled_pressures]
    result = {
        "model": "Literal Case et al. 2020 S5-S13 four-state closure, including positive minor-loss terms and both Q1-dependent vortex coefficients",
        "primary_source": "https://arxiv.org/pdf/2006.09400",
        "journal_doi": "10.1126/sciadv.aay6761",
        "state_order": ["Q1=x", "Q2=y", "Q3=z", "r"],
        "scope": "Algebra, conditional boundedness, and equilibrium linear analysis only; no trajectories or scientific novelty claim",
        "forward_domain": {"definition": "x>0, y>=0, 0<=z<=x, r>=0",
            "globally_forward_invariant": False,
            "counterexample_inlet_pressures": pressures.tolist(), "counterexample_state": boundary.tolist(),
            "counterexample_vector_field": field(boundary, pressures).tolist(),
            "z_zero_boundary_outward_above_y": float(boundary_y),
            "r_zero_boundary": "r'=epsilon*z; inward conditional on z>=0",
            "consequence": "Strict interior starts sufficiently close to the outward z=0 boundary also leave the forward domain; no clipping or sign repair is authorized by the source equations."},
        "conditional_energy_bound": conditional_bound_constants(pressures),
        "equilibrium_flow_bound_z_at_source": float(z_max(pressures[1])),
        "source_pressure_equilibria": [equilibrium_record(s, pressures) for s in equilibria(pressures)],
        "equal_pressure_linear_controls": linear_controls,
        "root_search_qualification": "1801 z samples on each admissible resistance branch, then Brent brackets. No certification against unbracketed multiple/tangential roots or bifurcations between pressure samples.",
        "vortex_xi_positive_only_above_Q1": float((QS/RS**3-P0)/Q0),
        "vortex_contraction_condition": "For Q1 below the stated threshold, xi<0 and partial r'/partial r=epsilon*(xi-3*eta*(r-rstar)^2)<0; this is not a four-state chaos exclusion.",
        "energy_identity_check": {"seed": 2026103017, "sample_count": 128,
                                   "maximum_scaled_residual": float(max(errors))},
        "global_chaos_exclusion": "Not established. The forward cone is not invariant; no global planar reduction or monotone theorem was proved.",
        "reproduction": "python research/audit_microfluidic_structure.py"}
    target = ROOT/'results/microfluidic/structural-audit.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({"output": str(target), "counterexample_vector": result['forward_domain']['counterexample_vector_field'],
                      "source_equilibria": result['source_pressure_equilibria'],
                      "energy_identity_check": result['energy_identity_check']}, indent=2))


if __name__ == '__main__':
    main()
