"""Algebraic controls for the exact, unsaturated single-cell HHS equation S2.

No trajectories or parameter survey are computed. Run from any working directory:
    python research/audit_hhs_structure.py
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
MR = float(brentq(lambda m: m**3 - m*m + .4*m - .1, .25, 1,
                 xtol=1e-15))
ER = 3*MR/(MR*MR + .1)


@dataclass(frozen=True)
class Parameters:
    a: float = 1.
    b: float = 3.
    c: float = 1.
    D0: float = .1
    Er: float = ER
    tauD: float = 1.1
    tauL: float = 1.2
    tauR: float = 1.


def field(z, mu, p=Parameters()):
    M, E, D, y, x = z
    return np.array([p.a-(p.b+1)*M+p.c*M*M*E,
                     p.b*M-p.c*M*M*E-D*E,
                     (p.D0-D+mu*D*x)/p.tauD,
                     (E-p.Er-y)/p.tauL, (y-x)/p.tauR])


def jacobian(z, mu, p=Parameters()):
    M, E, D, y, x = z
    return np.array([[-p.b-1+2*p.c*M*E, p.c*M*M, 0, 0, 0],
                     [p.b-2*p.c*M*E, -p.c*M*M-D, -E, 0, 0],
                     [0, 0, (-1+mu*x)/p.tauD, 0, mu*D/p.tauD],
                     [0, 1/p.tauL, 0, -1/p.tauL, 0],
                     [0, 0, 0, 1/p.tauR, -1/p.tauR]])


def equilibrium_polynomial(mu, p=Parameters()):
    C = 1+mu*p.Er
    return np.array([p.c*C, -p.c*C*p.a-mu*(p.b+1),
                     (p.b+1)*p.D0+mu*p.a*(p.b+2),
                     -p.a*p.D0-mu*p.a*p.a])


def additional_quadratic(mu):
    q, rem = np.polydiv(equilibrium_polynomial(mu), [1, -MR])
    return q, float(np.max(np.abs(rem)))


def discriminant(mu):
    q, _ = additional_quadratic(mu)
    return q[1]**2-4*q[0]*q[2]


def equilibrium_state(M, p=Parameters()):
    E = ((p.b+1)*M-p.a)/(p.c*M*M)
    D = (p.a-M)/E
    return np.array([M, E, D, E-p.Er, E-p.Er])


def equilibrium_states(mu):
    """Factor the persistent root before computing the other two roots."""
    q, _ = additional_quadratic(mu)
    ds = discriminant(mu)
    roots = [MR]
    if ds >= -1e-12:
        sqrt_ds = np.sqrt(max(0., ds))
        roots += [(-q[1]-sqrt_ds)/(2*q[0]),
                  (-q[1]+sqrt_ds)/(2*q[0])]
    return [equilibrium_state(M) for M in sorted(roots)
            if .25 < M < 1.]


def lower_branch_max_real(mu):
    return float(np.max(np.linalg.eigvals(
        jacobian(equilibrium_states(mu)[0], mu)).real))


def spectrum_record(z, mu):
    eig = sorted(np.linalg.eigvals(jacobian(z, mu)),
                 key=lambda v: (float(v.real), float(v.imag)))
    return {"state_order": ["M", "E", "D", "y", "x"],
            "state": z.tolist(),
            "field_residual_inf": float(np.max(np.abs(field(z, mu)))),
            "eigenvalues_real_imag": [[float(v.real), float(v.imag)] for v in eig],
            "maximum_real_part": float(max(v.real for v in eig))}


def bound_constants(z0, mu, p=Parameters()):
    """Return a conservative finite global bound, in logarithmic form for D."""
    M, E, D, y, x = map(float, z0)
    if min(M, E) < 0 or D <= 0 or mu < 0:
        raise ValueError("Require M,E >= 0, D > 0, mu >= 0")
    h = min(x, y, -p.Er)
    Ch = 1-mu*h
    delta = min(D, p.D0/Ch)
    B = max(M+E, p.a/min(1., delta))
    H = max(x, y, B-p.Er)
    C = 1+mu*p.Er
    rho = C/(8*p.tauD*p.b*B)
    Dstar = max(4*p.D0/C, mu/(rho*p.tauD))
    Zlo = mu/p.tauD*(p.tauL+p.tauR)*h
    Zhi = mu/p.tauD*(p.tauL+p.tauR)*H+rho*B
    V0 = np.log(D)+mu/p.tauD*(p.tauL*y+p.tauR*x)+rho*E
    log_bound = max(V0, np.log(Dstar)+Zhi)-Zlo
    return dict(h=h, Ch=Ch, delta=delta, chemical_sum_bound=B,
                filter_upper_bound=H, C=C, rho=rho, Dstar=Dstar,
                Zlo=Zlo, Zhi=Zhi, V0=V0,
                log_global_D_bound=log_bound,
                large_D_derivative_upper_bound=-C/(2*p.tauD))


def algebra_checks():
    rng = np.random.default_rng(20261030)
    p = Parameters()
    jac_errors, barrier_errors, poly_errors, determinant_errors = [], [], [], []
    for _ in range(64):
        z = np.r_[rng.uniform(.05, 6, 3), rng.uniform(-5, 5, 2)]
        mu = float(rng.uniform(0, 12))
        J = jacobian(z, mu)
        Jcs = np.column_stack([field(z.astype(complex)+1e-25j*np.eye(5)[k], mu).imag/1e-25
                              for k in range(5)])
        jac_errors.append(np.max(np.abs(J-Jcs)))
        M, E, D, y, x = z
        rho = .123
        grad = np.array([0, rho, 1/D, mu*p.tauL/p.tauD, mu*p.tauR/p.tauD])
        rhs = p.D0/(p.tauD*D)-(1+mu*p.Er)/p.tauD + \
            (mu/p.tauD-rho*D)*E+rho*p.b*M-rho*p.c*M*M*E
        barrier_errors.append(abs(grad@field(z, mu)-rhs))
        direct = ((p.b+1)*M-p.a)*(p.D0+mu*(p.a-M)) - \
            p.c*(1+mu*p.Er)*M*M*(p.a-M)
        poly_errors.append(abs(np.polyval(equilibrium_polynomial(mu), M)-direct))
        ref = np.array([MR, ER, p.D0, 0, 0])
        A, B, Q = -4+2*MR*ER, MR*MR, 3-2*MR*ER
        analytic_det = (A*(B+p.D0+mu*ER*p.D0)+B*Q)/(p.tauD*p.tauL*p.tauR)
        determinant_errors.append(abs(np.linalg.det(jacobian(ref, mu))-analytic_det))
    return {"samples": 64, "seed": 20261030,
            "complex_step_jacobian_max_absolute_error": float(max(jac_errors)),
            "barrier_derivative_identity_max_absolute_error": float(max(barrier_errors)),
            "equilibrium_polynomial_identity_max_absolute_error": float(max(poly_errors)),
            "reference_determinant_identity_max_absolute_error": float(max(determinant_errors))}


def main():
    sn = brentq(discriminant, 1., 1.1, xtol=1e-15)
    collision = brentq(lambda mu: np.polyval(np.polyder(equilibrium_polynomial(mu)), MR),
                       1., 1.1, xtol=1e-15)
    hops = [brentq(lower_branch_max_real, *bracket, xtol=1e-13)
            for bracket in [(1.06, 1.1), (2., 3.)]]
    ref = np.array([MR, ER, .1, 0., 0.])
    q, _ = additional_quadratic(sn)
    sn_state = equilibrium_state(-q[1]/(2*q[0]))
    A, B, Q = -4+2*MR*ER, MR*MR, 3-2*MR*ER
    source_literal_mu = (B*Q/A-B-.1)/(.1*ER)
    source = ROOT/'outputs/candidate_research/hbr-2026-v3.html'
    result = {
        "scope": "Independent algebra, global boundedness proof and linear equilibrium controls; no trajectories",
        "model": "Exact single-cell S2 with f=0; no spatial cubic terms",
        "source": {"url": "https://arxiv.org/html/2504.19655v3",
                   "version": "arXiv:2504.19655v3",
                   "journal_doi": "10.1103/p1bw-v87g",
                   "cached_html_sha256": hashlib.sha256(source.read_bytes()).hexdigest() if source.exists() else None},
        "parameters": asdict(Parameters()),
        "reference_M": MR, "reference_E": ER,
        "boundedness": {
            "status": "Analytic initial-dependent global bound for every finite admissible initial state",
            "assumptions": "Positive a,b,c,D0,Er,tauD,tauL,tauR; mu>=0; Minit,Einit>=0; Dinit>0; xinit,yinit finite real",
            "not_claimed": ["parameter-only absorbing bound", "positive physical lengths", "no chaos"],
            "example_initial_state": [1.05, 3.43, .1, 0., 0.], "example_mu": 1.5,
            "example_constants": bound_constants([1.05, 3.43, .1, 0., 0.], 1.5)},
        "equilibrium_events": {
            "additional_branch_saddle_node": {"mu": sn, **spectrum_record(sn_state, sn)},
            "persistent_reference_zero_eigenvalue": {"mu": collision,
                "interpretation": "Collision with an additional branch; not a Hopf event",
                **spectrum_record(ref, collision)},
            "compressed_branch_Hopf_candidates": [
                {"mu": mu, "max_real_slope_centered_difference":
                 (lower_branch_max_real(mu+1e-5)-lower_branch_max_real(mu-1e-5))/2e-5,
                 "qualification": "Linear crossing only; no first Lyapunov coefficient or cycle continuation",
                 **spectrum_record(equilibrium_states(mu)[0], mu)} for mu in hops]},
        "selected_equilibria": [
            {"mu": mu, "states": [spectrum_record(z, mu) for z in equilibrium_states(mu)]}
            for mu in [0., .8, 1.0564, 1.0566, 1.1, 1.5, 2., 3.]],
        "source_displayed_determinant_comparison": {
            "matrix_S6": "tauD*tauL*tauR*detJ=A*(B+D0+mu*Er*D0)+B*Q",
            "expanded_S7_S8": "Displayed source expression uses -B*Q instead of +B*Q",
            "definitions": "A=-(b+1)+2cMrEr; B=cMr^2; Q=b-2cMrEr",
            "literal_expanded_expression_zero_mu": float(source_literal_mu),
            "matrix_zero_mu": collision,
            "qualification": "A discrepancy between displayed formulas; no inference about unpublished implementation"},
        "figure_S1_reproduction_question": {
            "fixed_reference_equilibrium_count_at_mu_0_8": len(equilibrium_states(.8)),
            "fixed_reference_equilibrium_E_at_mu_0_8": ER,
            "interpretation": "Different exact constant asymptotes at this gain require another vector field or parameter convention; finite traces alone do not identify the implementation. Source uses E0 for both reference and sampled initial values."},
        "checks": algebra_checks(),
        "reproduction": "python research/audit_hhs_structure.py"}
    out = ROOT/'results/hhs/structural-audit.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({"output": str(out), "saddle_node_mu": sn,
                      "reference_zero_mu": collision, "Hopf_candidates_mu": hops,
                      "checks": result['checks']}, indent=2))


if __name__ == '__main__':
    main()
