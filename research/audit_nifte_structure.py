"""Algebraic NIFTE controls: no ODE trajectories or parameter survey.

Reproduce with: python research/audit_nifte_structure.py
Parameters are read from the independently transcribed primary-source record.
All eigenvalues and frequencies below use physical seconds, not t/tau.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as Poly
from scipy.linalg import solve_continuous_lyapunov
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]


def linear_parts(p: dict, loss: bool) -> tuple[np.ndarray, np.ndarray]:
    """Dimensionless state, physical time; exact published component equations."""
    a = p["Ld"] + p["Lf"]
    den = p["Ll"] * a + p["Lp"] * (a + p["Ll"])
    c, b = p["U0"] / p["P0"], p["P0"] / p["U0"]
    n = 9 if loss else 5
    mat = np.zeros((n, n))
    inp = np.zeros(n)
    inp[0] = 1 / (p["Rth"] * p["Cad"])
    mat[0, 0] = -inp[0]
    mat[0, 3:5] = c / p["Cad"]
    mat[1, 3] = c / p["Cd"]
    mat[2, 4] = c / p["Cp"]
    mat[3, :3] = b * np.array([-p["Lp"], -(p["Lp"] + p["Ll"]), p["Ll"]]) / den
    mat[3, 3] = -(p["Ll"] * p["Rf"] + p["Lp"] * (p["Rf"] + p["Rl"])) / den
    mat[3, 4] = -p["Lp"] * p["Rl"] / den
    mat[4, :3] = b * np.array([-a, p["Ll"], -(a + p["Ll"])]) / den
    mat[4, 3] = (p["Ll"] * p["Rf"] - a * p["Rl"]) / den
    mat[4, 4] = -a * p["Rl"] / den
    if loss:
        mat[0, 0] -= 1 / (p["Rtl0"] * p["Cad"])
        mat[0, 5:] = 1 / (p["Rtl0"] * p["Cad"])
        for i, (res, cap) in enumerate(zip(p["Rtl"], p["Ctl"])):
            mat[5 + i, 0] = 1 / (p["Rtl0"] * cap)
            mat[5 + i, 5:] = -1 / (p["Rtl0"] * cap)
            mat[5 + i, 5 + i] -= 1 / (res * cap)
    return mat, inp


def transfer_polynomials(p: dict, loss: bool) -> tuple[Poly, Poly, Poly, dict]:
    """Independent pressure/flow elimination: transfer = -numerator/denominator."""
    s = Poly([0, 1])
    a = p["Ld"] + p["Lf"]
    h11 = 1 + p["Cd"] * (p["Rf"] + p["Rl"]) * s + p["Cd"] * (a + p["Ll"]) * s**2
    h12 = p["Cp"] * p["Rl"] * s + p["Cp"] * p["Ll"] * s**2
    h21 = p["Cd"] * p["Rl"] * s + p["Cd"] * p["Ll"] * s**2
    h22 = 1 + p["Cp"] * p["Rl"] * s + p["Cp"] * (p["Lp"] + p["Ll"]) * s**2
    hydraulic = h11 * h22 - h12 * h21
    bpoly = 1 + p["Lp"] * p["Cp"] * s**2
    cpoly = 1 + p["Cd"] * p["Rf"] * s + p["Cd"] * a * s**2
    qpoly = (1 + p["Rth"] * p["Cad"] * s) * hydraulic + p["Rth"] * s * (p["Cd"] * bpoly + p["Cp"] * cpoly)
    wpoly, vpoly = Poly([1]), Poly([1])
    if loss:
        factors = [1 + res * cap * s for res, cap in zip(p["Rtl"], p["Ctl"])]
        for factor in factors:
            wpoly *= factor
        vpoly = p["Rtl0"] * wpoly
        for i, res in enumerate(p["Rtl"]):
            term = Poly([res])
            for j, factor in enumerate(factors):
                if i != j:
                    term *= factor
            vpoly += term
        den = vpoly * qpoly + p["Rth"] * wpoly * hydraulic
        num = vpoly * bpoly
    else:
        den, num = qpoly, bpoly
    extra = {"hydraulic_J": hydraulic.coef.tolist(), "B": bpoly.coef.tolist(), "C": cpoly.coef.tolist(), "Q": qpoly.coef.tolist(), "V": vpoly.coef.tolist(), "W": wpoly.coef.tolist()}
    return den, num, vpoly, extra


def complex_pairs(values) -> list[list[float]]:
    return [[float(z.real), float(z.imag)] for z in sorted(values, key=lambda z: (z.real, z.imag))]


def crossing_candidates(den: Poly, num: Poly, vpoly: Poly, mat: np.ndarray, inp: np.ndarray) -> tuple[list, list]:
    # B(i*w) is real. Im(den(i*w)*conj(V(i*w)))/w is a polynomial in w**2.
    prod = np.polynomial.polynomial.polymul(den.coef * 1j**np.arange(len(den.coef)), vpoly.coef * (-1j)**np.arange(len(vpoly.coef)))
    condition = Poly(prod.imag[1::2])
    candidates = []
    for root in condition.roots():
        if abs(root.imag) > 1e-7 * (1 + abs(root.real)) or root.real <= 0:
            continue
        omega = float(np.sqrt(root.real))
        z = 1j * omega
        ratio = -den(z) / num(z)
        if ratio.real <= 0 or abs(ratio.imag) > 1e-6 * (1 + abs(ratio.real)):
            continue
        gain = float(ratio.real)
        closed = mat.copy()
        closed[:, 1] += gain * inp
        eig = np.linalg.eigvals(closed)
        derivative = -num(z) / (den.deriv()(z) + gain * num.deriv()(z))
        candidates.append({"gain_FK": gain, "omega_rad_per_s": omega, "linear_period_seconds": float(2 * np.pi / omega), "eigenvalue_crossing_d_real_d_gain": float(derivative.real), "nearest_eigenvalue_error": float(np.min(abs(eig - z))), "other_unstable_eigenvalues": int(np.count_nonzero(eig.real > 1e-7)), "eigenvalues_per_s": complex_pairs(eig), "interpretation": "Imaginary-pair crossing only; nonlinear Hopf nondegeneracy not established."})
    return sorted(candidates, key=lambda r: r["gain_FK"]), condition.coef.tolist()


def audit_case(p: dict, loss: bool) -> dict:
    mat, inp = linear_parts(p, loss)
    den, num, vpoly, polys = transfer_polynomials(p, loss)
    gain = p["F"] * p["K"]
    closed = mat.copy()
    closed[:, 1] += gain * inp
    expected = (den + gain * num).coef[::-1]
    expected /= expected[0]
    observed = np.poly(closed)
    poly_error = float(np.max(abs(expected - observed) / (1 + abs(expected))))
    selector = np.eye(len(mat))[1]
    transfer_error = 0.
    for omega in np.geomspace(.001, 1000, 80):
        actual = selector @ np.linalg.solve(1j * omega * np.eye(len(mat)) - mat, inp)
        formula = -num(1j * omega) / den(1j * omega)
        transfer_error = max(transfer_error, float(abs(actual - formula) / (1 + abs(formula))))
    assert poly_error < 1e-8 and transfer_error < 1e-8
    gramian = solve_continuous_lyapunov(mat.T, -np.eye(len(mat)))
    gramian = (gramian + gramian.T) / 2
    w = gramian @ inp
    sector_mat = np.outer(w, selector) + np.outer(selector, w)
    gain_limit = 1 / np.linalg.eigvalsh(sector_mat)[-1]
    certified_control = gain_limit / 2
    gramian_residual = np.linalg.norm(mat.T @ gramian + gramian @ mat + np.eye(len(mat)), ord=2)
    contraction_margin = np.linalg.eigvalsh(np.eye(len(mat)) - certified_control * sector_mat)[0]
    assert np.linalg.eigvalsh(gramian)[0] > 0 and gramian_residual < 1e-7 and contraction_margin > .49
    rng = np.random.default_rng(20261030)
    scaling_error = odd_error = 0.
    for _ in range(32):
        x = rng.normal(0, .01, len(mat));factor = 2.7
        fx = mat @ x + inp * p["F"] * np.tanh(p["K"] * x[1])
        scaled_fx = mat @ (factor * x) + inp * (factor * p["F"]) * np.tanh(p["K"] / factor * (factor * x[1]))
        reverse_fx = mat @ (-x) + inp * p["F"] * np.tanh(-p["K"] * x[1])
        scaling_error = max(scaling_error, float(np.max(abs(scaled_fx - factor * fx))))
        odd_error = max(odd_error, float(np.max(abs(reverse_fx + fx))))
    crossings, crossing_poly = crossing_candidates(den, num, vpoly, mat, inp)
    return {"dimension": len(mat), "loss_model": loss, "Lf": p["Lf"], "gain_FK": gain, "time_unit": "physical seconds", "passive_eigenvalues_per_s": complex_pairs(np.linalg.eigvals(mat)), "origin_eigenvalues_per_s": complex_pairs(np.linalg.eigvals(closed)), "trace_per_s": float(np.trace(mat)), "characteristic_coefficients_descending_monic": expected.tolist(), "transfer_coefficients_ascending": {"denominator": den.coef.tolist(), "numerator_with_minus_sign_external": num.coef.tolist(), **polys}, "hopf_frequency_condition_coefficients_ascending_in_omega_squared": crossing_poly, "imaginary_pair_crossings": crossings, "quadratic_small_gain_control": {"P": gramian.tolist(), "P_min_eigenvalue": float(np.linalg.eigvalsh(gramian)[0]), "lyapunov_equation_residual_2norm": float(gramian_residual), "sufficient_strict_gain_bound": float(gain_limit), "control_gain": float(certified_control), "control_margin_min_eigenvalue": float(contraction_margin), "claim": "Numerically constructed global incremental contraction certificate in the stated coordinates at control_gain. Sufficient, not necessary; not interval-certified."}, "checks": {"characteristic_scaled_error": poly_error, "frequency_transfer_scaled_error": transfer_error, "scaling_absolute_error": scaling_error, "odd_symmetry_absolute_error": odd_error}}


def double_pair_candidate(p: dict) -> dict:
    """Scalar algebraic root and cubic local amplitude coefficients, no ODE solve."""
    def residual(lf):
        den, num, _, _ = transfer_polynomials(dict(p, Lf=lf), False)
        q = den.coef
        gain = q[4] * q[1] / q[5] - q[0]
        return q[2] + gain * num.coef[2] - q[4] * q[3] / q[5]

    lf = brentq(residual, 2e5, 3e5, xtol=1e-7)
    params = dict(p, Lf=lf)
    den, num, _, _ = transfer_polynomials(params, False)
    q = den.coef
    gain = q[4] * q[1] / q[5] - q[0]
    omega = np.sqrt(np.sort(np.roots([q[5], -q[3], q[1]])))
    gamma = q[4] / q[5]
    s = Poly([0, 1])
    factored = q[5] * (s + gamma) * (s**2 + omega[0]**2) * (s**2 + omega[1]**2)
    factor_error = float(np.max(abs(factored.coef - (den + gain * num).coef)))
    mat, inp = linear_parts(params, False)
    mat[:, 1] += gain * inp
    eig, vectors = np.linalg.eig(mat)
    vectors /= vectors[1:2, :]
    dual = np.linalg.inv(vectors)
    alphas = []
    alpha_errors = []
    for w in omega:
        idx = int(np.argmin(abs(eig - 1j * w)))
        alpha = dual[idx] @ inp
        derivative = -num(1j * w) / (den.deriv()(1j * w) + gain * num.deriv()(1j * w))
        alphas.append(alpha)
        alpha_errors.append(abs(alpha - derivative))
    cubic = -gain * np.real(alphas)[:, None] * np.array([[1., 2.], [2., 1.]])

    def real_parts(g, log_lf_shift):
        m, b = linear_parts(dict(p, Lf=lf * np.exp(log_lf_shift)), False)
        m[:, 1] += g * b
        ev = np.linalg.eigvals(m)
        return np.array([ev[np.argmin(abs(ev - 1j * w))].real for w in omega])

    unfolding = []
    for shrink in (1., .5):
        dg, dl = shrink * 1e-3, shrink * 1e-5
        jac = np.column_stack(((real_parts(gain + dg, 0) - real_parts(gain - dg, 0)) / (2 * dg), (real_parts(gain, dl) - real_parts(gain, -dl)) / (2 * dl)))
        unfolding.append({"gain_step": dg, "log_Lf_step": dl, "jacobian_real_eigenvalues_vs_gain_logLf": jac.tolist(), "determinant": float(np.linalg.det(jac))})
    assert factor_error < 1e-9 and max(alpha_errors) < 1e-10 and np.linalg.det(cubic) < 0
    return {"model": "five-state NTP only", "Lf": float(lf), "gain_FK": float(gain), "F_reference": p["F"], "K_at_reference_F": float(gain / p["F"]), "omega_rad_per_s": omega.tolist(), "frequency_ratio": float(omega[1] / omega[0]), "stable_real_eigenvalue": float(-gamma), "eigenvalues_per_s": complex_pairs(eig), "factorization_absolute_error": factor_error, "root_bracket_Lf": [2e5, 3e5], "scalar_residual": float(residual(lf)), "parameter_transversality_checks": unfolding, "normal_form": {"coordinates": "xi=K*X, critical right eigenvectors normalized to d component 1, modal expansion xi=sum_j(z_j*q_j+conj(z_j*q_j))+stable component", "gain_eigenvalue_derivatives_real_imag": [[float(a.real), float(a.imag)] for a in alphas], "independent_derivative_errors": [float(a) for a in alpha_errors], "radial_cubic_matrix": cubic.tolist(), "radial_cubic_determinant": float(np.linalg.det(cubic)), "derivation": "xi'=A_G*xi-(G/3)*b*xi_d^3+O(xi_d^5). No quadratic term. Cubic resonances project to -G*alpha_j*z_j*|z_j|^2 and -2*G*alpha_j*z_j*|z_k|^2. Frequencies have no cubic 1:1 or 1:3 resonance.", "implication": "Leading local cubic amplitude system has mode competition: any mixed positive-amplitude equilibrium has saddle radial stability. It does not predict an attracting mixed torus near this point at leading order.", "limits": "No interval certification, center-manifold remainder bound, periodic continuation, or global chaos exclusion. Away from this local expansion other dynamics remain possible."}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/nifte/structural-audit.json")
    args = parser.parse_args()
    source_file = ROOT / "results/candidate-audits/physical-candidate-screening.json"
    source = json.loads(source_file.read_text(encoding="utf-8"))
    candidate = source["retained_candidates"][0]
    p = candidate["nominal_parameters"]
    cases = {}
    for control in (False, True):
        params = dict(p)
        if control:
            params["Lf"] = 4.74e5
        for loss in (False, True):
            name = ("ntp_tl" if loss else "ntp") + ("_control" if control else "_nominal")
            cases[name] = audit_case(params, loss)
    result = {"audit_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "scope": "Algebraic controls only; no trajectories, parameter survey, attractor or novelty claim.", "primary_doi": candidate["source_doi"], "source_record": str(source_file.relative_to(ROOT)).replace("\\", "/"), "source_record_sha256": hashlib.sha256(source_file.read_bytes()).hexdigest(), "parameter_reference": p, "script": "research/audit_nifte_structure.py", "documentation": "docs/nifte-structural-audit.md", "exact_scaling": "X2=(F2/F1)*X1 when F1*K1=F2*K2, with all components fixed and initial states scaled consistently.", "cases": cases, "double_imaginary_pair_candidate": double_pair_candidate(p)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for name, case in cases.items():
        print(name, "gain=", case["gain_FK"], "crossings=", [(r["gain_FK"], r["omega_rad_per_s"]) for r in case["imaginary_pair_crossings"]], "small_gain=", case["quadratic_small_gain_control"]["control_gain"])
    print("Wrote", args.output)


if __name__ == "__main__":
    main()
