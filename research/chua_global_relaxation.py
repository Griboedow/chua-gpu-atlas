"""Reproduce the Chua global-convergence criterion and its prospective controls."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

import matplotlib
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FORECAST = ROOT / "results/chua-global-relaxation/forecast.json"


def coefficients(alpha, beta, gamma):
    return alpha + 1 + gamma**2 - 2 * beta, (beta + gamma) ** 2 - alpha * (
        beta - gamma**2
    )


def classify(alpha, beta, gamma):
    """Exact rational input. A witness uses the larger positive root of F."""
    if alpha <= 0 or beta <= 0 or gamma < 0:
        raise ValueError("Expected alpha > 0, beta > 0 and gamma >= 0")
    B, C = coefficients(alpha, beta, gamma)
    delta = B**2 - 4 * C
    guaranteed = bool(C >= 0 and (B >= 0 or delta < 0))
    result = {
        "B": str(B),
        "C": str(C),
        "all_diodes_bounded_orbits_approach_equilibria": guaranteed,
    }
    if not guaranteed:
        q = (-B + sp.sqrt(delta)) / 2
        a, b = 1 + gamma, beta + gamma
        # Re H(i sqrt(q)) = (gamma*b + q)/|D(i sqrt(q))|^2.
        m0 = sp.simplify((gamma * b + q) / ((b - q) ** 2 + a**2 * q) - 1)
        result["specified_diode"] = "Inconclusive without its characteristic."
        result["periodic_orbit_witness"] = {
            "omega_squared": str(q),
            "inner_diode_slope_m0": str(m0),
            "range": "Any 0 < x amplitude < 1, with any outer diode slope.",
        }
    return result


def symbolic_checks():
    alpha, beta, gamma, c = sp.symbols("alpha beta gamma c")
    w, v, dv = sp.symbols("w v dv")
    a, b, ell = 1 + gamma, beta + gamma, beta - gamma**2
    B, C = coefficients(alpha, beta, gamma)
    S = (
        (a + alpha * gamma / (2 * b)) * v**2
        + (2 * b - alpha * ell / b - 2 * c) * w * v
        + (a * b - alpha * ell * a / (2 * b) + alpha * gamma / 2) * w**2
    )
    u = dv + a * v + b * w
    dV = (alpha * (ell * w - gamma * v) / b - u) * u
    dV += sp.diff(S, w) * v + sp.diff(S, v) * dv
    Q = (dv + c * w) ** 2 + (B + 2 * c) * v**2 + (C - c**2) * w**2
    discrim = (alpha + 1 - gamma**2) ** 2 - 4 * a**2 * beta
    assert sp.cancel(dV + Q) == 0
    assert sp.expand(B**2 - 4 * C - discrim) == 0
    epsilon = sp.Rational(1, 10**30)
    cases = [
        (1, 1, 0, True),
        (sp.Rational(1, 2), sp.Rational(1, 2), 0, True),
        (9, 25 - epsilon, 0, False),
        (9, 25, 0, False),
        (9, 25 + epsilon, 0, True),
        (9, 15, sp.Rational(1, 4), False),
        (9, 19, sp.Rational(1, 4), True),
    ]
    for alpha, beta, gamma, expected in cases:
        answer = classify(alpha, beta, gamma)
        assert answer["all_diodes_bounded_orbits_approach_equilibria"] == expected
        if not expected:
            witness = answer["periodic_orbit_witness"]
            q, m0 = (
                sp.sympify(witness[k])
                for k in ("omega_squared", "inner_diode_slope_m0")
            )
            # The original central Jacobian must have the factor s^2+q.
            s = sp.Symbol("s")
            D = s**2 + (1 + gamma) * s + beta + gamma
            characteristic = (s + alpha * (1 + m0)) * D - alpha * (s + gamma)
            assert q > 0 and sp.simplify(sp.rem(characteristic, s**2 + q, s)) == 0


def integrate(alpha, beta, gamma, m0, m1, start, duration, method="LSODA"):
    def rhs(_, state):
        x, y, z = state
        h = m1 * x + (m0 - m1) * max(-1.0, min(1.0, x))
        return [alpha * (y - x - h), x - y + z, -beta * y - gamma * z]

    sol = solve_ivp(
        rhs,
        (0, duration),
        start,
        method=method,
        dense_output=True,
        rtol=2e-11,
        atol=2e-13,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def equilibrium(beta, gamma, m0, m1):
    b = beta + gamma
    return (m1 - m0) / (beta / b + m1) * np.array([1, gamma / b, -beta / b])


def well_controls(case):
    alpha = case["alpha"]
    m0, m1 = (float(Fraction(case[k])) for k in ("m0", "m1"))
    rows, traces = [], []
    for beta, gamma in case["parameter_pairs_beta_gamma"]:
        eq = equilibrium(beta, gamma, m0, m1)
        for x0 in case["initial_x"]:
            sol = integrate(
                alpha,
                beta,
                gamma,
                m0,
                m1,
                [x0, *case["initial_yz"]],
                case["test_duration"],
            )
            # The integral of clip(x,-1,1) from zero to x.
            sat_integral = (x0**2 - max(0, abs(x0) - 1) ** 2) / 2
            b = beta + gamma
            U = (1 + m1 - gamma / b) * x0**2 / 2 + (m0 - m1) * sat_integral
            V0 = alpha * U + (1 + gamma + alpha * gamma / (2 * b)) * x0**2
            error = float(np.max(abs(sol.y[:, -1] - np.sign(x0) * eq)))
            signed_min = float(np.min(np.sign(x0) * sol.y[0]))
            rows.append(
                {
                    "beta": beta,
                    "gamma": gamma,
                    "x0": x0,
                    "V0": V0,
                    "minimum_signed_x_at_solver_nodes": signed_min,
                    "terminal_max_error": error,
                    "passed": V0 < 0
                    and signed_min > 0
                    and error < case["test_terminal_max_norm"],
                }
            )
            if gamma == 0 and x0 == 0.9:
                times = np.linspace(0, 360, 18001)
                traces.append((times, sol.sol(times)[0]))
    return rows, traces[0]


def sharpness_control(case):
    alpha, beta, gamma = (case[k] for k in ("alpha", "beta", "gamma"))
    m0, m1 = (float(Fraction(case[k])) for k in ("m0", "m1"))
    omega = np.sqrt(20)
    s = 1j * omega
    D = s**2 + (1 + gamma) * s + beta + gamma
    amplitude = case["x_amplitude"] * np.array([1, (s + gamma) / D, -beta / D])
    duration = case["test_periods"] * 2 * np.pi / omega
    sol = integrate(alpha, beta, gamma, m0, m1, amplitude.real, duration, "DOP853")
    times = np.linspace(0, duration, 2401)
    exact = (amplitude[:, None] * np.exp(s * times)).real
    errors = np.max(abs(sol.sol(times) - exact), axis=1)
    return {
        "start": amplitude.real.tolist(),
        "period": 2 * np.pi / omega,
        "max_errors_xyz": errors.tolist(),
        "passed": bool(max(errors) < case["test_waveform_max_error"]),
    }


def scaling_controls(case):
    alpha, gamma = case["test_alpha"], case["gamma"]
    kappa = float(Fraction(case["kappa"]))
    beta = kappa * alpha**2
    m0, m1 = (float(Fraction(case["bistable_atlas_diode"][k])) for k in ("m0", "m1"))
    eq = equilibrium(beta, gamma, m0, m1)
    rows = []
    for start in case["test_initial_states"]:
        sol = integrate(alpha, beta, gamma, m0, m1, start, case["test_duration"])
        errors = [np.max(abs(sol.y[:, -1] - sign * eq)) for sign in (-1, 0, 1)]
        rows.append(
            {
                "start": start,
                "equilibrium_sign": int(np.argmin(errors)) - 1,
                "terminal_max_error": float(min(errors)),
                "passed": bool(min(errors) < case["test_terminal_max_norm"]),
            }
        )
    return {"cutoff": 1 / (2 * np.sqrt(kappa) - 1), "beta": beta, "runs": rows}


def scout_overlap():
    with np.load(ROOT / "results/chua-scaling/scout.npz") as archive:
        p = archive["parameters"]
        alpha, kappa = p[:, 0], p[:, 1]
        beta = kappa * alpha**2
        covered = np.where(alpha > 1, beta > (alpha + 1) ** 2 / 4, beta >= alpha)
        return {
            "total_runs": len(p),
            "covered_runs": int(covered.sum()),
            "covered_parameter_tuples": len(np.unique(p[covered], axis=0)),
            "covered_preliminary_candidates": int(
                sum(covered & archive["preliminary_candidate"])
            ),
        }


def figure(path, trace):
    fig, (plane, time) = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    alpha = np.linspace(0, 13, 500)
    boundary = np.where(alpha <= 1, alpha, (alpha + 1) ** 2 / 4)
    plane.fill_between(alpha, boundary, 46, color="#c4e3d6")
    plane.plot(alpha, boundary, color="#23614d", lw=1.8)
    plane.text(1, 35, "Every bounded orbit\napproaches equilibrium", color="#194b3d")
    plane.plot(9, 25, "o", color="#a43b32", ms=5)
    plane.annotate(
        "neutral ellipse", (9, 25), (9.2, 15), arrowprops={"arrowstyle": "-"}
    )
    plane.set(
        xlim=(0, 13),
        ylim=(0, 46),
        xlabel=r"$\alpha$",
        ylabel=r"$\beta$",
        title=r"All diode laws, $\gamma=0$",
    )
    times, voltage = trace
    time.axhline(2, color="#23614d", ls="--", lw=1)
    turns = np.diff(np.sign(np.diff(voltage)))
    for mask in (turns < 0, turns > 0):
        indices = np.flatnonzero(mask) + 1
        time.plot(times[indices], voltage[indices], color="#305486", lw=1.5)
    time.set(
        xlabel="Time",
        ylabel="Voltage at turning points",
        title=r"One well, $\alpha=9,\ \beta=26$",
        xlim=(0, times[-1]),
    )
    for ax in (plane, time):
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--classify",
        nargs=3,
        metavar=("ALPHA", "BETA", "GAMMA"),
        help="Classify exact decimals or fractions and print a periodic witness when available.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "outputs/chua-global-relaxation"
    )
    args = parser.parse_args()
    if args.classify:
        try:
            params = [sp.Rational(Fraction(value)) for value in args.classify]
            result = classify(*params)
        except (ValueError, ZeroDivisionError) as exc:
            parser.error(str(exc))
        print(json.dumps(result, indent=2))
        return
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    forecast = json.loads(FORECAST.read_text(encoding="utf-8"))
    symbolic_checks()
    wells, trace = well_controls(forecast["piecewise_linear"])
    ellipse = sharpness_control(forecast["sharpness"])
    scaling = scaling_controls(forecast["quadratic_path"])
    results = {
        "symbolic_identities_verified": True,
        "well_controls": wells,
        "ellipse": ellipse,
        "quadratic_path": scaling,
        "archived_scout": scout_overlap(),
    }
    (out / "checks.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    figure(out / "relaxation.png", trace)
    if not all(r["passed"] for r in [*wells, ellipse, *scaling["runs"]]):
        raise SystemExit("A forecast check failed. See checks.json.")


if __name__ == "__main__":
    main()
