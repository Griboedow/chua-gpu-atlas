"""Small independent controls for beta=kappa*alpha^2 in canonical Chua.

No parameter survey or finite-alpha chaos classification. Only exact identities,
linear eigenvalues, and short scalar/affine controls are computed.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[1]


def g(x, m0, m1):
    return (1+m1)*x+.5*(m0-m1)*(abs(x+1)-abs(x-1))


def original_field(q, alpha, kappa, m0, m1):
    x, y, z = q
    return np.array([alpha*(y-g(x, m0, m1)), x-y+z, -kappa*alpha*alpha*y])


def scaled_field(q, eps, kappa, m0, m1):
    x, y, Z = q
    return np.array([y-g(x, m0, m1), Z+eps*(x-y), -kappa*y])


def branch_matrix(a, eps, kappa):
    return np.array([[-a, 1., 0.], [eps, -eps, 1.], [0., -kappa, 0.]])


def branch_polynomial(a, eps, kappa):
    return np.array([1., a+eps, kappa+eps*(a-1), a*kappa])


def affine_fiber_control(slope, offset, amplitude, kappa, m0, m1, branch):
    omega = np.sqrt(kappa)
    period = 2*np.pi/omega
    equilibrium = -offset/slope
    initial = equilibrium+amplitude*slope/(slope*slope+kappa)

    def rhs(s, x):
        return [amplitude*np.cos(omega*s)-g(x[0], m0, m1)]

    sol = solve_ivp(rhs, (0, period), [initial], method='DOP853',
                    rtol=1e-12, atol=1e-14, max_step=.02, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    times = np.linspace(0, period, 8193)
    y = amplitude*np.cos(omega*times)
    exact = equilibrium+amplitude*(slope*np.cos(omega*times)+omega*np.sin(omega*times))/(slope*slope+kappa)
    actual = sol.sol(times)[0]
    yy = np.trapezoid(y*y, times)/period
    xy = np.trapezoid(actual*y, times)/period
    return {"branch": branch, "slope_g": slope, "offset_g": offset,
            "forcing_amplitude": amplitude, "kappa": kappa, "period_fast_s": float(period),
            "initial_x": float(initial), "x_min": float(actual.min()), "x_max": float(actual.max()),
            "max_trajectory_error_against_affine_solution": float(max(abs(actual-exact))),
            "one_period_closure": float(abs(actual[-1]-initial)),
            "measured_covariance_ratio": float(xy/yy), "exact_covariance_ratio": float(slope/(slope*slope+kappa)),
            "scalar_multiplier_exact": float(np.exp(-slope*period)),
            "scalar_exponent_fast_s_exact": -slope,
            "interpretation": "Positive exponent on an unstable periodic orbit is not a chaotic attractor" if slope < 0 else "Attracting scalar periodic response on one affine outer branch"}


def main():
    rng = np.random.default_rng(2026103018)
    scaling_errors, action_errors, polynomial_errors, boundary_factor_errors = [], [], [], []
    for _ in range(256):
        alpha = np.exp(rng.uniform(np.log(2), np.log(1e5)))
        eps = 1/alpha
        kappa = rng.uniform(.02, 2.)
        m0, m1 = rng.uniform(-8, -1.01), rng.uniform(-.95, -.05)
        q = rng.uniform(-5, 5, 3)
        T = np.array([1., 1., 1/alpha])
        transformed = T*original_field(q/T, alpha, kappa, m0, m1)/alpha
        actual = scaled_field(q, eps, kappa, m0, m1)
        scaling_errors.append(max(abs(transformed-actual)/(1+abs(actual))))
        gradient = np.array([0., kappa*q[1], q[2]])
        action_errors.append(abs(gradient@actual-eps*kappa*q[1]*(q[0]-q[1])))
        a = 1+m1
        polynomial_errors.append(max(abs(np.poly(branch_matrix(a, eps, kappa))-branch_polynomial(a, eps, kappa))))
        critical = (1-a)*(a+eps)
        expected = np.polymul([1., a+eps], [1., 0., a*(1-a)])
        boundary_factor_errors.append(max(abs(branch_polynomial(a, eps, critical)-expected)))

    m0, m1 = -8/7, -5/7
    a, inner = 1+m1, 1+m0
    kappa0 = a*(1-a)
    eigen_controls = []
    for alpha in [32, 128, 512, 2048, 8192, 32768]:
        eps = 1/alpha
        eig = np.linalg.eigvals(branch_matrix(a, eps, kappa0))
        pair = eig[np.abs(eig.imag) > 1e-8]
        growth = float(pair[0].real*alpha)
        eigen_controls.append({"alpha": alpha, "kappa_fixed": kappa0,
                               "eigenvalues_fast_s_real_imag": [[float(v.real), float(v.imag)] for v in eig],
                               "outer_pair_growth_original_t": growth,
                               "alpha_times_growth_original_t": alpha*growth,
                               "asymptotic_alpha_times_growth": (1-a)/(2*a),
                               "exact_kappa_stability_boundary": (1-a)*(a+eps)})

    alpha = 64.
    eps = 1/alpha
    kc = (1-a)*(a+eps)
    J = branch_matrix(a, eps, kc)
    omega = np.sqrt(kappa0)
    period = 2*np.pi/omega
    xe = (m1-m0)/a
    equilibrium = np.array([xe, 0., -eps*xe])
    perturbation = .05*np.array([1., a, -kc])
    initial = equilibrium+perturbation
    grid = np.linspace(0, period, 1025)
    exact = np.array([equilibrium+expm(J*s)@perturbation for s in grid])
    sol = solve_ivp(lambda s, q: scaled_field(q, eps, kc, m0, m1), (0, period), initial,
                    method='DOP853', rtol=1e-12, atol=1e-14, max_step=.02, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    boundary_control = {"alpha": alpha, "kappa": kc, "initial_state_scaled": initial.tolist(),
                        "equilibrium_scaled": equilibrium.tolist(), "period_fast_s": float(period),
                        "x_min": float(exact[:, 0].min()), "x_max": float(exact[:, 0].max()),
                        "full_piecewise_field_vs_matrix_exponential_max_error": float(np.max(abs(sol.sol(grid).T-exact))),
                        "full_piecewise_field_closure_inf": float(np.max(abs(sol.y[:, -1]-initial))),
                        "interpretation": "Exact local center ellipses: affine outer field makes this linear Hopf boundary degenerate, not a generic isolated-limit-cycle bifurcation."}

    fiber_controls = [affine_fiber_control(a, m0-m1, .1, kappa0, m0, m1, 'positive outer'),
                      affine_fiber_control(inner, 0., .1, kappa0, m0, m1, 'inner unstable')]
    starts = np.array([-4., -2., -.5, 0., .5, 2., 4.])
    returns = []
    for start in starts:
        sol = solve_ivp(lambda s, x: [.1*np.cos(omega*s)-g(x[0], m0, m1)],
                        (0, period), [start], method='DOP853', rtol=1e-12, atol=1e-14, max_step=.03)
        if not sol.success:
            raise RuntimeError(sol.message)
        returns.append(float(sol.y[0, -1]))

    result = {
        "scope": "Exact scaling and analytic structure; small scalar/affine controls only, no finite-alpha chaos survey",
        "conventions": {"m0_m1": "Slopes of f, so slopes of g=x+f are 1+m0 and 1+m1",
                        "epsilon": "1/alpha", "fast_time": "s=alpha*t", "scaled_state": "(x,y,Z=z/alpha)",
                        "beta_path": "beta=kappa*alpha^2 with kappa>0 and m0,m1 fixed"},
        "checks": {"random_seed": 2026103018, "random_cases": 256,
                   "scaled_field_max_relative_error": float(max(scaling_errors)),
                   "action_identity_max_absolute_error": float(max(action_errors)),
                   "characteristic_polynomial_max_absolute_error": float(max(polynomial_errors)),
                   "critical_polynomial_factorization_max_absolute_error": float(max(boundary_factor_errors))},
        "standard_slopes": {"m0": m0, "m1": m1, "outer_g_slope": a, "inner_g_slope": inner,
                            "limiting_critical_kappa": kappa0, "growth_coefficient": (1-a)/(2*a)},
        "fixed_critical_kappa_eigenvalue_controls": eigen_controls,
        "affine_boundary_center_control": boundary_control,
        "affine_fiber_controls": fiber_controls,
        "scalar_stroboscopic_order_control": {"forcing_amplitude": .1, "kappa": kappa0,
             "initial_x": starts.tolist(), "one_period_returns": returns,
             "smallest_adjacent_return_difference": float(min(np.diff(returns))),
             "qualification": "Finite numerical illustration; the order theorem follows analytically from uniqueness"},
        "fiber_inner_residence_thresholds": [
            {"m0": slope, "m1": m1, "critical_inner_fraction": (1+m1)/(m1-slope),
             "identity": "log(P_R_prime)/T = -[(1+m1)+(m0-m1)*inner_fraction]",
             "qualification": "Regular frozen scalar periodic response only; not a full finite-epsilon Lyapunov formula. Uniform attraction needs a strict uniform gap below the threshold."}
            for slope in [-2., -5., -10., -30.]],
        "exact_limit_result": "At epsilon=0, each bounded scalar stroboscopic orbit at fixed oscillator amplitude converges to a fixed point. No recurrent chaotic fiber attractor; unstable periodic responses can have positive exponents.",
        "finite_epsilon_result": "Unresolved globally. No theorem excluding compact chaotic attractors at arbitrarily large finite alpha is proved here.",
        "normal_hyperbolicity_caveat": "A smooth uniformly attracting periodic fiber family suggests a persistent two-dimensional amplitude-phase manifold locally. Regularity of the actual piecewise-linear return map, trapping and uniformity must be established before invoking a persistence theorem.",
        "literature": [
            {"authors": "Chua, Komuro and Matsumoto", "year": 1986, "title": "The Double Scroll Family",
             "doi": "10.1109/TCS.1986.1085869", "primary_url": "https://people.eecs.berkeley.edu/~chua/papers/Chua86.pdf",
             "inspected": "Primary full PDF available; introductory equations and Poincare/global-bifurcation scope inspected; bounded text search for asymptotic discussion",
             "relevance": "Classical piecewise-linear double-scroll equations and return-map analysis; no claim here of a new Chua mechanism"},
            {"author": "Bruno Rossetto", "year": 1993, "title": "Chua's circuit as a slow-fast autonomous dynamical system",
             "journal_reference": "Journal of Circuits, Systems and Computers 3, 483-496; reprinted in Chua's Circuit: A Paradigm for Chaos, pp.711-724",
             "chapter_doi": "10.1142/9789812798855_0037", "primary_abstract_url": "https://doi.org/10.1142/1997",
             "inspected": "Publisher abstract; full chapter retrieval failed",
             "relevance": "Explicit singular-perturbation prior art; exact beta=kappa*alpha^2 coverage not verified"},
            {"authors": "Euzebio and Llibre", "year": 2017, "title": "Zero-Hopf bifurcation in a Chua system",
             "doi": "10.1016/j.nonrwa.2017.02.002", "primary_url": "https://files.cercomp.ufg.br/weby/up/1076/o/Euzebio_Llibre_chua.pdf",
             "inspected": "Full primary author manuscript, model equation(1), rescaling and averaging hypotheses in section2",
             "relevance": "Averaging in a six-parameter cubic Chua system, not an established equivalence to this fixed-slope piecewise-linear limit"},
            {"authors": "Llibre, Mereu and Novaes", "year": 2015, "title": "Averaging theory for discontinuous piecewise differential systems",
             "doi": "10.1016/j.jde.2015.01.022", "primary_metadata_url": "https://portalrecerca.uab.cat/en/publications/averaging-theory-for-discontinuous-piecewise-differential-systems/",
             "inspected": "Author institution abstract/metadata; full text not recovered in this pass",
             "relevance": "Specialized nonsmooth averaging exists, but its hypotheses have not been checked for the proposed invariant-manifold claim"}],
        "novelty": "No originality or literature-absence claim; these are explicitly derived checks and a scoped research question",
        "reproduction": "python research/audit_chua_scaling.py"}
    target = ROOT/'results/chua-scaling/structure-audit.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({"output": str(target), "checks": result['checks'],
                      "growth_controls": eigen_controls, "boundary_control": boundary_control,
                      "fiber_controls": fiber_controls}, indent=2))


if __name__ == '__main__':
    main()
