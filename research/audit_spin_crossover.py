"""Structural audit of literal Orlov et al. arXiv:2606.25506v1 equations.

No model correction, parameter survey or chaos claim is performed here.
Run from any directory with NumPy and SciPy installed.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "spin-crossover"


@dataclass(frozen=True)
class Parameters:
    J_s: float = 112.0
    J_tau: float = 0.0
    Delta_s: float = 309.5
    g_tau: float = 3.0
    Gamma: float = 0.05
    alpha: float = 4.5
    T_reservoir: float = 100.0
    intensity: float = 0.0
    rho: float = 0.5
    enthalpy_over_capacity: float = 100.0


def gibbs(n, m, temperature, p):
    """Return HS probability and the unconditioned spin mean for S=2.

    Finite spin sums exactly equal theta and S*B_S in Eq. 26;
    log scaling avoids overflow and the removable coth singularity at m=0.
    """
    if temperature <= 0:
        raise ValueError("The printed model requires positive temperature")
    q = p.J_s * n * m / temperature
    spin = np.arange(-2.0, 3.0)
    log_ls = (p.Delta_s - p.J_tau * (n - 0.5)) / temperature - np.log(p.g_tau)
    z = np.r_[log_ls, q * spin]
    weights = np.exp(z - np.max(z))
    weights /= np.sum(weights)
    return float(np.sum(weights[1:])), float(weights[1:] @ spin)


def rhs(_, y, p):
    """Eqs. 16 and 26 in the literal (n_HS, m, T) variables, J_s > 0."""
    n, m, temperature = y
    if n <= 0:
        raise ValueError("This audit stays on the regular n_HS > 0 branch")
    probability, spin_mean = gibbs(n, m, temperature, p)
    mdot = p.Gamma * p.J_s * n * (spin_mean - n * m)
    ndot = p.Gamma * p.J_tau * (probability - n) + m / n * mdot
    tdot = (-p.alpha * (temperature - p.T_reservoir)
            + p.intensity * (1.0 + (p.rho - 1.0) * n)
            - p.enthalpy_over_capacity * mdot)
    return np.array([ndot, mdot, tdot])


def photothermal_rhs(y, J_tau=152.0, gamma=1/19, log_g=7.0,
                     t_eq=112.6, alpha=2.8, t_reservoir=75.0, rho=0.5):
    """Literal Eqs. 27 and 16, with zero latent-heat term.

    Delta_s = T_eq*log(g), so the chosen stationary point is exact;
    the paper's separately rounded Delta_s=788 differs by 0.2 K.
    """
    tau, temperature = y
    delta = t_eq * log_g - temperature * log_g - J_tau * tau
    tau_dot = -gamma * J_tau * (tau + 0.5 * np.tanh(delta / (2 * temperature)))
    intensity = 2 * alpha * (t_eq - t_reservoir) / (1 + rho)
    tdot = -alpha * (temperature - t_reservoir) + intensity * (1 + (rho-1)*(tau+0.5))
    return np.array([tau_dot, tdot])


def audit():
    rng = np.random.default_rng(2026103001)
    p = Parameters()
    residual = 0.0
    for _ in range(512):
        n = rng.uniform(0.1, 1.0)
        m = rng.uniform(-1.9*n, 1.9*n)
        y = np.array([n, m, rng.uniform(40, 240)])
        f = rhs(0.0, y, p)
        residual = max(residual, abs(2*n*f[0] - 2*m*f[1]))

    # The nonzero-J_tau identity explicitly shows why the exclusion does not
    # automatically apply when both interactions are retained.
    combined = Parameters(J_tau=152.0)
    combined_residual = 0.0
    for _ in range(512):
        n = rng.uniform(0.1, 1.0)
        m = rng.uniform(-1.9*n, 1.9*n)
        y = np.array([n, m, rng.uniform(40, 240)])
        f = rhs(0.0, y, combined)
        prob, _ = gibbs(*y, combined)
        exact = 2 * n * combined.Gamma * combined.J_tau * (prob-n)
        combined_residual = max(combined_residual, abs(2*n*f[0]-2*m*f[1]-exact))

    y0 = np.array([0.8, 0.1, 105.0])
    invariant = y0[0]**2-y0[1]**2
    def reduced(_, z):
        m, temperature = z
        n = np.sqrt(invariant + m*m)
        return rhs(0.0, [n, m, temperature], p)[1:]
    full = solve_ivp(lambda t,y: rhs(t,y,p), (0, 20), y0,
                     method="DOP853", rtol=1e-12, atol=1e-14, max_step=.02, dense_output=True)
    red = solve_ivp(reduced, (0, 20), y0[1:], method="Radau",
                    rtol=1e-11, atol=1e-13, max_step=.02, dense_output=True)
    if not full.success or not red.success:
        raise RuntimeError("Independent control integration failed")
    times = np.linspace(0, 20, 2001)
    states = full.sol(times).T
    reduced_states = red.sol(times).T
    reconstructed = np.column_stack([np.sqrt(invariant+reduced_states[:,0]**2), reduced_states])
    drift = float(np.max(abs(states[:,0]**2-states[:,1]**2-invariant)))
    agreement = float(np.max(abs(states-reconstructed)))

    # A physical-domain counterexample, not a proposed initial state for chaos.
    boundary_y = [1.0, 1.9, 20.0]
    boundary = []
    for params in [p, combined]:
        f = rhs(0, boundary_y, params)
        boundary.append(dict(parameters=asdict(params), state=boundary_y, rhs=f.tolist(),
                             outward_from_population_upper_boundary=bool(f[0]>0)))

    j, gamma, teq, lng, alpha, tr, rho = 152.0, 1/19, 112.6, 7.0, 2.8, 75.0, .5
    analytical = np.array([[gamma*j*(j/(4*teq)-1), gamma*j*lng/(4*teq)],
                           [-2*alpha*(1-rho)/(1+rho)*(teq-tr), -alpha]])
    eps = 1e-5
    center = np.array([0.0, teq])
    finite = np.column_stack([(photothermal_rhs(center+eps*d)-photothermal_rhs(center-eps*d))/(2*eps)
                              for d in np.eye(2)])
    zero_alpha_trace = gamma*j*(j/(4*teq)-1)
    # This arithmetic illustrates a possible normalization source of the
    # discrepancy; it is NOT a corrected model used in any calculation.
    factor_four_diagnostic = gamma*j*(j/teq-1)
    critical_reservoir = teq*(1-2*(j/(4*teq)-1)/(((1-rho)/(1+rho))*lng))
    report = dict(
        scope="Literal arXiv:2606.25506v1 model; analytical screening, not a chaos search or priority claim",
        source_equations=["16", "26a", "26b", "27", "31/32 in PDF; 35/36 in HTML"],
        first_integral=dict(formula="n_HS^2 - m^2", condition="J_tau = 0, n_HS > 0",
                            random_states=512, derivative_residual_max=residual,
                            full_vs_reduced_error=agreement, full_control_drift=drift,
                            methods=["DOP853 in 3D", "Radau on exact 2D level surface"],
                            initial_state=y0.tolist(), parameters=asdict(p),
                            state_min=states.min(axis=0).tolist(),state_max=states.max(axis=0).tolist()),
        combined_identity=dict(formula="d(n_HS^2-m^2)/dt = 2*Gamma*J_tau*n_HS*(p_HS-n_HS)",
                               random_states=512, derivative_residual_max=combined_residual,
                               exclusion_applied=False),
        population_boundary_counterexamples=boundary,
        photothermal_parameter_check=dict(J_tau=j,Gamma=gamma,T_eq=teq,log_g=lng,alpha=alpha,
            analytical_jacobian=analytical.tolist(),finite_difference_jacobian=finite.tolist(),
            jacobian_error=float(np.max(abs(analytical-finite))),
            trace_at_alpha_zero=zero_alpha_trace,trace_at_reported_alpha=float(np.trace(analytical)),
            reported_Hopf_alpha=2.8,literal_critical_reservoir_temperature=critical_reservoir,
            reported_critical_reservoir_temperature=78.8,
            no_quarter_factor_diagnostic=factor_four_diagnostic,
            interpretation="Unresolved coefficient/parameter convention discrepancy; no correction applied"),
        validation_passed=bool(residual<1e-12 and combined_residual<1e-12 and drift<1e-10 and agreement<1e-8
                               and np.max(abs(analytical-finite))<1e-6),
    )
    source=ROOT/'outputs/candidate_research/orlov_2606.25506v1.pdf'
    if source.exists():report['source_pdf_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    np.savez_compressed(OUT/'regular-control.npz',times=times,states=states,reconstructed_states=reconstructed)
    print(json.dumps(report,indent=2,allow_nan=False))
    if not report['validation_passed']:raise RuntimeError('Algebra/control audit failed')


if __name__ == '__main__':
    audit()
