"""Algebra and a few modal Floquet controls for the exact quadratic-drag pair.

No 4D attractor trajectories or parameter grid are integrated. Time is the
canonical sigma of docs/late-physical-candidates-2.md.
"""
from pathlib import Path
import datetime as dt
import json
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]


def mode_control(mu, kappa, symmetric):
    base2, transverse2 = ((1 + 2*kappa, 1.) if symmetric else (1., 1 + 2*kappa))

    def rhs(t, y):
        x, v = y
        return [v, -base2*x + mu*(1-abs(v))*v]

    def event(t, y):
        return y[0]

    event.direction = 1
    event.terminal = True

    def return_map(amplitude, full=False):
        skip = solve_ivp(rhs, (0, 1e-6), [0., amplitude], method="DOP853", rtol=1e-12, atol=1e-14)
        sol = solve_ivp(rhs, (1e-6, 200), skip.y[:, -1], method="DOP853", rtol=1e-12,
                        atol=1e-14, max_step=.04, events=event)
        if not sol.success or len(sol.t_events[0]) != 1:
            raise RuntimeError("Modal section did not return")
        return (float(sol.t_events[0][0]), sol.y_events[0][0]) if full else float(sol.y_events[0][0][1]-amplitude)

    amplitude = brentq(return_map, .2, 3., xtol=1e-12)
    period, returned = return_map(amplitude, full=True)

    def augmented(t, y):
        x, v = y[:2]
        damping = mu*(2*abs(v)-1)
        jl = np.array([[0., 1.], [-base2, -damping]])
        jt = np.array([[0., 1.], [-transverse2, -damping]])
        return np.r_[rhs(t, y[:2]), (jl @ y[2:6].reshape(2, 2)).ravel(),
                     (jt @ y[6:10].reshape(2, 2)).ravel(), damping]

    initial = np.r_[0., amplitude, np.eye(2).ravel(), np.eye(2).ravel(), 0.]
    outputs = {}
    for method in ("DOP853", "Radau"):
        sol = solve_ivp(augmented, (0, period), initial, method=method,
                        rtol=1e-11, atol=1e-13, max_step=.025)
        if not sol.success:
            raise RuntimeError(sol.message)
        end = sol.y[:, -1]
        ml = np.linalg.eigvals(end[2:6].reshape(2, 2))
        mt = np.linalg.eigvals(end[6:10].reshape(2, 2))
        outputs[method] = {"endpoint": end.tolist(), "longitudinal_multipliers": [[float(z.real), float(z.imag)] for z in ml],
                           "transverse_multipliers": [[float(z.real), float(z.imag)] for z in mt],
                           "transverse_exponents": (np.log(abs(mt))/period).tolist(),
                           "maximum_transverse_multiplier": float(max(abs(mt))),
                           "mean_damping": float(end[-1]/period),
                           "determinant_trace_error": float(abs(np.linalg.det(end[6:10].reshape(2, 2))-np.exp(-end[-1]))),
                           "flow_multiplier_error": float(min(abs(ml-1)))}
    return {"mu": mu, "kappa": kappa, "mode": "symmetric" if symmetric else "antisymmetric",
            "base_frequency_squared": base2, "transverse_frequency_squared": transverse2,
            "section_velocity": amplitude, "period_sigma": period,
            "section_closure": float(abs(returned[1]-amplitude)),
            "dop853_radau_full_augmented_difference": float(np.max(abs(np.array(outputs['DOP853']['endpoint'])-outputs['Radau']['endpoint']))),
            "integrators": outputs}


def algebra_checks():
    rng = np.random.default_rng(2026103013)
    maximum, violation = 0., -float("inf")
    for _ in range(256):
        mu, kappa = np.exp(rng.uniform(-3, 2, 2))
        eps = min(.1, .25/mu)
        stiffness = np.eye(2)+kappa*np.ones((2, 2))
        y, z = rng.normal(size=(2, 2))*10**rng.uniform(-2, 2)
        r = np.linalg.norm(y)
        h = mu*(abs(y)-1)*y
        phi = y/np.sqrt(1+r*r)
        derivative = np.eye(2)/np.sqrt(1+r*r)-np.outer(y,y)/(1+r*r)**1.5
        yp, zp = z-h, -stiffness@y
        direct = y@yp + z@np.linalg.solve(stiffness, zp)-eps*(yp@derivative@z+phi@zp)
        identity = -y@h-eps*z@derivative@z+eps*h@derivative@z+eps*phi@stiffness@y
        bound = (mu+eps*mu*mu/np.sqrt(2))*r*r-mu*r**3/(2*np.sqrt(2)) + eps*(1+2*kappa)*r-eps*np.dot(z,z)/(2*(1+r*r)**1.5)
        maximum = max(maximum, abs(direct-identity)/(1+abs(direct)+abs(identity)))
        violation = max(violation, (identity-bound)/(1+abs(identity)+abs(bound)))
    assert maximum < 1e-12 and violation <= 1e-12
    return {"seed": 2026103013, "samples": 256, "max_scaled_modified_energy_identity_error": float(maximum),
            "max_scaled_upper_bound_violation": float(violation), "purpose": "Algebra transcription checks, not a numerical proof of the analytic bound."}


def main():
    started = time.perf_counter()
    mu0, kappa0 = .00880947943, .557378006
    controls = [mode_control(mu0, kappa0, mode) for mode in (True, False)]
    controls += [mode_control(.2, 0., True)]
    controls += [mode_control(mu, 1.5, mode) for mu in (1., 3.) for mode in (True, False)]
    output = {"audit_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "equations": "x_i''=-x_i-kappa*(x1+x2)+mu*(1-abs(x_i'))*x_i'",
              "scope": "Seven modal orbit/Floquet controls only; no off-modal 4D attractor or parameter survey.",
              "primary_model_doi": "10.3390/polym15163349", "time_unit": "canonical sigma",
              "documentation": "docs/rayleigh-pair-structure.md", "script": "research/audit_rayleigh_pair.py",
              "source_control_scope": "Source-derived mu,kappa rounded as in the screening note; complete original physical initial positions were not supplied.",
              "algebra_checks": algebra_checks(), "modal_controls": controls,
              "wall_seconds": time.perf_counter()-started}
    path = ROOT/'results/candidate-audits/rayleigh-pair-structure.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2)+'\n', encoding='utf8')
    for c in controls:
        print(c['mu'], c['kappa'], c['mode'], 'T', c['period_sigma'],
              'rho', c['integrators']['DOP853']['maximum_transverse_multiplier'],
              'exponents', c['integrators']['DOP853']['transverse_exponents'], flush=True)
    print('Wrote', path, 'wall', output['wall_seconds'], flush=True)


if __name__ == '__main__':
    main()
