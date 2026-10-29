"""Linear stability and independent trajectories of the assumed photoresponse law."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from photofeedback import parameter_at_equilibrium, rhs, stability


def main():
    x, b, k, width = .8, .184, 2.36, .2
    p, eq = parameter_at_equilibrium(x, b, 1., k, width)
    _, jac = stability(p)
    stiffness, rate0, gain, spatial = -jac[1, 0], -jac[2, 2], -jac[1, 2], jac[2, 0]
    polynomial = [b * stiffness, b * b * rate0 - gain * spatial, b * rate0 * rate0]
    roots = np.sort(np.roots(polynomial).real)
    boundaries = []
    for tau in roots:
        pp, _ = parameter_at_equilibrium(x, b, tau, k, width)
        eig, j = stability(pp)
        slope = []
        for factor in [1 - 1e-5, 1 + 1e-5]:
            p1 = pp.copy()
            p1[2] *= factor
            slope.append(stability(p1)[0].real.max())
        boundaries.append(dict(tau=float(tau), eigenvalues=[[float(z.real), float(z.imag)] for z in eig],
                               crossing_derivative=float((slope[1] - slope[0]) / (2e-5 * tau))))
    tau_grid = np.geomspace(.01, 10., 300)
    eigmax = []
    for tau in tau_grid:
        pp, _ = parameter_at_equilibrium(x, b, tau, k, width)
        eigmax.append(stability(pp)[0].real.max())
    trace = {}
    for tau in [.05, .5, 3.]:
        pp, eq = parameter_at_equilibrium(x, b, tau, k, width)
        sol = solve_ivp(lambda t, y: rhs(t, y, pp), (0, 1500), eq + [.01, .001, 0],
                        method='DOP853', rtol=1e-10, atol=1e-12,
                        t_eval=np.linspace(1200, 1500, 6001), max_step=.05)
        if not sol.success:
            raise RuntimeError(sol.message)
        trace[str(tau)] = dict(height_span=float(np.ptp(sol.y[0])), final_state=sol.y[:, -1].tolist())
        trace[f'state_{tau}'] = sol.y.T
    path = Path('results/photofeedback/stability')
    path.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path / 'curves.npz', tau=tau_grid, eigenvalue_real=eigmax, time=sol.t,
                        **{key: val for key, val in trace.items() if isinstance(val, np.ndarray)})
    report = dict(equilibrium_height=x, b=b, k=k, width=width,
                  characteristic_coefficients=dict(K=float(stiffness), r0=float(rate0), c=float(gain), d0=float(spatial)),
                  tau_polynomial=polynomial, linear_hopf_boundaries=boundaries,
                  trajectories={key: val for key, val in trace.items() if not isinstance(val, np.ndarray)},
                  caveat='Linear imaginary-axis crossings, not a computed first Lyapunov coefficient or a calibrated device prediction.')
    (path / 'analysis.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
