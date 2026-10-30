"""Exact three-state collective-motor core, not a new oscillator model.

Oriola et al. (2014), EPL 107, 18002, equations (7) and (10),
DOI 10.1209/0295-5075/107/18002. Also the exact closed core of
Alouges et al. (2025), DOI 10.1142/S0218202525500423, equation (33).

State: (h, Y, Z). Parameters: (rho, r), both positive; r = mu*Omega.
Time is scaled by the constant molecular transition rate. The module provides
equations, coordinate checks and opt-in independent CPU periodic controls.
See docs/motor-pair-source-audit.md for physical and priority limitations.
"""

import numpy as np


SOURCE_CONTROLS = (
    {"name": "oriola_2014_figure_2", "rho": .05, "mu": 200.,
     "Omega": .2, "r": 40., "reported": "complex periodic orbit"},
    {"name": "oriola_2014_figure_3a", "rho": .7, "mu": 1000.,
     "Omega": .4, "r": 400., "reported": "periodic cusp-like orbit"},
    {"name": "oriola_2014_figure_3b", "rho": 4.5, "mu": 1000.,
     "Omega": .4, "r": 400., "reported": "periodic orbit"},
    {"name": "oriola_2014_figure_4a_left", "rho": 6., "mu": 63096.,
     "Omega": .425, "r": 26815.8, "reported": "subharmonic periodic orbit"},
    {"name": "oriola_2014_figure_4a_right", "rho": 3., "mu": 6500.,
     "Omega": .4, "r": 2600., "reported": "subharmonic periodic orbit"},
    {"name": "oriola_2014_figure_4d_left", "rho": 6., "mu": 19000.,
     "Omega": .4, "r": 7600., "reported": "subharmonic periodic orbit"},
)


def field_jac(state, parameters):
    """Return the normalized vector field and its analytic Jacobian."""
    h, y, z = state
    rho, r = parameters
    f = np.array([(-rho + r*(1-z))*h-r*y,
                  h*(1-z)-y, h*y-z])
    jac = np.array([[-rho+r*(1-z), -r, -r*h],
                    [1-z, -1, -h], [y, h, -1]])
    return f, jac


def rhs(t, state, parameters):
    """SciPy-compatible autonomous right-hand side."""
    return field_jac(state, parameters)[0]


def oriola_field(state, rho, mu, Omega):
    """Independent equation (7), ordered (a1_bar, b1_bar, X_bar)."""
    a, b, x = state
    return np.array([-a-2*np.pi*rho*x*b-mu*a*b,
                     -Omega-b+2*np.pi*rho*x*a+mu*a*a,
                     -rho*x-mu*a/(2*np.pi)])


def from_oriola(state, rho, mu, Omega):
    """Map equation (7) coordinates to (h,Y,Z); requires Omega>0."""
    if Omega <= 0:
        raise ValueError("Omega must be positive for this normalization")
    a, b, x = state
    return np.array([2*np.pi*rho*x+mu*a, a/Omega, (b+Omega)/Omega])


def to_oriola(state, rho, r, Omega):
    """Inverse map; requires rho>0, Omega>0. mu=r/Omega."""
    if rho <= 0 or Omega <= 0:
        raise ValueError("rho and Omega must be positive")
    h, y, z = state
    return np.array([Omega*y, Omega*(z-1), (h-r*y)/(2*np.pi*rho)])


def original_field(state, parameters, Omega=.2):
    """Reconstruct the normalized derivative independently from equation (7)."""
    rho, r = parameters
    mu = r/Omega
    a, b, x = to_oriola(state, rho, r, Omega)
    da, db, dx = oriola_field((a, b, x), rho, mu, Omega)
    return np.array([2*np.pi*rho*dx+mu*da, da/Omega, db/Omega])


def physical_diagnostics(state, parameters):
    """Return invariant-disk radius squared and displacement-memory coordinate.

    The set radius_squared<=1, |memory|<=r is positively invariant.
    A point outside it is not automatically an unphysical density: admissible
    source densities can have transients outside the smaller attracting disk.
    """
    h, y, z = state
    _, r = parameters
    return {"radius_squared": y*y+(1-z)*(1-z), "memory": h-r*y}


def source_phase_diagram_start(rho, mu, Omega):
    """Published Fig.5 state (a1_bar,b1_bar,X_bar)=(0,0,.074)."""
    return from_oriola((0., 0., .074), rho, mu, Omega)


def equilibrium_density_start(rho, displacement=.01):
    """Our explicit perturbation, not a source-reported Fig.2 initial state."""
    return np.array([2*np.pi*rho*displacement, 0., 0.])


def source_algebra_audit(samples=100):
    """Reproducible independent source-coordinate and variational checks.

    These are field identities, not simulated attractor or novelty evidence.
    Errors use |difference|/(1+|reference|), apart from stated matrix maxima.
    """
    rng = np.random.default_rng(70)
    jacobian_error = field_error = inverse_error = 0.
    for _ in range(samples):
        p = np.array([10**rng.uniform(-2, 2), 10**rng.uniform(0, 4)])
        y = rng.normal(size=3)
        f, jac = field_jac(y, p)
        complex_jac = np.column_stack([
            np.imag(field_jac(y.astype(complex)+1e-28j*np.eye(3)[i], p)[0])/1e-28
            for i in range(3)])
        jacobian_error = max(jacobian_error, float(np.max(
            np.abs(jac-complex_jac)/(1+np.abs(jac)))))
        field_error = max(field_error, float(np.max(
            np.abs(f-original_field(y, p))/(1+np.abs(f)))))
        back = from_oriola(to_oriola(y, *p, .2), p[0], p[1]/.2, .2)
        inverse_error = max(inverse_error, float(np.max(
            np.abs(back-y)/(1+np.abs(y)))))

    rng = np.random.default_rng(71)
    core_error = identity_error = 0.
    for _ in range(samples):
        a, ell, k, eta, potential = np.exp(rng.uniform(-1, 1, 5))
        a1, b1 = -.2*a, .1*a
        zeta = 2*np.pi*k/(eta*ell)
        lam = 2*np.pi**2*potential/(eta*ell)
        pc, ps, qc, qs, x = rng.normal(size=5)
        phase_velocity = zeta*x+lam*(ps-qs)/2
        # Literal five-component source Eq.(33), before any reduction.
        f5 = np.array([-a*pc+phase_velocity*ps+a1/ell,
                       -a*ps-phase_velocity*pc+b1/ell,
                       -a*qc-phase_velocity*qs+a1/ell,
                       -a*qs+phase_velocity*qc+b1/ell,
                       -ell*phase_velocity/(2*np.pi)])
        c, s = (pc+qc)/2, (ps-qs)/2
        dc, ds = (f5[0]+f5[2])/2, (f5[1]-f5[3])/2
        Omega, rho, r = -a1/a, k/(eta*a), -lam*a1/(a*a*ell)
        state = np.array([phase_velocity/a, ell*s/Omega, (ell*c+Omega)/Omega])
        mapped = np.array([(zeta*f5[4]+lam*ds)/(a*a),
                           ell*ds/(a*Omega), ell*dc/(a*Omega)])
        f, _ = field_jac(state, (rho, r))
        core_error = max(core_error, float(np.max(
            np.abs(f-mapped)/(1+np.abs(f)))))
        h, y, z = state
        q, dq = h-r*y, f[0]-r*f[1]
        calculated = q*dq/(rho*r)+y*f[1]+z*f[2]
        expected = -q*q/r+(r-1)*y*y-z*z
        identity_error = max(identity_error, float(
            abs(calculated-expected)/(1+abs(expected))))
    return {
        "samples": samples, "seed": 70, "alouges_and_identity_seed": 71,
        "error_normalization": "abs(actual-reference)/(1+abs(reference))",
        "complex_step": 1e-28,
        "oriola_parameters": "rho=10**U(-2,2), r=10**U(0,4), Omega=.2; normal states",
        "alouges_parameters": "positive scales exp(U(-1,1)); a1=-.2a,b1=.1a; normal states",
        "jacobian_complex_step_scaled_error": jacobian_error,
        "independent_eq7_field_scaled_error": field_error,
        "inverse_map_scaled_error": inverse_error,
        "alouges_full5_to_core_scaled_error": core_error,
        "lyapunov_identity_scaled_error": identity_error,
    }


def periodic_control(transient_time=2000.):
    """Independently integrate and close the published Fig.2 parameter case.

    The initial displacement is ours. Shooting tests the reached cycle, not
    every basin or a theorem that this parameter point has only one attractor.
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import root

    p = np.array([.05, 40.])
    initial = equilibrium_density_start(p[0])
    f = lambda t, y: rhs(t, y, p)
    j = lambda t, y: field_jac(y, p)[1]
    times = np.linspace(0., 20., 4001)
    short = {}
    trajectories = []
    for method in ("DOP853", "Radau"):
        kw = {"jac": j} if method == "Radau" else {}
        sol = solve_ivp(f, (0., 20.), initial, method=method, t_eval=times,
                        rtol=2e-11, atol=2e-12, max_step=.02, **kw)
        if not sol.success:
            raise RuntimeError(sol.message)
        short[method] = {"nfev": sol.nfev, "njev": sol.njev}
        trajectories.append(sol.y)
    short["maximum_coordinate_difference"] = np.max(
        np.abs(trajectories[0]-trajectories[1]), axis=1).tolist()

    def section(t, y):
        return y[0]
    section.direction = 1
    sol = solve_ivp(f, (0., transient_time), initial, method="DOP853",
                    rtol=2e-11, atol=2e-12, max_step=.05, events=section)
    if not sol.success or len(sol.t_events[0]) < 5:
        raise RuntimeError("Long control failed or did not provide section returns")
    hits = sol.y_events[0]
    periods = np.diff(sol.t_events[0])
    guess = np.array([hits[-1, 1], hits[-1, 2], periods[-1]])

    def augmented(t, w):
        fv, jac = field_jac(w[:3], p)
        return np.r_[fv, (jac @ w[3:].reshape(3, 3)).ravel()]

    def shoot(u, method="DOP853", max_step=.01):
        start = np.array([0., u[0], u[1]])
        w0 = np.r_[start, np.eye(3).ravel()]
        out = solve_ivp(augmented, (0., u[2]), w0, method=method,
                        rtol=2e-12, atol=2e-13, max_step=max_step)
        if not out.success:
            raise RuntimeError(out.message)
        end = out.y[:3, -1]
        monodromy = out.y[3:, -1].reshape(3, 3)
        residual = end-start
        derivative = np.column_stack([(monodromy-np.eye(3))[:, 1:], f(0., end)])
        return residual, derivative, monodromy, end

    closed = root(lambda u: shoot(u)[0], guess,
                  jac=lambda u: shoot(u)[1], tol=1e-10)
    residual, _, matrix, end = shoot(closed.x)
    independent, _, matrix_radau, _ = shoot(closed.x, "Radau", .005)
    start = np.r_[0., closed.x[:2]]
    velocity = f(0., end)
    projection = np.eye(3)-np.outer(velocity, np.array([1., 0., 0.]))/velocity[0]
    poincare = (projection @ matrix)[1:, 1:]
    disk = sol.y[1]**2+(1-sol.y[2])**2
    memory = sol.y[0]-p[1]*sol.y[1]
    return {
        "parameters": {"rho": p[0], "r": p[1]},
        "start": initial.tolist(),
        "start_provenance": "Our displacement 0.01, equilibrium motor density; not printed Fig.2 initial state",
        "short_control": short,
        "long_control": {"time": transient_time, "nfev": sol.nfev,
                         "section_hits": len(hits),
                         "last_ten_section_intervals": periods[-10:].tolist(),
                         "last_return_max_difference": float(np.max(np.abs(hits[-1]-hits[-2]))),
                         "max_disk_radius_squared_at_solver_nodes": float(np.max(disk)),
                         "max_abs_memory_at_solver_nodes": float(np.max(np.abs(memory)))},
        "shooting": {"optimizer_success": bool(closed.success),
                     "optimizer_message": str(closed.message),
                     "section_state": start.tolist(), "period": float(closed.x[2]),
                     "dop853_closure_max_abs": float(np.max(np.abs(residual))),
                     "radau_closure_max_abs": float(np.max(np.abs(independent))),
                     "monodromy_max_solver_difference": float(np.max(np.abs(matrix-matrix_radau))),
                     "monodromy": matrix.tolist(),
                     "poincare_matrix": poincare.tolist(),
                     "poincare_eigenvalues": [[float(v.real), float(v.imag)] for v in np.linalg.eigvals(poincare)],
                     "flow_identity_max_abs": float(np.max(np.abs(matrix @ f(0., start)-velocity))),
                     "section_velocity": float(velocity[0]),
                     "scope": "Local stability of the reached numerically closed cycle; not uniqueness across basins"},
    }


def refine_periodic_control(result):
    """Repeat shooting more tightly in the independent memory coordinate.

    The coordinate change q=h-rY is exact. It is used for a conditioning
    check, not as a model modification or a replacement of physical variables.
    """
    from scipy.integrate import solve_ivp
    rho, r = result["parameters"]["rho"], result["parameters"]["r"]
    transform = np.array([[1., r, 0.], [0., 1., 0.], [0., 0., 1.]])
    inverse = np.array([[1., -r, 0.], [0., 1., 0.], [0., 0., 1.]])
    scale = np.array([max(1., r), 1., 1.])
    guess = np.r_[result["shooting"]["section_state"][1:], result["shooting"]["period"]]

    def aug(t, w):
        q, y, z = w[:3]
        h = q+r*y
        f = np.array([-rho*h, h*(1-z)-y, h*y-z])
        j = np.array([[-rho, -rho*r, 0.],
                      [1-z, r*(1-z)-1, -h], [y, q+2*r*y, -1.]])
        return np.r_[f, (j @ w[3:].reshape(3, 3)).ravel()]

    def shoot(u, method="DOP853"):
        start = np.array([0., u[0], u[1]])
        out = solve_ivp(aug, (0., u[2]), np.r_[inverse @ start, np.eye(3).ravel()],
                        method=method, rtol=3e-14, atol=1e-15, max_step=.01)
        if not out.success:
            raise RuntimeError(out.message)
        end = transform @ out.y[:3, -1]
        monodromy = transform @ out.y[3:, -1].reshape(3, 3) @ inverse
        residual = end-start
        velocity = field_jac(end, (rho, r))[0]
        derivative = np.column_stack([(monodromy-np.eye(3))[:, 1:], velocity])
        return residual, derivative, monodromy, end

    # Bounded Newton iteration: avoid an unbounded solver chasing integration
    # roundoff on a long, highly conditioned relaxation orbit.
    solution = guess.copy()
    iteration_residuals = []
    for _ in range(2):
        residual, derivative, _, _ = shoot(solution)
        iteration_residuals.append(float(np.max(np.abs(residual))))
        correction = np.linalg.solve(derivative/scale[:, None], residual/scale)
        solution -= correction
    resid, _, matrix, end = shoot(solution)
    independent, _, matrix_radau, end_radau = shoot(solution, "Radau")
    start = np.r_[0., solution[:2]]
    def projected(matrix, end):
        velocity = field_jac(end, (rho, r))[0]
        projection = np.eye(3)-np.outer(velocity, [1., 0., 0.])/velocity[0]
        return (projection @ matrix)[1:, 1:]
    poincare = projected(matrix, end)
    poincare_radau = projected(matrix_radau, end_radau)
    def eigs(m):
        return [[float(v.real), float(v.imag)] for v in np.linalg.eigvals(m)]
    result["refined_shooting"] = {
        "coordinates": "q=h-rY, Y, Z; results transformed back to h,Y,Z",
        "rtol": 3e-14, "atol": 1e-15, "max_step": .01,
        "method": "Two bounded Newton corrections; no solver-success flag used as an accuracy claim",
        "iteration_residuals_max_abs": iteration_residuals,
        "section_state": start.tolist(), "period": float(solution[2]),
        "dop853_closure_max_abs": float(np.max(np.abs(resid))),
        "radau_closure_max_abs": float(np.max(np.abs(independent))),
        "monodromy_max_solver_difference": float(np.max(np.abs(matrix-matrix_radau))),
        "poincare_matrix": poincare.tolist(),
        "poincare_eigenvalues": eigs(poincare),
        "radau_poincare_eigenvalues": eigs(poincare_radau),
        "flow_identity_max_abs": float(np.max(np.abs(
            matrix @ field_jac(start, (rho, r))[0]-field_jac(end, (rho, r))[0]))),
        "scope": "Numerical local stability of this reached periodic branch; not uniqueness or absence of chaos elsewhere",
    }
    return result


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-output", type=Path,
                        help="Run the independent periodic control and save JSON")
    parser.add_argument("--transient-time", type=float, default=2000.)
    parser.add_argument("--refine-control", type=Path,
                        help="Refine an existing control JSON in the exact memory coordinate")
    parser.add_argument("--audit-output", type=Path,
                        help="Run only the reproducible source-algebra/Jacobian audit")
    args = parser.parse_args()
    if args.audit_output:
        result = source_algebra_audit()
        args.audit_output.parent.mkdir(parents=True, exist_ok=True)
        args.audit_output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
        print(json.dumps(result, indent=2))
    elif args.control_output or args.refine_control:
        if args.refine_control:
            result = refine_periodic_control(json.loads(args.refine_control.read_text(encoding="utf-8")))
            if not args.control_output:
                args.control_output = args.refine_control
        else:
            result = periodic_control(args.transient_time)
        args.control_output.parent.mkdir(parents=True, exist_ok=True)
        args.control_output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
        print(json.dumps(result, indent=2))
    else:
        parser.print_help()
