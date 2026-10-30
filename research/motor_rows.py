"""Exact first-harmonic N-row model of Alouges et al. (2025).

Source: arXiv:2409.03506v2, Eq.(15), constant total rate and cosine potential.
DOI: 10.1142/S0218202525500423. This is not a new model or chaos claim.
State is (c[0:N], s[0:N], X[0:N-1]); X[N-1]=-sum(X[0:N-1]).
Parameters are (rho,r,b). Source controls use N=8, b=1.
The command-line controls are CPU-only. No source author code is used.
"""

import numpy as np


def dimension(n):
    if n < 2:
        raise ValueError("The source N-row model requires N>=2")
    return 3*n-1


def unpack(state, n):
    if len(state) != dimension(n):
        raise ValueError("State must have 3*N-1 components")
    state = np.asarray(state)
    c, s = state[:n], state[n:2*n]
    independent = state[2*n:]
    x = np.r_[independent, -np.sum(independent)]
    return c, s, x


def field_jac(state, parameters, n):
    """Faithful constrained vector field and analytic reduced Jacobian."""
    rho, r, b = parameters
    c, s, x = unpack(state, n)
    v = -rho*x+r*(s-np.mean(s))
    field = np.r_[-c-v*s+1, -s+v*c+b, v[:-1]]
    eye = np.eye(n)
    projection = eye-np.ones((n, n))/n
    embedding = np.vstack((np.eye(n-1), -np.ones((1, n-1))))
    jac = np.zeros((3*n-1, 3*n-1), dtype=np.result_type(field))
    jac[:n, :n] = -eye
    jac[:n, n:2*n] = -np.diag(v)-r*s[:, None]*projection
    jac[:n, 2*n:] = rho*s[:, None]*embedding
    jac[n:2*n, :n] = np.diag(v)
    jac[n:2*n, n:2*n] = -eye+r*c[:, None]*projection
    jac[n:2*n, 2*n:] = -rho*c[:, None]*embedding
    jac[2*n:, n:2*n] = r*projection[:-1]
    jac[2*n:, 2*n:] = -rho*np.eye(n-1)
    return field, jac


def rhs(t, state, parameters, n):
    # Avoid constructing a matrix for every explicit solver evaluation.
    rho, r, b = parameters
    c, s, x = unpack(state, n)
    v = -rho*x+r*(s-np.mean(s))
    return np.r_[-c-v*s+1, -s+v*c+b, v[:-1]]


def source_force_constraint_field(state, parameters, n):
    """Independent force-difference solve and unscaled Fourier equations.

    A is chosen to keep transition rates nonnegative. Its magnitude cancels
    from the normalized field; B=b*A and gamma=-r/A.
    """
    rho, r, b = parameters
    c, s, x = unpack(state, n)
    a = -.1/np.sqrt(1+b*b)
    big_b, gamma = b*a, -r/a
    cosine, sine = a*c, a*s
    differences = np.eye(n)[:-1]-np.eye(n)[1:]
    constraints = np.vstack((differences, np.ones(n)))
    forcing = np.r_[-rho*(differences @ x)-gamma*(differences @ sine), 0.]
    velocity = np.linalg.solve(constraints, forcing)
    dc = -cosine-velocity*sine+a
    ds = -sine+velocity*cosine+big_b
    return np.r_[dc/a, ds/a, velocity[:-1]]


def source_parameters(delta):
    """Source Table1 convention, reproducing its printed d=-56.4588.

    Chemical activity is the numerical energy in nN*nm, not the bare number
    15. This convention is recorded explicitly instead of silently fitting.
    """
    ell, eta, stiffness = 1e-8, 1e-7, 9.5e-5
    kbt, potential, prefactor = .0042668, .042668e-18, 1000.
    activity0 = 15*kbt
    lam = 2*np.pi**2*potential/(eta*ell)
    zeta = 2*np.pi*stiffness/(eta*ell)
    d = -prefactor*ell/(2*np.pi*lam)*(2*np.pi*prefactor+zeta*ell/np.sqrt(activity0))
    activity = activity0*(1+delta)
    total_rate = prefactor*np.sqrt(activity)
    first_rate = d*activity
    a = first_rate/total_rate
    gamma = lam/(total_rate*ell)
    rho, r, b = stiffness/(eta*total_rate), -gamma*a, 1.
    return {"delta": delta, "rho": rho, "r": r, "b": b,
            "activity_nN_nm": activity, "a0_per_second": total_rate,
            "a1_per_second": first_rate, "A": a, "B": a,
            "gamma": gamma, "d_source_convention": d,
            "ell_nm": 10., "source": "Table1, Eqs.16/65, Fig9"}


def source_start(n=8, seed=20251030, amplitude_nm=.01):
    """Source random-start prescription with our declared reproducible seed.

    Draw absolute positions x1..x_(N-1), set x0=xN=0, then take differences.
    Drawing independent relative displacements would be a different start.
    """
    rng = np.random.default_rng(seed)
    positions = np.r_[0., rng.uniform(-amplitude_nm, amplitude_nm, n-1), 0.]
    x = 2*np.pi*np.diff(positions)/10.
    return np.r_[np.ones(n), np.ones(n), x[:-1]], positions


def source_algebra_audit():
    rng = np.random.default_rng(72)
    plane_rng = np.random.default_rng(73)
    errors = {"jacobian_complex_step_scaled_error": 0.,
              "source_force_constraint_scaled_error": 0.,
              "n2_core_equivalence_scaled_error": 0.,
              "constraint_velocity_sum_max_abs": 0.,
              "stable_plane_field_identity_scaled_error": 0.,
              "stable_eigenspace_identity_max_abs": 0.}
    for n in (2, 3, 8):
        for _ in range(30):
            p = (10**rng.uniform(-1, 1), 10**rng.uniform(-1, 1), rng.uniform(-2, 2))
            state = rng.normal(size=dimension(n))
            field, jac = field_jac(state, p, n)
            cj = np.column_stack([np.imag(rhs(0., state.astype(complex)
                 +1e-28j*np.eye(dimension(n))[i], p, n))/1e-28 for i in range(dimension(n))])
            errors["jacobian_complex_step_scaled_error"] = max(
                errors["jacobian_complex_step_scaled_error"], float(np.max(np.abs(jac-cj)/(1+np.abs(jac)))))
            independently = source_force_constraint_field(state, p, n)
            errors["source_force_constraint_scaled_error"] = max(
                errors["source_force_constraint_scaled_error"], float(np.max(np.abs(field-independently)/(1+np.abs(field)))))
            c, s, x = unpack(state, n)
            rho, r, b = p
            v = -rho*x+r*(s-np.mean(s))
            errors["constraint_velocity_sum_max_abs"] = max(
                errors["constraint_velocity_sum_max_abs"], float(abs(np.sum(v))))
            plane = np.r_[plane_rng.normal(size=n), np.full(n, plane_rng.normal()), np.zeros(n-1)]
            expected_plane = np.r_[1-plane[:n], b-plane[n:2*n], np.zeros(n-1)]
            errors["stable_plane_field_identity_scaled_error"] = max(
                errors["stable_plane_field_identity_scaled_error"],
                float(np.max(np.abs(rhs(0., plane, p, n)-expected_plane)/(1+np.abs(expected_plane)))))
            equilibrium = np.r_[np.ones(n), np.full(n, b), np.zeros(n-1)]
            stable_basis = np.zeros((dimension(n), n+1))
            stable_basis[:n, :n] = np.eye(n)
            stable_basis[n:2*n, n] = 1.
            errors["stable_eigenspace_identity_max_abs"] = max(
                errors["stable_eigenspace_identity_max_abs"],
                float(np.max(np.abs(field_jac(equilibrium, p, n)[1]@stable_basis+stable_basis))))
            if n == 2:
                y = -(s[0]-s[1])/2
                z = 1-(c[0]+c[1])/2
                h = rho*x[0]+r*y
                dy, dz = -(field[2]-field[3])/2, -(field[0]+field[1])/2
                mapped = np.array([rho*field[4]+r*dy, dy, dz])
                core = np.array([(-rho+r*(1-z))*h-r*y, h*(1-z)-y, h*y-z])
                errors["n2_core_equivalence_scaled_error"] = max(
                    errors["n2_core_equivalence_scaled_error"], float(np.max(np.abs(mapped-core)/(1+np.abs(core)))))
    return {"seed": 72, "stable_plane_seed": 73, "samples_per_N": 30, "N_values": [2, 3, 8],
            "complex_step": 1e-28, "normalization": "abs(error)/(1+abs(reference))",
            "sample_ranges": "rho,r=10**U(-1,1); b=U(-2,2); standard-normal states",
            **errors}


def control_metrics(times, states, n, source):
    c, s = states[:n], states[n:2*n]
    x = np.vstack((states[2*n:], -np.sum(states[2*n:], axis=0)))
    last = times >= times[-1]-min(100., .4*times[-1])
    radius = c*c+s*s
    dimensional_radius = abs(source["A"])*np.sqrt(radius)
    return {"final_state": states[:, -1].tolist(),
            "last_window_start": float(times[last][0]),
            "relative_position_amplitudes_X": np.ptp(x[:, last], axis=1).tolist(),
            "max_abs_relative_position_X": float(np.max(np.abs(x))),
            "max_motor_radius_squared": float(np.max(radius)),
            "invariant_radius_squared": 1+source["b"]**2,
            "max_constraint_sum_abs": float(np.max(np.abs(np.sum(x, axis=0)))),
            "ell_times_density_lower_bound": float(.5-np.max(dimensional_radius)),
            "ell_times_density_upper_bound": float(.5+np.max(dimensional_radius)),
            "final_derivative_max_abs": float(np.max(np.abs(rhs(0., states[:, -1],
                (source['rho'], source['r'], source['b']), n))))}


def run_source_controls(long_time=1000., seed=20251030):
    """Independent source-regime controls; no parameter survey or chaos test."""
    from scipy.integrate import solve_ivp
    import time
    started = time.monotonic()
    n = 8
    initial, positions = source_start(n, seed)
    cases = []
    for delta in (-.1, .1):
        source = source_parameters(delta)
        p = (source["rho"], source["r"], source["b"])
        short_times = np.linspace(0., 20., 2001)
        short = {}
        trajectories = []
        for method in ("DOP853", "Radau"):
            kwargs = {"jac": lambda t, y: field_jac(y, p, n)[1]} if method == "Radau" else {}
            t0 = time.monotonic()
            sol = solve_ivp(lambda t, y: rhs(t, y, p, n), (0., 20.), initial,
                            method=method, t_eval=short_times, rtol=2e-10,
                            atol=2e-12, max_step=.1, **kwargs)
            short[method] = {"success": bool(sol.success), "message": sol.message,
                             "nfev": sol.nfev, "njev": sol.njev,
                             "wall_seconds": time.monotonic()-t0}
            if not sol.success:
                raise RuntimeError(sol.message)
            trajectories.append(sol.y)
        short["max_coordinate_error"] = np.max(np.abs(trajectories[0]-trajectories[1]), axis=1).tolist()
        short["max_abs_error"] = float(np.max(np.abs(trajectories[0]-trajectories[1])))

        # A long explicit run and a separate implicit run to a shared shorter
        # horizon. Retain actual solver diagnostics and no silent exclusions.
        long_results = {}
        long_solutions = {}
        for method, end in (("DOP853", long_time), ("Radau", min(long_time, 300.))):
            eval_times = np.linspace(0., end, int(round(end/.025))+1)
            kwargs = {"jac": lambda t, y: field_jac(y, p, n)[1]} if method == "Radau" else {}
            t0 = time.monotonic()
            sol = solve_ivp(lambda t, y: rhs(t, y, p, n), (0., end), initial,
                            method=method, t_eval=eval_times, rtol=2e-10,
                            atol=2e-12, max_step=.1, **kwargs)
            record = {"success": bool(sol.success), "message": sol.message,
                      "end_time": float(sol.t[-1]), "nfev": sol.nfev,
                      "njev": sol.njev, "wall_seconds": time.monotonic()-t0,
                      **control_metrics(sol.t, sol.y, n, source)}
            long_results[method] = record
            long_solutions[method] = sol
        explicit, implicit = long_solutions["DOP853"], long_solutions["Radau"]
        overlap = implicit.y.shape[1]
        if np.allclose(explicit.t[:overlap], implicit.t, rtol=0, atol=1e-12):
            long_results["shared_horizon_max_abs_difference"] = float(np.max(
                np.abs(explicit.y[:, :overlap]-implicit.y)))
        # Rising zeros of the first relative displacement, linearly interpolated
        # from saved samples. These are diagnostics, not precise event solves.
        x = explicit.y[2*n]
        idx = np.flatnonzero((x[:-1]<0)&(x[1:]>=0)&(explicit.t[:-1]>.5*long_time))
        crossings = explicit.t[idx]-x[idx]*(explicit.t[idx+1]-explicit.t[idx])/(x[idx+1]-x[idx])
        long_results["sampled_section_periods_last_ten"] = np.diff(crossings)[-10:].tolist()
        long_results["sampled_section_method"] = "X1 rising zero, linear interpolation at spacing.025; not shooting"
        cases.append({"source_parameters": source, "short": short, "long": long_results})
    return {"N": n, "dimension": dimension(n), "seed": seed,
            "initial_absolute_positions_nm": positions.tolist(),
            "initial_state": initial.tolist(),
            "start_provenance": "Source uniform-displacement prescription; our seed, not published exact draws",
            "cases": cases, "total_wall_seconds": time.monotonic()-started,
            "scope": "Source-regime qualitative controls only; no Lyapunov survey or new-attractor claim"}


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-output", type=Path)
    parser.add_argument("--control-output", type=Path)
    parser.add_argument("--long-time", type=float, default=1000.)
    parser.add_argument("--seed", type=int, default=20251030)
    args = parser.parse_args()
    if args.audit_output:
        result, path = source_algebra_audit(), args.audit_output
    elif args.control_output:
        result, path = run_source_controls(args.long_time, args.seed), args.control_output
    else:
        parser.print_help()
        raise SystemExit(0)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
