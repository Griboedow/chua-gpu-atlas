"""Instantaneous algebra/geometry checks for the pinned DCR author-code variant.

No trajectories, parameter survey, event policy, or Lyapunov calculation.
This is an independent implementation of the equations described in the note.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "outputs/candidate_research/warmington2026-github-pinned.jl"
OUT = ROOT / "results/dcr/structure-audit.json"
N, C, D = 100, 0.14, 0.02
LR, LT = 8.6*C, 7.8*C
NT = int(np.ceil(N*LR/LT))
L = NT*LT
TWO_PI = 2*np.pi


def force(omega):
    """Algebraic form of the source switch=0 force, evaluated in FP64."""
    omega = np.asarray(omega)
    return 1-omega/(4*np.pi)*(1+np.tanh(40*omega))


def force_derivative(omega):
    u = 40*np.asarray(omega)
    th = np.tanh(u)
    return -(1+th+u*(1-th*th))/(4*np.pi)


def geometry(x, theta):
    delta = np.mod(x, LT)
    r = np.arange(N)*LR-np.sin(theta)
    index = np.clip(np.rint((r-delta)/LT).astype(int), 0, NT-1)
    tooth = delta+index*LT
    angle = np.mod(theta, TWO_PI)
    active_angle = (angle > np.pi-np.arcsin(C)) & (angle < np.pi+np.arcsin(C))
    contact = active_angle & (np.abs(r-tooth) < D)
    return contact, active_angle, r, tooth


def velocity(x, theta, alpha=1., beta=1., mask=None):
    contact, active, _, _ = geometry(x, theta)
    if mask is not None:
        contact = np.asarray(mask, dtype=bool)
    cosines = np.abs(np.cos(theta[contact]))
    # In the pinned FP64 code F(2*pi)=0. Retain this expression explicitly.
    free_force = float(force(TWO_PI)*np.abs(np.cos(theta[active & ~contact])).sum())
    A = float(cosines.sum()+free_force-beta*x)

    def residual(v):
        return float(alpha*v+beta*x-(cosines*force(v/cosines)).sum()-free_force)

    if A == 0:
        v = 0.
    else:
        lo, hi = sorted([A/alpha, 0.])
        # An endpoint can be the exact sharp/backward limit in FP64.
        if residual(lo) == 0:
            v = lo
        elif residual(hi) == 0:
            v = hi
        else:
            v = brentq(residual, lo, hi, xtol=5e-15)
    derivative = float(alpha-force_derivative(v/cosines).sum())
    return v, residual(v), derivative, contact


def main():
    started = time.perf_counter()
    rng = np.random.default_rng(20261030)
    w = brentq(lambda z:z*np.tanh(z)-1, 1., 2., xtol=5e-15)
    maximum_force_derivative = (w-1)/(4*np.pi)
    worst_omega = -w/40
    derivative_identity_error = abs(float(force_derivative(worst_omega))-maximum_force_derivative)
    step = 1e-6
    derivative_fd_error = abs(float((force(worst_omega+step)-force(worst_omega-step))/(2*step))-maximum_force_derivative)
    multiplicities = np.bincount((4*np.arange(N)) % 39, minlength=39)
    block_counts = np.array([sum(multiplicities[(j+np.arange(12)) % 39]) for j in range(39)])
    eta = float(np.exp(-160*np.pi)/(1+np.exp(-160*np.pi)))
    sorted_error = nearest_error = 0.
    for _ in range(256):
        x = rng.uniform(-1000, 1000)
        theta = rng.uniform(-20*np.pi, 20*np.pi, N)
        literal_teeth = np.sort(np.mod(np.arange(NT)*LT+x, L))
        compact_teeth = np.mod(x, LT)+np.arange(NT)*LT
        sorted_error = max(sorted_error, float(np.max(np.abs(literal_teeth-compact_teeth))))
        _, _, rotor, tooth = geometry(x, theta)
        literal_distance = np.min(np.abs(rotor[:, None]-literal_teeth[None, :]), axis=1)
        nearest_error = max(nearest_error, float(np.max(np.abs(np.abs(rotor-tooth)-literal_distance))))

    counts = []
    force_residual = 0.
    derivative_min = np.inf
    sign_errors = 0
    invariance_checks = []
    # Deliberately recruit every geometrically eligible rotor at each state.
    # These are isolated algebraic states, not an integration or a survey.
    for j in range(128):
        alpha, beta = [(1.,1.), (1.,4.), (4.,1.), (4.,4.)][j % 4]
        x = rng.uniform(0, 36/beta)
        delta = np.mod(x, LT)
        centers = np.arange(N)*LR
        nearest = np.clip(np.rint((centers-delta)/LT), 0, NT-1)*LT+delta
        difference = centers-nearest
        eligible = np.abs(difference) < C+D-1e-10
        theta = np.zeros(N)
        sine = np.clip(difference[eligible], -C+1e-10, C-1e-10)
        theta[eligible] = np.pi-np.arcsin(sine)
        v, residual, derivative, contact = velocity(x, theta, alpha, beta)
        counts.append(int(contact.sum()))
        force_residual = max(force_residual, abs(residual))
        derivative_min = min(derivative_min, derivative)
        A = np.abs(np.cos(theta[contact])).sum()-beta*x
        sign_errors += int(np.sign(v) != np.sign(A))
        if j < 8:
            low = velocity(0., theta, alpha, beta)[0]
            high = velocity(36/beta, theta, alpha, beta)[0]
            invariance_checks.append({"alpha": alpha, "beta": beta,
                                      "v_at_zero": low, "v_at_36_over_beta": high})

    # Algebraic contact kinematics: theta is on the negative-cosine branch.
    theta = np.pi-rng.uniform(-np.arcsin(C), np.arcsin(C), 128)
    v = rng.uniform(-3, 3, 128)
    theta_dot = v/np.abs(np.cos(theta))
    separation_rate_error = float(np.max(np.abs(-np.cos(theta)*theta_dot-v)))
    sine_rate_error = float(np.max(np.abs(np.cos(theta)*theta_dot+v)))
    sine_a = rng.uniform(-C/2, C/2, 128)
    dx = rng.uniform(-C/4, C/4, 128)
    theta_a = np.pi-np.arcsin(sine_a)
    theta_b = np.pi-np.arcsin(sine_a-dx)
    branch_displacement_error = float(np.max(np.abs(-np.sin(theta_b)+np.sin(theta_a)-dx)))
    edge_theta = np.zeros(N); edge_theta[0] = np.pi-np.arcsin(.01)
    edge_mask, _, edge_r, edge_tooth = geometry(LT-.01, edge_theta)
    edge_ordinary = float(abs(edge_r[0]-edge_tooth[0]))
    edge_teeth = np.mod(np.arange(NT)*LT+LT-.01, L)
    edge_periodic = float(np.min(np.abs(np.mod(edge_r[0]-edge_teeth+L/2,L)-L/2)))

    # Distance-entry example: rotor 0, tooth 0, no other active angles.
    theta = np.zeros(N); theta[0] = np.pi
    x = D
    free_v, _, _, free_mask = velocity(x, theta)
    included = free_mask.copy(); included[0] = True
    contact_v, _, _, _ = velocity(x, theta, mask=included)
    # sin(pi) is not exactly zero in FP64, so the analytic geometry is also
    # recorded explicitly rather than using rounding to decide the theorem.
    entry = {"alpha":1., "beta":1., "gamma":0., "x":x,
             "rotor_index_zero_based":0, "theta":"pi", "all_other_angles":0.,
             "analytic_gap":-D, "literal_strict_contact":False,
             "free_velocity":free_v, "free_gap_derivative":TWO_PI-free_v,
             "included_contact_velocity":contact_v, "included_contact_gap_derivative":0.,
             "interpretation":"The boundary is excluded by the literal strict test. Contact preserves gap; exact capture needs a stated boundary convention."}

    # Lower-angle example, away from the tooth-distance and wrap boundaries.
    theta = np.zeros(N); theta[1] = np.pi-np.arcsin(C)
    x = LR-C
    free_v, _, _, free_mask = velocity(x, theta)
    included = free_mask.copy(); included[1] = True
    contact_v, _, _, _ = velocity(x, theta, mask=included)
    angular = {"alpha":1., "beta":1., "gamma":0., "x":x,
               "rotor_index_zero_based":1, "theta":"pi-asin(c)",
               "all_other_angles":0., "analytic_gap":0.,
               "free_velocity":free_v, "free_theta_dot":TWO_PI,
               "included_contact_velocity":contact_v,
               "included_contact_theta_dot":float(contact_v/np.sqrt(1-C*C)),
               "interpretation":"Opposing normal directions at the lower angular boundary require explicit treatment of possible sliding/chattering."}

    assert sorted_error < 1e-10 and nearest_error < 1e-10
    assert max(counts) <= int(block_counts.max()) <= 36
    assert force_residual < 1e-10 and derivative_min > 0 and sign_errors == 0
    assert all(p['v_at_zero'] >= 0 and p['v_at_36_over_beta'] <= 0 for p in invariance_checks)
    assert entry['free_gap_derivative'] > 0 and angular['included_contact_theta_dot'] < 0
    assert separation_rate_error < 2e-14 and sine_rate_error < 2e-14
    assert derivative_identity_error < 1e-14 and derivative_fd_error < 1e-9
    assert branch_displacement_error < 2e-14 and not edge_mask[0] and edge_periodic < 1e-12
    report = {
        "scope":"Source-code geometry and switch=0, gamma=0, N=100; instantaneous algebra only. No trajectory or event continuation implemented.",
        "source_pinned_commit":"796c5e39468b06fd6fc2c378437d618593948e98",
        "source_file":"outputs/candidate_research/warmington2026-github-pinned.jl",
        "source_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "parameters":{"N":N,"c":C,"d":D,"rotor_spacing":LR,"tooth_spacing":LT,"tooth_count":NT,"tooth_domain_length":L},
        "geometry_bound":{"spacing_ratio":"43/39","residue_site_spacing":LT/39,
                          "allowed_residue_arc_length":2*(C+D),"sites_per_arc_upper_bound":12,
                          "max_site_multiplicity":int(multiplicities.max()),"simple_contact_count_bound":36,
                          "exact_integer_residue_count_bound":int(block_counts.max()),
                          "residue_multiplicities":multiplicities.tolist(),"cyclic_12_site_counts":block_counts.tolist(),
                          "qualification":"Finite ordinary-distance tooth boundaries can only reduce the periodic residue eligibility count."},
        "force_root":{"w_solving_w_tanh_w_equals_one":float(w),"maximum_positive_F_prime":float(maximum_force_derivative),
                      "worst_omega":float(-w/40),"G_prime_lower_bound_alpha1_n36":float(1-36*maximum_force_derivative),
                      "G_prime_lower_bound_alpha1_n33":float(1-33*maximum_force_derivative),
                      "proof_scope":"alpha>=1 with the fixed source geometry; not arbitrary N/geometry or arbitrary low drag"},
        "free_force":{"analytic_F_at_2pi":eta,"fp64_F_at_2pi":float(force(TWO_PI)),
                      "analytic_formula":"eta=1/(1+exp(160*pi))", "analytic_safe_x_upper_bound":"N/beta",
                      "analytic_sharper_upper_bound":"(36+64*eta)/beta", "fp64_sharper_upper_bound":"36/beta",
                      "qualification":"Bounds assume beta>0, gamma=0, unique algebraic roots and a well-defined forward event continuation respecting the branch sign."},
        "checks":{"seed":20261030,"geometry_states":256,"recruited_instantaneous_states":128,
                  "max_sorted_tooth_error":sorted_error,"max_nearest_distance_error":nearest_error,
                  "max_instantaneous_contact_count":max(counts),"max_force_root_residual":force_residual,
                  "min_tested_force_residual_derivative":float(derivative_min),"velocity_sign_failures":sign_errors,
                  "contact_kinematic_gap_rate_error":separation_rate_error,"contact_sine_rate_error":sine_rate_error,
                  "fixed_contact_reconstruction_error":branch_displacement_error,
                  "maximum_F_prime_identity_error":derivative_identity_error,
                  "maximum_F_prime_finite_difference_error":derivative_fd_error,
                  "nonperiodic_boundary_example":{"x":LT-.01,"rotor0_theta":"pi-asin(0.01)",
                      "ordinary_nearest_distance":edge_ordinary,"incorrect_periodic_distance":edge_periodic,
                      "literal_contact":bool(edge_mask[0])},
                  "sampled_inward_boundaries":invariance_checks,"passed":True},
        "distance_entry_example":entry,"lower_angular_boundary_example":angular,
        "event_gate":"The exact fixed-branch reduction does not determine capture, grazing, simultaneous-contact or angular-boundary continuation. No global well-posedness or chaos claim is made.",
        "references":["https://doi.org/10.1098/rsif.2025.0438","https://doi.org/10.5281/zenodo.18196853",
                      "https://github.com/polymaths-lab/actomyosin-kuramoto-rotor-model/blob/796c5e39468b06fd6fc2c378437d618593948e98/DCR_Model_available.jl"],
        "wall_seconds":time.perf_counter()-started,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({k:report[k] for k in ["force_root","checks","distance_entry_example","lower_angular_boundary_example","wall_seconds"]},indent=2))


if __name__ == "__main__":
    main()
