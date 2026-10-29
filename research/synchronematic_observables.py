"""Physical observables on the synchronematic phase-angle quotient.

Each particle has the equivalence (theta, chi) ~ (theta + pi, chi + pi).
The functions below also remove the unobservable uniform phase shift. They are
projections, not an embedding of the seven-dimensional four-particle state.
"""
import numpy as np


def collective_order(state):
    """Return synchronematic coherence C and complex nematic order Q.

    C^2 = (|mean exp(i(theta+chi))|^2 +
           |mean exp(i(theta-chi))|^2) / 2.
    Equivalently C^2 is the all-pairs mean of cos(delta theta)*cos(delta chi).
    """
    state = np.asarray(state)
    n = state.shape[-1] // 2
    if n < 1 or state.shape[-1] != 2 * n:
        raise ValueError('Expected an even final dimension: all phases, then orientations')
    theta, chi = state[..., :n], state[..., n:]
    u = np.exp(1j * (theta + chi)).mean(axis=-1)
    v = np.exp(1j * (theta - chi)).mean(axis=-1)
    coherence = np.sqrt((abs(u)**2 + abs(v)**2) / 2)
    nematic = np.exp(2j * chi).mean(axis=-1)
    return coherence, nematic


def pair_signal(state, i, j, *, odd=False):
    """Return cos(delta theta)*cos(delta chi), or its phase-odd counterpart."""
    state = np.asarray(state)
    n = state.shape[-1] // 2
    delta = state[..., j] - state[..., i]
    phase = np.sin(delta) if odd else np.cos(delta)
    return phase * np.cos(state[..., n + j] - state[..., n + i])


def audit_saved_data():
    """Check all plotted observables and the saved doubled orbit; no integration."""
    import itertools
    import json
    from pathlib import Path
    from synchronematic import flow_and_jac, geometry

    root = Path(__file__).resolve().parents[1] / 'results/synchronematic'
    trace = np.load(root / 'n4-cpu-extended/trace.npz')
    state = trace['state']
    pairs = geometry(trace['positions'])
    n = state.shape[-1] // 2
    sample_indices = [0, len(state) // 3, 2 * len(state) // 3, len(state) - 1]

    def values(y):
        c, q = collective_order(y)
        return np.column_stack([c, q.real, q.imag, pair_signal(y, 0, 1),
                                pair_signal(y, 0, 1, odd=True),
                                pair_signal(y, 0, 2, odd=True), np.sin(2*y[:, n])])

    reference = values(state)
    ranges = []
    error = field_error = jacobian_error = 0.
    for mask in itertools.product([0, 1], repeat=n):
        shift = np.pi * np.r_[mask, mask]
        shifted = state + shift
        error = max(error, float(np.max(abs(values(shifted) - reference))))
        ranges.append(abs(np.exp(1j * shifted[:, :n]).mean(axis=1)))
        for index in sample_indices:
            f, j = flow_and_jac(state[index], pairs)
            ff, jj = flow_and_jac(shifted[index], pairs)
            field_error = max(field_error, float(np.max(abs(f - ff))))
            jacobian_error = max(jacobian_error, float(np.max(abs(j - jj))))
    spread = np.ptp(ranges, axis=0)
    shifted = state.copy()
    shifted[:, :n] += 1.234
    uniform_error = float(np.max(abs(values(shifted) - reference)))
    phase = state[:, :n, None] - state[:, None, :n]
    angle = state[:, n:, None] - state[:, None, n:]
    identity_error = float(np.max(abs(reference[:, 0]**2 -
                                     (np.cos(phase)*np.cos(angle)).mean(axis=(1, 2)))))

    orbit = np.load(root / 'n4-periodic-branch/orbit-1.078.npz')['state']
    mid = (len(orbit) - 1) // 2
    difference = orbit[mid:] - orbit[:mid+1]
    residuals = []
    for mask in itertools.product([0, 1], repeat=n):
        mask = np.asarray(mask)
        shift = np.pi * np.r_[mask[1:] - mask[0], mask]
        residuals.append(float(np.max(abs(np.angle(np.exp(1j*(difference-shift)))))))
    full = np.column_stack([np.zeros(len(orbit)), orbit[:, :3], orbit[:, 3:]])
    pair = pair_signal(full, 0, 1)
    axis = np.sin(2*full[:, n])
    result = dict(trace_samples=len(state), gauge_masks=2**n,
                  trace_indices_for_field_and_jacobian=sample_indices,
                  max_field_difference=field_error, max_jacobian_difference=jacobian_error,
                  max_all_render_observables_gauge_error=error,
                  max_uniform_phase_error=uniform_error,
                  max_coherence_identity_error=identity_error,
                  phase_order_R_gauge_spread_max=float(spread.max()),
                  phase_order_R_gauge_spread_median=float(np.median(spread)),
                  phase_order_R_gauge_spread_min=float(spread.min()),
                  doubled_half_best_uniform_gauge_residual=min(residuals),
                  doubled_half_pair_correlation_difference=float(np.max(abs(pair[mid:]-pair[:mid+1]))),
                  doubled_half_nematic_axis_difference=float(np.max(abs(axis[mid:]-axis[:mid+1]))),
                  new_trajectories_integrated=0)
    assert max(error, uniform_error, identity_error, field_error, jacobian_error) < 1e-12
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    audit_saved_data()
