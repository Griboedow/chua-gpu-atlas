"""Algebraic checks for the exact motor-row memory and energy identities.

No new trajectories are integrated. The proof and finite-past qualifications
are in docs/motor-rows-memory-structure.md. This is an analysis of the existing
published model, not a new physical extension or a priority claim.
"""
import argparse
import datetime
import json
from pathlib import Path

import numpy as np

from motor_rows import dimension, source_force_constraint_field, unpack


def audit():
    rng = np.random.default_rng(2026103007)
    errors = dict(complex_memory_field=0., projected_velocity_equation=0.,
                  energy_derivative=0., integrating_factor_derivative=0.)
    count = 0
    for n in (2, 3, 4, 8):
        for _ in range(64):
            rho, gain = 10**rng.uniform(-1, 1, 2)
            b = rng.uniform(-2, 2)
            radius = np.sqrt(1+b*b)
            phase = rng.uniform(-np.pi, np.pi, n)
            radial = radius*np.sqrt(rng.random(n))
            state = np.r_[radial*np.cos(phase), radial*np.sin(phase),
                          rng.normal(size=n-1)]
            assert len(state) == dimension(n)
            c, s, x = unpack(state, n)
            force_field = source_force_constraint_field(state, (rho, gain, b), n)
            dc, ds = force_field[:n], force_field[n:2*n]
            velocity = np.r_[force_field[2*n:], -force_field[2*n:].sum()]
            vdot = -rho*velocity + gain*(ds-ds.mean())
            predicted = -(rho+1)*velocity-rho*x+gain*(velocity*c-(velocity*c).mean())
            z = c+1j*s
            zdot = dc+1j*ds
            complex_rhs = 1+1j*b+(-1+1j*velocity)*z
            energy_dot = velocity@vdot+rho*x@velocity
            balance = np.sum((gain*c-rho-1)*velocity**2)
            # d/dt(exp(t-iX)z) = exp(t-iX)(1+ib), evaluated at t=0.
            factor_derivative = np.exp(-1j*x)*((1-1j*velocity)*z+zdot)
            factor_expected = np.exp(-1j*x)*(1+1j*b)
            pairs = dict(complex_memory_field=(zdot, complex_rhs),
                         projected_velocity_equation=(vdot, predicted),
                         energy_derivative=(energy_dot, balance),
                         integrating_factor_derivative=(factor_derivative, factor_expected))
            for key, (left, right) in pairs.items():
                err = float(np.max(np.abs(left-right))/(1+np.max(np.abs(left))+np.max(np.abs(right))))
                errors[key] = max(errors[key], err)
            count += 1
    assert max(errors.values()) < 1e-11, errors
    return dict(checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                seed=2026103007, random_states=count, row_counts=[2, 3, 4, 8],
                source_field='Independent source force-difference reconstruction',
                maximum_scaled_identity_errors=errors,
                new_trajectories_integrated=0,
                proof_scope='Existing normalized model; rho,r>0 and invariant row disks c_j^2+s_j^2<=1+b^2. Strict sufficient exclusion r*sqrt(1+b^2)<rho+1. No exclusion is asserted outside it.',
                priority='These are reproducible derivations, not a claim that the identities or memory formulation were previously unknown.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('results/motor-rows/memory-structure.json'))
    args = parser.parse_args()
    result = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
