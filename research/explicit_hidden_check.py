"""Check exported enclosures and cross-check nominal periods with two ODE solvers."""

import json
from decimal import Decimal

import numpy as np
from flint import arb, ctx
from scipy.integrate import solve_ivp

from explicit_hidden_cycles import OUT, PARAMETERS


def read(name):
    return json.loads((OUT / name).read_text())


def endpoints(values):
    values = list(values)
    return (
        min(arb(pair[0]).lower() for pair in values),
        max(arb(pair[1]).upper() for pair in values),
    )


def check_cover():
    cover = read("beta-cover.json")
    lo, hi = map(Decimal, cover["beta_interval"])
    r = Decimal(cover["cell_radius"])
    edge = lo
    for row in cover["rows"]:
        center = Decimal(row["center"])
        assert center - r == edge
        edge = center + r
        for name in ("right_lobe", "outer"):
            assert arb(row["contraction_upper"][name]) < 1
            assert all(arb(m) > 0 for m in row["schur_lower"][name])
        assert arb(row["origin_radius_lower"]) > arb("2.32e-7")
        assert arb(row["outer_radius_lower"]) > arb("0.00588")
    assert edge == hi and len(cover["rows"]) == cover["cell_count"]
    intervals = {
        k: endpoints(row["periods"][k] for row in cover["rows"])
        for k in ("right_lobe", "outer")
    }
    assert intervals["outer"][1] < intervals["right_lobe"][0]
    return {k: [str(v) for v in pair] for k, pair in intervals.items()}


def integrate(name, source, method):
    a, b, m0, m1, g = map(float, PARAMETERS)
    central = source["shooting_center"]
    q0 = np.array([1.0, float(central[0]["midpoint"]), float(central[1]["midpoint"])])
    period = sum(float(item["midpoint"]) for item in central[2::3])

    def rhs(t, state):
        x, y, z = state
        h = m1 * x + (m0 - m1) * min(1.0, max(-1.0, x))
        return [a * (y - x - h), x - y + z, -b * y - g * z]

    def section(t, state):
        return state[0] - 1

    section.direction = 1
    sol = solve_ivp(
        rhs,
        (0, period * 1.01),
        q0,
        method=method,
        atol=2e-13,
        rtol=2e-12,
        max_step=0.015,
        events=section,
        dense_output=True,
    )
    assert sol.success
    returns = [t for t in sol.t_events[0] if t > 1e-6]
    assert len(returns) == 1
    dt = abs(returns[0] - period)
    closure = np.linalg.norm(sol.sol(returns[0]) - q0, ord=np.inf)
    assert dt < 2e-8 and closure < 5e-7
    return dict(
        period=period,
        computed_return=returns[0],
        time_error=dt,
        section_closure_inf=closure,
        solver=method,
    )


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled")
    ctx.dps = 80
    periods = check_cover()
    base = read("certificate-70.json")
    other = read("certificate-90.json")
    starts = read("starts-70.json")
    other_starts = read("starts-90.json")
    points = {}
    for k in ("right_lobe", "outer"):
        assert base["displaced_box_rejected"] and other["displaced_box_rejected"]
        for axis in (0, 1):
            for field in ("lower", "upper"):
                x = arb(base["cycles"][k]["section_box"][axis][field])
                y = arb(other["cycles"][k]["section_box"][axis][field])
                assert abs(x - y) < arb("1e-10")
        for record in (starts, other_starts):
            start = record["starts"][k]
            assert start["inside_uniform_shooting_box"]
            assert arb(start["derivative_norm_upper"]["upper"]) < 1
            assert arb(start["final_distance_upper"]["upper"]) < arb(
                start["trapping_radius"]["lower"]
            )
        points[k] = [
            integrate(k, base["cycles"][k], method) for method in ("DOP853", "Radau")
        ]
    summary = dict(
        beta_cover_period_enclosures=periods,
        independent_integrations=points,
        cover_closed=True,
        shifted_roots_rejected=True,
        precision_checks_passed=True,
    )
    (OUT / "independent-check.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(points, indent=2))


if __name__ == "__main__":
    main()
