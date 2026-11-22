"""Verified bounds for the period of the Fig. 10 saddle on the Table 2 parameter box."""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
from itertools import product
from pathlib import Path

from certify_chua_coexistence_box import MultipleShooting
from chua_interval_map import IntervalHalfMap, col, eye, midpoint
from flint import arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parents[1]
NAMES = ("alpha", "beta", "m0", "m1", "gamma")
RANGES = (
    ("8.3886", "8.4222"),
    ("12.1812", "12.2790"),
    ("-1.3689", "-1.3634"),
    ("-0.1715", "-0.1700"),
    ("0.04328", "0.04363"),
)
DIVISIONS = (8, 8, 8, 4, 1)
GUESS = (".69457", "3.40794", ".30457", ".94129", ".19495", ".899315")
# The mathematical column's m1=-0.17 is the upper boundary of Table 2.
NOMINAL = ("8.41", "12.23", "-1.366", "-.17", ".0435")


class SymmetricShooting(MultipleShooting):
    def __init__(self):
        self.original = IntervalHalfMap.__new__(IntervalHalfMap)
        self.n = 2
        self.schedule = [(0, 1, 0.30457, 1), (1, 1, 0.899315, 1)]
        self.start_x = [-1, 1]

    def evaluate(self, q, parameters):
        # Match the final +1 crossing to the central inversion of the initial
        # -1 crossing. Differentiation changes the two endpoint coefficients.
        f, j, arcs = super().evaluate(q, parameters)
        for k in (0, 1):
            f[4 + k, 0] += 2 * q[k, 0]
            j[4 + k, k] += 2
        return f, j, arcs


def grid(index):
    endpoints = []
    for i, n, (lo, hi) in zip(index, DIVISIONS, RANGES):
        a, b = Decimal(lo), Decimal(hi)
        h = (b - a) / n
        endpoints.append((a + i * h, a + (i + 1) * h))
    return endpoints


def newton(system, parameters):
    q = col([arb(s) for s in GUESS])
    for _ in range(6):
        f, j, _ = system.evaluate(q, parameters)
        q = midpoint(q - j.inv() * f)
    return q


def itinerary(system, arc):
    A, c, start, duration = (arc[k] for k in ("A", "c", "start", "duration"))
    region = arc["region"]
    flow = system.original.flow
    assert duration > 0
    pending = [(arb(0), arb(1), 0)]
    accepted = 0
    xmax = arb(0)
    while pending:
        u, v, depth = pending.pop()
        state, _ = flow(A, c, start, u.union(v) * duration)
        x = state[0, 0]
        if region == 0:
            # At both ends the shooting equations prescribe x=-1,+1.
            valid = (A * (state - c))[0, 0] > 0
        else:
            valid = x > 1
            if not valid and v <= arb("0.5"):
                collar, _ = flow(A, c, start, arb(0).union(v) * duration)
                valid = (A * (collar - c))[0, 0] > 0
            if not valid and u >= arb("0.5"):
                collar, _ = flow(A, c, start, u.union(arb(1)) * duration)
                valid = (A * (collar - c))[0, 0] < 0
        if valid and abs(x) < 4:
            accepted += 1
            xmax = max(xmax, x.abs_upper())
        else:
            if depth >= 16:
                raise ArithmeticError("Flight guard failed")
            mid = (u + v) / 2
            pending.extend(((u, mid, depth + 1), (mid, v, depth + 1)))
    return accepted, xmax


def point_witness(system, par):
    q = newton(system, par)
    f, j, _ = system.evaluate(q, par)
    C = midpoint(j.inv())
    box = col([arb(v, arb("1e-13")) for v in q.entries()])
    _, jb, _ = system.evaluate(box, par)
    remainder = eye(6) - C * jb
    assert all(sum(abs(remainder[i, j]) for j in range(6)) < 1 for i in range(6))
    image = q - C * f + remainder * (box - q)
    assert all(b.contains_interior(y) for b, y in zip(box.entries(), image.entries()))
    return box


def certify(task):
    index, digits = task
    ctx.dps = digits
    system = SymmetricShooting()
    endpoints = grid(index)
    param = [arb(str(a)).union(arb(str(b))) for a, b in endpoints]
    nominal = [p.mid() for p in param]
    q = newton(system, nominal)
    f, j, _ = system.evaluate(q, nominal)
    C = midpoint(j.inv())
    fp = system.parameter_jacobian(q, param)
    forcing = -C * f - (C * fp) * col([p - p.mid() for p in param])
    widths = [2 * max(v.abs_upper(), arb("1e-18")) for v in forcing.entries()]
    box = col([arb(v, r) for v, r in zip(q.entries(), widths)])
    _, jb, _ = system.evaluate(box, param)
    rem = eye(6) - C * jb
    image = q + forcing + rem * (box - q)
    norm = max(
        (sum(abs(rem[i, j]) * widths[j] for j in range(6)) / widths[i]).abs_upper()
        for i in range(6)
    )
    if not (
        norm < 1
        and all(b.contains_interior(y) for b, y in zip(box.entries(), image.entries()))
    ):
        raise ArithmeticError(f"No shooting inclusion at {index}: norm {norm}")

    # Work with the sharpened root enclosure for flight and saddle tests.
    _, _, arcs = system.evaluate(image, param)
    d = eye(2)
    integral = arb(0)
    count, xmax = 0, arb(0)
    for arc in arcs:
        cells, ampl = itinerary(system, arc)
        count += cells
        xmax = max(xmax, ampl)
        velocity, E = arc["v"], arc["E"]
        if arc["region"] == 0:
            assert velocity[0, 0] > 0
        else:
            assert velocity[0, 0] < 0
        event = arb_mat(
            [
                [E[i, j] - velocity[i, 0] * E[0, j] / velocity[0, 0] for j in (1, 2)]
                for i in (1, 2)
            ]
        )
        d = event * d
        integral += arc["duration"] * arc["A"].trace()
    # The two section velocities at the antisymmetric endpoints have equal
    # magnitude and opposite sign. Liouville gives the exact negative sign.
    trace = -d.trace()
    det = -integral.exp()
    if not (det < 0 and 1 - trace + det < 0 and 1 + trace + det > 0):
        raise ArithmeticError(
            f"Saddle check failed at {index}: trace {trace}, det {det}"
        )
    period = 2 * (image[2, 0] + image[5, 0])
    # Bind neighboring roots. On a common face its midpoint has an enclosed
    # zero in each adjoining uniqueness box. Uniqueness and invertibility then
    # extend equality from this midpoint to the entire connected closed face.
    faces = 0
    for axis, (i, n) in enumerate(zip(index, DIVISIONS)):
        for side in (0, 1):
            if (side == 0 and i == 0) or (side == 1 and i == n - 1):
                continue
            face = [arb(str(((a + b) / 2).normalize())) for a, b in endpoints]
            face[axis] = arb(str(endpoints[axis][side].normalize()))
            common = point_witness(system, face)
            assert all(
                b.contains_interior(v) for b, v in zip(box.entries(), common.entries())
            )
            faces += 1
    nominal_link = all(a <= Decimal(v) <= b for v, (a, b) in zip(NOMINAL, endpoints))
    if nominal_link:
        common = point_witness(system, [arb(v) for v in NOMINAL])
        assert all(
            b.contains_interior(v) for b, v in zip(box.entries(), common.entries())
        )
        radius = arb("1e-6")
        assert arb(".69457429", radius).contains_interior(common[0, 0])
        assert arb("3.40794063", radius).contains_interior(common[1, 0])
    return (
        index,
        period.lower().str(22),
        period.upper().str(22),
        norm.upper().str(12),
        count,
        xmax.upper().str(12),
        faces,
        nominal_link,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--digits", type=int, default=50)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/chua-saddle-period/interval.json"
    )
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError("The certificate requires assertions.")
    ctx.dps = args.digits
    keys = list(product(*(range(n) for n in DIVISIONS)))
    tasks = [(idx, args.digits) for idx in keys]
    if args.workers == 1:
        rec = list(map(certify, tasks))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            rec = list(pool.map(certify, tasks, chunksize=4))
    assert [r[0] for r in rec] == keys and any(r[7] for r in rec)
    T = min(arb(r[1]).lower() for r in rec).union(max(arb(r[2]).upper() for r in rec))
    # R in ohms, C2 in nF. Product /1000 gives microseconds.
    rc = (
        (arb(788) * arb("0.99")).union(arb(788) * arb("1.01"))
        * ((arb(311) * arb("0.999")).union(arb(311) * arb("1.001")))
        / 1000
    )
    microseconds = T * rc
    # beta = C2*R**2/L couples the clock to the location of each
    # parameter cell. sqrt(beta*L_mH*C2_nF) is in microseconds.
    lcr_factor = arb("0.999").union(arb("1.001"))
    lc = arb("15.79") * arb(311) * lcr_factor**2
    physical_periods = []
    for row in rec:
        beta_lo, beta_hi = grid(row[0])[1]
        beta = arb(str(beta_lo)).union(arb(str(beta_hi)))
        cell_rc = (beta * lc).sqrt().intersection(rc)
        physical_periods.append(arb(row[1]).union(arb(row[2])) * cell_rc)
    coupled = min(v.lower() for v in physical_periods).union(
        max(v.upper() for v in physical_periods)
    )
    # With all diode resistors independently inside ±1%, the smallest
    # second breakpoint still lies beyond the entire certified orbit tube.
    k1_min = 1 + (12020 * arb(".99")) / (1008 * arb("1.01"))
    k2_max = 1 + (1003 * arb("1.01")) / (24410 * arb(".99"))
    break_lower = (k1_min / k2_max).lower()
    max_x = max(arb(r[5]).upper() for r in rec)
    assert max_x < break_lower and coupled > arb(218) * arb("2.5")
    print(
        "T",
        T.lower().str(14),
        T.upper().str(14),
        "us",
        microseconds.lower().str(14),
        microseconds.upper().str(14),
        "coupled us",
        coupled.lower().str(14),
        coupled.upper().str(14),
    )

    def endpoints(v):
        return [v.lower().str(20), v.upper().str(20)]

    face_count = sum((n - 1) * len(keys) // n for n in DIVISIONS)
    assert sum(r[6] for r in rec) == 2 * face_count
    info = {
        "source": "https://doi.org/10.1007/s11071-022-08078-y",
        "parameter_order": NAMES,
        "parameter_ranges": RANGES,
        "grid_dimensions": DIVISIONS,
        "digits": args.digits,
        "cells": len(keys),
        "face_inclusions": sum(r[6] for r in rec),
        "distinct_internal_faces": face_count,
        "nominal_seed_matched": True,
        "period": endpoints(T),
        "rc_microseconds": endpoints(rc),
        "time_scale": {
            "R_ohms": 788,
            "relative_R_radius": "0.01",
            "C2_nanofarads": 311,
            "relative_C2_radius": "0.001",
            "L_millihenries": "15.79",
            "relative_L_radius": "0.001",
        },
        "physical_period_microseconds": endpoints(microseconds),
        "samples_at_2_5_microseconds": endpoints(microseconds / arb("2.5")),
        "coupled_period_microseconds": endpoints(coupled),
        "coupled_samples_at_2_5_microseconds": endpoints(coupled / arb("2.5")),
        "maximum_abs_x_upper": max_x.str(12),
        "second_breakpoint_lower_at_1_percent_resistors": break_lower.str(15),
        "maximum_contraction_upper": max(arb(r[3]).upper() for r in rec).str(12),
        "verified_flight_subintervals": sum(r[4] for r in rec),
        "records_columns": [
            "grid_index",
            "period_lower",
            "period_upper",
            "contraction_upper",
            "flight_subintervals",
            "abs_x_upper",
            "matched_faces",
            "contains_nominal",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Keep the methodology and summary readable and allow one line per cell.
    body = json.dumps(info, indent=2)[:-2] + ',\n  "records": [\n'
    body += ",\n".join("    " + json.dumps(row, separators=(",", ":")) for row in rec)
    args.output.write_text(body + "\n  ]\n}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
