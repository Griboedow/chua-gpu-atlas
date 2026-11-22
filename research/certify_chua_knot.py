"""Certify that the coexisting attracting S02 and S08 Chua cycles form a non-split link."""

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from certify_chua_link_witnesses import FullReturnMap, full_record, rational
from chua_interval_map import bounds, col, eye, midpoint
from chua_link_bracket import frontier_jones
from flint import arb, ctx, fmpz_mat
from spherogram import Link

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/chua-link-search"
GRID = 2**32


def upper_integer(value):
    return int(value.upper().ceil().unique_fmpz())


def certified_arcs(name, digits):
    reference = json.loads((SOURCE / f"overlaps/S02-S08/{name}-d90-n384.json").read_text())
    certificate = json.loads((SOURCE / "certificate-S02-S08-d140.json").read_text())
    record = next(r for r in certificate["records"] if r["id"] == name)
    ctx.dps = certificate["digits"]
    assert ctx.dps == 140 and certificate["name"] == "S02-S08"
    pars = [rational(certificate["parameters"][k]) for k in ("alpha", "beta", "m0", "m1")]
    seed = full_record(reference)
    curve = FullReturnMap(seed, ctx.dps, pars)
    # Reproduce the actual Arb box before serialization. The archived short
    # ball strings alone lose too much precision for the S02 return.
    center = col([arb(str(v)) for v in seed["cycle"]["q"][:2]])
    for _ in range(8):
        value, jac, _ = curve.evaluate(center.entries())
        center = midpoint(center - (jac - eye(2)).inv() * (value - center))
    root = col([arb(v, arb(certificate["root_radius"])) for v in center.entries()])
    assert [bounds(v) for v in root.entries()] == record["root"]
    value, jac, _ = curve.evaluate(center.entries())
    _, derivative, _ = curve.evaluate(root.entries())
    C = midpoint((jac - eye(2)).inv())
    remainder = eye(2) - C * (derivative - eye(2))
    inclusion = center - C * (value - center) + remainder * (root - center)
    assert all(b.contains_interior(k) for b, k in zip(root.entries(), inclusion.entries()))
    assert all(sum(abs(remainder[i, j]) for j in range(2)) < 1 for i in range(2))
    assert record["primitive_period"] and record["passed"]
    ctx.dps = digits
    _, _, arcs = curve.evaluate(root.entries(), verify_path=True)
    assert [{k: a[k] for k in ("region", "target", "sign")} for a in arcs] == record[
        "section_itinerary"
    ]
    return curve, arcs


def enclose_polygon(curve, arcs, steps_per_time):
    vertices, endpoint_errors, second_order, directions = [], [], [], []
    for arc in arcs:
        A, center, _ = curve.matrices(None, arc["region"])
        duration = arc["duration"]
        count = max(1, upper_integer(duration * steps_per_time))
        dt = duration / count
        # Re-anchor each narrow time interval at its individual midpoint.
        # Exponentiating one wide interval over a long flight would wrap.
        window = (A * arb(0, dt.upper() / 2)).exp()
        for k in range(count):
            start = curve.flow(A, center, arc["start"], k * dt)[0]
            mid = curve.flow(A, center, arc["start"], (2 * k + 1) * dt / 2)[0]
            segment = center + window * (mid - center)
            velocity = A * (segment - center)
            acceleration = A * velocity
            # Python's rounding only chooses a convenient integer vertex.
            # Subsequent Arb comparisons bound its actual distance to the orbit.
            point = tuple(round(float(v.mid()) * GRID) for v in start.entries())
            vertices.append(point)
            endpoint_errors.append(
                [upper_integer(abs(v * GRID - q)) for v, q in zip(start.entries(), point)]
            )
            second_order.append(
                [upper_integer(dt**2 * abs(v) * GRID / 8) for v in acceleration.entries()]
            )
            directions.append([1 if v > 0 else -1 if v < 0 else 0 for v in velocity.entries()])
    assert all(abs(v) < 2**50 for point in vertices for v in point)
    assert all(
        0 <= r < 2**48 for rows in (endpoint_errors, second_order) for row in rows for r in row
    )
    coords = np.asarray(vertices, dtype=np.int64)
    differences = np.roll(coords, -1, axis=0) - coords
    signs = np.asarray(directions) * (np.sign(differences) == np.asarray(directions))
    e = np.asarray(endpoint_errors)
    tube = np.maximum(e, np.roll(e, -1, axis=0)) + np.asarray(second_order)
    return coords, tube, signs


def tube_boxes(points, errors):
    # All operations here stay within 2**52 on the accepted inputs.
    assert int(np.max(np.abs(points))) < 2**50
    assert int(np.max(errors)) < 2**50 and np.all(errors >= 0)
    ends = np.roll(points, -1, axis=0)
    return np.minimum(points, ends) - errors, np.maximum(points, ends) + errors


def isotopy(points, errors, signs):
    lo, hi = tube_boxes(points, errors)
    assert np.all(np.any(signs != 0, axis=1)), "A homotopy segment may fold"
    adjacent = 0
    for i in range(len(points) - 1):
        candidates = np.flatnonzero(np.all((hi[i + 1 :] >= lo[i]) & (hi[i] >= lo[i + 1 :]), axis=1))
        for j in candidates + i + 1:
            if j - i not in (1, len(points) - 1) or not np.any(
                (signs[i] != 0) & (signs[i] == signs[j])
            ):
                raise ValueError(f"Unresolved homotopy segments {i}, {j} of {len(points)}")
            adjacent += 1
    assert adjacent == len(points)
    return {
        "coordinate_separation": math.comb(len(points), 2) - adjacent,
        "adjacent_monotonicity": adjacent,
    }


def joint_separation(first, second):
    lo_a, hi_a = tube_boxes(*first[:2])
    lo_b, hi_b = tube_boxes(*second[:2])
    for lo, hi in zip(lo_a, hi_a):
        if np.any(np.all((hi >= lo_b) & (hi_b >= lo), axis=1)):
            raise ValueError("Intercomponent homotopy boxes overlap")
    return len(lo_a) * len(lo_b)


def exact_crossings(curves, view):
    # Transform integer grid coordinates to two screen coordinates and depth.
    # All products after the bounding-box filter use unbounded Python integers.
    M = np.array(view, dtype=np.int64)
    points = np.vstack(curves)
    assert fmpz_mat(M.tolist()).det() != 0 and int(np.max(np.abs(points))) < 2**50
    assert np.max(np.sum(np.abs(M), axis=1)) <= 16
    P = points @ M.T
    start = 0
    ends = []
    components = []
    for c, curve in enumerate(curves):
        ends.extend(start + (np.arange(len(curve)) + 1) % len(curve))
        components.extend([c] * len(curve))
        start += len(curve)
    next_p = P[ends]
    low = np.minimum(P[:, :2], next_p[:, :2])
    high = np.maximum(P[:, :2], next_p[:, :2])
    vec = next_p - P

    def det(x, y):
        return int(x[0]) * int(y[1]) - int(x[1]) * int(y[0])

    assert all(det(vec[i], vec[ends[i]]) != 0 for i in range(len(P))), "Nongeneric projected vertex"
    events = []
    for i in range(len(P) - 1):
        js = (
            np.flatnonzero(np.all((high[i + 1 :] >= low[i]) & (high[i] >= low[i + 1 :]), axis=1))
            + i
            + 1
        )
        for j in js:
            if ends[i] == j or ends[j] == i:
                continue
            u, v, w = vec[i], vec[j], P[j] - P[i]
            D, n1, n2 = det(u, v), det(w, v), det(w, u)
            if D == 0:
                assert n1 != 0, "Collinear projected edges"
                continue
            orientation = 1 if D > 0 else -1
            if D < 0:
                D, n1, n2 = -D, -n1, -n2
            if not (0 <= n1 <= D and 0 <= n2 <= D):
                continue
            assert n1 not in (0, D) and n2 not in (0, D), "Crossing at a vertex"
            z = -int(w[2]) * D + int(u[2]) * n1 - int(v[2]) * n2
            assert z != 0
            events.append(
                (
                    (components[i], int(i), Fraction(n1, D), z > 0),
                    (components[j], int(j), Fraction(n2, D), z < 0),
                    orientation,
                )
            )
    return events


def coloring(events):
    ordered = sorted((e, i) for i, (a, b, _) in enumerate(events) for e in (a, b))
    assert all(a[0][:3] != b[0][:3] for a, b in itertools.pairwise(ordered)), "Triple crossing"
    n = len(events)
    if n == 0:
        return 1, 0
    matrix = [[0] * n for _ in range(n)]
    arc = 0
    for (_, _, _, over), crossing in ordered:
        if over:
            matrix[crossing][arc] += 2
        else:
            matrix[crossing][arc] -= 1
            arc = (arc + 1) % n
            matrix[crossing][arc] -= 1
    assert arc == 0
    determinant = abs(int(fmpz_mat([r[:-1] for r in matrix[:-1]]).det()))
    return determinant, n


def jones(events, count):
    if not events:
        from chua_link_diagrams import unlink_jones

        return unlink_jones(count), 0
    labels = {}
    before = {}
    offset = 0
    free = 0
    for c in range(count):
        ordered = sorted(
            (p[1], p[2], i, side)
            for i, (a, b, _) in enumerate(events)
            for side, p in enumerate((a, b))
            if p[0] == c
        )
        if not ordered:
            free += 1
            continue
        assert all(a[:2] != b[:2] for a, b in itertools.pairwise(ordered)), "Triple crossing"
        for k, (_, _, i, side) in enumerate(ordered):
            labels[i, side] = offset + k
            before[offset + k] = offset + (k - 1) % len(ordered)
        offset += len(ordered)
    pd = []
    for i, (a, b, orientation) in enumerate(events):
        under, over = (1, 0) if a[3] else (0, 1)
        u, o = labels[i, under], labels[i, over]
        if (orientation if under == 0 else -orientation) < 0:
            pd.append([before[u], o, u, before[o]])
        else:
            pd.append([before[u], before[o], u, o])
    link = Link(pd)
    link.unlinked_unknot_components = free
    assert len(link.link_components) + free == count and link.is_planar()
    value, stats = frontier_jones(link)
    return value, stats["maximum_states"]


def span(polynomial):
    return max(map(int, polynomial)) - min(map(int, polynomial))


def plot_comparison(diagram, target):
    from collections import Counter

    import matplotlib.pyplot as plt

    first, second = (c["jones"] for c in diagram["components"])
    split = Counter()
    for e, a in first.items():
        for f, b in second.items():
            for shift in (-2, 2):
                split[int(e) + int(f) + shift] -= a * b
    joint = {int(e): c for e, c in diagram["joint_jones"].items()}
    exponents = np.array(sorted(set(joint) | set(split)))
    fig, ax = plt.subplots(figsize=(7.6, 2.9), layout="constrained")
    ax.axhline(0, color="0.65", linewidth=0.8)
    ax.bar(
        exponents - 0.65,
        [joint.get(e, 0) for e in exponents],
        width=1.2,
        color="#23556f",
        label="Coexisting S02 and S08",
    )
    ax.bar(
        exponents + 0.65,
        [split[e] for e in exponents],
        width=1.2,
        color="#c77847",
        label="Same knots separated",
    )
    ax.set(xlabel="Exponent of A", ylabel="Jones coefficient", xticks=np.arange(-30, 31, 4))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(ncols=2, frameon=False, loc="upper center")
    ax.margins(y=0.27)
    fig.savefig(target, dpi=180)
    plt.close(fig)


VIEWS = [[[3, 5, 1], [1, -2, 2], [0, 1, 0]], [[2, 1, -3], [1, 4, 2], [0, 0, 1]]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps-per-time", type=int, default=128)
    parser.add_argument("--digits", type=int, default=140)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/chua-link-certificate/certificate.json"
    )
    parser.add_argument(
        "--figure", action="store_true", help="Write the Jones comparison beside the certificate"
    )
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError("The certificate requires assertions")
    polygons, component_info = [], []
    for name in ("S02", "S08"):
        curve, arcs = certified_arcs(name, args.digits)
        points, tube, signs = enclose_polygon(curve, arcs, args.steps_per_time)
        checks = isotopy(points, tube, signs)
        polygons.append((points, tube, signs))
        component_info.append(
            {
                "id": name,
                "vertices": len(points),
                "max_coordinate_tube_units": int(tube.max()),
                "polygon_sha256": hashlib.sha256(points.astype("<i8").tobytes()).hexdigest(),
                "segment_pairs": checks,
            }
        )
    joint_checks = joint_separation(*polygons)
    coords = [p[0] for p in polygons]
    diagrams = []
    for view in VIEWS:
        components = []
        for poly in coords:
            events = exact_crossings([poly], view)
            determinant, crossings = coloring(events)
            polynomial, states = jones(events, 1)
            components.append(
                {
                    "crossings": crossings,
                    "determinant": determinant,
                    "jones": polynomial,
                    "max_bracket_states": states,
                }
            )
        events = exact_crossings(coords, view)
        polynomial, states = jones(events, 2)
        signed_crossings = sum(
            orientation * (1 if a[3] else -1) for a, b, orientation in events if a[0] != b[0]
        )
        assert signed_crossings == 0
        width = span(polynomial)
        split_width = 4 + sum(span(c["jones"]) for c in components)
        assert width != split_width and components[0]["determinant"] == 25
        diagrams.append(
            {
                "projection": view,
                "components": components,
                "joint_crossings": len(events),
                "joint_jones": polynomial,
                "span": width,
                "split_span": split_width,
                "linking_number": signed_crossings // 2,
                "max_bracket_states": states,
            }
        )
    assert diagrams[0]["span"] == diagrams[1]["span"]
    reference_polynomial = json.loads(
        (SOURCE / "overlaps/S02-S08/exhaustive-jones.json").read_text()
    )["polynomial"]
    assert all(d["joint_jones"] == reference_polynomial for d in diagrams)
    for i in range(2):
        assert diagrams[0]["components"][i]["jones"] == diagrams[1]["components"][i]["jones"]
    record = {
        "parameter_source": "results/chua-link-search/certificate-S02-S08-d140.json",
        "digits": args.digits,
        "steps_per_unit_time": args.steps_per_time,
        "dyadic_grid": GRID,
        "components": component_info,
        "separated_cross_component_pairs": joint_checks,
        "diagrams": diagrams,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8", newline="\n")
    if args.figure:
        plot_comparison(diagrams[0], args.output.with_name("jones-comparison.png"))
    print(
        f"{len(coords[0])} + {len(coords[1])} edges. Joint span {width}, split span {split_width}."
    )


if __name__ == "__main__":
    main()
