"""Validated local ingredients for the proposed layer organizer.

Each reported inclusion is checked with Arb balls, never inferred from a
floating-point residual. Failure leaves an explicit incomplete certificate.
"""

import json
from pathlib import Path
import numpy as np
from flint import arb, arb_mat, ctx
from chua_interval_map import IntervalHalfMap, col, eye, midpoint, upper, bounds
from predict_chua_layers import ROOT, parent_at

OUT = ROOT / "results/chua-layer-infinity"


def submatrix(J):
    return arb_mat([[J[i, j] for j in range(2)] for i in range(2)])


def saddle_center(f, beta):
    prior = parent_at(float(beta.mid()))
    p = col([arb(str(x)) for x in prior["q"][:2]])
    for _ in range(7):
        value, J, H, arcs = f.evaluate(list(p.entries()), beta)
        residual = value - p
        if all(abs(x) < arb(10) ** (-(ctx.dps - 25)) for x in residual.entries()):
            break
        p = midpoint(p - (submatrix(J) - eye(2)).inv() * residual)
    return p


def enclose_saddle(f, center, beta, radius, inflation_floor="1e-60"):
    box = col([arb(x, arb(radius)) for x in center.entries()])
    value0, J0, _, _ = f.evaluate(list(center.entries()), beta.mid())
    C = midpoint((submatrix(J0) - eye(2)).inv())
    for iteration in range(6):
        value, J, H, arcs = f.evaluate(list(box.entries()), beta)
        at_center = f.evaluate(list(center.entries()), beta, order=2)[0]
        R = eye(2) - C * (submatrix(J) - eye(2))
        K = center - C * (at_center - center) + R * (box - center)
        inside = all(box[i, 0].contains_interior(K[i, 0]) for i in range(2))
        norm = arb(0)
        for i in range(2):
            norm = norm.max(sum(abs(R[i, j]) for j in range(2)))
        if not inside or not norm < 1:
            raise ValueError("Saddle Krawczyk inclusion or contraction failed")
        if iteration == 0:
            initial = dict(
                box=[bounds(x) for x in box.entries()],
                image=[bounds(x) for x in K.entries()],
                contraction=bounds(norm),
            )
        # Intersect with a slightly inflated inclusion, retaining strictness.
        box = col(
            [
                box[i, 0].intersection(
                    arb(K[i, 0].mid(), K[i, 0].rad() * arb("1.01") + arb(inflation_floor))
                )
                for i in range(2)
            ]
        )
    value, J, H, arcs = f.evaluate(list(box.entries()), beta, verify_path=True)
    M = eye(2) - submatrix(J)
    pb = M.inv() * col([J[i, 2] for i in range(2)])
    v = col([pb[0, 0], pb[1, 0], 1])
    pbb = M.inv() * col([(v.transpose() * h * v)[0, 0] for h in H])
    tr = J[0, 0] + J[1, 1]
    det = submatrix(J).det()
    disc = (tr * tr - 4 * det).sqrt()
    eig = [(tr - disc) / 2, (tr + disc) / 2]
    return (
        box,
        pb,
        pbb,
        dict(
            krawczyk=initial,
            root=[bounds(x) for x in box.entries()],
            beta_derivative=[bounds(x) for x in pb.entries()],
            beta_second_derivative=[bounds(x) for x in pbb.entries()],
            half_eigenvalues=[bounds(x) for x in eig],
            half_determinant=bounds(det),
            flight_times=[bounds(a["duration"]) for a in arcs],
            path_cells=[a["path_certificate"]["accepted_cells"] for a in arcs],
        ),
    )


def local_coordinates(f, p, beta, V, W, pb, delta):
    z = col([arb(0, delta), arb(0, delta)])
    state = p + V * z
    value, J, H, arcs = f.evaluate(list(state.entries()), beta, verify_path=True)
    A = W * submatrix(J) * V
    # Coordinates are z = V^-1 (state - p(beta)), with V held fixed.
    S = arb_mat([[V[0, 0], V[0, 1], pb[0, 0]], [V[1, 0], V[1, 1], pb[1, 0]], [0, 0, 1]])
    local_hess = [S.transpose() * h * S for h in H]
    transformed = [
        sum_matrices([local_hess[i] * W[k, i] for i in range(2)]) for k in range(2)
    ]
    second = arb(0)
    mixed = arb(0)
    for k in range(2):
        for i in range(2):
            mixed = mixed.max(transformed[k][i, 2].abs_upper())
            for j in range(2):
                second = second.max(transformed[k][i, j].abs_upper())
    # This bounds derivatives of the centered map with respect to beta since
    # the centered map vanishes at z=0 for every beta.
    parameter_bound = 2 * mixed * delta
    return (
        A,
        transformed,
        dict(
            domain_half_width=str(delta),
            derivative=[[bounds(A[i, j]) for j in range(2)] for i in range(2)],
            second_derivative_upper=bounds(second),
            mixed_parameter_derivative_upper=bounds(mixed),
            centered_parameter_derivative_bound=bounds(parameter_bound),
            path_cells=[a["path_certificate"]["accepted_cells"] for a in arcs],
        ),
    )


def sum_matrices(items):
    result = items[0]
    for x in items[1:]:
        result += x
    return result


def main():
    OUT.mkdir(exist_ok=True)
    f = IntervalHalfMap(80)
    beta0 = arb("14.3063590999613")
    beta = arb(beta0, arb("1e-10"))
    center = saddle_center(f, beta0)
    p0, pb0, pbb0, cert0 = enclose_saddle(f, center, beta0, "1e-50")
    p, pb, pbb, cert = enclose_saddle(f, center, beta, "1e-8")
    _, J, _, _ = f.evaluate(list(p0.entries()), beta0)
    values, Vnp = np.linalg.eig(
        np.array([[float(J[i, j].mid()) for j in range(2)] for i in range(2)])
    )
    order = np.argsort(-abs(values))
    Vnp = np.real_if_close(Vnp[:, order])
    if Vnp[0, 0] < 0:
        Vnp[:, 0] *= -1
    V = arb_mat([[arb(str(v)) for v in row] for row in Vnp])
    W = V.inv()
    A, H, local = local_coordinates(f, p, beta, V, W, pb, arb("1e-7"))
    result = dict(
        status="local_ingredients_only",
        precision_bits=ctx.prec,
        beta=bounds(beta),
        fixed_parameter_saddle=cert0,
        parameter_family_saddle=cert,
        coordinates=[[str(V[i, j]) for j in range(2)] for i in range(2)],
        local=local,
        scope="Local saddle and itinerary enclosures. These alone do not certify a homoclinic tangency or infinitely many windows.",
    )
    (OUT / "local-certificate.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(dict(saddle=cert["half_eigenvalues"], local=local), indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
