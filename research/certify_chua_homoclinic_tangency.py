"""Check sufficient interval inequalities for a generically unfolding tangency.

The certificate uses graph-transform bounds on the two exact local invariant
manifolds and Taylor bounds on an actual return iterate. A tangent-line root
supplies only the center of the tested box, never a substitute for a manifold.
See docs/chua-infinity-certificate.md for the mathematical implication.
"""

import argparse, json, hashlib, time, platform
from importlib.metadata import version
from pathlib import Path
from flint import arb, arb_mat, ctx
from chua_interval_map import IntervalHalfMap, col, bounds
from certify_chua_layer_organizer import (
    OUT,
    saddle_center,
    enclose_saddle,
    local_coordinates,
)
from chua_interval_connection import connect


def all_bound(matrix, bound):
    return all(abs(x) < bound for x in matrix.entries())


def certify(digits=100, beta_shift="0", configuration=None):
    started = time.monotonic()
    config=configuration or {}
    f = config['map_factory'](digits) if configuration else IntervalHalfMap(digits)
    source = config.get('source',OUT / "connection-center.json")
    center = json.loads(source.read_text())
    local_source = center if configuration else json.loads((OUT / "local-certificate.json").read_text())
    a0 = arb(center["a"])
    b0 = arb(center["beta"]) + arb(beta_shift)
    v0 = arb(0)
    V = arb_mat([[arb(x) for x in row] for row in local_source["coordinates"]])
    W = V.inv()
    ra = arb(config.get('ra',"1e-14"))
    rb = arb(config.get('rb',"5e-14"))
    rv = arb(config.get('rv',"8e-11"))
    delta = arb(config.get('delta',"1e-7"))
    A0 = arb(a0, ra)
    B0 = arb(b0, rb)
    V0 = arb(0, rv)
    pcenter = config['saddle_center'](f,b0,center['parent_guess']) if configuration else saddle_center(f, b0)
    floor=config.get('inflation_floor','1e-60')
    p0, pb0, pbb0, point_cert = enclose_saddle(f, pcenter, b0, config.get('point_radius',"1e-55"),inflation_floor=floor)
    pB, pbB, pbbB, branch_cert = enclose_saddle(f, pcenter, B0, "1e-8",inflation_floor=floor)
    D, H, local = local_coordinates(f, pB, B0, V, W, pbB, delta)
    if config.get('centered_local_derivative',False):
        # Mean-value enclosure at each fixed beta. H bounds the complete
        # segment from z=0 to z in the local square; pB encloses p(beta).
        _, Jcenter, _, _ = f.evaluate(list(pB.entries()), B0)
        Dcenter = W * arb_mat([[Jcenter[i,j] for j in range(2)] for i in range(2)]) * V
        for i in range(2):
            for j in range(2):
                enclosure=Dcenter[i,j]+sum(H[i][j,k]*arb(0,delta) for k in range(2))
                D[i,j]=D[i,j].intersection(enclosure)
        local['centered_derivative']=[[bounds(D[i,j]) for j in range(2)] for i in range(2)]
    # Exact decimal guard constants; every inequality is checked in Arb.
    guards=config.get('guards',{})
    a = arb(guards.get('a',"1.6469"))
    b = arb(guards.get('b',".00035"))
    c = arb(guards.get('c',".00024"))
    d = arb(guards.get('d',".00032"))
    M = arb(guards.get('M','1500'))
    P = arb(guards.get('P',".000041"))
    L = arb(guards.get('L',".001"))
    K = arb(guards.get('K','2000'))
    Gamma = arb(guards.get('Gamma',".0001"))
    local_gates = dict(
        expansion=D[0, 0] < -a,
        cross_u=abs(D[0, 1]) < b,
        cross_s=abs(D[1, 0]) < c,
        contraction=abs(D[1, 1]) < d,
        second_derivatives=all(
            abs(H[k][i, j]) < M for k in range(2) for i in range(2) for j in range(2)
        ),
        parameter_derivatives=all(
            2 * abs(H[k][i, 2]) * delta < P for k in range(2) for i in range(2)
        ),
    )
    graph_bounds = dict(
        unstable_slope=(c + d * L) / (a - b * L),
        stable_slope=(b + L * d) / (a - L * c),
        unstable_graph_contraction=d + L * b,
        stable_graph_contraction=1 / (a - L * c),
        unstable_curvature=(M * (1 + L) ** 3 + (d + L * b) * K) / (a - b * L) ** 2,
        stable_curvature=(M * (1 + L) ** 3 + K * (c * L + d) ** 2) / (a - L * c),
        unstable_parameter=P * (1 + L) + (d + L * b) * Gamma,
        stable_parameter=(Gamma + P * (1 + L)) / (a - L * c),
    )
    graph_gates = dict(
        unstable_projection_covers=a - b * L > 1,
        stable_image_inside=c * L + d < 1,
        unstable_slope=graph_bounds["unstable_slope"] < L,
        stable_slope=graph_bounds["stable_slope"] < L,
        unstable_contraction=graph_bounds["unstable_graph_contraction"] < 1,
        stable_contraction=graph_bounds["stable_graph_contraction"] < 1,
        unstable_curvature=graph_bounds["unstable_curvature"] < K,
        stable_curvature=graph_bounds["stable_curvature"] < K,
        unstable_parameter=graph_bounds["unstable_parameter"] < Gamma,
        stable_parameter=graph_bounds["stable_parameter"] < Gamma,
    )
    if not all(local_gates.values()) or not all(graph_gates.values()):
        raise ValueError("Local invariant-graph certificate failed: "+str((local_gates,graph_gates,local)))
    g0, J0, H0, refs = connect(
        f, a0, v0, b0, b0, p0, pb0, pb0, pbb0, V, W, center["iterations"]
    )
    g, J, H, records = connect(
        f,
        A0,
        V0,
        B0,
        b0,
        p0,
        pb0,
        pbB,
        pbbB,
        V,
        W,
        center["iterations"],
        refs,
        verify_path=True,
    )
    da = A0 - a0
    dv = V0
    db = B0 - b0
    offsets = [da, dv, db]
    # For the exact unstable graph v=psi(a,beta) and stable graph u=chi(v,beta).
    psia = arb(0, L)
    psiaa = arb(0, K)
    psib = arb(0, Gamma)
    chis = arb(0, L)
    chiss = arb(0, K)
    chib = arb(0, Gamma)
    stable_graph_value = arb(0, L * g[1, 0].abs_upper())
    face_values = {}
    for label, sign in [("beta_lower", -1), ("beta_upper", 1)]:
        dd = [da, dv, sign * rb]
        enclosure = g0[0, 0] - stable_graph_value
        for i in range(3):
            enclosure += J0[0, i] * dd[i]
            for j in range(3):
                enclosure += H[0][i, j] * dd[i] * dd[j] / 2
        face_values[label] = enclosure
    # A derivative Taylor enclosure along the segment from the center to a face.
    for label, sign in [("a_lower", -1), ("a_upper", 1)]:
        Ua = J0[0, 0] + H[0][0, 0] * (sign * ra) + H[0][0, 1] * dv + H[0][0, 2] * db
        face_values[label] = Ua + J[0, 1] * psia - chis * (J[1, 0] + J[1, 1] * psia)
    tangent_s = J[1, 0] + J[1, 1] * psia
    Uaa = (
        H[0][0, 0] + 2 * H[0][0, 1] * psia + H[0][1, 1] * psia * psia + J[0, 1] * psiaa
    )
    Saa = (
        H[1][0, 0] + 2 * H[1][0, 1] * psia + H[1][1, 1] * psia * psia + J[1, 1] * psiaa
    )
    Gaa = Uaa - chis * Saa - chiss * tangent_s * tangent_s
    Gb = J[0, 2] + J[0, 1] * psib - chis * (J[1, 2] + J[1, 1] * psib) - chib
    departure = []
    for k, r in enumerate(records):
        z = W * (r["state"] - pB)
        if z[0, 0].abs_lower() > delta:
            departure.append(dict(iterate=k + 1, unstable_coordinate=bounds(z[0, 0])))
    eig = branch_cert["half_eigenvalues"]
    # Use full interval evaluation of the derivative, not rounded report text.
    _, JP, _, _ = f.evaluate(list(pB.entries()), B0)
    tr = JP[0, 0] + JP[1, 1]
    det = JP[0, 0] * JP[1, 1] - JP[0, 1] * JP[1, 0]
    root = (tr * tr - 4 * det).sqrt()
    mu = (tr - root) / 2
    nu = (tr + root) / 2
    parameter_sign=arb(config.get('parameter_sign','1'))
    curvature_sign=arb(config.get('curvature_sign','1'))
    gates = dict(
        local=all(local_gates.values()),
        graphs=all(graph_gates.values()),
        initial_on_graph_box=L * A0.abs_upper() < rv,
        initial_inside_local_box=A0.abs_upper() < delta,
        final_inside_stable_graph_domain=g[1, 0].abs_upper() < delta,
        negative_lower_parameter=parameter_sign*face_values["beta_lower"] < 0,
        positive_upper_parameter=parameter_sign*face_values["beta_upper"] > 0,
        negative_left_derivative=curvature_sign*face_values["a_lower"] < 0,
        positive_right_derivative=curvature_sign*face_values["a_upper"] > 0,
        quadratic=curvature_sign*Gaa > 0,
        transverse_parameter=parameter_sign*Gb > 0,
        nontrivial_excursion=bool(departure),
        saddle_unstable=mu < -1,
        saddle_stable=(abs(nu)<1 and not nu.contains(0)) if configuration else (nu > 0 and nu < 1),
        dissipative=abs(det) < 1,
    )
    result = dict(
        status="passed" if all(gates.values()) else "failed",
        decimal_digits=digits,
        precision_bits=ctx.prec,
        beta_shift=beta_shift,
        software=dict(
            python=platform.python_version(),
            python_flint=version("python-flint"),
            numpy=version("numpy"),
            scipy=version("scipy"),
        ),
        parameter_interval=bounds(B0),
        initial_coordinate_interval=bounds(A0),
        initial_stable_enclosure=bounds(V0),
        center_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        return_iterations=center["iterations"],
        fixed_parameter_saddle=point_cert,
        parameter_family_saddle=branch_cert,
        local_map=local,
        local_gates=local_gates,
        graph_constants=dict(
            domain=str(delta), slope=str(L), curvature=str(K), parameter_derivative=str(Gamma)
        ),
        graph_transform_bounds={k: bounds(v) for k, v in graph_bounds.items()},
        graph_gates=graph_gates,
        return_value=[bounds(x) for x in g.entries()],
        return_jacobian=[[bounds(J[i, j]) for j in range(3)] for i in range(2)],
        center_return_value=[bounds(x) for x in g0.entries()],
        center_return_jacobian=[[bounds(J0[i, j]) for j in range(3)] for i in range(2)],
        center_return_hessians=[
            [[bounds(h[i, j]) for j in range(3)] for i in range(3)] for h in H0
        ],
        return_hessians=[
            [[bounds(h[i, j]) for j in range(3)] for i in range(3)] for h in H
        ],
        face_signs={k: bounds(v) for k, v in face_values.items()},
        normal_curvature=bounds(Gaa),
        parameter_unfolding=bounds(Gb),
        full_return_multipliers=[bounds(mu * mu), bounds(nu * nu)],
        full_return_determinant=bounds(det * det),
        departure_witness=departure[0] if departure else None,
        itinerary_cells=[r["path_cells"] for r in records],
        gates=gates,
        elapsed_seconds=time.monotonic() - started,
        claim=(
            "The inequalities certify a quadratic homoclinic tangency of the local half-return map, transverse in beta, subject to the documented graph-transform and event-map implementation."
            if all(gates.values())
            else "At least one sufficient inequality failed; no tangency claim."
        ),
    )
    if configuration:
        result['moving_frame']=getattr(f,'moving_frame',False)
        result['centered_local_derivative']=bool(config.get('centered_local_derivative',False))
        result['claim']='The inequalities certify a quadratic homoclinic tangency of the configured return map, transverse in beta, subject to the documented graph-transform and event-map implementation.' if all(gates.values()) else 'At least one sufficient inequality failed; no tangency claim.'
        result['map_description']=config.get('map_description','Configured transverse return')
        result['sign_conventions']=dict(parameter=str(parameter_sign),curvature=str(curvature_sign),
            note='Face gates use the displayed orientation signs. Recorded face bounds and derivatives retain their raw coordinate signs.')
        result['parameters']=center.get('parameters')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--digits", type=int, default=100)
    p.add_argument(
        "--beta-shift", default="0", help="Shift the tested box for a negative control."
    )
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    report = certify(args.digits, args.beta_shift)
    path = args.output or OUT / f"tangency-certificate-{args.digits}.json"
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in [
                    "status",
                    "parameter_interval",
                    "face_signs",
                    "normal_curvature",
                    "parameter_unfolding",
                    "gates",
                    "elapsed_seconds",
                ]
            },
            indent=2,
        ),
        flush=True,
    )
    if report["status"] != "passed":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
