"""Locate a tangent-line root to seed the separate interval certificate.

The Newton iteration is numerical. Its residual is never used as evidence
that the exact stable and unstable manifolds touch.
"""

import json
from flint import arb, arb_mat
from chua_interval_map import IntervalHalfMap, col, midpoint, bounds
from certify_chua_saddle import OUT, saddle_center, enclose_saddle
from chua_interval_connection import connect


def main():
    f = IntervalHalfMap(100)
    source = json.loads((OUT / "local-certificate.json").read_text())
    V = arb_mat([[arb(x) for x in row] for row in source["coordinates"]])
    W = V.inv()
    beta = arb("14.3063590999613")
    a = -arb(".0006075330794225644") / arb("2.7134317") ** 9
    history = []
    for iteration in range(10):
        center = saddle_center(f, beta)
        p, pb, pbb, _ = enclose_saddle(f, center, beta, "1e-55")
        g, J, H, _ = connect(f, a, arb(0), beta, beta, p, pb, pb, pbb, V, W, 27)
        residual = col([g[0, 0], J[0, 0]])
        history.append(
            dict(
                iteration=iteration,
                a=bounds(a),
                beta=bounds(beta),
                gap=bounds(residual[0, 0]),
                slope=bounds(residual[1, 0]),
            )
        )
        print(iteration, "gap", residual[0, 0], "slope", residual[1, 0], flush=True)
        if abs(residual[0, 0]) < arb("1e-35") and abs(residual[1, 0]) < arb("1e-25"):
            break
        jac = arb_mat([[J[0, 0], J[0, 2]], [H[0][0, 0], H[0][0, 2]]])
        step = midpoint(jac.inv() * residual)
        a = (a - step[0, 0]).mid()
        beta = (beta - step[1, 0]).mid()
    else:
        (OUT / "center-reconstruction-failure.json").write_text(
            json.dumps(history, indent=2) + "\n", encoding="utf-8"
        )
        raise ValueError("Tangent-line Newton iteration did not converge")
    # Keep the recorded certificate center fixed; report whether it is reproduced.
    frozen = json.loads((OUT / "connection-center.json").read_text())
    agreement = abs(a - arb(frozen["a"])) < arb("1e-35") and abs(
        beta - arb(frozen["beta"])
    ) < arb("1e-35")
    result = dict(
        passed=agreement,
        a=bounds(a),
        beta=bounds(beta),
        history=history,
        scope="Numerical center reconstruction; not a manifold existence proof.",
    )
    (OUT / "center-reconstruction.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print("Center reconstruction:", "passed" if agreement else "failed", flush=True)
    if not agreement:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
