"""Independent mpmath controls of the 27-fold composition and its derivatives."""

import json
import time
import mpmath as mp
from flint import arb, arb_mat
from chua_interval_map import IntervalHalfMap
from chua_interval_connection import connect
from certify_chua_saddle import OUT, saddle_center, enclose_saddle
from check_chua_interval_derivatives import mp_return


def number(x):
    return mp.mpf(x.mid().str(95, radius=False))


def main():
    started = time.monotonic()
    mp.mp.dps = 100
    f = IntervalHalfMap(100)
    center = json.loads((OUT / "connection-center.json").read_text())
    local = json.loads((OUT / "local-certificate.json").read_text())
    a, beta = arb(center["a"]), arb(center["beta"])
    V = arb_mat([[arb(x) for x in row] for row in local["coordinates"]])
    p, pb, pbb, _ = enclose_saddle(f, saddle_center(f, beta), beta, "1e-55")
    value, J, H, records = connect(
        f, a, arb(0), beta, beta, p, pb, pb, pbb, V, V.inv(), 27
    )
    starts = [p + V * arb_mat([[a], [0]])] + [r["state"] for r in records[:-1]]
    guesses = [f.schedule(list(s.entries()), beta) for s in starts]
    Vm = mp.matrix([[number(V[i, j]) for j in range(2)] for i in range(2)])
    Wm = Vm**-1
    pm, pbm, pbbm = [mp.matrix([number(x) for x in q.entries()]) for q in (p, pb, pbb)]
    am, bm = number(a), number(beta)

    def normal_gap(da, db):
        # This Taylor branch has the same first two parameter derivatives at 0.
        saddle = pm + pbm * db + pbbm * (db * db / 2)
        state = saddle + Vm * mp.matrix([am + da, 0])
        for ts in guesses:
            state = mp.matrix(mp_return([state[0], state[1], bm + db], ts))
        return (Wm * (state - saddle))[0]

    base = normal_gap(mp.mpf(0), mp.mpf(0))
    rows = []
    for hs in ["1e-18", "3e-19"]:
        h = mp.mpf(hs)
        left, right = normal_gap(-h, 0), normal_gap(h, 0)
        lo_beta, hi_beta = normal_gap(0, -h), normal_gap(0, h)
        curvature = (right - 2 * base + left) / (h * h)
        row = dict(
            step=hs,
            value_error=mp.nstr(abs(base - number(value[0, 0])), 16),
            a_derivative_error=mp.nstr(
                abs((right - left) / (2 * h) - number(J[0, 0])), 16
            ),
            beta_derivative_error=mp.nstr(
                abs((hi_beta - lo_beta) / (2 * h) - number(J[0, 2])), 16
            ),
            relative_curvature_error=mp.nstr(
                abs(curvature / number(H[0][0, 0]) - 1), 16
            ),
        )
        row["passed"] = all(
            mp.mpf(row[k]) < mp.mpf("1e-12")
            for k in [
                "a_derivative_error",
                "beta_derivative_error",
                "relative_curvature_error",
            ]
        )
        row["passed"] = row["passed"] and mp.mpf(row["value_error"]) < mp.mpf("1e-40")
        rows.append(row)
        print(row, flush=True)
    report = dict(
        passed=all(r["passed"] for r in rows),
        checks=rows,
        elapsed_seconds=time.monotonic() - started,
        scope="Independent numerical controls of the full composition; not interval proofs.",
    )
    (OUT / "global-derivative-controls.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    if not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
