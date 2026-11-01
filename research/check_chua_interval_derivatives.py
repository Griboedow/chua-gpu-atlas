"""Independent mpmath event solves check the enclosed return derivatives.

These are numerical implementation controls, separate from interval inclusion
tests. Centered differences at two scales check the analytic first and second
event derivatives, including beta derivatives.
"""

import json, time
import mpmath as mp
from flint import arb
from chua_interval_map import IntervalHalfMap
from certify_chua_saddle import OUT


def mp_return(w, guesses):
    y, z, beta = w
    state = mp.matrix([0, y, z])
    for region, target, guess in zip([0, 1, 0, 1, 0], [1, 1, 1, 1, 0], guesses):
        slope = -mp.mpf(8) / 7 if region == 0 else -mp.mpf(5) / 7
        A = mp.matrix([[-9 * (1 + slope), 9, 0], [1, -1, 1], [0, -beta, 0]])
        c = mp.matrix([mp.mpf(3) * region / 2, 0, -mp.mpf(3) * region / 2])
        start = state - c

        def orbit(t):
            return c + mp.expm(A * t) * start

        t = mp.findroot(
            lambda t: orbit(t)[0] - target,
            (mp.mpf(str(guess)) * 0.9999, mp.mpf(str(guess)) * 1.0001),
            solver="secant",
            tol=mp.mpf("1e-85"),
        )
        state = orbit(t)
        state[0] = target
    return -state[1, 0], -state[2, 0]


def main():
    started = time.monotonic()
    mp.mp.dps = 100
    f = IntervalHalfMap(100)
    rows = []
    # Both a local point and a point on the outgoing excursion.
    for name, wstr in [
        (
            "saddle_neighborhood",
            ["0.23909563991385", "-0.4390043251667", "14.30635909996146"],
        ),
        ("outgoing_neighborhood", ["0.2386", "-0.4395", "14.30635909996146"]),
    ]:
        w = [mp.mpf(s) for s in wstr]
        a = [arb(s) for s in wstr]
        value, J, H, arcs = f.evaluate(a[:2], a[2])
        guesses = [float(z["duration"].mid()) for z in arcs]
        base = mp_return(w, guesses)
        scales = []
        for hs in ["1e-15", "3e-16"]:
            h = mp.mpf(hs)
            first = []
            second = []
            cache = {}

            def ev(shift):
                key = tuple(shift)
                if key not in cache:
                    cache[key] = mp_return(
                        [w[i] + h * shift[i] for i in range(3)], guesses
                    )
                return cache[key]

            for i in range(3):
                e = [0, 0, 0]
                e[i] = 1
                left = ev([-x for x in e])
                right = ev(e)
                for k in range(2):
                    jac = (right[k] - left[k]) / (2 * h)
                    first.append(
                        abs(jac - mp.mpf(J[k, i].mid().str(95, radius=False)))
                        / (1 + abs(jac))
                    )
                for j in range(i, 3):
                    if i == j:
                        d2 = [
                            (right[k] - 2 * base[k] + left[k]) / (h * h)
                            for k in range(2)
                        ]
                    else:
                        shifts = []
                        for si, sj in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                            sh = [0, 0, 0]
                            sh[i] = si
                            sh[j] = sj
                            shifts.append(ev(sh))
                        d2 = [
                            (shifts[0][k] - shifts[1][k] - shifts[2][k] + shifts[3][k])
                            / (4 * h * h)
                            for k in range(2)
                        ]
                    for k in range(2):
                        second.append(
                            abs(d2[k] - mp.mpf(H[k][i, j].mid().str(95, radius=False)))
                            / (1 + abs(d2[k]))
                        )
            scales.append(
                dict(
                    step=hs,
                    maximum_scaled_jacobian_error=mp.nstr(max(first), 15),
                    maximum_scaled_hessian_error=mp.nstr(max(second), 15),
                )
            )
            print(name, scales[-1], flush=True)
        value_error = max(
            abs(base[k] - mp.mpf(value[k, 0].mid().str(95, radius=False)))
            for k in range(2)
        )
        rows.append(
            dict(
                point=name,
                coordinates=wstr,
                value_error=mp.nstr(value_error, 15),
                differences=scales,
                passed=bool(
                    value_error < mp.mpf("1e-65")
                    and all(
                        mp.mpf(s["maximum_scaled_jacobian_error"]) < mp.mpf("1e-20")
                        and mp.mpf(s["maximum_scaled_hessian_error"]) < mp.mpf("1e-20")
                        for s in scales
                    )
                ),
            )
        )
    report = dict(
        passed=all(r["passed"] for r in rows),
        checks=rows,
        elapsed_seconds=time.monotonic() - started,
        scope="Independent numerical derivative controls, not themselves interval proofs.",
    )
    (OUT / "derivative-controls.json").write_text(json.dumps(report, indent=2) + "\n")
    if not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
