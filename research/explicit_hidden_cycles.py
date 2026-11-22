"""Interval shooting and equilibrium exclusion on the Fig. 6 positive-loss Chua slice."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
import json
from pathlib import Path

import numpy as np
from flint import acb, arb, arb_mat, ctx

from certify_chua_coexistence_box import MultipleShooting
from certify_chua_hidden_cycle import lyapunov, outer_path
from chua_interval_map import IntervalHalfMap, bounds, col, eye, midpoint

OUT = Path(__file__).resolve().parents[1] / "results/chua-explicit-hidden"
PARAMETERS = ("8.41", "12.23", "-1.366", "-0.12", "0.0435")
# Floating-point guesses from direct integration. The subsequent enclosure
# validates the closed orbit without relying on the accuracy of these guesses.
SEEDS = {
    "right_lobe": ([0.2095, 0.50766], [(1, 1, 1.275, 1), (0, 1, 0.810, 1)]),
    "outer": (
        [5.62, 30.7],
        [
            (1, 1, 0.924, 1),
            (0, 0, 0.026, 1),
            (0, -1, 0.022, -1),
            (-1, -1, 0.924, -1),
            (0, 0, 0.026, -1),
            (0, 1, 0.022, 1),
        ],
    ),
}


class Shooting(MultipleShooting):
    def __init__(self, schedule):
        self.schedule = schedule
        self.n = len(schedule)
        self.start_x = [1] + [row[1] for row in schedule[:-1]]
        self.original = IntervalHalfMap.__new__(IntervalHalfMap)


def cycle(label, parameters, displaced=False):
    guess, schedule = SEEDS[label]
    germ = Shooting(schedule)
    nominal = [p.mid() for p in parameters]
    state = col([1, *guess])
    raw = []
    for region, target, guess_time, sign in germ.schedule:
        duration = arb(str(guess_time))
        raw.extend([state[1, 0], state[2, 0], duration])
        matrix, center = germ.matrices(nominal, region)
        state = germ.original.flow(matrix, center, state, duration)[0]
        state[0, 0] = target
    center = col(raw)
    for _ in range(9):
        residual, jacobian, _ = germ.evaluate(center, nominal)
        center = midpoint(center - jacobian.inv() * residual)
    if displaced:
        center[0, 0] += arb("0.01")
    residual, jacobian, _ = germ.evaluate(center, nominal)
    inverse = midpoint(jacobian.inv())
    beta_derivative = germ.parameter_jacobian(center, parameters)
    forcing = -inverse * residual - (inverse * beta_derivative) * col(
        [p - p.mid() for p in parameters]
    )
    scale = max(
        arb("1e-24"),
        max(v.abs_upper() for v in forcing.entries())
        / (20 if label == "right_lobe" else 4000),
    )
    widths = [max(scale, v.abs_upper()) * 3 for v in forcing.entries()]
    if displaced:
        widths = [arb("1e-5")] * len(widths)
    box = col([arb(v, r) for v, r in zip(center.entries(), widths)])
    _, jacobian, arcs = germ.evaluate(box, parameters)
    remainder = eye(box.nrows()) - inverse * jacobian
    image = center + forcing + remainder * (box - center)
    contraction = max(
        (
            sum(abs(remainder[i, j]) * widths[j] for j in range(box.nrows()))
            / widths[i]
        ).abs_upper()
        for i in range(box.nrows())
    ).abs_upper()
    if not (
        contraction < 1
        and all(v.contains_interior(w) for v, w in zip(box.entries(), image.entries()))
    ):
        raise ValueError(
            f"{label} inclusion failed, contraction upper bound {float(contraction):.6g}"
        )
    if label == "outer":
        # Inversion followed by a half-period phase shift maps the orbit
        # into the same uniqueness box, proving its central symmetry.
        half = germ.n // 2
        for i in range(germ.n):
            reflected = (i + half) % germ.n
            for component, sign in enumerate((-1, -1, 1)):
                assert box[3 * i + component, 0].contains_interior(
                    sign * image[3 * reflected + component, 0]
                )
    derivative = eye(2)
    log_det = arb(0)
    cells = []
    for arc in arcs:
        speed = arc["v"][0, 0]
        assert arc["duration"] > 0 and not speed.contains(0)
        matrix, velocity = arc["E"], arc["v"]
        event = arb_mat(
            [
                [matrix[i, j] - velocity[i, 0] * matrix[0, j] / speed for j in (1, 2)]
                for i in (1, 2)
            ]
        )
        derivative = event * derivative
        log_det += arc["A"].trace() * arc["duration"]
        sign = arc["sign"]
        record = dict(
            region=abs(arc["region"]),
            target=sign * arc["target"],
            start=sign * arc["start"],
            duration=arc["duration"],
        )
        cells.append(
            germ.original.verify_arc(arc["A"], sign * arc["c"], record)[
                "accepted_cells"
            ]
        )
    determinant, trace = log_det.exp(), derivative.trace()
    margins = [1 - trace + determinant, 1 + trace + determinant, 1 - determinant]
    assert all(v > 0 for v in margins), (label, [str(v) for v in margins])
    # The first and last events are respectively inward and outward at x=1.
    assert arcs[0]["v"][0, 0] < 0 and arcs[-1]["v"][0, 0] > 0
    assert sum(a["target"] == 1 and a["v"][0, 0] > 0 for a in arcs) == 1
    return dict(
        inversion_symmetric=(label == "outer"),
        period=bounds(sum(arc["duration"] for arc in arcs)),
        section_box=[bounds(box[i, 0]) for i in (0, 1)],
        itinerary=germ.schedule,
        contraction=bounds(contraction),
        schur_margins=[bounds(v) for v in margins],
        period_trace=bounds(trace),
        period_determinant=bounds(determinant),
        path_cells=cells,
        shooting_center=[bounds(v) for v in center.entries()],
        shooting_box=[bounds(v) for v in box.entries()],
        krawczyk_image=[bounds(v) for v in image.entries()],
    )


def spectrum(matrix):
    a1 = -matrix.trace()
    a2 = sum(
        matrix[i, i] * matrix[j, j] - matrix[i, j] * matrix[j, i]
        for i, j in ((0, 1), (0, 2), (1, 2))
    )
    a3 = -matrix.det()
    coefficients = [arb(1), a1, a2, a3]
    lambdas = np.roots([float(v.mid()) for v in coefficients])
    real = next(v.real for v in lambdas if abs(v.imag) < 1e-10)
    t = arb(str(real))
    polynomial = lambda x: ((x + a1) * x + a2) * x + a3
    derivative = lambda x: (3 * x + 2 * a1) * x + a2
    for _ in range(6):
        t = (t - polynomial(t) / derivative(t)).mid()
    for radius in ("1e-22", "1e-12", "1e-7", "1e-5", ".001", ".1"):
        box = arb(t, arb(radius))
        if derivative(box).contains(0):
            continue
        root = t - polynomial(t) / derivative(box)
        if box.contains_interior(root):
            break
    else:
        raise ValueError("Real eigenvalue did not enclose")
    real = root
    sigma = (-a1 - real) / 2
    omega = (-a3 / real - sigma * sigma).sqrt()
    assert omega > 0 and sigma < 0
    return real, sigma, omega


class SpectralFlow:
    """Outward-rounded exponential of a three-dimensional real focus matrix."""

    def __init__(self, matrix):
        self.real, self.sigma, self.omega = spectrum(matrix)
        I = eye(3)
        a = matrix - self.sigma * I
        self.P = (a * a + self.omega**2 * I) / (
            (self.real - self.sigma) ** 2 + self.omega**2
        )
        self.Q = I - self.P
        self.B = a * self.Q / self.omega

    def flow(self, A, c, state, time):
        angle = self.omega * time
        E = (self.real * time).exp() * self.P + (self.sigma * time).exp() * (
            angle.cos() * self.Q + angle.sin() * self.B
        )
        return c + E * (state - c), E


def exclude_equilibria(parameters):
    a, b, m0, m1, g = parameters
    A0, _ = MultipleShooting.matrices(parameters, 0)
    A, c = MultipleShooting.matrices(parameters, 1)
    assert m0 < -1 and m1 > -1 and c[0, 0] > 1
    lam, sig, w = spectrum(A0)
    assert lam > 0

    def vector(root):
        y = 1 + m0 + root / a
        return [1, y, -b * y / (root + g)]

    u = col(vector(lam))
    v = vector(acb(sig, w))
    S = arb_mat(
        [[1, 1, 0], [u[1, 0], v[1].real, v[1].imag], [u[2, 0], v[2].real, v[2].imag]]
    )
    delta = arb("1e-6")
    inv = S.inv()
    r_inner = (delta / (2 * sum(x * x for x in inv.entries()).sqrt())).lower()
    # At first positive exit, u_coord = 1 - s_1. The stable plane
    # contracts in Euclidean modal coordinates up to either exit.
    e = arb(0, delta)
    start = col(
        [
            1,
            u[1, 0] + (S[1, 1] - u[1, 0]) * e + S[1, 2] * e,
            u[2, 0] + (S[2, 1] - u[2, 0]) * e + S[2, 2] * e,
        ]
    )
    P, lyap = lyapunov(A)
    level = (c[0, 0] - 1) ** 2 / (2 * P.inv()[0, 0])
    flow = SpectralFlow(A)
    assert flow.real < 0
    for duration in (10, 20, 40, 80, 160, 320, 640):
        state, _ = flow.flow(A, c, start, arb(duration))
        d = state - c
        energy = (d.transpose() * P * d)[0, 0]
        if energy < level:
            break
    else:
        raise ValueError("No capture")
    cells = outer_path(flow, A, c, start, arb(duration))
    rowsum = max(sum(abs(P[i, j]) for j in range(3)).abs_upper() for i in range(3))
    return dict(
        origin_radius=bounds(r_inner),
        outer_radius=bounds((level / (2 * rowsum)).sqrt().lower()),
        eigenvalues=[bounds(v) for v in (lam, sig, w)],
        modal_basis=[[bounds(S[i, j]) for j in range(3)] for i in range(3)],
        exit_state=[bounds(v) for v in start.entries()],
        outer_path_cells=cells,
        outer_time=duration,
        capture_energy=bounds(energy),
        ellipsoid_level=bounds(level),
        outer_lyapunov=lyap,
    )


def tile(task):
    center, radius, digits = task
    ctx.dps = digits
    par = [arb(v) for v in PARAMETERS]
    par[1] = arb(arb(center), arb(radius))
    cycles = {k: cycle(k, par) for k in SEEDS}
    hidden = exclude_equilibria(par)

    def interval(record):
        return [record["lower"], record["upper"]]

    return dict(
        center=center,
        periods={k: interval(v["period"]) for k, v in cycles.items()},
        contraction_upper={k: v["contraction"]["upper"] for k, v in cycles.items()},
        schur_lower={
            k: [m["lower"] for m in v["schur_margins"]] for k, v in cycles.items()
        },
        origin_radius_lower=hidden["origin_radius"]["lower"],
        outer_radius_lower=hidden["outer_radius"]["lower"],
        outer_capture_time=hidden["outer_time"],
    )


def cover(digits, workers, output):
    lo, hi, step = map(Decimal, ("12.21775", "12.24225", "0.00005"))
    radius = step / 2
    count = int((hi - lo) / step)
    assert lo + count * step == hi
    jobs = [(str(lo + i * step + radius), str(radius), digits) for i in range(count)]
    if workers == 1:
        rows = list(map(tile, jobs))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(tile, jobs, chunksize=4))
    assert len(rows) == count
    out = output or OUT / "beta-cover.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    record = dict(
        beta_interval=[str(lo), str(hi)],
        cell_radius=str(radius),
        digits=digits,
        nominal=PARAMETERS,
        cell_count=count,
        rows=rows,
        passed=True,
    )
    out.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print("Certified", count, "touching cells over", lo, hi)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--digits", type=int, default=70)
    parser.add_argument("--beta-radius", default="0.000025")
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--cover", action="store_true", help="Certify the entire tiled beta interval"
    )
    parser.add_argument("--workers", type=int, default=1)
    options = parser.parse_args()
    if not __debug__:
        raise RuntimeError("Assertions are required for the certificate")
    if options.cover:
        cover(options.digits, options.workers, options.output)
        return
    ctx.dps = options.digits
    par = [arb(value) for value in PARAMETERS]
    par[1] = arb(par[1], arb(options.beta_radius))
    results = {k: cycle(k, par) for k in SEEDS}
    hidden = exclude_equilibria(par)
    try:
        cycle("right_lobe", par, displaced=True)
    except ValueError as error:
        if "inclusion failed" not in str(error):
            raise
    else:
        raise AssertionError("Displaced shooting box accepted")
    record = dict(
        displaced_box_rejected=True,
        digits=ctx.dps,
        parameter_order=["alpha", "beta", "m0", "m1", "gamma"],
        nominal=PARAMETERS,
        beta_radius=options.beta_radius,
        cycles=results,
        equilibrium_exclusion=hidden,
        passed=True,
    )
    out = options.output or OUT / f"certificate-{ctx.dps}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    for k, r in results.items():
        print(k, r["period"]["ball"], r["contraction"]["ball"], flush=True)
    print("Capture", hidden["origin_radius"]["ball"], hidden["outer_time"], flush=True)


if __name__ == "__main__":
    main()
