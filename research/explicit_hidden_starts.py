"""Validate voltage-only starts by interval flights and a contracting return box."""

import argparse
import json

import numpy as np
from flint import arb, arb_mat, ctx

from explicit_hidden_cycles import OUT, PARAMETERS, SEEDS, Shooting, SpectralFlow, cycle
from certify_chua_coexistence_box import enclosure
from chua_interval_map import IntervalHalfMap, bounds, col, midpoint


class Capture:
    def __init__(self, name, parameters, source):
        self.field = IntervalHalfMap.__new__(IntervalHalfMap)
        self.schedule = SEEDS[name][1]
        self.matrices = {r: Shooting.matrices(parameters, r) for r in (-1, 0, 1)}
        self.flows = {r: SpectralFlow(ac[0]) for r, ac in self.matrices.items()}
        self.center = col([arb(v["midpoint"]) for v in source["shooting_center"][:2]])
        _, matrix = self.map(self.center, self.schedule)
        B = np.linalg.eig(
            np.array([[float(matrix[i, j].mid()) for j in range(2)] for i in range(2)])
        )[1]
        B = np.real_if_close(B)
        if np.iscomplexobj(B):
            B = np.column_stack((B[:, 0].real, B[:, 0].imag))
        if name == "right_lobe":
            # Balance the shear of the section derivative in a well-conditioned
            # chart with one eigenvector and a scaled perpendicular.
            v = B[:, 0] / np.linalg.norm(B[:, 0])
            B = np.column_stack((v, 0.08 * np.array([-v[1], v[0]])))
        self.B = arb_mat(
            [[arb(format(B[i, j], ".15g")) for j in range(2)] for i in range(2)]
        )
        self.inverse = self.B.inv()

    def path(self, A, c, start, time, target, region, sign):
        """Include endpoint collars extending across several subdivision cells."""
        start, c = sign * start, sign * c
        target = sign * target
        lo = 1 if region else 0
        hi = None if region else 1
        work = [(arb(0), arb(1), 0)]
        while work:
            u, v, depth = work.pop()
            x = self.field.flow(A, c, start, u.union(v) * time)[0][0, 0]
            low = x > lo
            high = hi is None or x < hi
            if not (low and high) and u < arb("0.4"):
                span = self.field.flow(A, c, start, arb(0).union(v) * time)[0]
                dx = (A * (span - c))[0, 0]
                low = low or (start[0, 0] >= lo and dx > 0)
                high = high or (hi is not None and start[0, 0] <= hi and dx < 0)
            if not (low and high) and v > arb("0.6"):
                span = self.field.flow(A, c, start, u.union(arb(1)) * time)[0]
                dx = (A * (span - c))[0, 0]
                low = low or (target >= lo and dx < 0)
                high = high or (hi is not None and target <= hi and dx > 0)
            if low and high:
                continue
            if depth >= 14:
                raise ValueError(f"Unresolved full-time flight in region {region}")
            mid = (u + v) / 2
            work.extend([(u, mid, depth + 1), (mid, v, depth + 1)])

    def flights(self, state, schedule, initial_derivative=None):
        jac = (
            initial_derivative
            if initial_derivative is not None
            else arb_mat([[0, 0], [1, 0], [0, 1]])
        )
        for region, target, guess, sign in schedule:
            A, c = self.matrices[region]
            self.field.flow = self.flows[region].flow
            guess_time = arb(str(guess))
            nominal_state = midpoint(state)
            for _ in range(6):
                hit = self.field.flow(A, c, nominal_state, guess_time)[0]
                guess_time = (
                    guess_time - (hit[0, 0] - target) / (A * (hit - c))[0, 0]
                ).mid()
            time = self.field.crossing(A, c, state, target, float(guess_time))
            self.path(A, c, state, time, target, region, sign)
            end, exponential = self.field.flow(A, c, state, time)
            velocity = A * (end - c)
            assert not velocity[0, 0].contains(0)
            propagated = exponential * jac
            jac = arb_mat(
                [
                    [
                        propagated[i, j]
                        - velocity[i, 0] * propagated[0, j] / velocity[0, 0]
                        for j in range(2)
                    ]
                    for i in range(3)
                ]
            )
            # A centered enclosure retains time/state cancellation at a crossing.
            if any(v.rad() > arb("1e-30") for v in state.entries()):
                tcenter = self.field.crossing(
                    A, c, nominal_state, target, float(guess_time)
                )
                at_center = self.field.flow(A, c, nominal_state, tcenter)[0]
                jump = arb_mat(
                    [
                        [
                            exponential[i, j]
                            - velocity[i, 0] * exponential[0, j] / velocity[0, 0]
                            for j in range(3)
                        ]
                        for i in range(3)
                    ]
                )
                sharp = at_center + jump * (state - nominal_state)
                for j in (1, 2):
                    end[j, 0] = end[j, 0].intersection(sharp[j, 0])
            end[0, 0] = target
            state = end
        return col([state[1, 0], state[2, 0]]), arb_mat(
            [[jac[i, j] for j in range(2)] for i in (1, 2)]
        )

    def map(self, point, schedule):
        return self.flights(col([1, *point.entries()]), schedule)

    def boundary_box(self, radius):
        offset = col([arb(0, radius)] * 2)
        region = self.center + self.B * offset
        value, _ = self.map(self.center, self.schedule)
        _, jac = self.map(region, self.schedule)
        jac = self.inverse * jac * self.B
        norm = max(sum(abs(jac[i, j]) for j in range(2)).abs_upper() for i in range(2))
        error = max(
            v.abs_upper() for v in (self.inverse * (value - self.center)).entries()
        )
        if not error + norm * radius < radius:
            raise ValueError(f"Trapping contraction failed, norm={float(norm)}")
        return norm, error

    def transport(self, point):
        center = midpoint(point)
        val, _ = self.map(self.center + self.B * center, self.schedule)
        _, jac = self.map(self.center + self.B * point, self.schedule)
        return self.inverse * (val - self.center) + (self.inverse * jac * self.B) * (
            point - center
        )


def certify(name, parameters, source):
    capture = Capture(name, parameters, source)
    if name == "right_lobe":
        x, y = ".55", "-.36"
        tail = [(0, 1, 0.49, 1)]
        initial_radius = "0.0000004"
    else:
        x, y = "-10.65", "-10.05"
        tail = [(-1, -1, 0.629, -1), *capture.schedule[-2:]]
        initial_radius = "0.001"
    radius = arb("0.0000001")
    norm, error = capture.boundary_box(radius)
    root = cycle(name, parameters)
    assert all(
        arb(a["lower"]).upper() < arb(b["lower"]).lower()
        and arb(b["upper"]).upper() < arb(a["upper"]).lower()
        for a, b in zip(source["shooting_box"], root["shooting_box"])
    )
    qroot = col([enclosure(q) for q in root["section_box"]])
    root_distance = max(
        v.abs_upper() for v in (capture.inverse * (qroot - capture.center)).entries()
    )
    assert root_distance < radius
    r = arb(initial_radius)
    start = col([arb(arb(x), r), arb(arb(y), r), 0])
    dv = arb_mat([[1, 0], [0, 1], [0, 0]])
    center_state = col([arb(x), arb(y), 0])
    dest, _ = capture.flights(center_state, tail, dv)
    _, J = capture.flights(start, tail, dv)
    point = capture.inverse * (dest - capture.center) + (capture.inverse * J) * col(
        [arb(0, r), arb(0, r)]
    )
    for i in range(350):
        displacement = max(v.abs_upper() for v in point.entries())
        if displacement < radius:
            break
        point = capture.transport(point)
    else:
        raise ValueError(f"{name} did not enter trapping box")
    print(
        name, "returns", i, "bound", displacement.str(8), "q", norm.str(5), flush=True
    )
    return dict(
        inside_uniform_shooting_box=True,
        x=x,
        y=y,
        z="0",
        voltage_radius=initial_radius,
        section_basis=[[bounds(capture.B[i, j]) for j in range(2)] for i in range(2)],
        section_center=[bounds(v) for v in capture.center.entries()],
        trapping_radius=bounds(radius),
        derivative_norm_upper=bounds(norm),
        center_defect=bounds(error),
        root_distance_upper=bounds(root_distance),
        full_returns=i,
        final_distance_upper=bounds(displacement),
        terminal_enclosure=[bounds(v) for v in point.entries()],
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--digits", type=int, default=70)
    args = p.parse_args()
    if not __debug__:
        raise RuntimeError("Assertions must be enabled")
    ctx.dps = args.digits
    source = json.loads((OUT / "certificate-70.json").read_text())
    par = [arb(v) for v in PARAMETERS]
    report = dict(
        nominal=PARAMETERS,
        digits=ctx.dps,
        starts={
            name: certify(name, par, rec) for name, rec in source["cycles"].items()
        },
    )
    (OUT / f"starts-{ctx.dps}.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
