"""Ball-arithmetic affine return map for the classical Chua transect.

Arb encloses each matrix exponential. An interval-Newton inclusion encloses
each crossing time. Derivatives follow implicit differentiation of the event
equation, through order two in (initial y, initial z, beta).
Itinerary admissibility is a separate, explicit all-time check.
"""

from flint import arb, arb_mat, ctx
import numpy as np
from estimate_chua_homoclinic_tangency import HalfMap


def col(v):
    return arb_mat([[x] for x in v])


def eye(n):
    return arb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def midpoint(M):
    return arb_mat(
        [[M[i, j].mid() for j in range(M.ncols())] for i in range(M.nrows())]
    )


def upper(x):
    return float(x.abs_upper())


def bounds(x):
    return dict(
        ball=x.str(30),
        lower=x.lower().str(30),
        upper=x.upper().str(30),
        midpoint=x.mid().str(ctx.dps, radius=False),
        radius=x.rad().str(20),
    )


class IntervalHalfMap:
    regions = (0, 1, 0, 1, 0)
    targets = (1, 1, 1, 1, 0)
    section_x = 0
    return_sign = -1

    def __init__(self, digits=70):
        ctx.dps = digits
        self.float_map = HalfMap()

    def matrices(self, beta, region):
        slope = -arb(8) / 7 if region == 0 else -arb(5) / 7
        A = arb_mat([[-9 * (1 + slope), 9, 0], [1, -1, 1], [0, -beta, 0]])
        c = col([arb(3) * region / 2, 0, -arb(3) * region / 2])
        B = arb_mat([[0, 0, 0], [0, 0, 0], [0, -1, 0]])
        return A, c, B

    def schedule(self, point, beta):
        state = np.array([0.0, float(point[0].mid()), float(point[1].mid())])
        b = float(beta.mid())
        mats = self.float_map.matrices(b)
        times = []
        for region, target in zip(self.regions, self.targets):
            dt, actual = self.float_map.next_event(state, region, mats[region])
            if actual != target:
                raise ValueError("Unexpected midpoint itinerary")
            A, c, ev, V, W = mats[region]
            state = c + ((V * np.exp(ev * dt)) @ W).real @ (state - c)
            state[0] = target
            times.append(dt)
        return times

    @staticmethod
    def flow(A, c, state, t):
        E = (A * t).exp()
        return c + E * (state - c), E

    def crossing(self, A, c, state, target, guess):
        t0 = arb(str(guess))
        f0 = self.flow(A, c, state, t0)[0][0, 0] - target
        v0 = (A * (self.flow(A, c, state, t0)[0] - c))[0, 0]
        if v0.contains(0):
            raise ValueError("Crossing midpoint speed enclosure contains zero")
        radius = max(1e-13, 4 * upper(f0 / v0))
        certified = False
        for _ in range(8):
            T = arb(t0, arb(str(radius)))
            stateT, E = self.flow(A, c, state, T)
            velocity = A * (stateT - c)
            if not velocity[0, 0].contains(0):
                N = t0 - f0 / velocity[0, 0]
                if T.contains_interior(N):
                    certified = True
                    T = N
                    break
            radius *= 2
        if not certified:
            raise ValueError("Event interval-Newton inclusion failed")
        for _ in range(5):
            mid = T.mid()
            end, E = self.flow(A, c, state, T)
            v = A * (end - c)
            if v[0, 0].contains(0):
                break
            N = mid - (self.flow(A, c, state, mid)[0][0, 0] - target) / v[0, 0]
            new = T.intersection(N)
            if upper(new.rad()) >= 0.99 * upper(T.rad()):
                break
            T = new
        if not T > 0:
            raise ValueError("Nonpositive flight time enclosure")
        return T

    def exponential_derivatives(self, A, B, t, order):
        if order < 1:
            return (A * t).exp(), None, None
        count = 3 if order > 1 else 2
        block = arb_mat(3 * count, 3 * count)
        for k in range(count):
            for i in range(3):
                for j in range(3):
                    block[3 * k + i, 3 * k + j] = A[i, j]
                    if k < count - 1:
                        block[3 * k + i, 3 * (k + 1) + j] = B[i, j]
        E = (block * t).exp()
        part = lambda k: arb_mat(
            [[E[i, 3 * k + j] for j in range(3)] for i in range(3)]
        )
        return part(0), part(1), 2 * part(2) if order > 1 else None

    def evaluate(self, point, beta, order=2, verify_path=False, centered=True):
        center_arcs = None
        delta = [x - x.mid() for x in [*point, beta]]
        if centered and order and any(x.rad() > arb("1e-45") for x in [*point, beta]):
            center_arcs = self.evaluate(
                [x.mid() for x in point], beta.mid(), order, False, False
            )[3]
        state = col([self.section_x, *point])
        D = [col([0, 1, 0]), col([0, 0, 1]), col([0, 0, 0])]
        DD = [[col([0, 0, 0]) for j in range(3)] for i in range(3)]
        times = self.schedule(point, beta)
        arcs = []
        for arc_index, (reg, target, guess) in enumerate(
            zip(self.regions, self.targets, times)
        ):
            A, c, B = self.matrices(beta, reg)
            start = state
            T = self.crossing(A, c, start, target, guess)
            E, Eb, Ebb = self.exponential_derivatives(A, B, T, order)
            z = start - c
            end = c + E * z
            v = A * (end - c)
            if v[0, 0].contains(0):
                raise ValueError("Nontransverse return enclosure")
            if order:
                F = [
                    E * D[i] + (Eb * z if i == 2 else col([0, 0, 0])) for i in range(3)
                ]
                dt = [-f[0, 0] / v[0, 0] for f in F]
                newD = [F[i] + v * dt[i] for i in range(3)]
                if order > 1:
                    Ft = [
                        A * F[i] + (B * (end - c) if i == 2 else col([0, 0, 0]))
                        for i in range(3)
                    ]
                    Ftt = A * v
                    newDD = []
                    for i in range(3):
                        row = []
                        for j in range(3):
                            mixed = E * DD[i][j]
                            if i == 2:
                                mixed += Eb * D[j]
                            if j == 2:
                                mixed += Eb * D[i]
                            if i == 2 and j == 2:
                                mixed += Ebb * z
                            mixed += Ft[i] * dt[j] + Ft[j] * dt[i] + Ftt * dt[i] * dt[j]
                            ddt = -mixed[0, 0] / v[0, 0]
                            out = mixed + v * ddt
                            out[0, 0] = 0
                            row.append(out)
                        newDD.append(row)
                    DD = newDD
                D = newD
                for d in D:
                    d[0, 0] = 0
            end[0, 0] = target
            state = end
            if center_arcs is not None:
                ca = center_arcs[arc_index]
                enclosure = ca["end"]
                if order > 1:
                    tightened = []
                    for i in range(3):
                        derivative = ca["D"][i]
                        for j in range(3):
                            derivative += DD[i][j] * delta[j]
                        tightened.append(
                            col(
                                [
                                    D[i][k, 0].intersection(derivative[k, 0])
                                    for k in range(3)
                                ]
                            )
                        )
                    D = tightened
                    for i in range(3):
                        enclosure += ca["D"][i] * delta[i]
                    for i in range(3):
                        for j in range(3):
                            enclosure += DD[i][j] * delta[i] * delta[j] / 2
                else:
                    for i in range(3):
                        enclosure += D[i] * delta[i]
                state = col(
                    [state[i, 0].intersection(enclosure[i, 0]) for i in range(3)]
                )
                state[0, 0] = target
            arc = dict(
                region=reg,
                target=target,
                duration=T,
                start=start,
                end=state,
                velocity=v,
                D=D,
                DD=DD,
            )
            if verify_path:
                arc["path_certificate"] = self.verify_arc(A, c, arc)
            arcs.append(arc)
        value = self.return_sign*col([state[1, 0], state[2, 0]])
        jac = (
            arb_mat([[self.return_sign*D[j][i, 0] for j in range(3)] for i in [1, 2]])
            if order
            else None
        )
        hess = (
            [
                arb_mat([[self.return_sign*DD[j][k][i, 0] for k in range(3)] for j in range(3)])
                for i in [1, 2]
            ]
            if order > 1
            else None
        )
        return value, jac, hess, arcs

    def verify_arc(self, A, c, arc):
        """Cover every normalized time in [0,1], including endpoint collars."""
        T = arc["duration"]
        start = arc["start"]
        reg = arc["region"]
        target = arc["target"]
        work = [(arb(0), arb(1), 0)]
        cells = 0
        min_clearance = None
        while work:
            left, right, depth = work.pop()
            tau = left.union(right) * T
            state, E = self.flow(A, c, start, tau)
            x = state[0, 0]
            vx = (A * (state - c))[0, 0]
            if reg == 1:
                ok = x > 1 or (left.is_zero() and vx > 0) or (right == 1 and vx < 0)
            else:
                ok = x > 0 and x < 1
                if left.is_zero():
                    ok = (
                        ok
                        or (start[0, 0].is_zero() and vx > 0 and x < 1)
                        or (start[0, 0] == 1 and vx < 0 and x > 0)
                    )
                if right == 1:
                    ok = (
                        ok
                        or (target == 0 and vx < 0 and x < 1)
                        or (target == 1 and vx > 0 and x > 0)
                    )
            if ok:
                cells += 1
                continue
            if depth >= 20:
                raise ValueError(
                    f"Itinerary interval unresolved: region {reg}, time {tau}, x {x}"
                )
            mid = (left + right) / 2
            work.extend([(left, mid, depth + 1), (mid, right, depth + 1)])
        return dict(covered_normalized_time="[0,1]", accepted_cells=cells)
