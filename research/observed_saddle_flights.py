"""Affine flights and positive-velocity x=-1 returns for the lossy Chua circuit."""
import numpy as np
from scipy.optimize import brentq, root

from study_chua_circuit_constraints import LossyChua, LossyFlights


NOMINAL = dict(alpha=8.41, beta=12.23, gamma=0.0435, m0=-1.366, m1=-0.17)


class Section:
    def __init__(self, parameters=None, step=0.015):
        p = dict(NOMINAL if parameters is None else parameters)
        self.beta = p.pop("beta")
        self.model = LossyChua(**p)
        self.matrices = LossyFlights(self.model, [-1, 0, 1], [], False).matrices(self.beta)
        self.grid = np.arange(step, 8 + step, step)
        self.exponentials = {
            reg: np.exp(self.grid[:, None] * mat[2]) for reg, mat in self.matrices.items()
        }

    def flight(self, state, region):
        matrix, center, eig, V, W = self.matrices[region]
        coeff = V[0] * (W @ (state - center))
        x = center[0] + np.einsum("ij,j->i", self.exponentials[region], coeff).real
        outside = abs(x) > 1 if region == 0 else region * x < 1
        indices = np.flatnonzero(outside)
        if not len(indices):
            raise ValueError("No breakpoint detected within the eight-unit flight horizon")
        i = indices[0]
        wall = float(np.sign(x[i])) if region == 0 else float(region)
        left = self.grid[i-1] if i else 1e-10
        right = self.grid[i]
        def gap(t):
            return center[0] + np.dot(np.exp(eig*t), coeff).real - wall
        dt = brentq(gap, left, right, xtol=8e-16)
        flow = ((V * np.exp(eig*dt)) @ W).real
        end = center + flow @ (state-center)
        end[0] = wall
        return dt, end, flow

    def next_section(self, state, region, arcs=None):
        time = 0.
        derivative = np.eye(3)
        for _ in range(80):
            dt, end, flow = self.flight(state, region)
            if arcs is not None:
                arcs.append((region, state.copy(), dt))
            time += dt
            derivative = flow @ derivative
            if region == -1:
                speed = self.model.field(self.beta, end)
                projection = np.eye(3) - np.outer(speed, [1, 0, 0])/speed[0]
                return time, end, projection @ derivative
            state, region = end, int(end[0]) if region == 0 else 0
        raise ValueError("Positive x=-1 crossing not reached in 80 flights")

    def launch(self, x, arcs=None):
        return self.next_section(np.array([x, 0., 0.]), int(np.sign(x)) if abs(x)>1 else 0, arcs)

    def full_return(self, yz, arcs=None):
        time, end, jacobian = self.next_section(np.r_[-1., yz], 0, arcs)
        return time, end[1:], jacobian[1:, 1:]

    def absolute_x_peak(self, arcs):
        """Numerical within-flight extrema, checking both endpoints of every arc."""
        peak = 0.
        for region, state, duration in arcs:
            _, center, eig, V, W = self.matrices[region]
            coeff = V[0]*(W @ (state-center))
            def coordinate(t):
                return float(center[0] + (np.exp(eig*t) @ coeff).real)
            def speed(t):
                return float((np.exp(eig*t) @ (eig*coeff)).real)
            grid = np.linspace(0., duration, max(4, int(np.ceil(duration/.025))+1))
            rates = np.einsum("ij,j->i", np.exp(grid[:,None]*eig), eig*coeff).real
            indices = np.flatnonzero(rates[:-1]*rates[1:] < 0.)
            times = [0., duration] + [brentq(speed, grid[k], grid[k+1]) for k in indices]
            peak = max(peak, max(abs(coordinate(t)) for t in times))
        return peak

    def saddle(self):
        # The article's rounded point first travels from the left outer region.
        _, first, _ = self.next_section(np.array([-2.321, -0.787, 3.129]), -1)
        def shooting(yz):
            _, end, derivative = self.full_return(yz)
            return end-yz, derivative-np.eye(2)
        answer = root(shooting, first[1:], jac=True, tol=1e-12)
        period, end, derivative = self.full_return(answer.x)
        if np.linalg.norm(end-answer.x) > 1e-10:
            raise ArithmeticError(answer.message)
        values, vectors = np.linalg.eig(derivative)
        unstable = int(np.argmax(abs(values)))
        left = np.linalg.inv(vectors)[unstable].real
        left /= np.linalg.norm(left)
        if left[0] < 0:
            left *= -1
        return answer.x, period, sorted(values.real), left
