"""Continue the atlas's classical periodic window using full-field shooting.

The ODE always evaluates the actual continuous piecewise-linear field. Events
only measure crossings. The variational flow is a product of matrix exponentials
over measured affine residence intervals; no saltation jump is needed because
the vector field itself is continuous. Results are numerical, not interval proofs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import root


class Chua:
    def __init__(self, alpha=9., m0=-8/7, m1=-5/7, rtol=2e-11, step=.04, method="DOP853"):
        self.alpha, self.m0, self.m1 = alpha, m0, m1
        self.rtol, self.step, self.method = rtol, step, method

    def field(self, beta, u):
        x, y, z = u
        g = (1+self.m1)*x+(self.m0-self.m1)*max(-1., min(1., x))
        return np.array([self.alpha*(y-g), x-y+z, -beta*y])

    def jacobian(self, beta, inner):
        m = self.m0 if inner else self.m1
        return np.array([[-self.alpha*(1+m), self.alpha, 0.], [1., -1., 1.], [0., -beta, 0.]])

    def integrate(self, beta, u, duration, dense=True, method=None):
        def lower(t, q): return q[0]+1
        def upper(t, q): return q[0]-1
        def section(t, q): return q[0]
        section.direction = 1
        def maximum(t, q): return self.field(beta, q)[0]
        maximum.direction = -1
        sol = solve_ivp(lambda t, q: self.field(beta, q), [0., duration], u,
                        rtol=self.rtol, atol=self.rtol/10, max_step=self.step,
                        method=method or self.method, dense_output=dense,
                        events=[lower, upper, section, maximum])
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    def flow(self, beta, v):
        y, z, period = v
        if not 1 < period < 200:
            raise ValueError(f"Invalid period {period}")
        sol = self.integrate(beta, [0., y, z], period)
        times = sorted([0., period]+[float(t) for e in sol.t_events[:2] for t in e if 1e-10<t<period-1e-10])
        mono = np.eye(3)
        inner_time = 0.
        for left, right in zip(times[:-1], times[1:]):
            inner = abs(sol.sol((left+right)/2)[0]) < 1.
            if inner: inner_time += right-left
            mono = expm(self.jacobian(beta, inner)*(right-left)) @ mono
        end = sol.y[:, -1]
        residual = end-np.array([0., y, z])
        jac = np.column_stack(((mono-np.eye(3))[:, 1:], self.field(beta, end)))
        return residual, jac, sol, mono, inner_time

    def shoot(self, beta, guess):
        cache = {}
        def evaluate(v):
            k = tuple(v)
            if k not in cache:
                cache[k] = self.flow(beta, v)
            return cache[k]
        ans = root(lambda v: evaluate(v)[0], guess, jac=lambda v: evaluate(v)[1], tol=2e-9)
        res = np.max(np.abs(evaluate(ans.x)[0]))
        if res > 2e-8:
            raise RuntimeError(f"shoot failed beta={beta}, residual={res}: {ans.message}")
        return ans.x

    def record(self, beta, v):
        residual, _, sol, mono, inner_time = self.flow(beta, v)
        end = sol.y[:, -1]
        f = self.field(beta, end)
        section_matrix = ((np.eye(3)-np.outer(f, [1., 0., 0.])/f[0])@mono)[1:, 1:]
        multipliers = np.linalg.eigvals(section_matrix)
        multipliers = sorted(multipliers, key=abs, reverse=True)
        mesh = np.linspace(0., v[2], 2001)
        trace = sol.sol(mesh)
        symmetry = np.max(np.abs(trace+sol.sol((mesh+v[2]/2)%v[2])))
        crossings = np.concatenate(sol.t_events[:2])
        speed = min(abs(self.field(beta, sol.sol(t))[0]) for t in crossings)
        return dict(beta=float(beta), alpha=self.alpha, m0=self.m0, m1=self.m1,
                    shooting_vector=np.asarray(v).tolist(), period=float(v[2]),
                    residual=float(np.max(np.abs(residual))),
                    multipliers=[[float(q.real),float(q.imag)] for q in multipliers],
                    poincare_matrix=section_matrix.tolist(),
                    monodromy=mono.tolist(),
                    neutral_defect=float(np.linalg.norm(mono@self.field(beta, [0.,v[0],v[1]])-f)),
                    determinant_liouville=float(np.exp(-self.alpha*((1+self.m0)*inner_time+(1+self.m1)*(v[2]-inner_time))-v[2])),
                    determinant_poincare=float(np.linalg.det(section_matrix)),
                    inner_fraction=float(inner_time/v[2]),
                    minimum_breakpoint_speed=float(speed),
                    breakpoint_crossings=len(crossings),
                    positive_section_crossings=int(np.count_nonzero((sol.t_events[2]>1e-7)&(sol.t_events[2]<v[2]-1e-7))),
                    maxima=sol.y_events[3][:,0].tolist(),
                    half_period_symmetry_error=float(symmetry),
                    x_range=[float(trace[0].min()),float(trace[0].max())])


def initialize(model, beta):
    warm = model.integrate(beta, [.7, 0., 0.], 700., dense=False)
    late = model.integrate(beta, warm.y[:,-1], 100.)
    times, states = late.t_events[2], late.y_events[2]
    for lag in range(1, min(20,len(times)//2)):
        if np.max(np.abs(states[-1]-states[-1-lag])) < 1e-6:
            guess = [*states[-1,1:], times[-1]-times[-1-lag]]
            return model.shoot(beta, guess)
    raise RuntimeError("No recurrent section found")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1]/'results/chua-layers')
    p.add_argument('--initialize-only', action='store_true')
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    model=Chua()
    initial=initialize(model,14.32)
    record=model.record(14.32,initial)
    (args.output/'initial.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    if args.initialize_only: return
    records=[record]
    failures=[]
    for direction in [-1,1]:
        v=initial.copy()
        for beta in 14.32+direction*np.arange(1,51)*.0005:
            try:
                v=model.shoot(float(beta),v)
                r=model.record(float(beta),v)
                records.append(r)
                print(beta,v[2],r['multipliers'],flush=True)
            except (RuntimeError, ValueError) as ex:
                failures.append(dict(beta=float(beta),message=str(ex)))
                print('STOP',failures[-1],flush=True)
                break
    report=dict(scope='One classical five-maximum periodic window',
                settings=dict(rtol=model.rtol,atol=model.rtol/10,max_step=model.step,method=model.method),
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                elapsed_seconds=time.monotonic()-start,records=records,failures=failures)
    (args.output/'continuation.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__': main()
