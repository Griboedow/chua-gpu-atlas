"""Exact three-state normalization of Li, Li and Zhou (2025), Polymers 17:617.

No small-thermal-time approximation. The paper assumes heating Q=Q0-beta*u;
it is not replaced with a newly invented resistance law. Parameters c,r,g.
"""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.stats import qmc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from chua_gpu import GPU


def field_jacobian(state, p):
    x, v, y = state
    c, r, g = p
    return (np.array([v, -x-c*v-abs(v)*v+y, -r*(y+g*x)]),
            np.array([[0., 1., 0.], [-1., -c-2*abs(v), 1.], [-r*g, 0., -r]]))


def normalize(raw, state):
    """Raw dimensionless parameters j,alpha,beta,Q0,C1,C2,tau; state u,u_dot,T."""
    j, a, beta, Q0, C1, C2, tau = raw
    K = 1+j
    u0, T0 = a*Q0/(K+a*beta), K*Q0/(K+a*beta)
    p = np.array([C1/np.sqrt(K), 1/(tau*np.sqrt(K)), a*beta/K])
    scale = np.array([C2, C2/np.sqrt(K), a*C2/K])
    return p, (state-np.array([u0, 0., T0]))*scale, scale, np.sqrt(K)


def critical_gain(c, r):
    return c/r+c*c+c*r


class RheostatGPU(GPU):
    def __init__(self):
        super().__init__('', 'float64')
        self.kernel = cl.Kernel(cl.Program(self.context, Path(__file__).with_suffix('.cl').read_text()).build(), 'advance')
        self.info.update(model='Li, Li and Zhou 2025 full rheostat ODE, exact normalization',
                         integrator='Dormand-Prince 5(4), analytic tangent, per-accepted-step renormalization')

    def scan(self, p, starts, transient=500., duration=1000., rtol=1e-8, step_factor=1.,
             max_attempts=3000000, max_wall_seconds=600, progress=None):
        p, starts = np.asarray(p, float), np.asarray(starts, float)
        n = len(p)
        if p.shape != (n, 3) or starts.shape != (n, 3) or n == 0:
            raise ValueError('Expected nonempty (N,3) arrays')
        if not np.isfinite(p).all() or not np.isfinite(starts).all() or (p[:, 0]<0).any() or (p[:, 1:]<=0).any():
            raise ValueError('Require c>=0, r,g>0 and finite inputs')
        if transient < 0 or min(duration, rtol, step_factor, max_attempts, max_wall_seconds) <= 0:
            raise ValueError('Invalid integration settings')
        params = np.zeros((n, 4)); params[:, :3] = p
        state = np.zeros((n, 8)); state[:, :3] = starts; state[:, 3:6] = np.array([1,2,3])/np.sqrt(14)
        clocks = np.zeros((n, 4)); clocks[:, 1] = .001
        logs = np.zeros((n, 4)); stats = np.zeros((n, 8))
        stats[:, [0, 2, 6]] = np.inf; stats[:, [1, 3, 7]] = -np.inf
        status = np.zeros(n, np.int32); counts = np.zeros((n, 2), np.int32)
        arrays = [state, params, clocks, logs, stats, status, counts]
        buffers = [self.buffer(a) for a in arrays]
        wall = time.perf_counter(); gpu_seconds = 0.; calls = 0; last = wall
        while True:
            event = self.kernel(self.queue, (n,), None, *buffers, np.int32(n), np.int32(128),
                                np.float64(transient), np.float64(duration), np.float64(rtol), np.float64(rtol/1000),
                                np.float64(step_factor), np.int32(max_attempts))
            event.wait(); gpu_seconds += (event.profile.end-event.profile.start)*1e-9; calls += 1
            if calls % 8 == 0:
                cl.enqueue_copy(self.queue, clocks, buffers[2]).wait()
                cl.enqueue_copy(self.queue, status, buffers[5]).wait()
                active = (status == 0) & (clocks[:, 0] < transient+duration)
                if progress and time.perf_counter()-last > 10:
                    progress(np.minimum(clocks[:, 0], transient+duration).mean()/(transient+duration), time.perf_counter()-wall)
                    last = time.perf_counter()
                if not active.any() or time.perf_counter()-wall > max_wall_seconds:
                    break
        for a, b in zip(arrays, buffers): cl.enqueue_copy(self.queue, a, b).wait()
        status[(status == 0) & (clocks[:, 0] < transient+duration)] = 4
        windows = logs/(duration/4); windows[status != 0] = np.nan
        return dict(parameters=p, starts=starts, final_states=state[:, :3], lle=windows.mean(axis=1),
                    lle_windows=windows, status=status, stats=stats, step_counts=counts, elapsed=clocks[:, 0],
                    backend=self.info, settings=dict(transient=transient, duration=duration, rtol=rtol,
                    atol=rtol/1000, step_factor=step_factor, max_attempts=max_attempts, max_wall_seconds=max_wall_seconds),
                    timing=dict(wall_seconds=time.perf_counter()-wall, gpu_seconds=gpu_seconds),
                    parameter_columns=['c','r','g'], state_columns=['X','V','Y'],
                    stats_columns=['min_X','max_X','min_Y','max_Y','max_abs_V','max_abs_state','late_min_X','late_max_X'],
                    status_legend={'0':'completed','1':'amplitude cutoff 1e5','2':'numerical failure/step underflow',
                                   '3':'attempt budget exhausted; unresolved','4':'wall budget exhausted; unresolved'})


def validate(gpu, out):
    rng = np.random.default_rng(2026102954)
    algebra, derivative = [], []
    for _ in range(64):
        raw = 10**rng.uniform(-1, 1, 7); state = rng.uniform(-1, 1, 3)
        j,a,beta,Q0,C1,C2,tau = raw; u,v,T = state
        original = np.array([v, -(1+j)*u+a*T-C1*v-C2*abs(v)*v, (Q0-beta*u-T)/tau])
        p, scaled, scale, rate = normalize(raw, state)
        f, jac = field_jacobian(scaled, p)
        algebra.append(np.max(abs(f-original*scale/rate)))
        h = 1e-6
        fd = np.column_stack([(field_jacobian(scaled+h*e,p)[0]-field_jacobian(scaled-h*e,p)[0])/(2*h) for e in np.eye(3)])
        derivative.append(np.max(abs(fd-jac)))
    p = np.array([[.1,1.,1.], [.02,20.,2.], [1.,.1,30.]])
    starts = np.array([[.1,0,-.1], [.02,.03,.01], [1.,-.5,.1]])
    result = gpu.scan(p, starts, transient=0, duration=3, rtol=1e-10)
    errors, stretches = [], []
    for pp, s, end, lle in zip(p, starts, result['final_states'], result['lle']):
        def rhs(t, z):
            f,j = field_jacobian(z[:3], pp)
            return np.r_[f, j@z[3:]]
        sol = solve_ivp(rhs, (0,3), np.r_[s,np.array([1,2,3])/np.sqrt(14)], method='DOP853', rtol=1e-12, atol=1e-14)
        assert sol.success
        errors.append(float(max(abs(sol.y[:3,-1]-end))))
        stretches.append(float(abs(np.log(np.linalg.norm(sol.y[3:,-1]))/3-lle)))
    check = dict(normalization_max_error=float(max(algebra)), jacobian_fd_max_error=float(max(derivative)),
                 gpu_dop853_state_errors=errors, gpu_dop853_lle_errors=stretches,
                 passed=bool(max(algebra)<1e-9 and max(derivative)<1e-5 and max(errors)<1e-7 and max(stretches)<1e-7 and (result['status']==0).all()))
    out.mkdir(parents=True, exist_ok=True); (out/'validation.json').write_text(json.dumps(check, indent=2))
    print(json.dumps(check), flush=True)
    if not check['passed']: raise RuntimeError('Validation failed')


def save(result, out, extra):
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out/'scan.npz', **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    candidate = (result['status']==0)&(result['lle']>.005)&(result['lle_windows'].min(axis=1)>.002)
    meta = {k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    meta.update(extra, trajectories=len(candidate), candidate_indices=np.flatnonzero(candidate).tolist(),
                status_counts={str(i):int((result['status']==i).sum()) for i in np.unique(result['status'])},
                criterion='Mean exponent > .005 and all four windows > .002; finite-time scout only',
                physical_caveat='No full-device calibration; Q>=0, positive temperature/tension and small strain must be checked after choosing original scaling.')
    (out/'summary.json').write_text(json.dumps(meta, indent=2))
    print({k:meta[k] for k in ['trajectories','candidate_indices','status_counts','timing']}, flush=True)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--validate',action='store_true')
    parser.add_argument('--power',type=int,default=11); parser.add_argument('--transient',type=float,default=500)
    parser.add_argument('--duration',type=float,default=1000); parser.add_argument('--rtol',type=float,default=1e-8)
    parser.add_argument('--step-factor',type=float,default=1); parser.add_argument('--refine',type=Path)
    parser.add_argument('--indices',type=int,nargs='+')
    parser.add_argument('--max-attempts',type=int,default=3000000)
    parser.add_argument('--top',type=int,default=24); parser.add_argument('--output',type=Path,default=Path('results/rheostat/scout'))
    args=parser.parse_args(); gpu=RheostatGPU()
    if args.validate: return validate(gpu,Path('results/rheostat'))
    if args.refine:
        old=np.load(args.refine/'scan.npz'); valid=np.flatnonzero(old['status']==0)
        ids=np.array(args.indices,int) if args.indices is not None else valid[np.argsort(-old['lle'][valid])[:args.top]]
        if len(ids)==0 or len(np.unique(ids))!=len(ids) or not np.isin(ids,valid).all():
            raise ValueError('Follow-up indices must be unique completed source trajectories')
        p,starts=old['parameters'][ids],old['starts'][ids]
        extra=dict(parent_source=str(args.refine),parent_indices=ids.tolist())
    else:
        rows=qmc.Sobol(3,scramble=True,seed=2026102954).random_base2(args.power)
        c=10**(-2.3+2.6*rows[:,0]);r=10**(-1.5+3*rows[:,1])
        gain=critical_gain(c,r)*10**(-.3+1.8*rows[:,2])
        p=np.repeat(np.column_stack([c,r,gain]),2,axis=0)
        starts=np.tile([[.01,0.,0.],[1.,.2,-.1]],(len(rows),1))
        extra=dict(seed=2026102954,sampling='Sobol log10c[-2.3,.3],log10r[-1.5,1.5],log10(g/g_H)[-.3,1.5]; two common starts',
                   exact_normalization='s=sqrt(K)t,X=C2(u-u0),V=C2*u_dot/sqrt(K),Y=alpha*C2*(T-T0)/K, K=1+j',
                   source_doi='10.3390/polym17050617')
    print('Scanning',len(p),'rheostat trajectories',flush=True)
    result=gpu.scan(p,starts,args.transient,args.duration,args.rtol,args.step_factor,
                    max_attempts=args.max_attempts,
                    progress=lambda f,t:print(f'{f:.1%}, {t:.1f}s',flush=True))
    save(result,args.output,extra)


if __name__=='__main__':main()
