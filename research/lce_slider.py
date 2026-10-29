"""Event-resolved model of Ge et al. (2024), DOI 10.3390/polym16233298.

Equations 9, 11--13 reduce to u'=v, v'=k*(d-u-phi)-b*v,
phi'=I*H(u)*(1-phi)-phi, where u=(x-h)/C and d=(g/k-h)/C.
Only underdamped mechanical parameters are implemented here. Between switches
the affine flow is evaluated analytically; no light switching is time-sampled.
"""
import json
from pathlib import Path

import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def flow(state, parameters, on, normalized=False):
    k, b, intensity, d = parameters
    u, v, phi = state
    photo=(1+intensity)*on-(1+intensity*on)*phi if normalized else intensity*on*(1-phi)-phi
    return np.array([v, k * (d-u-phi)-b*v, photo])


def affine_flow(state, parameters, on, time, homogeneous=False, normalized=False):
    k, b, intensity, d = parameters
    r = 1 + intensity * on
    q = float(on) if normalized else intensity * on / r
    if homogeneous:
        q, d = 0., 0.
    omega = np.sqrt(k - b*b/4)
    u, v, phi = state
    amplitude = -k * (phi-q) / (r*r-b*r+k)
    equilibrium = d-q
    c = u-equilibrium-amplitude
    s = (v+r*amplitude+b*c/2)/omega
    er, eb = np.exp(-r*time), np.exp(-b*time/2)
    co, si = np.cos(omega*time), np.sin(omega*time)
    return np.array([
        equilibrium+amplitude*er+eb*(c*co+s*si),
        -r*amplitude*er+eb*((-b*c/2+omega*s)*co+(-b*s/2-omega*c)*si),
        q+(phi-q)*er])


def next_crossing(section, parameters, on, subdivisions=64, normalized=False):
    state = np.array([0., *section])
    k, b, intensity, d = parameters
    direction = 1 if on else -1
    if direction * state[1] <= 1e-12:
        raise ValueError('Nontransverse or incorrectly directed initial crossing')
    step = 2*np.pi/(np.sqrt(k-b*b/4)*subdivisions)
    lo = 1e-10 / max(1., np.sqrt(k), 1+intensity)
    for index in range(1, 8193):
        hi = index*step
        if direction * affine_flow(state, parameters, on, hi, normalized=normalized)[0] < 0:
            root = brentq(lambda t: affine_flow(state, parameters, on, t, normalized=normalized)[0], lo, hi,
                          xtol=5e-14, rtol=1e-14)
            end = affine_flow(state, parameters, on, root, normalized=normalized)
            if abs(end[1]) < 1e-10:
                raise ValueError('Near grazing crossing')
            matrix = np.column_stack([affine_flow(e, parameters, on, root, True, normalized) for e in np.eye(3)])
            f = flow(end, parameters, on, normalized)
            projection = np.eye(3) - np.outer(f, [1., 0., 0.])/f[0]
            jac = (projection @ matrix)[1:, 1:]
            return end[1:], root, jac
        lo = hi
    raise RuntimeError('No return within the event-bracketing horizon')


def return_map(section, parameters, subdivisions=64, normalized=False):
    middle, t1, j1 = next_crossing(section, parameters, True, subdivisions, normalized)
    end, t2, j2 = next_crossing(middle, parameters, False, subdivisions, normalized)
    return end, t1+t2, j2@j1


class SliderGPU:
    def __init__(self, normalized=False):
        devices = [d for p in cl.get_platforms() for d in p.get_devices()
                   if d.type & cl.device_type.GPU and d.double_fp_config]
        if not devices:
            raise RuntimeError('FP64 GPU required')
        self.device = devices[0]
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context)
        self.normalized=normalized
        options=['-DNORMALIZED_PHOTO'] if normalized else []
        self.kernel = cl.Kernel(cl.Program(self.context, Path(__file__).with_suffix('.cl').read_text()).build(options=options), 'advance')

    def scan(self, parameters, starts, transient=200, cycles=800, subdivisions=64, progress=None):
        import time
        p = np.ascontiguousarray(parameters, dtype=float)
        starts = np.ascontiguousarray(starts, dtype=float)
        lanes = len(p)
        if p.shape != (lanes, 4) or starts.shape != (lanes, 2):
            raise ValueError('Input shape mismatch')
        if cycles % 4 or cycles < 4 or transient < 0 or subdivisions < 16:
            raise ValueError('Invalid observation or bracketing settings')
        k,b,I,d = p.T
        valid_photo=np.all(I>=0) and np.all(d<1) if self.normalized else np.all(I>0) and np.all(d<I/(1+I))
        if not (np.isfinite(p).all() and np.isfinite(starts).all() and np.all(k>b*b/4) and np.all(b>0) and valid_photo
                and np.all(d>0) and np.all(starts[:,0]>1e-10)):
            raise ValueError('Require underdamping, positive parameters, interior threshold and transverse start')
        state = np.zeros((lanes, 6)); state[:,:2] = starts
        state[:,2:4] = np.array([1., .3])/np.sqrt(1.09)
        logs = np.zeros((lanes, 4)); durations = np.zeros((lanes, 4)); status = np.zeros(lanes, np.int32)
        arrays = [state, p, logs, durations, status]
        buffers = [cl.Buffer(self.context, cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR, hostbuf=a) for a in arrays]
        wall = time.perf_counter(); last=wall; done=0
        for window, count in [(-1, transient)]+[(i, cycles//4) for i in range(4)]:
            for offset in range(0, count, 4):
                length=min(4, count-offset)
                self.kernel(self.queue, (lanes,), None, *buffers, np.int32(lanes), np.int32(length),
                            np.int32(window), np.int32(subdivisions)).wait()
                done+=length
                if progress and time.perf_counter()-last>10:
                    progress(done/(cycles+transient), time.perf_counter()-wall);last=time.perf_counter()
        for a, buffer in zip(arrays,buffers): cl.enqueue_copy(self.queue,a,buffer).wait()
        with np.errstate(divide='ignore', invalid='ignore'):
            windows=logs/durations; lle=logs.sum(axis=1)/durations.sum(axis=1)
        windows[status!=0]=np.nan;lle[status!=0]=np.nan
        return dict(parameters=p,starts=starts,final_sections=state[:,:2],lle=lle,lle_windows=windows,
                    status=status,durations=durations,log_stretch=logs,
                    minimum_crossing_speed=state[:,4],maximum_abs_u=state[:,5],
                    settings=dict(transient_cycles=transient,observation_cycles=cycles,subdivisions=subdivisions),
                    method='Analytic affine segments; bisection at each crossing; derivative of event-defined return map',
                    normalized_photo=self.normalized,
                    state_coordinates='U=u/q,V=v/q,Y=phi/q,delta=d/q; continuous rescaled I=0 limit' if self.normalized else 'u,v,phi,d',
                    device=self.device.name,wall_seconds=time.perf_counter()-wall)


def validate():
    params=np.array([5.8,.1,.1,(1.2/5.8-.2)/.2])
    s=np.array([.3,.03]);state=np.r_[0.,s]
    ref=solve_ivp(lambda t,y:flow(y,params,True),(0,.4),state,method='DOP853',rtol=1e-12,atol=1e-14)
    exact=affine_flow(state,params,True,.4)
    end, period, jac=return_map(s,params)
    eps=1e-7
    numerical=np.column_stack([(return_map(s+eps*e,params)[0]-return_map(s-eps*e,params)[0])/(2*eps) for e in np.eye(2)])
    initial=s.copy()
    for _ in range(4):s=return_map(s,params)[0]
    gpu=SliderGPU().scan(params[None],initial[None],transient=0,cycles=4)
    result=dict(affine_vs_DOP853=float(np.max(abs(exact-ref.y[:,-1]))),
                return_jacobian_error=float(np.max(abs(jac-numerical))),
                gpu_cpu_four_return_error=float(np.max(abs(s-gpu['final_sections'][0]))),
                gpu_status=int(gpu['status'][0]),first_period=float(period),baseline_parameters=params.tolist())
    result['passed']=(result['affine_vs_DOP853']<1e-10 and result['return_jacobian_error']<1e-6
                      and result['gpu_cpu_four_return_error']<1e-9 and result['gpu_status']==0)
    path=Path('results/lce-slider');path.mkdir(parents=True,exist_ok=True)
    (path/'validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2),flush=True)
    if not result['passed']:raise RuntimeError('Validation failed')


if __name__=='__main__':validate()
