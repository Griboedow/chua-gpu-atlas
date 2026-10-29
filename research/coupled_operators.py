"""Continuous-heat, six-state optical operators of Guo et al. (2025).

Source: 10.1038/s41467-025-66395-3, SI Eqs. 10, 17--20. Each
displacement is normalized by its own illuminated thermal deflection;
time is t/tau_heat. Temperatures are never reset at a relay transition.
Only underdamped mechanics is implemented by the analytic segment solver.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def flow(y, p, modes):
    k, b, dm, dp = p
    xm, vm, ym, xp, vp, yp = y
    return np.array([vm, k*(ym-xm)-b*vm, modes[0]-ym,
                     vp, k*(yp-xp)-b*vp, modes[1]-yp])


def initial_modes(y, p):
    return np.array([int(y[3] > p[3] or (y[3] == p[3] and y[4] > 0)),
                     int(y[0] < p[2] or (y[0] == p[2] and y[1] < 0))])


def segment3(y, k, b, target, t):
    x, v, h = y
    omega = np.sqrt(k-b*b/4)
    a = k*(h-target)/(1-b+k)
    c = x-target-a
    s = (v+a+b*c/2)/omega
    eh, eb = np.exp(-t), np.exp(-b*t/2)
    co, si = np.cos(omega*t), np.sin(omega*t)
    return np.array([target+a*eh+eb*(c*co+s*si),
                     -a*eh+eb*((-b*c/2+omega*s)*co+(-b*s/2-omega*c)*si),
                     target+(h-target)*eh])


def segment(y, p, modes, t):
    return np.r_[segment3(y[:3], *p[:2], modes[0], t),
                 segment3(y[3:], *p[:2], modes[1], t)]


def segment_matrix(p, t):
    a = np.column_stack([segment3(e, *p[:2], 0., t) for e in np.eye(3)])
    m = np.zeros((6, 6)); m[:3, :3] = a; m[3:, 3:] = a
    return m


def next_event(y, p, modes, subdivisions=64):
    """First sign-changing root on a bounded bracketing mesh.

    Mesh refinement is necessary: this routine does not rigorously rule out
    two crossings within a single bracket. Grazing/simultaneous events fail.
    """
    dt = min(.1, 2*np.pi/(np.sqrt(p[0]-p[1]**2/4)*subdivisions))
    sides = np.array([1-2*modes[1], 2*modes[0]-1])
    indices = [0, 3]; thresholds = p[2:]
    lo = 0.
    for index in range(1, 16385):
        hi = index*dt
        end = segment(y, p, modes, hi)
        crossed = sides*(end[indices]-thresholds) < 0
        if crossed.any():
            roots = [np.inf, np.inf]
            for j in np.flatnonzero(crossed):
                lower = lo if lo > 0 else 1e-13
                roots[j] = brentq(lambda t: segment(y, p, modes, t)[indices[j]]-thresholds[j],
                                 lower, hi, xtol=5e-14, rtol=1e-14)
            which = int(np.argmin(roots)); elapsed = roots[which]
            if abs(roots[0]-roots[1]) < 1e-10:
                raise RuntimeError('Simultaneous events need separate analysis')
            end = segment(y, p, modes, elapsed)
            normal_index = indices[which]; speed = end[normal_index+1]
            if abs(speed) < 1e-9 or elapsed < 1e-11:
                raise RuntimeError('Grazing or unresolved short event interval')
            end[normal_index] = thresholds[which]
            changed = modes.copy(); changed[1-which] = 1-changed[1-which]
            jump = flow(end, p, changed)-flow(end, p, modes)
            saltation = np.eye(6); saltation[:, normal_index] += jump/speed
            return end, changed, elapsed, saltation, which
        lo = hi
    raise RuntimeError('No crossing in finite bracketing horizon')


def fixed_time(y, p, duration, subdivisions=128):
    """Analytic segments and saltation, returning a fixed-time flow matrix."""
    state = y.copy(); modes = initial_modes(state, p); matrix = np.eye(6)
    remaining = duration; events = 0
    while remaining > 1e-13:
        end, changed, elapsed, saltation, which = next_event(state, p, modes, subdivisions)
        if elapsed >= remaining:
            matrix = segment_matrix(p, remaining) @ matrix
            state = segment(state, p, modes, remaining); break
        matrix = saltation @ segment_matrix(p, elapsed) @ matrix
        state, modes = end, changed; remaining -= elapsed; events += 1
    return state, matrix, events


class OperatorsGPU:
    def __init__(self):
        devices = [d for platform in cl.get_platforms() for d in platform.get_devices()
                   if d.type & cl.device_type.GPU and d.double_fp_config]
        if not devices: raise RuntimeError('An FP64 GPU is required')
        self.device = devices[0]; self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context)
        self.kernel = cl.Kernel(cl.Program(self.context, Path(__file__).with_suffix('.cl').read_text()).build(), 'advance')

    def scan(self, parameters, starts, transient=1000, events=4000, subdivisions=64):
        p = np.ascontiguousarray(parameters, dtype=float)
        starts = np.ascontiguousarray(starts, dtype=float); lanes = len(p)
        if p.shape != (lanes, 4) or starts.shape != (lanes, 6): raise ValueError('Input shape mismatch')
        if events < 4 or events % 4 or transient < 0 or subdivisions < 16: raise ValueError('Invalid integration settings')
        k, b, dm, dp = p.T
        if not (np.isfinite(p).all() and np.isfinite(starts).all() and np.all(k>b*b/4)
                and np.all(b>0) and np.all((p[:,2:]>0)&(p[:,2:]<1))): raise ValueError('Invalid parameters')
        if np.any((starts[:,[2,5]]<0)|(starts[:,[2,5]]>1)): raise ValueError('Thermal fractions must be in [0,1]')
        state = np.zeros((lanes, 20)); state[:,:6] = starts
        state[:,6:12] = np.array([1.,.3,.4,-.2,.5,-.1])/np.sqrt(1.55)
        state[:,[7,10]] *= np.sqrt(k)[:,None]
        state[:,12] = (starts[:,3]>dp)|((starts[:,3]==dp)&(starts[:,4]>0))
        state[:,13] = (starts[:,0]<dm)|((starts[:,0]==dm)&(starts[:,1]<0))
        state[:,14:16] = np.inf
        state[:,16] = np.max(abs(starts[:,[0,3]]),axis=1)
        state[:,17] = np.min(starts[:,[2,5]],axis=1); state[:,18] = np.max(starts[:,[2,5]],axis=1)
        logs = np.zeros((lanes,4)); durations = np.zeros((lanes,4)); status = np.zeros(lanes,np.int32)
        arrays=[state,p,logs,durations,status]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter(); last=wall; done=0
        for window,count in [(-1,transient)]+[(i,events//4) for i in range(4)]:
            for offset in range(0,count,4):
                length=min(4,count-offset)
                self.kernel(self.queue,(lanes,),None,*buffers,np.int32(lanes),np.int32(length),np.int32(window),np.int32(subdivisions)).wait()
                done+=length
                if time.perf_counter()-last>20:
                    print(f'{done}/{transient+events} events, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
        for a,buffer in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,buffer).wait()
        with np.errstate(divide='ignore',invalid='ignore'):
            lle=logs.sum(axis=1)/durations.sum(axis=1); windows=logs/durations
        lle[status!=0]=np.nan;windows[status!=0]=np.nan
        return dict(parameters=p,starts=starts,final_states=state[:,:6],lle=lle,lle_windows=windows,
                    log_stretch=logs,durations=durations,status=status,minimum_crossing_speed=state[:,14],
                    minimum_event_interval=state[:,15],maximum_abs_displacement_bracket_evaluations=state[:,16],
                    minimum_temperature=state[:,17],maximum_temperature=state[:,18],event_minus_fraction=state[:,19]/(transient+events),
                    settings=dict(transient_events=transient,observation_events=events,subdivisions=subdivisions),
                    method='Analytic affine flow, first root on bracketing mesh, saltation tangent; Euclidean metric on (x,v/sqrt(k),y) blocks',
                    status_codes={'0':'completed','1':'event horizon exhausted','2':'grazing/short interval','3':'nonfinite state','4':'zero tangent','5':'simultaneous events'},
                    device=self.device.name,wall_seconds=time.perf_counter()-wall)


def save_result(path, data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,**{key:value for key,value in data.items() if isinstance(value,np.ndarray)})
    info={key:value for key,value in data.items() if not isinstance(value,np.ndarray)}
    good=data['status']==0; candidate=good&(data['lle']>.005)&np.all(data['lle_windows']>.002,axis=1)
    info.update(lanes=len(good),completed=int(good.sum()),preliminary_candidates=int(candidate.sum()),
                candidate_indices=np.flatnonzero(candidate).tolist(),maximum_lle=float(np.max(data['lle'][good])) if good.any() else None)
    path.with_suffix('.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2),flush=True)


def validate():
    p=np.array([400.,16.,1.5/(25*.0017*176),1/(25*.001*200)])
    start=np.array([.1,.3,.25,.3,-.2,.4]);modes=initial_modes(start,p)
    exact=segment(start,p,modes,.07)
    numerical=solve_ivp(lambda t,y:flow(y,p,modes),(0,.07),start,method='DOP853',rtol=1e-12,atol=1e-14).y[:,-1]
    end,matrix,event_count=fixed_time(start,p,1.2)
    eps=2e-7
    fd=np.column_stack([(fixed_time(start+eps*e,p,1.2)[0]-fixed_time(start-eps*e,p,1.2)[0])/(2*eps) for e in np.eye(6)])
    state=start.copy(); tangent=np.array([1.,.3*np.sqrt(p[0]),.4,-.2,.5*np.sqrt(p[0]),-.1])/np.sqrt(1.55)
    elapsed=0.;logs=0.
    for _ in range(12):
        state,changed,t,salt,which=next_event(state,p,modes)
        tangent=salt@segment_matrix(p,t)@tangent
        norm=np.linalg.norm(tangent/np.array([1,np.sqrt(p[0]),1,1,np.sqrt(p[0]),1]))
        tangent/=norm;logs+=np.log(norm);elapsed+=t;modes=changed
    gpu=OperatorsGPU().scan(p[None],start[None],transient=0,events=12)
    boundary_start=state.copy()
    for _ in range(4):state,modes,_,_,_=next_event(state,p,modes)
    resumed=OperatorsGPU().scan(p[None],boundary_start[None],transient=0,events=4)
    out=dict(affine_vs_DOP853=float(np.max(abs(exact-numerical))),
             fixed_time_saltation_jacobian_error=float(np.max(abs(matrix-fd))),fixed_time_events=event_count,
             determinant_log_error=float(abs(np.linalg.slogdet(matrix)[1]+2*(p[1]+1)*1.2)),
             gpu_cpu_12_event_state_error=float(np.max(abs(boundary_start-gpu['final_states'][0]))),
             gpu_cpu_12_event_lle_error=float(abs(logs/elapsed-gpu['lle'][0])),gpu_status=int(gpu['status'][0]),parameters=p.tolist())
    out['boundary_restart_error']=float(np.max(abs(state-resumed['final_states'][0])))
    out['passed']=out['affine_vs_DOP853']<1e-10 and out['fixed_time_saltation_jacobian_error']<1e-5 and out['gpu_cpu_12_event_state_error']<1e-8 and out['gpu_cpu_12_event_lle_error']<1e-8 and out['gpu_status']==0 and out['boundary_restart_error']<1e-8 and out['determinant_log_error']<1e-7
    Path('results/light-operators').mkdir(parents=True,exist_ok=True)
    Path('results/light-operators/coupled-validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
    if not out['passed']:raise RuntimeError('Validation failed')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--scout',action='store_true');parser.add_argument('--refine');parser.add_argument('--output');parser.add_argument('--transient',type=int,default=1000);parser.add_argument('--events',type=int,default=4000);parser.add_argument('--subdivisions',type=int,default=64)
    args=parser.parse_args()
    if args.scout:
        from scipy.stats import qmc
        samples=qmc.Sobol(4,scramble=True,seed=2026102963).random_base2(11)
        k=10**(1+3*samples[:,0]); damping_bar=10**(-1.3+1.5*samples[:,1]);b=damping_bar*np.sqrt(k)
        # The exact source reference is inserted, then logarithmic mechanics and interior thresholds are explored.
        p=np.column_stack([k,b,.005+.99*samples[:,2],.005+.99*samples[:,3]])
        p[0]=[400,16,1.5/(25*.0017*176),1/(25*.001*200)]
        p=np.repeat(p,2,axis=0);rng=np.random.default_rng(2026102964)
        starts=rng.uniform(.05,.95,(len(p),6));starts[:,[1,4]]=0
        data=OperatorsGPU().scan(p,starts,args.transient,args.events,args.subdivisions)
        save_result(args.output or 'results/light-operators/coupled-scout.npz',data)
    elif args.refine:
        source=np.load(args.refine); good=source['status']==0
        candidates=np.flatnonzero(good&(source['lle']>.005)&np.all(source['lle_windows']>.002,axis=1))
        if not len(candidates): candidates=np.flatnonzero(good)[np.argsort(source['lle'][good])[-24:]]
        data=OperatorsGPU().scan(source['parameters'][candidates],source['final_states'][candidates],args.transient,args.events,args.subdivisions)
        data['source_indices']=candidates;data['source_file']=str(args.refine)
        save_result(args.output or 'results/light-operators/coupled-refined.npz',data)
    else:validate()
