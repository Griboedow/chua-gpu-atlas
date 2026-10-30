"""Chua along beta=kappa*alpha**2, integrated in exact scaled coordinates.

State is (x,y,Z=z/alpha); fast time s=alpha*t. Requested durations and
reported exponents are in ORIGINAL dimensionless time t. No altered diode.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results'/'chua-scaling'
CLOCK_MODE = 'compensated'


def field_jac(y, p):
    alpha, kappa, m0, m1 = p
    x, v, Z = y
    diode = m1*x+(m0-m1)*np.clip(x, -1, 1)
    slope = m0 if abs(x) < 1 else m1
    return (np.array([v-x-diode, Z+(x-v)/alpha, -kappa*v]),
            np.array([[-1-slope, 1, 0], [1/alpha, -1/alpha, 1], [0, -kappa, 0]]))


def save(path, data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    metadata={k:v for k,v in data.items() if k not in arrays}
    np.savez_compressed(path.with_suffix('.npz'),**arrays)
    metadata['arrays']={k:dict(shape=list(v.shape),dtype=str(v.dtype)) for k,v in arrays.items()}
    path.with_suffix('.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n',encoding='utf-8')


class ScalingGPU:
    def __init__(self,clock_mode=None):
        devices=[d for p in cl.get_platforms() for d in p.get_devices()
                 if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 GPU required')
        self.device=devices[0];self.context=cl.Context([self.device])
        self.queue=cl.CommandQueue(self.context)
        self.clock_mode=clock_mode or CLOCK_MODE
        source=Path(__file__).with_name('chua_scaling_gpu_legacy.cl' if self.clock_mode=='legacy' else 'chua_scaling_gpu.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(),'adaptive_advance')

    def scan(self,p,starts,transient=100.,duration=200.,rtol=1e-8,max_step=.25,
             max_attempts=12000000,max_wall_seconds=600.,escape=1e6,resume_arrays=None):
        p=np.ascontiguousarray(p,dtype=float);starts=np.asarray(starts,dtype=float);n=len(p)
        if p.shape!=(n,4) or starts.shape!=(n,3) or n<1 or not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Invalid arrays')
        if np.min(p[:,:2])<=0 or transient<0 or min(duration,rtol,max_step,max_wall_seconds,escape,max_attempts)<=0:raise ValueError('Invalid settings')
        state=np.zeros((n,8));state[:,:3]=starts
        tangent=np.array([1.,2.,3.]);tangent/=np.linalg.norm(tangent);state[:,3:6]=tangent
        clocks=np.zeros((n,4));clocks[:,1]=.001
        logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,[0,2,4,6]]=np.inf;stats[:,[1,3,5,7]]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        boundaries=np.full((n,5,3),np.nan);boundary_times=np.full((n,5),np.nan)
        resume_time=np.zeros(n)
        if resume_arrays is not None:
            state[:,:3]=resume_arrays['final_states'];state[:,3:6]=resume_arrays['final_tangents']
            clocks[:]=resume_arrays['raw_fast_clocks'];logs[:]=resume_arrays['raw_window_log_sums']
            stats[:]=resume_arrays['stats'];counts[:]=resume_arrays['step_counts']
            boundaries[:]=resume_arrays['boundary_states'];boundary_times[:]=resume_arrays['boundary_original_times']
            resume_time=clocks[:,0]/p[:,0]
            # A rejected step at the old absolute-clock guard did not advance
            # state/time/logs. Retry that step in the compensated-clock backend.
            clocks[:,1]=np.maximum(clocks[:,1],1e-15)
        arrays=[state,p,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();last=wall;target=0.;calls=0
        for b,(window,length) in enumerate([(-1,transient)]+[(i,duration/4) for i in range(4)]):
            target+=length
            if length:
                while True:
                    self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(128),
                                np.float64(target),np.int32(window),np.float64(rtol),np.float64(rtol*.001),
                                np.float64(escape),np.float64(max_step),np.int32(max_attempts)).wait();calls+=1
                    if calls%8==0:
                        cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                        active=(status==0)&(clocks[:,0]<target*p[:,0]*(1-1e-13))
                        if time.perf_counter()-last>25:
                            print(f'window {window}: {active.sum()} active, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
                        if not active.any():break
                        if time.perf_counter()-wall>max_wall_seconds:
                            status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
            cl.enqueue_copy(self.queue,state,buffers[0]).wait();cl.enqueue_copy(self.queue,clocks,buffers[2]).wait()
            update=resume_time<=target
            boundaries[update,b]=state[update,:3];boundary_times[update,b]=clocks[update,0]/p[update,0]
        for a,buff in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,buff).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,initial_states=starts,initial_tangent=tangent,final_states=state[:,:3].copy(),
                    final_tangents=state[:,3:6].copy(),lle=windows.mean(axis=1),lle_windows=windows,
                    raw_window_log_sums=logs,raw_fast_clocks=clocks,elapsed_original_time=clocks[:,0]/p[:,0],
                    status=status,step_counts=counts,stats=stats,boundary_states=boundaries,boundary_original_times=boundary_times,
                    settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol*.001,max_fast_step=max_step,
                                  max_attempts=max_attempts,max_wall_seconds=max_wall_seconds,escape_scaled_norm_max=escape),
                    method='FP64 DP5(4), full piecewise field and tangent, state and tangent error control; no explicit kink location',
                    clock_mode=self.clock_mode,
                    device=self.device.name.strip(),wall_seconds=time.perf_counter()-wall,
                    parameter_columns=['alpha','kappa=beta/alpha^2','m0','m1'],state_columns=['x','y','Z=z/alpha'],
                    time_unit='Requested times and LLE are ORIGINAL dimensionless t; raw clocks are fast s=alpha*t',
                    tangent_norm='Euclidean norm in scaled (x,y,Z), not the historical (x,y,z) finite-time norm',
                    stats_columns=['xmin','xmax','ymin','ymax','Zmin','Zmax','last_window_xmin','last_window_xmax'],
                    status_codes={'0':'completed','1':'scaled magnitude cutoff (not proof of divergence)','2':'numerical failure',
                                  '3':'step budget unresolved','4':'wall budget unresolved'})


def controls():
    gpu=ScalingGPU();records=[]
    for p,y,S in [(np.array([9.,100/567,-8/7,-5/7]),np.array([.7,0,0]),90.),
                  (np.array([4096.,.2,-5.,-5/7]),np.array([.7,.5,-.1]),25.),
                  (np.array([128.,.24,-8/7,-5/7]),np.array([1.5,0,-1.5/128]),20.)]:
        out=gpu.scan(p[None],y[None],transient=0,duration=S/p[0],rtol=1e-11,max_step=.02)
        def rhs(s,z):
            f,j=field_jac(z[:3],p);return np.r_[f,j@z[3:]]
        sol=solve_ivp(rhs,(0,S),np.r_[y,out['initial_tangent']],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.01)
        end=sol.y[:,-1];norm=np.linalg.norm(end[3:])
        # Independent original coordinates and original time: no call to scaled field.
        a,k,m0,m1=p
        def original(t,v):
            x,y,z=v;h=m1*x+.5*(m0-m1)*(abs(x+1)-abs(x-1))
            return [a*(y-x-h),x-y+z,-k*a*a*y]
        orig=solve_ivp(original,(0,S/a),y*np.array([1,1,a]),method='Radau',rtol=2e-12,atol=2e-14,max_step=.02/a)
        errors=dict(state_scaled=float(np.max(abs(out['final_states'][0]-end[:3]))),
                    normalized_tangent=float(np.max(abs(out['final_tangents'][0]-end[3:]/norm))),
                    fast_log_growth=float(abs(out['raw_window_log_sums'][0].sum()-np.log(norm))),
                    independent_original_coordinate_state=float(np.max(abs(out['final_states'][0]-orig.y[:,-1]/[1,1,a]))))
        passed=bool(sol.success and orig.success and out['status'][0]==0 and max(errors.values())<2e-7)
        records.append(dict(parameters=p.tolist(),start=y.tolist(),fast_duration=S,errors=errors,passed=passed))
        print(records[-1],flush=True)
        if not passed:raise RuntimeError('Scaling control failed')
    p=np.array([4096.,.24,-8/7,-5/7]);xeq=(p[3]-p[2])/(1+p[3]);y=np.array([xeq,0,-xeq/p[0]])
    out=gpu.scan(p[None],y[None],transient=0,duration=10/p[0],rtol=1e-12,max_step=.01)
    expected=expm(field_jac(y,p)[1]*10)@out['initial_tangent'];norm=np.linalg.norm(expected)
    errors=dict(state=float(np.max(abs(out['final_states'][0]-y))),
                tangent=float(np.max(abs(out['final_tangents'][0]-expected/norm))),
                log_growth=float(abs(out['raw_window_log_sums'][0].sum()-np.log(norm))))
    passed=bool(out['status'][0]==0 and max(errors.values())<1e-9)
    records.append(dict(kind='exact outer equilibrium and matrix exponential',errors=errors,passed=passed))
    if not passed:raise RuntimeError('Linear control failed')
    save(OUT/'gpu-validation',dict(records=records,passed=True))


def survey():
    alphas=np.array([9.,16.,32.,64.,128.,256.,512.,1024.,2048.])
    kappas=np.unique(np.r_[np.linspace(.03,.28,41),10/49])
    slopes=np.array([-8/7,-2.,-5.,-10.,-30.])
    tuples=np.array([[a,k,m,-5/7] for m in slopes for a in alphas for k in kappas])
    # Fixed scaled finite oscillator action, plus conventional z=0 start.
    starts=np.tile([[.7,0,0],[.7,.5,0],[-.7,2.,0]],(len(tuples),1))
    p=np.repeat(tuples,3,axis=0)
    out=ScalingGPU().scan(p,starts)
    candidate=(out['status']==0)&(out['lle']>.01)&np.all(out['lle_windows']>.005,axis=1)
    out.update(preliminary_candidate=candidate,alphas=alphas,kappas=kappas,slopes=slopes,
               lanes=len(p),completed=int((out['status']==0).sum()),preliminary_positive=int(candidate.sum()),
               start_scope='Three fixed scaled-coordinate starts, one zero-action, two finite-action; not a basin census')
    save(OUT/'scout',out)
    print({k:out[k] for k in ['lanes','completed','preliminary_positive','wall_seconds']},flush=True)


def refine():
    old=np.load(OUT/'scout.npz')
    indices=np.flatnonzero(old['preliminary_candidate'])
    out=ScalingGPU().scan(old['parameters'][indices],old['final_states'][indices],
                          transient=100,duration=400,rtol=1e-11,max_step=.125,
                          max_attempts=40000000,max_wall_seconds=900.)
    candidate=(out['status']==0)&(out['lle']>.01)&np.all(out['lle_windows']>.005,axis=1)
    out.update(source_indices=indices,retained_candidate=candidate,lanes=len(indices),
               completed=int((out['status']==0).sum()),retained_positive=int(candidate.sum()),
               source_scout_sha256=hashlib.sha256((OUT/'scout.npz').read_bytes()).hexdigest(),
               selection='Every scout preliminary positive, continued from its saved endpoint with a fresh tangent')
    save(OUT/'refined',out)
    print({k:out[k] for k in ['lanes','completed','retained_positive','wall_seconds']},flush=True)


def resume_clock():
    old=np.load(OUT/'refined.npz');indices=np.flatnonzero(old['status']==2)
    if not len(indices):
        print('No numerical-failure lanes to resume. Any older resume dataset is not applicable.');return
    resume={k:old[k][indices] for k in ['final_states','final_tangents','raw_fast_clocks','raw_window_log_sums',
                                      'stats','step_counts','boundary_states','boundary_original_times']}
    out=ScalingGPU().scan(old['parameters'][indices],old['initial_states'][indices],
                          transient=100,duration=400,rtol=1e-11,max_step=.125,max_attempts=40000000,
                          max_wall_seconds=600.,resume_arrays=resume)
    candidate=(out['status']==0)&(out['lle']>.01)&np.all(out['lle_windows']>.005,axis=1)
    out.update(refined_indices=indices,source_indices=old['source_indices'][indices],retained_candidate=candidate,
               lanes=len(indices),completed=int((out['status']==0).sum()),retained_positive=int(candidate.sum()),
               source_refined_sha256=hashlib.sha256((OUT/'refined.npz').read_bytes()).hexdigest(),
               selection=f'Only {len(indices)} previous refined numerical failures, resumed from exact saved state/tangent/logs with compensated clock; original results retained')
    save(OUT/'refined-clock-resume',out)
    print({k:out[k] for k in ['lanes','completed','retained_positive','wall_seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['controls','survey','refine','resume_clock'])
    parser.add_argument('--clock-mode',choices=['compensated','legacy'],default='compensated',help='legacy reproduces the historical large-clock guard failures')
    args=parser.parse_args();CLOCK_MODE=args.clock_mode;globals()[args.mode]()
