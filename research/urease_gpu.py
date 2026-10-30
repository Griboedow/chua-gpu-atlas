"""Validated implicit OpenCL integration of the full four-species urease model."""
import argparse
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from urease import ROOT,SCALE,COLUMNS,TRANSFORM,INVERSE,parameters,transport_rhs_jac,save

TANGENT=np.array([.3,.2,-.5,.8]);TANGENT/=np.linalg.norm(TANGENT)


class UreaseGPU:
    def __init__(self):
        devices=[d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 OpenCL GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        self.kernel=cl.Kernel(cl.Program(self.context,Path(__file__).with_name('urease.cl').read_text()).build(),'advance')

    def scan(self,p,starts,transient=20000.,duration=40000.,rtol=1e-7,max_step=10.,max_attempts=1000000,max_wall=600.):
        p=np.ascontiguousarray(p,dtype=float);starts=np.ascontiguousarray(starts,dtype=float);n=len(p)
        if p.shape!=(n,11) or starts.shape!=(n,4) or n<1 or not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Invalid inputs')
        if np.min(p)<=0 or np.min(starts)<0 or np.min(starts[:,1])<=0:raise ValueError('Nonphysical inputs')
        if min(duration,rtol,max_step,max_wall)<=0 or transient<0:raise ValueError('Invalid settings')
        states=np.zeros((n,8));states[:,:4]=starts@TRANSFORM.T;states[:,4:]=TANGENT
        clocks=np.zeros((n,2));clocks[:,1]=.001;logs=np.zeros((n,4));ranges=np.empty((n,8));ranges[:,::2]=np.inf;ranges[:,1::2]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[states,p,clocks,logs,ranges,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();last=wall;calls=0
        while True:
            self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(32),np.float64(transient),np.float64(duration),np.float64(rtol),np.float64(rtol*.001),np.float64(max_step),np.int32(max_attempts)).wait();calls+=1
            if calls%16==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<transient+duration)
                if time.perf_counter()-last>20:print(f'{active.sum()} active, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-wall>max_wall:status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,initial_states=starts,final_states=states[:,:4]@INVERSE.T,transport_final_states=states[:,:4].copy(),
            final_tangents=states[:,4:].copy(),lle=windows.mean(axis=1),lle_windows=windows,status=status,ranges=ranges,elapsed=clocks[:,0],step_counts=counts,
            settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol*.001,max_step=max_step,max_attempts=max_attempts,max_wall=max_wall),
            device=self.device.name,wall_seconds=time.perf_counter()-wall,state_scale_molar=SCALE,time_unit='s',parameter_columns=COLUMNS,
            state_columns=['S','H','P','Q'],tangent_coordinates=['S','P-H','H+Q','H'],
            method='SDIRK2 base, step doubling, Richardson extrapolation; analytic discrete tangent through converged implicit stages',
            stability_scope='The SDIRK2 base is L-stable, but its accepted Richardson extrapolate is not unconditionally A-stable; error control and independent CPU comparisons are essential',
            status_codes={'0':'completed','2':'step underflow or numerical failure','3':'step budget unresolved','4':'wall budget unresolved'})


def validate():
    p=parameters();gpu=UreaseGPU();reference=np.load(ROOT/'reference.npz')
    starts=np.array([[.5,.1,0,0],reference['states'][12000]])
    rows=[]
    for i,y in enumerate(starts):
        duration=30. if i==0 else 1200.
        def aug(t,z):
            f,j=transport_rhs_jac(z[:4],p);return np.r_[f,j@z[4:]]
        initial=np.r_[TRANSFORM@y,TANGENT]
        cpu=solve_ivp(aug,(0,duration),initial,method='Radau',rtol=1e-10,atol=1e-13,max_step=1.)
        if not cpu.success:raise RuntimeError(cpu.message)
        target=cpu.y[:4,-1];v=cpu.y[4:,-1];lle=np.log(np.linalg.norm(v))/duration
        for tol in [1e-6,1e-8]:
            result=gpu.scan(p[None],y[None],transient=0,duration=duration,rtol=tol,max_step=5.,max_wall=120.)
            row=dict(control=i,rtol=tol,duration=duration,status=int(result['status'][0]),
                     state_error=float(np.max(abs(result['transport_final_states'][0]-target)/(1+abs(target)))),
                     lle_error=float(abs(result['lle'][0]-lle)),cpu_lle=float(lle),gpu_lle=float(result['lle'][0]),
                     tangent_direction_error=float(np.max(abs(result['final_tangents'][0]-v/np.linalg.norm(v)))),
                     accepted_rejected=result['step_counts'][0].tolist(),wall_seconds=result['wall_seconds'])
            rows.append(row);print(row,flush=True)
    passed=all(r['status']==0 and r['state_error']<2e-4 and r['lle_error']<1e-5 and r['tangent_direction_error']<1e-3 for r in rows)
    out=dict(rows=rows,passed=passed,device=gpu.device.name,comparison='Independent implicit Radau in same exact linear coordinates; full state and tangent',
             scope='One initial fast relaxation and one long cycle-scale control; not asymptotic classification')
    save(ROOT/'gpu-validation',out)
    if not passed:raise RuntimeError('Implicit GPU validation failed')


def survey(tuples=512):
    rng=np.random.default_rng(2026103006);p=np.tile(parameters(),(tuples,1))
    # Fixed acid/base chemistry; vary open-reservoir and transport conditions.
    for col,lo,hi in [(5,2e-4,.02),(6,.001,.08),(8,2e-5,.002),(9,5e-5,.003),(10,1e-5,.002)]:
        p[:,col]=np.exp(rng.uniform(np.log(lo),np.log(hi),tuples))
    p[:,7]=p[:,5]
    p[tuples//2:,7]*=np.exp(rng.uniform(np.log(.1),np.log(10),tuples-tuples//2))
    p[0]=parameters();p[1]=parameters();p[1,6]=p[1,5]
    p=np.repeat(p,2,axis=0)
    starts=[]
    for i,row in enumerate(p):
        se,he=row[9:11]/SCALE
        starts.append([.3*se,.2*he,0.,0.] if i%2==0 else [.1*se,.02*he,.1*se,1.7*se])
    result=UreaseGPU().scan(p,np.array(starts))
    positive=(result['status']==0)&(result['lle']>2e-5)&np.all(result['lle_windows']>1e-5,axis=1)
    result.update(preliminary_candidate=positive,preliminary_count=int(positive.sum()),seed=2026103006,
                  scope='Original four-state source reaction network; varied transport and reservoir settings, not fitted experiments',
                  cohort_rule='First half tuples use k=kS; second half vary product outflow. Tuple0 source, tuple1 equal-transport no-cycle control.')
    save(ROOT/'scout',result);print({'preliminary':int(positive.sum()),'completed':int((result['status']==0).sum()),'wall_seconds':result['wall_seconds']},flush=True)


def parameter_controls():
    gpu=UreaseGPU();base=parameters();rows=[]
    variants=[base.copy() for _ in range(4)]
    variants[0][6]=base[5]
    variants[1][7]*=.1
    variants[2][7]*=10
    variants[3][5:8]*=5;variants[3][8]*=3;variants[3][9]*=2
    for i,p in enumerate(variants):
        y=np.array([.3*p[9]/SCALE,.2*p[10]/SCALE,0.,0.]);initial=np.r_[TRANSFORM@y,TANGENT]
        def aug(t,z):
            f,j=transport_rhs_jac(z[:4],p);return np.r_[f,j@z[4:]]
        cpu=solve_ivp(aug,(0,3000),initial,method='Radau',rtol=1e-10,atol=1e-13,max_step=1.)
        if not cpu.success:raise RuntimeError(cpu.message)
        out=gpu.scan(p[None],y[None],transient=0,duration=3000,rtol=1e-7,max_step=5.,max_wall=120.)
        v=cpu.y[4:,-1];ref=np.log(np.linalg.norm(v))/3000
        row=dict(parameters=p.tolist(),control=i,status=int(out['status'][0]),state_error=float(np.max(abs(out['transport_final_states'][0]-cpu.y[:4,-1])/(1+abs(cpu.y[:4,-1])))),
                 lle_error=float(abs(out['lle'][0]-ref)),cpu_lle=float(ref),gpu_lle=float(out['lle'][0]),wall_seconds=out['wall_seconds'])
        rows.append(row);print(row,flush=True)
    passed=all(r['status']==0 and r['state_error']<1e-4 and r['lle_error']<1e-6 for r in rows)
    save(ROOT/'gpu-parameter-controls',dict(rows=rows,passed=passed,method='Independent Radau state+tangent, 3000 seconds at four modified rate sets'))
    if not passed:raise RuntimeError('Parameter controls failed')


def refine():
    old=np.load(ROOT/'scout.npz');ok=np.flatnonzero(old['status']==0)
    ids=np.union1d(np.flatnonzero(old['lle']>0),ok[np.argsort(old['lle'][ok])[-16:]])
    result=UreaseGPU().scan(old['parameters'][ids],old['final_states'][ids],transient=200000.,duration=400000.,rtol=1e-9,max_step=5.,max_attempts=3000000,max_wall=600.)
    positive=(result['status']==0)&(result['lle']>1e-6)&np.all(result['lle_windows']>5e-7,axis=1)
    result.update(source_indices=ids,retained_candidate=positive,retained_count=int(positive.sum()),source='scout.npz',
        selection='All positive mean scouts plus the largest 16 means; long criterion lowered to mean>1e-6 and every window>5e-7 per second')
    save(ROOT/'refined',result);print({'retained':int(positive.sum()),'completed':int((result['status']==0).sum()),'lanes':len(ids),'wall_seconds':result['wall_seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--survey',action='store_true');parser.add_argument('--refine',action='store_true');parser.add_argument('--parameter-controls',action='store_true');parser.add_argument('--tuples',type=int,default=512);a=parser.parse_args()
    if a.survey:survey(a.tuples)
    elif a.refine:refine()
    elif a.parameter_controls:parameter_controls()
    else:validate()
