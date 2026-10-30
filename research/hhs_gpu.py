"""Adaptive FP64 search of the literal autonomous five-state HHS single cell."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.stats import qmc
from hhs import ROOT,MREF,EREF,field_jac,save


class HHSGPU:
    def __init__(self):
        devices=[d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        source=Path(__file__).with_name('hhs.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(),'advance')

    def scan(self,p,starts,transient=500.,duration=1000.,rtol=1e-8,max_step=.05,max_wall_seconds=600.,max_attempts=6000000,resume_arrays=None):
        p=np.ascontiguousarray(p,dtype=float);starts=np.asarray(starts,dtype=float);n=len(p)
        if p.shape!=(n,3) or starts.shape!=(n,5) or n<1 or not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Invalid inputs')
        if np.min(p[:,0])<0 or np.min(p[:,1:])<=0 or np.min(starts[:,:3])<0 or transient<0 or min(duration,rtol,max_step,max_wall_seconds)<=0:raise ValueError('Invalid settings')
        pars=np.zeros((n,4));pars[:,:3]=p
        states=np.zeros((n,16));states[:,:5]=starts;tangent=np.arange(1,6,dtype=float);tangent/=np.linalg.norm(tangent);states[:,8:13]=tangent
        clocks=np.zeros((n,4));clocks[:,1]=.001;logs=np.zeros((n,4));stats=np.zeros((n,16));stats[:,:10:2]=np.inf;stats[:,1:10:2]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        if resume_arrays is not None:
            states[:,:5]=resume_arrays['final_states'];states[:,8:13]=resume_arrays['final_tangents'];clocks[:]=resume_arrays['raw_clocks']
            logs[:]=resume_arrays['raw_window_log_sums'];stats[:]=resume_arrays['stats'];counts[:]=resume_arrays['step_counts']
        arrays=[states,pars,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();last=wall;calls=0
        while True:
            self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(128),np.float64(transient),np.float64(duration),np.float64(rtol),np.float64(rtol*.001),np.float64(max_step),np.int32(max_attempts)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<transient+duration)
                if time.perf_counter()-last>25:print(f'{active.sum()} active, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-wall>max_wall_seconds:status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        maximum_filter=np.maximum(stats[:,7],stats[:,9])
        alpha_positive_bound=np.full(n,np.inf);positive=maximum_filter>0;alpha_positive_bound[positive]=1/maximum_filter[positive]
        return dict(parameters=p,initial_states=starts,initial_tangent=tangent,final_states=states[:,:5].copy(),final_tangents=states[:,8:13].copy(),
            lle=windows.mean(axis=1),lle_windows=windows,raw_window_log_sums=logs,raw_clocks=clocks,status=status,stats=stats,elapsed=clocks[:,0],step_counts=counts,
            alpha_upper_bound_observed=alpha_positive_bound,minimum_observed_length_at_alpha_point1=1-.1*maximum_filter,
            settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol*.001,max_step=max_step,max_attempts=max_attempts),
            method='Adaptive DP5(4), exact polynomial variational equations, state and tangent error control',device=self.device.name,wall_seconds=time.perf_counter()-wall,
            parameter_columns=['mu','tauD','tauL'],fixed_parameters=dict(a=1,b=3,c=1,D0=.1,tauR=1,Eref=EREF),state_columns=['M','E','D','y','x'],
            time_unit='source time in units of tauR, not calibrated seconds',
            physical_scope='Observed alpha bound from accepted observation states only; L=1-alpha*x,L0=1-alpha*y,beta=mu/alpha. No clipping or physical calibration.',
            status_codes={'0':'completed','1':'numerical magnitude cutoff 1e6','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved','5':'negative concentration beyond 1e-8, invalid numerical trajectory'})


def controls():
    gpu=HHSGPU();records=[]
    for p in [np.array([0.,1.1,1.2]),np.array([.8,1.1,1.2]),np.array([1.5,1.1,1.2]),np.array([4.,.1,10.])]:
        for equilibrium in [False,True]:
            y=np.array([MREF,EREF,.1,0,0]) if equilibrium else np.array([.44,2.71,.1,0,0])
            out=gpu.scan(p[None],y[None],transient=0,duration=2,rtol=1e-10,max_step=.005)
            def augmented(t,z):
                f,j=field_jac(z[:5],p);return np.r_[f,j@z[5:]]
            sol=solve_ivp(augmented,(0,2),np.r_[y,out['initial_tangent']],method='Radau',rtol=1e-11,atol=1e-13,max_step=.002)
            end=sol.y[:,-1];norm=np.linalg.norm(end[5:])
            errors=dict(state=float(np.max(abs(out['final_states'][0]-end[:5]))),tangent=float(np.max(abs(out['final_tangents'][0]-end[5:]/norm))),lle=float(abs(out['lle'][0]-np.log(norm)/2)))
            if equilibrium:
                ref=expm(field_jac(y,p)[1]*2)@out['initial_tangent'];errors['linear_variational_expm']=float(np.max(abs(out['final_tangents'][0]-ref/np.linalg.norm(ref))))
            passed=bool(sol.success and out['status'][0]==0 and max(errors.values())<1e-7)
            records.append(dict(parameters=p.tolist(),at_reference_equilibrium=equilibrium,errors=errors,passed=passed));print(records[-1],flush=True)
            if not passed:raise RuntimeError('GPU control failed')
    save(ROOT/'gpu-validation',dict(records=records,passed=True,scope='Independent implicit Radau state/tangent, and exact equilibrium matrix-exponential variational controls'))


def survey():
    samples=qmc.Sobol(3,scramble=True,seed=2026103015).random_base2(10)
    lower=np.log([.001,.05,.05]);upper=np.log([30.,30.,30.]);broad=np.exp(lower+samples*(upper-lower))
    section=np.column_stack([np.linspace(0,3,65),np.full(65,1.1),np.full(65,1.2)])
    tuples=np.vstack([section,broad]);p=np.repeat(tuples,3,axis=0)
    seeds=np.array([[1.05,3.43,.1,0,0],[.44,2.71,.1,0,0],[.49,4.21,.1,0,0]]);starts=np.tile(seeds,(len(tuples),1))
    out=HHSGPU().scan(p,starts)
    candidate=(out['status']==0)&(out['lle']>.002)&np.all(out['lle_windows']>.001,axis=1)
    out.update(preliminary_candidate=candidate,source_section_count=len(section),broad_tuple_count=len(broad),lanes=len(p),completed=int((out['status']==0).sum()),
        preliminary_positive=int(candidate.sum()),selection_seed=2026103015,
        start_scope='Three rounded chemical starts printed in source Fig.S1, D=.1,x=y=0; fixed main-text Eref convention, not an asserted reproduction of ambiguous source plots.')
    save(ROOT/'scout',out);print({k:out[k] for k in ['lanes','completed','preliminary_positive','wall_seconds']},flush=True)


def refine(file):
    old=np.load(file);ok=np.flatnonzero(old['status']==0)
    weak=(old['status']==0)&np.all(old['lle_windows']>0,axis=1)
    indices=np.union1d(np.flatnonzero(weak|old['preliminary_candidate']),ok[np.argsort(old['lle'][ok])[-24:]])
    out=HHSGPU().scan(old['parameters'][indices],old['final_states'][indices],transient=2000,duration=8000,rtol=1e-10,max_step=.025,max_wall_seconds=900,max_attempts=12000000)
    candidate=(out['status']==0)&(out['lle']>.0005)&np.all(out['lle_windows']>.00025,axis=1)
    out.update(source_indices=indices,retained_candidate=candidate,retained=int(candidate.sum()),source=Path(file).name,selection='All four positive windows, union preliminary passes and top24 means')
    save(ROOT/(Path(file).stem+'-refined'),out);print({'lanes':len(indices),'completed':int((out['status']==0).sum()),'retained':out['retained'],'wall_seconds':out['wall_seconds']},flush=True)


def long_control():
    source=np.load(ROOT/'uncoupled-cycle.npz');p=source['parameters'];y=source['initial_state'];gpu=HHSGPU();records=[]
    tangent=np.arange(1,6,dtype=float);tangent/=np.linalg.norm(tangent)
    def augmented(t,z):
        f,j=field_jac(z[:5],p);return np.r_[f,j@z[5:]]
    sol=solve_ivp(augmented,(0,250),np.r_[y,tangent],method='Radau',rtol=1e-11,atol=1e-13,max_step=.02)
    end=sol.y[:,-1];norm=np.linalg.norm(end[5:])
    for rtol,step in [(1e-8,.05),(1e-10,.025)]:
        out=gpu.scan(p[None],y[None],transient=0,duration=250,rtol=rtol,max_step=step)
        errors=dict(state=float(np.max(abs(out['final_states'][0]-end[:5]))),tangent=float(np.max(abs(out['final_tangents'][0]-end[5:]/norm))),lle=float(abs(out['lle'][0]-np.log(norm)/250)))
        passed=bool(sol.success and out['status'][0]==0 and max(errors.values())<1e-7)
        records.append(dict(rtol=rtol,max_step=step,errors=errors,passed=passed));print(records[-1],flush=True)
        if not passed:raise RuntimeError('Long cycle control failed')
    save(ROOT/'long-gpu-validation',dict(records=records,passed=True,parameters=p,initial_state=y,duration=250,
        scope='250 source-time units on independently closed mu=0 cycle, two GPU tolerances/steps vs implicit Radau state+tangent'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true');parser.add_argument('--long-control',action='store_true');parser.add_argument('--refine');a=parser.parse_args()
    if a.controls:controls()
    elif a.long_control:long_control()
    elif a.refine:refine(a.refine)
    else:survey()
