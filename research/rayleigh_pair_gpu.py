"""FP64 searches of the exact canonical quadratic-drag thermal spring pair."""
import argparse
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from rayleigh_pair import ROOT,SOURCE,field_jac,modal_energy,save


class RayleighGPU:
    def __init__(self):
        devices=[d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        source=Path(__file__).with_name('rayleigh_pair.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(),'advance')

    def scan(self,p,starts,transient=50.,duration=100.,rtol=1e-8,max_step=.05,max_wall_seconds=600.,max_attempts=6000000,time_scales=None,resume_arrays=None):
        p=np.ascontiguousarray(p,dtype=float);starts=np.asarray(starts,dtype=float);n=len(p)
        if p.shape!=(n,2) or starts.shape!=(n,4) or n<1 or not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Invalid input arrays')
        if np.min(p)<0 or transient<0 or min(duration,rtol,max_step,max_wall_seconds)<=0:raise ValueError('Invalid settings')
        scale=np.ones(n)
        if time_scales is None:
            positive=p[:,0]>0;scale[positive]=np.maximum(1.,np.maximum(p[positive,0],1/p[positive,0]))
        else:scale[:]=time_scales
        if not np.isfinite(scale).all() or np.min(scale)<=0:raise ValueError('Invalid time scales')
        pars=np.zeros((n,4));pars[:,:2]=p;pars[:,2]=scale
        states=np.zeros((n,8));states[:,:4]=starts;tangent=np.array([.2,.5,.8,1.]);tangent/=np.linalg.norm(tangent);states[:,4:]=tangent
        clocks=np.zeros((n,4));clocks[:,1]=.001/scale;logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,::2]=np.inf;stats[:,1::2]=-np.inf
        modal=np.zeros((n,4));status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        if resume_arrays is not None:
            states[:,:4]=resume_arrays['final_states'];states[:,4:]=resume_arrays['final_tangents']
            clocks[:]=resume_arrays['raw_clocks'];logs[:]=resume_arrays['raw_window_log_sums']
            stats[:]=resume_arrays['stats'];modal[:]=resume_arrays['modal_statistics'];counts[:]=resume_arrays['step_counts']
        arrays=[states,pars,clocks,logs,stats,modal,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();last=wall;calls=0
        while True:
            self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(128),np.float64(transient),np.float64(duration),np.float64(rtol),np.float64(rtol*.001),np.float64(max_step),np.int32(max_attempts)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[6]).wait()
                active=(status==0)&(clocks[:,0]<transient+duration)
                if time.perf_counter()-last>25:print(f'{active.sum()} active, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-wall>max_wall_seconds:status[active]=4;cl.enqueue_copy(self.queue,buffers[6],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        normalized=logs/(duration/4);normalized[status!=0]=np.nan;windows=normalized/scale[:,None]
        modal_mean=modal[:,:3]/duration;modal_mean[status!=0]=np.nan
        return dict(parameters=p,initial_states=starts,initial_tangent=tangent,final_states=states[:,:4].copy(),final_tangents=states[:,4:].copy(),
            time_scales=scale,lle=windows.mean(axis=1),lle_windows=windows,normalized_lle=normalized.mean(axis=1),normalized_lle_windows=normalized,
            raw_window_log_sums=logs,raw_clocks=clocks,status=status,stats=stats,modal_statistics=modal,modal_means=modal_mean,elapsed_sigma=clocks[:,0]*scale,step_counts=counts,
            settings=dict(transient_scaled=transient,duration_scaled=duration,rtol=rtol,atol=rtol*.001,max_step_sigma=max_step,max_attempts=max_attempts),
            method='Adaptive DP5(4), exact C1 variational equations, computational time tau=sigma/scale',device=self.device.name,wall_seconds=time.perf_counter()-wall,
            parameter_columns=['mu','kappa'],state_columns=['x1','v1','x2','v2'],time_unit='canonical sigma (not physical seconds)',
            scaling_rule='scale=max(1,mu,1/mu) for mu>0, or explicit control override; normalized_lle=lle*scale',
            modal_statistics_columns=['integral_symmetric_energy_dtau','integral_antisymmetric_energy_dtau','integral_minimum_energy_fraction_dtau','max_minimum_energy_fraction'],
            physical_scope='Exact reduced mathematical source family; no unprovided fiber/tension limits imposed',
            status_codes={'0':'completed','1':'numerical magnitude cutoff 1e6','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved'})


def controls():
    gpu=RayleighGPU();records=[];y=np.array([.1,.7,-.3,.2])
    for p in [np.array([0.,2.]),SOURCE,np.array([3.,1.5]),np.array([10.,.01])]:
        for scale in [1.,7.]:
            out=gpu.scan(p[None],y[None],transient=0,duration=2/scale,rtol=1e-10,max_step=.005,time_scales=[scale])
            def augmented(t,z):
                f,j=field_jac(z[:4],p);return np.r_[f,j@z[4:]]
            sol=solve_ivp(augmented,(0,2),np.r_[y,out['initial_tangent']],method='Radau',rtol=1e-11,atol=1e-13,max_step=.002)
            end=sol.y[:,-1];norm=np.linalg.norm(end[4:])
            errors=dict(state=float(np.max(abs(out['final_states'][0]-end[:4]))),tangent=float(np.max(abs(out['final_tangents'][0]-end[4:]/norm))),lle=float(abs(out['lle'][0]-np.log(norm)/2)))
            if p[0]==0:errors['linear_expm']=float(np.max(abs(out['final_states'][0]-expm(field_jac(y,p)[1]*2)@y)))
            passed=bool(sol.success and out['status'][0]==0 and max(errors.values())<1e-7)
            records.append(dict(parameters=p.tolist(),computational_scale=scale,errors=errors,passed=passed));print(records[-1],flush=True)
            if not passed:raise RuntimeError('GPU control failed')
    save(ROOT/'gpu-validation',dict(records=records,passed=True,scope='Independent implicit Radau state/tangent, linear exact control and two computational-time scalings'))


def resume_control():
    gpu=RayleighGPU();p=np.array([[3.,1.5]]);y=np.array([[.1,.7,-.3,.2]])
    settings=dict(transient=1.,duration=5.,rtol=1e-9,max_step=.01,time_scales=[1.])
    partial=gpu.scan(p,y,**settings,max_attempts=100)
    if partial['status'][0]!=3:raise RuntimeError('Control did not reach its forced step budget')
    continued=gpu.scan(p,y,**settings,max_attempts=100000,resume_arrays=partial)
    whole=gpu.scan(p,y,**settings,max_attempts=100000)
    keys=['final_states','final_tangents','raw_clocks','raw_window_log_sums','stats','modal_statistics','step_counts']
    errors={k:float(np.max(abs(continued[k]-whole[k]))) for k in keys}
    passed=bool(continued['status'][0]==0 and whole['status'][0]==0 and max(errors.values())<1e-12)
    save(ROOT/'resume-validation',dict(passed=passed,errors=errors,scope='Forced 100-attempt interruption followed by exact clock/state/tangent/log/statistics continuation versus uninterrupted run'))
    print(dict(passed=passed,errors=errors),flush=True)
    if not passed:raise RuntimeError('Continuation control failed')


def resume_partial(file):
    import json
    old=np.load(file);meta=json.loads(Path(file).with_suffix('.json').read_text());indices=np.flatnonzero(old['status']==3)
    if not len(indices):raise ValueError('No step-budget cases to resume')
    selected={k:old[k][indices] for k in ['final_states','final_tangents','raw_clocks','raw_window_log_sums','stats','modal_statistics','step_counts']}
    s=meta['settings'];out=RayleighGPU().scan(old['parameters'][indices],old['initial_states'][indices],
        transient=s['transient_scaled'],duration=s['duration_scaled'],rtol=s['rtol'],max_step=s['max_step_sigma'],
        time_scales=old['time_scales'][indices],resume_arrays=selected,max_attempts=24000000,max_wall_seconds=900)
    candidate=(out['status']==0)&(out['normalized_lle']>.01)&np.all(out['normalized_lle_windows']>.005,axis=1)
    out.update(source_indices=old['source_indices'][indices],resumed_rows=indices,previous_elapsed_sigma=old['elapsed_sigma'][indices],previous_step_counts=old['step_counts'][indices],
        retained_candidate=candidate,retained=int(candidate.sum()),source=Path(file).name,
        continuation_scope='Same trajectory, tangent, clocks, four raw log sums and modal statistics continued; no transient restart or replacement of previous data')
    save(ROOT/(Path(file).stem+'-resumed'),out)
    print({'lanes':len(indices),'completed':int((out['status']==0).sum()),'retained':out['retained'],'wall_seconds':out['wall_seconds']},flush=True)


def survey(resolution=32):
    mu=np.geomspace(.02,50,resolution);kappa=np.geomspace(.002,50,resolution)
    grid=np.array(np.meshgrid(mu,kappa,indexing='ij')).reshape(2,-1).T
    # Append the published parameter conversion and two independent modal-instability controls.
    grid=np.vstack([grid,SOURCE,[10,.01],[30,.1]])
    p=np.repeat(grid,3,axis=0)
    seeds=np.array([[.1,.7,-.3,.2],[1.7,-.4,.2,.95],[-.7,.9,1.1,-.6]])
    starts=np.tile(seeds,(len(grid),1));out=RayleighGPU().scan(p,starts)
    candidate=(out['status']==0)&(out['normalized_lle']>.02)&np.all(out['normalized_lle_windows']>.01,axis=1)
    out.update(preliminary_candidate=candidate,mu_grid=mu,kappa_grid=kappa,lanes=len(p),completed=int((out['status']==0).sum()),preliminary_positive=int(candidate.sum()),
        start_scope='Three chosen canonical states per tuple, strictly off both exact modal planes; not author-supplied physical initial conditions')
    save(ROOT/'scout',out);print({k:out[k] for k in ['lanes','completed','preliminary_positive','wall_seconds']},flush=True)


def refine(file):
    old=np.load(file);ok=np.flatnonzero(old['status']==0)
    weak=(old['status']==0)&np.all(old['normalized_lle_windows']>0,axis=1)
    indices=np.union1d(np.flatnonzero(weak|old['preliminary_candidate']),ok[np.argsort(old['normalized_lle'][ok])[-24:]])
    out=RayleighGPU().scan(old['parameters'][indices],old['final_states'][indices],transient=200,duration=400,rtol=1e-10,max_step=.025,max_wall_seconds=900,max_attempts=12000000)
    candidate=(out['status']==0)&(out['normalized_lle']>.01)&np.all(out['normalized_lle_windows']>.005,axis=1)
    out.update(source_indices=indices,retained_candidate=candidate,retained=int(candidate.sum()),source=Path(file).name,selection='All four positive normalized windows, union preliminary passes and top24 means')
    save(ROOT/(Path(file).stem+'-refined'),out);print({'lanes':len(indices),'completed':int((out['status']==0).sum()),'retained':out['retained'],'wall_seconds':out['wall_seconds']},flush=True)


def weak_followup():
    scout=np.load(ROOT/'scout.npz');refined=np.load(ROOT/'scout-refined.npz');ok=np.flatnonzero(refined['status']==0)
    row=ok[np.argmax(refined['normalized_lle_windows'][ok].min(axis=1))];source=int(refined['source_indices'][row]);indices=np.arange(3*(source//3),3*(source//3)+3)
    out=RayleighGPU().scan(scout['parameters'][indices],scout['initial_states'][indices],transient=1000,duration=2000,rtol=1e-10,max_step=.025,max_attempts=24000000,max_wall_seconds=900)
    candidate=(out['status']==0)&(out['normalized_lle']>.001)&np.all(out['normalized_lle_windows']>.0005,axis=1)
    out.update(source_indices=indices,selected_refined_row=int(row),selected_source_row=source,retained_candidate=candidate,retained=int(candidate.sum()),
        selection='Fresh three original starts at the completed refined tuple with largest minimum normalized window; not a parameter search',
        retained_thresholds=dict(normalized_mean=.001,each_normalized_window=.0005))
    save(ROOT/'weak-coupling-long',out);print({'source':source,'completed':int((out['status']==0).sum()),'retained':out['retained'],'wall_seconds':out['wall_seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true');parser.add_argument('--resume-control',action='store_true');parser.add_argument('--resume');parser.add_argument('--weak-followup',action='store_true');parser.add_argument('--resolution',type=int,default=32);parser.add_argument('--refine');a=parser.parse_args()
    if a.controls:controls()
    elif a.resume_control:resume_control()
    elif a.resume:resume_partial(a.resume)
    elif a.weak_followup:weak_followup()
    elif a.refine:refine(a.refine)
    else:survey(a.resolution)
