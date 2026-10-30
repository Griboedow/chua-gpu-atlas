"""FP64 OpenCL surveys of the exact five-/nine-state NIFTE models."""
import argparse
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from nifte import ROOT,COLUMNS,parameters,matrix,rhs_jac,save


class NifteGPU:
    def __init__(self,loss=False):
        devices=[d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        self.loss=loss;self.n=9 if loss else 5
        source=Path(__file__).with_name('nifte.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,f'#define N {self.n}\n'+source).build(),'advance')

    def scan(self,p,starts,transient=500.,duration=1000.,rtol=1e-8,max_step=.02,max_wall_seconds=600.,max_attempts=4000000):
        p=np.ascontiguousarray(p,dtype=float);starts=np.asarray(starts,dtype=float);lanes=len(p);n=self.n
        if p.shape!=(lanes,23) or starts.shape!=(lanes,n) or lanes<1 or not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Invalid input arrays')
        if np.any(p[:,1:]<=0) or np.any(p[:,0]<0) or transient<0 or min(duration,rtol,max_step,max_wall_seconds)<=0:raise ValueError('Invalid settings')
        states=np.zeros((lanes,2*n));states[:,:n]=starts
        tangent=np.linspace(.2,1.,n);tangent/=np.linalg.norm(tangent);states[:,n:]=tangent
        mats=np.array([matrix(row,self.loss)[0] for row in p]);pars=np.column_stack([p[:,0],p[:,1],1/(p[:,4]*p[:,9])])
        clocks=np.zeros((lanes,4));clocks[:,1]=.001;logs=np.zeros((lanes,4));stats=np.zeros((lanes,2*n));stats[:,::2]=np.inf;stats[:,1::2]=-np.inf
        status=np.zeros(lanes,np.int32);counts=np.zeros((lanes,2),np.int32)
        arrays=[states,mats,pars,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=np.ascontiguousarray(a)) for a in arrays]
        wall=time.perf_counter();last=wall;calls=0
        while True:
            self.kernel(self.queue,(lanes,),None,*buffers,np.int32(lanes),np.int32(128),np.float64(transient),np.float64(duration),np.float64(rtol),np.float64(rtol*.001),np.float64(max_step),np.int32(max_attempts)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[3]).wait();cl.enqueue_copy(self.queue,status,buffers[6]).wait()
                active=(status==0)&(clocks[:,0]<transient+duration)
                if time.perf_counter()-last>25:print(f'{n}D: {active.sum()} active, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-wall>max_wall_seconds:status[active]=4;cl.enqueue_copy(self.queue,buffers[6],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,initial_states=starts,initial_tangent=tangent,final_states=states[:,:n].copy(),final_tangents=states[:,n:].copy(),lle=windows.mean(axis=1),lle_windows=windows,status=status,stats=stats,elapsed=clocks[:,0],step_counts=counts,
            settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol*.001,max_step=max_step,max_attempts=max_attempts),
            method='Adaptive DP5(4), analytic variational equations, physical seconds',device=self.device.name,wall_seconds=time.perf_counter()-wall,parameter_columns=COLUMNS,
            model='NTP-TL 9D' if self.loss else 'NTP 5D',state_columns=['p','d','q','u','v']+(['z1','z2','z3','z4'] if self.loss else []),
            stats_layout='Interleaved minima/maxima over observation only',physical_scope='No vessel/interface calibration imposed; mathematical equations on R^n',
            status_codes={'0':'completed','1':'numerical magnitude cutoff 1e6','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved'})


def controls():
    records=[]
    for loss in [False,True]:
        gpu=NifteGPU(loss);n=gpu.n
        for gain in [0.,452.1,302.871676391795]:
            p=parameters(True)
            if gain==0:p[0]=0
            else:p[1]=gain/p[0]
            if gain<400 and gain>0:p[6]=286255.183178673
            y=np.zeros(n);y[:5]=[.01,.005,-.003,.02,-.01]
            out=gpu.scan(p[None],y[None],transient=0,duration=2,rtol=1e-10,max_step=.005)
            lin=matrix(p,loss)
            def augmented(t,z):
                f,j=rhs_jac(z[:n],p,loss,lin);return np.r_[f,j@z[n:]]
            sol=solve_ivp(augmented,(0,2),np.r_[y,out['initial_tangent']],method='Radau',rtol=1e-11,atol=1e-13,max_step=.002)
            end=sol.y[:,-1];norm=np.linalg.norm(end[n:])
            errors=dict(state=float(np.max(abs(out['final_states'][0]-end[:n]))),tangent=float(np.max(abs(out['final_tangents'][0]-end[n:]/norm))),lle=float(abs(out['lle'][0]-np.log(norm)/2)))
            if gain==0:errors['linear_expm']=float(np.max(abs(out['final_states'][0]-expm(lin[0]*2)@y)))
            passed=bool(sol.success and out['status'][0]==0 and max(errors.values())<1e-7)
            records.append(dict(dimension=n,gain=gain,errors=errors,passed=passed));print(records[-1],flush=True)
            if not passed:raise RuntimeError('GPU control failed')
    save(ROOT/'gpu-validation',dict(records=records,passed=True,scope='Independent Radau state/tangent, zero-input matrix exponential, both published source cycles and double-Hopf linear parameters'))


def survey(loss=False,resolution=32,focus=False,random_tuples=0):
    gpu=NifteGPU(loss)
    if random_tuples:
        rng=np.random.default_rng(2026103008);p=np.tile(parameters(),(random_tuples,1))
        for column,lo,hi in [(2,.1,10),(3,.1,10),(4,.2,5),(6,.02,4),(8,.25,4),(9,.5,2),(10,.5,2),(11,.5,2)]:
            p[:,column]*=np.exp(rng.uniform(np.log(lo),np.log(hi),random_tuples))
        p[:,1]=np.exp(rng.uniform(np.log(50),np.log(5000),random_tuples))/p[:,0]
        p[0]=parameters(True);lf=np.array([]);g=np.array([])
    elif focus:
        lf=np.linspace(2.0e5,4.4e5,resolution);g=np.linspace(230.,650.,resolution)
    else:
        lf=np.geomspace(2.0e5,1.5e7,resolution);g=np.geomspace(50.,4000.,resolution)
    if not random_tuples:
        grid=np.array(np.meshgrid(lf,g,indexing='ij')).reshape(2,-1).T
        p=np.tile(parameters(),(len(grid),1));p[:,6]=grid[:,0];p[:,1]=grid[:,1]/p[:,0]
    p=np.repeat(p,2,axis=0)
    y=np.zeros((len(p),gpu.n));y[::2,:5]=[.01,.005,-.003,.02,-.01];y[1::2,:5]=[-.2,.02,.01,.3,-.2]
    out=gpu.scan(p,y)
    positive=(out['status']==0)&(out['lle']>.01)&np.all(out['lle_windows']>.005,axis=1)
    out.update(preliminary_candidate=positive,lanes=len(p),preliminary_positive=int(positive.sum()),completed=int((out['status']==0).sum()),
        Lf_grid=lf,gain_grid=g,seed=2026103008 if random_tuples else None,
        scope='Fixed F=1.37; K=G/F. Native model components varied as mathematical parameters; no geometric design calibration.')
    name=f'{"random" if random_tuples else "focus" if focus else "scout"}-{gpu.n}d';save(ROOT/name,out)
    print({k:out[k] for k in ['lanes','completed','preliminary_positive','wall_seconds']},flush=True)


def refine(file):
    old=np.load(file);loss=old['initial_states'].shape[1]==9;gpu=NifteGPU(loss)
    ok=np.flatnonzero(old['status']==0)
    weak=(old['status']==0)&np.all(old['lle_windows']>0,axis=1)
    indices=np.union1d(np.flatnonzero(weak|old['preliminary_candidate']),ok[np.argsort(old['lle'][ok])[-24:]])
    out=gpu.scan(old['parameters'][indices],old['final_states'][indices],transient=2000,duration=8000,rtol=1e-10,max_step=.01,max_wall_seconds=900)
    positive=(out['status']==0)&(out['lle']>.001)&np.all(out['lle_windows']>.0005,axis=1)
    out.update(source_indices=indices,retained_candidate=positive,source=Path(file).name,retained=int(positive.sum()),lanes=len(indices),
        selection='All four scout windows positive, or preliminary threshold pass, union top 24 means',
        retained_thresholds=dict(mean=.001,each_window=.0005))
    save(ROOT/(Path(file).stem+'-refined'),out);print({'lanes':len(indices),'retained':out['retained'],'wall_seconds':out['wall_seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true');parser.add_argument('--loss',action='store_true');parser.add_argument('--resolution',type=int,default=32);parser.add_argument('--focus',action='store_true');parser.add_argument('--random-tuples',type=int,default=0);parser.add_argument('--refine');a=parser.parse_args()
    if a.controls:controls()
    elif a.refine:refine(a.refine)
    else:survey(a.loss,a.resolution,a.focus,a.random_tuples)
