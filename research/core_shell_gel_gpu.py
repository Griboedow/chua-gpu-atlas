"""FP64 adaptive GPU survey of the unreduced four-state core-shell gel."""
import argparse
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from core_shell_gel import ROOT,COLUMNS,parameters,start,save


class GelGPU:
    def __init__(self):
        devices=[d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        source=Path(__file__).with_name('core_shell_gel.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(),'advance')

    def scan(self,p,starts,transient=500.,duration=1000.,rtol=1e-8,max_step=.1,max_wall_seconds=600.,max_attempts=3000000):
        p=np.ascontiguousarray(p,dtype=float);starts=np.asarray(starts,dtype=float);n=len(p)
        if p.shape!=(n,8) or starts.shape!=(n,4) or n<1 or not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Invalid input arrays')
        if np.min(p)<=0 or np.any(p[:,7]<=1) or transient<0 or min(duration,rtol,max_step,max_wall_seconds)<=0:raise ValueError('Invalid settings')
        states=np.zeros((n,8));states[:,:4]=starts;states[:,4:]=np.array([1.,.3,-.2,.4])/np.sqrt(1.29)
        clocks=np.zeros((n,4));clocks[:,1]=.001;logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,::2]=np.inf;stats[:,1::2]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[states,p,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();last=wall;calls=0
        while True:
            self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(256),np.float64(transient),np.float64(duration),np.float64(rtol),np.float64(rtol*.001),np.float64(max_step),np.int32(max_attempts)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<transient+duration)
                if time.perf_counter()-last>20:print(f'{active.sum()} active, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-wall>max_wall_seconds:status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,initial_states=starts,final_states=states[:,:4].copy(),lle=windows.mean(axis=1),lle_windows=windows,status=status,stats=stats,elapsed=clocks[:,0],step_counts=counts,
            settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol*.001,max_step=max_step,max_attempts=max_attempts),
            method='Adaptive DP5(4), analytic variational equations, constant scaled concentrations',device=self.device.name,wall_seconds=time.perf_counter()-wall,parameter_columns=COLUMNS,state_columns=['G','U','x','phi'],
            stats_columns=['min_G','max_G','min_U','max_U','min_x','max_x','min_phi','max_phi'],status_codes={'0':'completed','1':'physical domain exit','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved'})


def survey(tuples=2048,broad=False):
    rng=np.random.default_rng(2026103003 if broad else 2026103002)
    p=np.tile(parameters(),(tuples,1))
    def loguniform(lo,hi):return np.exp(rng.uniform(np.log(lo),np.log(hi),tuples))
    p[:,0]=loguniform(1.,12.);p[:,1]=loguniform(30.,600.);p[:,4]=rng.uniform(3.,14.,tuples)
    p[:,5]=loguniform(.03,3.);p[:,6]=loguniform(.03,100.)
    if broad:p[:,2]=loguniform(.03,50.);p[:,3]=loguniform(.1,50.)
    p[0]=parameters();p=np.repeat(p,2,axis=0)
    y=np.array([start(param,x,phi) for i,param in enumerate(p) for x,phi in [((1.,.15) if i%2==0 else (.9,.5))]])
    out=GelGPU().scan(p,y)
    positive=(out['status']==0)&(out['lle']>.002)&np.all(out['lle_windows']>.001,axis=1)
    out.update(preliminary_candidate=positive,lanes=len(p),preliminary_positive=int(positive.sum()),completed=int((out['status']==0).sum()),seed=2026103003 if broad else 2026103002,
        scope='Published full equations; finite substrate rates scanned beyond the quasi-steady reduction, without experimental calibration')
    save(ROOT/('broad-scout' if broad else 'scout'),out)
    print({k:out[k] for k in ['lanes','completed','preliminary_positive','wall_seconds']},flush=True)


def refine(file):
    old=np.load(file);ok=np.flatnonzero(old['status']==0)
    indices=np.union1d(np.flatnonzero(old['preliminary_candidate']),ok[np.argsort(old['lle'][ok])[-24:]])
    out=GelGPU().scan(old['parameters'][indices],old['final_states'][indices],transient=2000,duration=4000,rtol=1e-10,max_step=.05)
    positive=(out['status']==0)&(out['lle']>.002)&np.all(out['lle_windows']>.001,axis=1)
    out.update(source_indices=indices,retained_candidate=positive,source=str(Path(file).name),retained=int(positive.sum()),lanes=len(indices))
    save(ROOT/(Path(file).stem+'-refined'),out)
    print({'lanes':len(indices),'retained':out['retained'],'wall_seconds':out['wall_seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--tuples',type=int,default=2048);parser.add_argument('--broad',action='store_true');parser.add_argument('--refine');args=parser.parse_args()
    if args.refine:refine(args.refine)
    else:survey(args.tuples,args.broad)
