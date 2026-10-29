"""FP64 adaptive GPU survey of the published thermal sheet ODE family."""
from pathlib import Path
import sys
import time

import numpy as np
import pyopencl as cl

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from chua_gpu import GPU


class ThermalSheetGPU(GPU):
    def __init__(self):
        super().__init__('', 'float64')
        source=Path(__file__).with_name('thermal_sheet.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(),'advance')
        self.info.update(model='Shen 2025 sheet equations with the reported cubic velocity damping',
                         integrator='Dormand-Prince 5(4), analytic tangent, per-accepted-step renormalization')

    def scan(self,parameters,starts,transient=300,duration=600,rtol=1e-7,atol=1e-10,
             step_factor=1.,progress=None,max_attempts=3000000,max_wall_seconds=600):
        parameters=np.asarray(parameters,float);starts=np.asarray(starts,float);n=len(parameters)
        if parameters.shape!=(n,5) or starts.shape!=(n,3):raise ValueError('Expected (N,5), (N,3)')
        if not np.isfinite(parameters).all() or not np.isfinite(starts).all():raise ValueError('Nonfinite inputs')
        if (parameters[:,:3]<=0).any() or (parameters[:,3:]<0).any():raise ValueError('Invalid parameters')
        if transient<0 or min(duration,rtol,atol,step_factor)<=0:raise ValueError('Invalid integration settings')
        if n==0 or not 0<max_attempts<2000000000 or not np.isfinite(max_wall_seconds) or max_wall_seconds<=0:
            raise ValueError('Invalid sample count or computational budget')
        p=np.zeros((n,8));p[:,:5]=parameters
        state=np.zeros((n,8));state[:,:3]=starts;state[:,3:6]=np.array([1,2,3])/np.sqrt(14)
        clocks=np.zeros((n,4));clocks[:,1]=np.minimum(.001,.05/np.sqrt(p[:,0]+p[:,1]))
        logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,[0,2,5,6]]=np.inf;stats[:,[1,3,7]]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32);domain=np.zeros((n,2),np.int32)
        arrays=[state,p,clocks,logs,stats,status,counts,domain];buffers=[self.buffer(a) for a in arrays]
        wall=time.perf_counter();gpu=0.;calls=0;last=wall;timed_out=False
        while True:
            ev=self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(128),
                           np.float64(transient),np.float64(duration),np.float64(rtol),np.float64(atol),
                           np.float64(step_factor),np.int32(max_attempts))
            ev.wait();gpu+=(ev.profile.end-ev.profile.start)*1e-9;calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<transient+duration)
                if progress and time.perf_counter()-last>10:
                    progress(float(np.minimum(clocks[:,0],transient+duration).mean()/(transient+duration)),time.perf_counter()-wall)
                    last=time.perf_counter()
                if not np.any(active):break
                if time.perf_counter()-wall>=max_wall_seconds:
                    timed_out=True;break
        for array,buffer in zip(arrays,buffers):cl.enqueue_copy(self.queue,array,buffer).wait()
        if timed_out:status[(status==0)&(clocks[:,0]<transient+duration)]=4
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=parameters,starts=starts,final_states=state[:,:3],lle_windows=windows,
                    lle=windows.mean(axis=1),stats=stats,status=status,step_counts=counts,domain_flags=domain,
                    elapsed=clocks[:,0],settings=dict(transient=transient,duration=duration,rtol=rtol,atol=atol,step_factor=step_factor,
                                                      max_attempts=max_attempts,max_wall_seconds=max_wall_seconds),
                    backend=self.info,timing=dict(wall_seconds=time.perf_counter()-wall,gpu_seconds=gpu),
                    parameter_columns=['K','G','B','S','D'],
                    stats_columns=['min_X','max_X','min_Y','max_Y','max_abs_V','min_stiffness_ratio','late_min_X','late_max_X'],
                    status_legend={'0':'completed','1':'state magnitude exceeded 1e4','2':'numerical failure or step underflow',
                                   '3':'per-trajectory attempt budget exhausted; unresolved','4':'batch wall-time budget exhausted; unresolved'},
                    domain_legend={'columns':['transient','observation'],'bits':{'1':'X<0','2':'Y<0','4':'Y>1','8':'1-SY<0'},'tolerance':1e-8})
