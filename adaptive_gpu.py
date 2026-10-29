"""Adaptive FP64 GPU scanner for the much wider, faster Chua regime.

Dormand-Prince 5(4) controls all six state/tangent components. Unlike an
event-resolved reference, this solver does not explicitly locate diode kinks;
independent checks and tolerance refinement remain necessary.
"""
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from chua_gpu import GPU, DEFAULT_SEEDS


class AdaptiveGPU(GPU):
    def __init__(self, device_filter=""):
        super().__init__(device_filter, "float64")
        root=Path(__file__).resolve().parent
        source=(root/'kernels.cl').read_text()+'\n'+(root/'adaptive_kernels.cl').read_text()
        program=cl.Program(self.context,source).build(options=['-DUSE_DOUBLE'])
        self.adaptive=cl.Kernel(program,'adaptive_advance')
        self.info={**self.info,'integrator':'Adaptive Dormand-Prince 5(4)'}

    def scan(self, parameters, initial_states=None, transient=100., duration=200.,
             rtol=1e-6, atol=1e-9, escape=1e6, threshold=.01, progress=None,
             attempts=128):
        p=np.asarray(parameters,dtype=np.float64)
        if p.ndim!=2 or p.shape[1]!=4 or not len(p) or not np.isfinite(p).all() or (p[:,:2]<=0).any():
            raise ValueError('Expected finite (N,4) parameters with positive alpha and beta')
        if initial_states is None: initial_states=np.tile(DEFAULT_SEEDS,(len(p),1,1))
        starts=np.asarray(initial_states,dtype=np.float64)
        if starts.ndim!=3 or starts.shape[0]!=len(p) or starts.shape[1]<1 or starts.shape[2]!=3 or not np.isfinite(starts).all():
            raise ValueError('Initial states must have shape (N,S,3)')
        if not np.isfinite([transient,duration,rtol,atol,escape,threshold]).all() or transient<0 or min(duration,rtol,atol,escape,threshold)<=0 or attempts<1:
            raise ValueError('Invalid integration settings')
        ns=starts.shape[1];n=len(p)*ns;params=np.ascontiguousarray(np.repeat(p,ns,axis=0))
        state=np.zeros((n,8));state[:,:3]=starts.reshape(n,3);state[:,3:6]=np.array([1,2,3])/np.sqrt(14)
        clocks=np.zeros((n,4));clocks[:,1]=np.minimum(.01,.1/np.maximum(params[:,0]*abs(1+params[:,2]),np.sqrt(params[:,1])))
        logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,[0,6]]=np.inf;stats[:,[1,7]]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[state,params,clocks,logs,stats,status,counts]
        buffers=[self.buffer(a,readonly=i==1) for i,a in enumerate(arrays)]
        wall=time.perf_counter();gpu_seconds=0.;max_kernel=0.;target=0.;calls=0
        for window,length in [(-1,transient)]+[(i,duration/4) for i in range(4)]:
            target+=length
            if length==0:continue
            while True:
                event=self.adaptive(self.queue,(n,),None,*buffers,np.int32(n),np.int32(attempts),
                    np.float64(target),np.int32(window),np.float64(rtol),np.float64(atol),np.float64(escape))
                event.wait();elapsed=(event.profile.end-event.profile.start)*1e-9
                gpu_seconds+=elapsed;max_kernel=max(max_kernel,elapsed);calls+=1
                if calls%8==0:
                    cl.enqueue_copy(self.queue,clocks,buffers[2]).wait()
                    cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                    active=status==0
                    if progress:progress(float(np.minimum(np.where(active,clocks[:,0],target),target).mean()/(transient+duration)),time.perf_counter()-wall)
                    if not np.any(active&(clocks[:,0]<target)):break
                    if calls>2000000:raise RuntimeError('Adaptive scan exceeded submission budget')
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);lle=windows.mean(axis=1);valid=status==0
        windows[~valid]=np.nan;lle[~valid]=np.nan
        late_range=stats[:,7]-stats[:,6];classes=np.full(n,3,np.uint8)
        classes[valid&(lle<=threshold)&(windows[:,-1]<=threshold)]=2
        classes[valid&(late_range<1e-5)&(lle<0)]=1
        classes[valid&(lle>threshold)&(windows.min(axis=1)>threshold/2)&(late_range>1e-4)]=4
        classes[status==1]=0;classes[status==2]=5
        def shaped(a):return a.reshape((len(p),ns)+a.shape[1:])
        return {**{k:shaped(a) for k,a in [('classes',classes),('lle',lle),('lle_windows',windows),('status',status),
                   ('stats',stats),('final_states',state[:,:3]),('step_counts',counts)]},
            'timing':{'wall_seconds':time.perf_counter()-wall,'gpu_seconds':gpu_seconds,'max_kernel_seconds':max_kernel},
            'settings':{'transient':transient,'duration':duration,'rtol':rtol,'atol':atol,'escape':escape,'threshold':threshold,
                        'renormalization':'Every accepted step','max_step':'min(0.1, 0.5/sqrt(beta))',
                        'error_norm':'Maximum scaled local error across all 6 state/tangent coordinates',
                        'stats_note':'x sums use adaptive step endpoints, not equally spaced times; do not interpret as time averages.'}}
