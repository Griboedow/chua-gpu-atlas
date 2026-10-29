"""Phenomenological optical-feedback extension inspired by Dinavahi-Harris.

This is NOT a fitted reconstruction of their device. A finite photoresponse
time is introduced explicitly to investigate the missing feedback mechanism.
State: height x>0, velocity v, illumination-memory z in [0,1].
Parameters: a,b,tau,k,width,chosen_equilibrium_height,unused,unused.
"""
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from chua_gpu import GPU

def optical(x,width):
    s=np.exp(-((x-1)/width)**2)
    return s,-2*(x-1)/width**2*s

def permeability(x):
    sigmoid=1/(1+np.exp(np.clip(-8.75*(x-6/7),-700,700)))
    return 8-7*sigmoid,-61.25*sigmoid*(1-sigmoid)

def rhs(t,state,p):
    a,b,tau,k,width,_,_,_=p;x,v,z=state[:3]
    mu,mux=permeability(x);light,lx=optical(x,width);j=1/(1.04+k*z)
    f=[v,a*mu*j/x**4-1-b*v,(j*light-z)/tau]
    if len(state)==6:
        ux,uv,uz=state[3:];fx=a*j*(mux/x**4-4*mu/x**5);fz=-a*mu*k*j*j/x**4
        qx=j*lx/tau;qz=(-k*j*j*light-1)/tau
        f += [uv,fx*ux-b*uv+fz*uz,qx*ux+qz*uz]
    return np.array(f)

def parameter_at_equilibrium(x,b,tau,k,width):
    light,_=optical(x,width);z=2*light/(1.04+np.sqrt(1.04**2+4*k*light))
    mu,_=permeability(x);a=x**4*(1.04+k*z)/mu
    return np.array([a,b,tau,k,width,x,0,0]),np.array([x,0,z])

def stability(p):
    _,s=parameter_at_equilibrium(p[5],p[1],p[2],p[3],p[4])
    jac=np.column_stack([rhs(0,np.r_[s,e],p)[3:] for e in np.eye(3)])
    return np.linalg.eigvals(jac),jac

class PhotoFeedbackGPU(GPU):
    def __init__(self):
        super().__init__('','float64')
        source=Path(__file__).with_suffix('.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(),'advance')
        self.info.update(model='Phenomenological finite-photoresponse extension; not a fitted hardware model',integrator='Dormand-Prince 5(4)')

    def scan(self,parameters,starts,transient=200,duration=400,rtol=1e-7,atol=1e-10,progress=None):
        p=np.ascontiguousarray(parameters,dtype=float);starts=np.asarray(starts,float);n=len(p)
        if p.shape!=(n,8) or starts.shape!=(n,3):raise ValueError('Shapes must be (N,8),(N,3)')
        if (p[:,:3]<=0).any() or (p[:,3]<0).any() or (p[:,4]<=0).any():raise ValueError('Invalid parameters')
        state=np.zeros((n,8));state[:,:3]=starts;state[:,3:6]=np.array([1,2,3])/np.sqrt(14)
        clocks=np.zeros((n,4));clocks[:,1]=np.minimum(.001,.1*p[:,2])
        logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,[0,6]]=np.inf;stats[:,[1,7]]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[state,p,clocks,logs,stats,status,counts];buffers=[self.buffer(a) for a in arrays]
        wall=time.perf_counter();gpu=0;target=0;calls=0;last=wall
        for window,length in [(-1,transient)]+[(i,duration/4) for i in range(4)]:
            target+=length
            if length==0:continue
            while True:
                ev=self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(128),np.float64(target),np.int32(window),np.float64(rtol),np.float64(atol));ev.wait()
                gpu+=(ev.profile.end-ev.profile.start)*1e-9;calls+=1
                if calls%8==0:
                    cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                    active=status==0
                    if progress and time.perf_counter()-last>10:
                        progress(float(np.minimum(np.where(active,clocks[:,0],target),target).mean()/(transient+duration)),time.perf_counter()-wall);last=time.perf_counter()
                    if not np.any(active&(clocks[:,0]<target)):break
                    if calls>1000000:raise RuntimeError('Submission limit')
        for a,bf in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,bf).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,starts=starts,final_states=state[:,:3],lle_windows=windows,lle=windows.mean(axis=1),stats=stats,status=status,step_counts=counts,
                    settings=dict(transient=transient,duration=duration,rtol=rtol,atol=atol),backend=self.info,timing=dict(wall_seconds=time.perf_counter()-wall,gpu_seconds=gpu))
