"""Independent implementation of Potter et al., PRA 110, 042604 (2024), Eq. 10.

States: X,Y,Z,V. Parameters: Vn,Rn,Gamma,Omega,alpha,ZT,x0,unused.
Geometry l=0.5 and lambda=0.13 is fixed. No stochastic forcing or reset.
The model is a nonlinear semiclassical feedback ODE, not a simulation of a
many-body quantum wavefunction. All times and exponents are dimensionless.
"""
from pathlib import Path
import time
import numpy as np
from scipy.special import roots_hermite
from scipy.interpolate import CubicHermiteSpline
import pyopencl as cl

ROOT = Path(__file__).resolve().parent
LENGTH, TUNNEL = .5, .13
PARAMETERS = ['Vn','Rn','Gamma','Omega','alpha','ZT','x0','unused']


def conductance(x, order=512):
    """F00,F01,F11 and their derivatives with respect to xV; Gauss-Hermite."""
    if not hasattr(conductance, 'cache'): conductance.cache = {}
    if order not in conductance.cache:
        u,w=roots_hermite(order)
        conductance.cache[order]=(u,np.array([w,w*np.sqrt(2)*u,w*2*u*u])/np.sqrt(np.pi))
    u,weights=conductance.cache[order]
    a=(np.atleast_1d(x)[...,None]+LENGTH*u)/TUNNEL
    # Stable sech; no clipping of the physical state.
    e=np.exp(-np.abs(a)); sech=2*e/(1+e*e)
    f=sech @ weights.T
    derivative=(-sech*np.tanh(a)/TUNNEL) @ weights.T
    return f,derivative


def rhs(t, s, p, order=512):
    vn,r,gamma,omega,alpha,zt,x0,_=p
    x,y,z,v=s[:4]
    f,df=conductance(x0-np.sqrt(2)*LENGTH*v,order)
    a,b,c=f[0]; av,bv,cv=df[0]*(-np.sqrt(2)*LENGTH)
    g=.5*((1+z)*a+2*x*b+(1-z)*c)
    gv=.5*((1+z)*av+2*x*bv+(1-z)*cv)
    dv=vn-(1+r*g)*v
    flow=[omega*y+2*dv*z-gamma*x, -omega*x-gamma*y,
          -2*dv*x-alpha*gamma*(z-zt),dv]
    if len(s)==8:
        ux,uy,uz,uv=s[4:]
        duv=-r*v*b*ux-r*v*.5*(a-c)*uz+(-1-r*g-r*v*gv)*uv
        flow += [omega*uy+2*duv*z+2*dv*uz-gamma*ux,
                 -omega*ux-gamma*uy,-2*duv*x-2*dv*ux-alpha*gamma*uz,duv]
    return np.asarray(flow)


class QuantumGPU:
    def __init__(self, spacing=.002, order=512):
        ds=[d for platform in cl.get_platforms() for d in platform.get_devices()
            if d.type & cl.device_type.GPU and d.double_fp_config]
        if not ds: raise RuntimeError('An OpenCL GPU with FP64 is required.')
        self.device=ds[0]; self.context=cl.Context([self.device])
        self.queue=cl.CommandQueue(self.context,properties=cl.command_queue_properties.PROFILING_ENABLE)
        self.grid=np.linspace(-32,32,round(64/spacing)+1)
        spacing=self.grid[1]-self.grid[0]
        blocks=[conductance(block,order) for block in np.array_split(self.grid,64)]
        values=np.concatenate([b[0] for b in blocks]); deriv=np.concatenate([b[1] for b in blocks])
        table=np.zeros((len(self.grid),8)); table[:,:3]=values;table[:,4:7]=deriv
        self.spline=CubicHermiteSpline(self.grid,values,deriv)
        flags=cl.mem_flags
        self.table=cl.Buffer(self.context,flags.READ_ONLY|flags.COPY_HOST_PTR,hostbuf=table)
        source=(ROOT/'quantum_neuron.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(options=[
            f'-DTABLE_N={len(table)}',f'-DTABLE_DX={spacing:.17g}']), 'advance')
        self.info={'device':self.device.name.strip(),'driver':self.device.driver_version,
                   'precision':'float64','integrator':'Dormand-Prince 5(4)',
                   'quadrature_order':order,'table_spacing':spacing,'table_extent':[-32,32],
                   'fast_math':False}

    def scan(self, parameters, starts, transient=200., duration=400.,rtol=1e-7,atol=1e-10,progress=None):
        p=np.ascontiguousarray(parameters,dtype=np.float64)
        starts=np.asarray(starts,dtype=np.float64)
        if p.ndim!=2 or p.shape[1]!=8 or starts.shape!=(len(p),4):raise ValueError('Expected (N,8) parameters and (N,4) starts')
        if not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Nonfinite input')
        if (p[:,:5]<=0).any() or (p[:,4]>2).any() or (p[:,5]<0).any() or (p[:,5]>1).any():raise ValueError('Outside physical parameter constraints')
        if (np.sum(starts[:,:3]**2,axis=1)>1+1e-12).any() or (starts[:,3]<0).any() or (starts[:,3]>p[:,0]).any():raise ValueError('Unphysical initial condition')
        if transient<0 or min(duration,rtol,atol)<=0:raise ValueError('Invalid solver settings')
        n=len(p);state=np.zeros((n,8));state[:,:4]=starts;state[:,4:]=np.array([1,2,3,4])/np.sqrt(30)
        clocks=np.zeros((n,4));clocks[:,1]=.001
        logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,[0,6]]=np.inf;stats[:,[1,7]]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[state,p,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();gpu_seconds=0.;calls=0;target=0.;last_print=wall
        for window,length in [(-1,transient)]+[(i,duration/4) for i in range(4)]:
            target+=length
            if length==0:continue
            while True:
                event=self.kernel(self.queue,(n,),None,*buffers,self.table,np.int32(n),np.int32(128),
                                  np.float64(target),np.int32(window),np.float64(rtol),np.float64(atol))
                event.wait();gpu_seconds+=(event.profile.end-event.profile.start)*1e-9;calls+=1
                if calls%8==0:
                    cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                    active=status==0
                    if progress and time.perf_counter()-last_print>10:
                        progress(float(np.minimum(np.where(active,clocks[:,0],target),target).mean()/(transient+duration)),time.perf_counter()-wall)
                        last_print=time.perf_counter()
                    if not np.any(active&(clocks[:,0]<target)):break
                    if calls>1000000:raise RuntimeError('Submission limit exceeded')
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return {'parameters':p,'starts':starts,'final_states':state[:,:4], 'lle_windows':windows,
                'lle':windows.mean(axis=1),'status':status,'stats':stats,'step_counts':counts,
                'settings':{'transient':transient,'duration':duration,'rtol':rtol,'atol':atol},
                'timing':{'wall_seconds':time.perf_counter()-wall,'gpu_seconds':gpu_seconds},'backend':self.info}
