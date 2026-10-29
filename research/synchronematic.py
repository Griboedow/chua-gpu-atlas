"""Independent deterministic implementation of Leyva et al. (2026), Eqs. 10-12.

The two noise strengths are set to zero. Geometry is fixed. We remove the
uniform phase rotation by projecting both phase velocities and phase tangents
onto the zero-sum subspace. Coefficients are the paper's Fourier approximation,
not a newly fitted physical model. Authors' software is not copied here.
"""
from pathlib import Path
import time
import numpy as np
import pyopencl as cl

COEFFICIENTS=np.array([[2.921,2.070,-.781,.312],[-5.351,.223,3.009,-1.359]])
# Mathematical coefficients implied by the NumPy branch at author commit
# 016614ad6debe11e0353f0a01b1fad937fc77461. See source-audit.json; sign differs.
REPOSITORY_NUMPY_COEFFICIENTS=np.array([[2.92076,2.06954,-.78064,.312294],[-5.35098,.222,3.00884,1.359334]])

def gradient_coefficients(coefficients=COEFFICIENTS):
    """Gradient part in the Euclidean angular metric; residual has zero divergence."""
    c=np.asarray(coefficients);result=np.zeros((2,4))
    for harmonic,even,odd in [(1,0,1),(3,2,3)]:
        potential=(c[1,even]-harmonic*c[0,odd])/(harmonic**2+1)
        result[0,odd]=-harmonic*potential
        result[1,even]=potential
    return result

def potential(state,pairs,coefficients=COEFFICIENTS):
    """Scalar potential of the gradient component, summed over unordered pairs."""
    n=len(state)//2;theta,chi=state[:n],state[n:];g=gradient_coefficients(coefficients);value=0.
    for i in range(n):
        for j in range(i):
            w,co,si,_=pairs[i,j];ca=np.cos(chi[i])*co+np.sin(chi[i])*si
            cb=np.cos(chi[j])*co+np.sin(chi[j])*si;delta=theta[j]-theta[i]
            value+=w*(g[1,0]*np.cos(delta)+g[1,2]*np.cos(3*delta))*ca*cb
    return value

def geometry(positions):
    p=np.asarray(positions,float);delta=p[None,:,:]-p[:,None,:]
    distance=np.linalg.norm(delta,axis=2);np.fill_diagonal(distance,np.inf)
    result=np.zeros((len(p),len(p),4));result[:,:,0]=distance**-3
    result[:,:,1]=np.divide(delta[:,:,0],distance);result[:,:,2]=np.divide(delta[:,:,1],distance)
    return result

def flow_and_jac(state,pairs,coefficients=COEFFICIENTS):
    n=len(state)//2;theta,chi=state[:n],state[n:]
    delta=theta[None,:]-theta[:,None];co=np.cos(delta);si=np.sin(delta)
    co3=4*co**3-3*co;si3=3*si-4*si**3
    ca=np.cos(chi)[:,None]*pairs[:,:,1]+np.sin(chi)[:,None]*pairs[:,:,2]
    sa=np.sin(chi)[:,None]*pairs[:,:,1]-np.cos(chi)[:,None]*pairs[:,:,2]
    cb=np.cos(chi)[None,:]*pairs[:,:,1]+np.sin(chi)[None,:]*pairs[:,:,2]
    sb=np.sin(chi)[None,:]*pairs[:,:,1]-np.cos(chi)[None,:]*pairs[:,:,2]
    functions=[];derivatives=[]
    for a,b,c,d in coefficients:
        functions.append(a*co+b*si+c*co3+d*si3)
        derivatives.append(-a*si+b*co-3*c*si3+3*d*co3)
    ft,fc=functions;dt,dc=derivatives;w=pairs[:,:,0]
    flow=np.r_[(w*ft*ca*cb).sum(axis=1),(w*fc*sa*cb).sum(axis=1)]
    jac=np.zeros((2*n,2*n))
    for i in range(n):
        for j in range(n):
            if i==j:continue
            v=w[i,j];f=ft[i,j];h=fc[i,j]
            d=v*dt[i,j]*ca[i,j]*cb[i,j]
            jac[i,j]+=d;jac[i,i]-=d
            jac[i,n+i]-=v*f*sa[i,j]*cb[i,j];jac[i,n+j]-=v*f*ca[i,j]*sb[i,j]
            d=v*dc[i,j]*sa[i,j]*cb[i,j]
            jac[n+i,j]+=d;jac[n+i,i]-=d
            jac[n+i,n+i]+=v*h*ca[i,j]*cb[i,j];jac[n+i,n+j]-=v*h*sa[i,j]*sb[i,j]
    flow[:n]-=flow[:n].mean();jac[:n]-=jac[:n].mean(axis=0)
    return flow,jac

class SynchronematicGPU:
    def __init__(self,n):
        self.n=n
        devices=[d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 OpenCL GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context,properties=cl.command_queue_properties.PROFILING_ENABLE)
        source=Path(__file__).with_suffix('.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(options=[f'-DN={n}',f'-DD={2*n}',f'-DAUG={4*n}']),'advance')
        self.info={'device':self.device.name.strip(),'precision':'float64','integrator':'RK4','fast_math':False,'oscillators':n,
                   'model':'Leyva et al. 2026 Eqs 10-12, deterministic Fourier approximation, global phase removed',
                   'coefficients':COEFFICIENTS.tolist(),'time':'s=t/L^3; actual positions L*normalized_positions; L=10 gives far-field separation if normalized minimum>=1'}

    def scan(self,positions,starts,dt=.02,transient=200.,duration=400.,progress=None,coefficients=None):
        positions=np.asarray(positions,float);starts=np.asarray(starts,float);n=self.n;lanes=len(starts)
        if positions.shape!=(lanes,n,2) or starts.shape!=(lanes,2*n):raise ValueError('Position/start shapes mismatch')
        if not np.isfinite(positions).all() or not np.isfinite(starts).all():raise ValueError('Nonfinite input')
        if lanes==0 or dt<=0 or duration<=0 or transient<0:raise ValueError('Invalid ensemble or integration interval')
        coefficients=np.broadcast_to(COEFFICIENTS if coefficients is None else coefficients,(lanes,2,4)).copy()
        if not np.isfinite(coefficients).all():raise ValueError('Nonfinite coefficients')
        nt=round(transient/dt);nw=round(duration/(4*dt))
        if not np.isclose(nt*dt,transient) or not np.isclose(4*nw*dt,duration):raise ValueError('Times must be multiples of dt')
        pairs=np.ascontiguousarray([geometry(p) for p in positions]);state=np.zeros((lanes,4*n));state[:,:2*n]=starts
        tangent=np.arange(1,2*n+1,dtype=float);tangent[:n]-=tangent[:n].mean();tangent/=np.linalg.norm(tangent);state[:,2*n:]=tangent
        logs=np.zeros((lanes,4));stats=np.zeros((lanes,4));status=np.zeros(lanes,np.int32)
        arrays=[state,pairs,coefficients,logs,stats,status];buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();gpu=0;done=0;total=nt+4*nw;last=wall
        for window,length in [(-1,nt)]+[(i,nw) for i in range(4)]:
            for offset in range(0,length,128):
                steps=min(128,length-offset)
                event=self.kernel(self.queue,(lanes,),None,*buffers,np.int32(lanes),np.int32(steps),np.float64(dt),np.int32(window));event.wait()
                gpu+=(event.profile.end-event.profile.start)*1e-9;done+=steps
                if progress and time.perf_counter()-last>10:progress(done/total,time.perf_counter()-wall);last=time.perf_counter()
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return {'positions':positions,'starts':starts,'coefficients':coefficients,'final_states':state[:,:2*n],'lle_windows':windows,'lle':windows.mean(axis=1),
                'status':status,'stats':stats,'settings':{'dt':dt,'transient':transient,'duration':duration},
                'timing':{'wall_seconds':time.perf_counter()-wall,'gpu_seconds':gpu},'backend':self.info}
