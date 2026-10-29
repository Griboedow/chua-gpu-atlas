"""Adaptive FP64 GPU study of an explicitly assumed finite optical width.

The logistic heating law modifies the ideal relay model. It is not a measured
profile for the source apparatus. Physical velocities are normalized internally
by sqrt(k); inputs and saved final states use the original six coordinates.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.special import expit


def rhs_and_jac(y,p):
    k,b,dm,dp,w=p;xm,vm,ym,xp,vp,yp=y
    sm=expit((xp-dp)/w);sp=expit((dm-xm)/w)
    f=np.array([vm,k*(ym-xm)-b*vm,sm-ym,vp,k*(yp-xp)-b*vp,sp-yp])
    a=np.array([[0,1,0],[-k,-b,k],[0,0,-1.]])
    jac=np.zeros((6,6));jac[:3,:3]=a;jac[3:,3:]=a;jac[2,3]=sm*(1-sm)/w;jac[5,0]=-sp*(1-sp)/w
    return f,jac


class SmoothOperatorsGPU:
    def __init__(self):
        devices=[d for platform in cl.get_platforms() for d in platform.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        self.kernel=cl.Kernel(cl.Program(self.context,Path(__file__).with_suffix('.cl').read_text()).build(),'advance')

    def scan(self,parameters,starts,transient=200.,duration=800.,rtol=1e-7,step_factor=1.,max_attempts=10000000,max_wall_seconds=900.,balanced_tangents=None):
        parameters=np.ascontiguousarray(parameters,dtype=float);starts=np.ascontiguousarray(starts,dtype=float);n=len(parameters)
        if parameters.shape!=(n,5) or starts.shape!=(n,6) or n<1:raise ValueError('Input shape mismatch')
        if not np.isfinite(parameters).all() or not np.isfinite(starts).all():raise ValueError('Nonfinite input')
        if np.any(parameters[:,:2]<=0) or np.any(parameters[:,4]<=0) or np.any((parameters[:,2:4]<=0)|(parameters[:,2:4]>=1)):raise ValueError('Invalid parameters')
        if np.any((starts[:,[2,5]]<0)|(starts[:,[2,5]]>1)) or transient<0 or min(duration,rtol,step_factor,max_wall_seconds)<=0:raise ValueError('Invalid settings or temperatures')
        p=np.zeros((n,8));p[:,:5]=parameters
        state=np.zeros((n,16));state[:,:6]=starts;state[:,[1,4]]/=np.sqrt(parameters[:,0,None]);state[:,6:12]=np.array([1,.3,.4,-.2,.5,-.1])/np.sqrt(1.55)
        if balanced_tangents is not None:
            tangent=np.asarray(balanced_tangents,dtype=float)
            if tangent.shape!=(n,6) or not np.isfinite(tangent).all() or np.any(np.linalg.norm(tangent,axis=1)==0):raise ValueError('Invalid balanced tangents')
            state[:,6:12]=tangent/np.linalg.norm(tangent,axis=1)[:,None]
        clocks=np.zeros((n,4));clocks[:,1]=np.minimum(.001,.1/np.sqrt(p[:,0]));logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,[0,2,4,6]]=np.inf;stats[:,[1,3,5,7]]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[state,p,clocks,logs,stats,status,counts];buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();last=wall;calls=0
        while True:
            self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(128),np.float64(transient),np.float64(duration),np.float64(rtol),np.float64(rtol*.001),np.float64(step_factor),np.int32(max_attempts)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<transient+duration)
                if time.perf_counter()-last>20:
                    print(f'{clocks[:,0].mean()/(transient+duration):.1%}, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-wall>max_wall_seconds:
                    status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for a,buf in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,buf).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        final=state[:,:6].copy();final[:,[1,4]]*=np.sqrt(parameters[:,0,None])
        return dict(parameters=parameters,starts=starts,final_states=final,balanced_tangents=state[:,6:12].copy(),lle=windows.mean(axis=1),lle_windows=windows,
                    elapsed=clocks[:,0],status=status,step_counts=counts,stats=stats,
                    settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol*.001,step_factor=step_factor,max_attempts=max_attempts,max_wall_seconds=max_wall_seconds),
                    method='Dormand-Prince 5(4) with analytic tangent and extra optical-transition step cap; explicit logistic model modification',
                    parameter_columns=['k','b','delta_minus','delta_plus','width'],
                    stats_columns=['min_x_minus','max_x_minus','min_x_plus','max_x_plus','min_temperature','max_temperature','late_min_x_minus','late_max_x_minus'],
                    status_codes={'0':'completed','1':'thermal invariant-domain failure','2':'numerical failure or underflow','3':'step budget unresolved','4':'wall budget unresolved'},
                    device=self.device.name,wall_seconds=time.perf_counter()-wall)


def save(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,**{k:v for k,v in data.items() if isinstance(v,np.ndarray)})
    info={k:v for k,v in data.items() if not isinstance(v,np.ndarray)}
    ok=data['status']==0;positive=ok&(data['lle']>.005)&np.all(data['lle_windows']>.002,axis=1)
    info.update(lanes=len(ok),completed=int(ok.sum()),retained=int(positive.sum()))
    path.with_suffix('.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2),flush=True)


def validate():
    p=np.array([74.1567393414624,1.784180421282182,.4139517339970917,.25683679452165964,.005]);y=np.array([.3,.4,.2,.4,-.2,.5])
    tangent=np.array([1,.3*np.sqrt(p[0]),.4,-.2,.5*np.sqrt(p[0]),-.1])/np.sqrt(1.55)
    def augmented(t,z):
        f,j=rhs_and_jac(z[:6],p);return np.r_[f,j@z[6:]]
    cpu=solve_ivp(augmented,(0,2),np.r_[y,tangent],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.005).y[:,-1]
    norm=np.linalg.norm(cpu[6:]/np.array([1,np.sqrt(p[0]),1,1,np.sqrt(p[0]),1]))
    gpu=SmoothOperatorsGPU().scan(p[None],y[None],transient=0,duration=2,rtol=1e-10,step_factor=.5)
    eps=1e-6;fd=np.column_stack([(rhs_and_jac(y+eps*e,p)[0]-rhs_and_jac(y-eps*e,p)[0])/(2*eps) for e in np.eye(6)])
    out=dict(cpu_gpu_state_error=float(max(abs(gpu['final_states'][0]-cpu[:6]))),cpu_gpu_lle_error=float(abs(gpu['lle'][0]-np.log(norm)/2)),jacobian_fd_error=float(np.max(abs(fd-rhs_and_jac(y,p)[1]))),status=int(gpu['status'][0]))
    out['passed']=out['cpu_gpu_state_error']<1e-7 and out['cpu_gpu_lle_error']<1e-7 and out['jacobian_fd_error']<1e-6 and out['status']==0
    Path('results/light-operators/smooth-gpu-validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
    if not out['passed']:raise RuntimeError('Validation failed')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--scan',action='store_true');parser.add_argument('--refine');parser.add_argument('--reentrant-starts',action='store_true');parser.add_argument('--output');parser.add_argument('--transient',type=float,default=200);parser.add_argument('--duration',type=float,default=800);parser.add_argument('--rtol',type=float,default=1e-7);parser.add_argument('--step-factor',type=float,default=1.)
    args=parser.parse_args()
    if args.scan:
        source=np.load('results/light-operators/coupled-scout.npz');rng=np.random.default_rng(2026102967)
        widths=[.0005,.001,.002,.003,.005,.007,.01,.015,.02,.03,.04,.05,.06,.065,.08,.1]
        parameters=[];starts=[];refs=[]
        for index in [1026,3054]:
            for width in widths:
                for _ in range(8):
                    parameters.append(np.r_[source['parameters'][index],width]);y=rng.uniform(.05,.95,6);y[[1,4]]=0;starts.append(y);refs.append(index)
        data=SmoothOperatorsGPU().scan(parameters,starts,args.transient,args.duration,args.rtol,args.step_factor);data['reference_index']=np.array(refs)
        save(args.output or 'results/light-operators/smooth-width-scout.npz',data)
    elif args.reentrant_starts:
        filename='results/light-operators/smooth-width-scout.npz';source=np.load(filename)
        ids=np.flatnonzero((source['reference_index']==1026)&(source['parameters'][:,4]==.04))
        data=SmoothOperatorsGPU().scan(source['parameters'][ids],source['final_states'][ids],args.transient,args.duration,args.rtol,args.step_factor)
        data['source_indices']=ids;data['source_file']=filename;data['reference_index']=source['reference_index'][ids]
        save(args.output or 'results/light-operators/smooth-reentrant-starts.npz',data)
    elif args.refine:
        source=np.load(args.refine);ids=[]
        for ref in np.unique(source['reference_index']):
            for width in np.unique(source['parameters'][:,4]):
                group=np.flatnonzero((source['reference_index']==ref)&(source['parameters'][:,4]==width)&(source['status']==0))
                if len(group):ids.extend(group[np.argsort(source['lle'][group])[-2:]])
        ids=np.array(ids);data=SmoothOperatorsGPU().scan(source['parameters'][ids],source['final_states'][ids],args.transient,args.duration,args.rtol,args.step_factor)
        data['source_indices']=ids;data['reference_index']=source['reference_index'][ids];data['source_file']=args.refine
        save(args.output or 'results/light-operators/smooth-width-refined.npz',data)
    else:validate()
