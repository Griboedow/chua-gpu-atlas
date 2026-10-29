"""Batched FP64 scans of the literal five-state two-strain model."""
import argparse
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from two_strain import parameters, initial, encode, decode, save, ROOT, COLUMNS


class TwoStrainGPU:
    def __init__(self):
        devices=[d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        self.kernel=cl.Kernel(cl.Program(self.context,Path(__file__).with_suffix('.cl').read_text()).build(),'advance')

    def scan(self,p,starts,transient=20000.,duration=40000.,rtol=1e-8,max_step=1.,max_wall_seconds=900.,max_attempts=5000000):
        p=np.ascontiguousarray(p,dtype=float);starts=np.asarray(starts,dtype=float);n=len(p)
        if p.shape!=(n,8) or starts.shape!=(n,5) or n<1 or not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Invalid input arrays')
        if np.min(p)<=0 or transient<0 or min(duration,rtol,max_step,max_wall_seconds)<=0:raise ValueError('Invalid settings')
        states=np.zeros((n,16));states[:,:5]=starts;states[:,5:10]=np.array([1.,.3,.4,-.2,.5])/np.sqrt(1.54)
        clocks=np.zeros((n,4));clocks[:,1]=.001;logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,[0,2,4,6]]=np.inf;stats[:,[1,3,5,7]]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[states,p,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();last=wall;calls=0
        while True:
            self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(256),np.float64(transient),np.float64(duration),np.float64(rtol),np.float64(rtol*.001),np.float64(max_step),np.int32(max_attempts)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<transient+duration)
                if time.perf_counter()-last>20:
                    print(f'{clocks[:,0].mean()/(transient+duration):.1%}, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-wall>max_wall_seconds:
                    status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,initial_log_states=starts,final_log_states=states[:,:5].copy(),
            final_states=np.array([decode(q,param) for q,param in zip(states[:,:5],p)]),
            lle=windows.mean(axis=1),lle_windows=windows,status=status,stats=stats,elapsed=clocks[:,0],step_counts=counts,
            settings=dict(transient=transient,duration=duration,rtol=rtol,atol=rtol*.001,max_step=max_step,max_attempts=max_attempts),
            method='Adaptive DP5(4); log prevalence for state, original compartment variational equations with fixed diagonal scaling',
            device=self.device.name,wall_seconds=time.perf_counter()-wall,parameter_columns=COLUMNS,
            log_state_columns=['S','log(I1/mu)','log((I2+I12)/mu)','I12/(I2+I12)','R1'],
            stats_columns=['min_log_I1_over_mu','max_log_I1_over_mu','min_log_total2_over_mu','max_log_total2_over_mu','min_S','max_S','min_recovered','max_population_sum'],
            status_codes={'0':'completed','1':'population domain failure','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved'})


def survey(levels=129):
    ps=[];starts=[];kinds=[]
    for r2 in np.linspace(1.24,2.04,levels):
        p=parameters(r2)
        for kind in ['coexistence','invasion']:
            ps.append(p);starts.append(encode(initial(p,kind),p));kinds.append(kind)
    data=TwoStrainGPU().scan(np.array(ps),np.array(starts))
    data['start_kinds']=np.array(kinds)
    ok=data['status']==0;candidate=ok&(data['lle']>1e-5)&np.all(data['lle_windows']>5e-6,axis=1)
    data['preliminary_candidate']=candidate
    data.update(lanes=len(ps),completed=int(ok.sum()),preliminary_positive=int(candidate.sum()),
        scope='Targeted source-parameter slice with two starts; positive finite-time exponents do not certify attractors or novelty')
    save(ROOT/'slice-scout',data)
    print({k:data[k] for k in ['lanes','completed','preliminary_positive','wall_seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--levels',type=int,default=129);args=parser.parse_args();survey(args.levels)
