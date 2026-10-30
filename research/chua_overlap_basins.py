"""Resolve attraction basins of four independently verified periodic cycles."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pyopencl as cl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from chua_gpu import GPU

KERNEL=r'''
__kernel void section_returns(__global real8 *state,__global const real4 *params,
 __global real4 *returns,__global int *counts,__global int *status,int n,int steps,real dt){
 int id=get_global_id(0);if(id>=n||status[id])return;
 real8 s=state[id];real4 p=params[id];int count=counts[id];
 for(int k=0;k<steps;k++){
  real8 v=rk4(s,p,dt);
  if(!isfinite(v.s0)||!isfinite(v.s1)||!isfinite(v.s2)||fmax(fabs(v.s0),fmax(fabs(v.s1),fabs(v.s2)))>100){status[id]=1;break;}
  if(s.s0<=0 && v.s0>0){real u=-s.s0/(v.s0-s.s0);returns[id*6+count%6]=(real4)(s.s1+u*(v.s1-s.s1),s.s2+u*(v.s2-s.s2),0,0);count++;}
  s=v;
 }
 state[id]=s;counts[id]=count;
}
'''


def classify(gpu,params,starts,dt,duration,reference):
    params=np.ascontiguousarray(params,dtype=gpu.dtype)
    n=len(starts);state=np.zeros((n,8));state[:,:3]=starts;returns=np.full((n,6,4),np.nan);counts=np.zeros(n,np.int32);status=np.zeros(n,np.int32)
    buffers=[gpu.buffer(x,i==1) for i,x in enumerate([state,params,returns,counts,status])]
    program=cl.Program(gpu.context,(ROOT/'kernels.cl').read_text()+KERNEL).build(options=['-DUSE_DOUBLE']);kernel=cl.Kernel(program,'section_returns')
    steps=round(duration/dt);last=-1
    for offset in range(0,steps,200):
        kernel(gpu.queue,(n,),None,*buffers,np.int32(n),np.int32(min(200,steps-offset)),gpu.scalar(dt)).wait()
        done=int((offset+200)*10/steps)
        if done>last:print('basins',done*10,'percent',flush=True);last=done
    for x,b in zip([state,params,returns,counts,status],buffers):cl.enqueue_copy(gpu.queue,x,b).wait()
    distances=np.max(np.linalg.norm(returns[:,:,:2,None]-reference.T[None,None,:,:],axis=2),axis=1)
    labels=np.argmin(distances,axis=1)+1;labels[(np.min(distances,axis=1)>2e-5)|(counts<6)|(status!=0)]=0
    return dict(labels=labels,distances=distances,counts=counts,status=status,final_states=state[:,:3],returns=returns[:,:,:2])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--nodes',type=int,default=201);p.add_argument('--tight',action='store_true')
    args=p.parse_args();out=args.output;data=json.loads((out/'four-cycles.json').read_text());gpu=GPU()
    if args.tight:
        raw=np.load(out/'overlap-basins.npz');initial=raw['initial'][::4];params=raw['parameters'][::4]
        result=classify(gpu,params,initial,.000625,4000,np.array([r['cycle']['q'][:2] for r in data['records']]))
        np.savez_compressed(out/'overlap-basins-tight.npz',initial=initial,**result)
        same=result['labels']==raw['labels'][::4];report=dict(nodes=len(same),matching_labels=int(same.sum()),fraction_matching=float(same.mean()),
            baseline_unresolved=int((raw['labels'][::4]==0).sum()),tight_unresolved=int((result['labels']==0).sum()),
            transition_counts=np.array([[int(((raw['labels'][::4]==a)&(result['labels']==b)).sum()) for b in range(5)] for a in range(5)]).tolist())
        (out/'overlap-basins-tight.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(report,flush=True);return
    y=np.linspace(.19,.30,args.nodes);z=np.linspace(-.54,-.30,args.nodes)
    initial=np.array([[0.,a,b] for a in y for b in z]);params=np.tile([data[k] for k in ['alpha','beta','m0','m1']],(len(initial),1))
    result=classify(gpu,params,initial,.00125,2000,np.array([r['cycle']['q'][:2] for r in data['records']]))
    np.savez_compressed(out/'overlap-basins.npz',y=y,z=z,initial=initial,parameters=params,**result)
    report=dict(shape=[len(y),len(z)],dt=.00125,duration=2000,section_tolerance=2e-5,device=gpu.info,
        counts={str(k):int((result['labels']==k).sum()) for k in range(5)},
        labels=['unresolved']+[r['id'] for r in data['records']],scope='Sampled area fractions of this initial-condition rectangle, not full basin volumes.')
    (out/'overlap-basins.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(report,flush=True)


if __name__=='__main__':main()
