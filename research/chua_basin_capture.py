"""Measure capture times and sensitivity of four-attractor basin outcomes."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pyopencl as cl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from chua_gpu import GPU

KERNEL=r'''
__kernel void capture(__global real8 *state,__global const real4 *params,
 __global const real4 *refs,__global int *labels,__global int *lastmatch,
 __global int *streak,__global real *times,int n,int steps,real dt,real offset){
 int id=get_global_id(0);if(id>=n||labels[id]!=0)return;
 real8 s=state[id];real4 p=params[id];int previous=lastmatch[id],run=streak[id];
 for(int k=0;k<steps;k++){
  real8 v=rk4(s,p,dt);
  if(!isfinite(v.s0)||!isfinite(v.s1)||!isfinite(v.s2)||fmax(fabs(v.s0),fmax(fabs(v.s1),fabs(v.s2)))>100){labels[id]=-1;break;}
  if(s.s0<=0&&v.s0>0){
   real u=-s.s0/(v.s0-s.s0);real y=s.s1+u*(v.s1-s.s1),z=s.s2+u*(v.s2-s.s2);
   real best=1e30;int match=0;
   for(int j=0;j<4;j++){real dy=y-refs[j].s0,dz=z-refs[j].s1,d=dy*dy+dz*dz;if(d<best){best=d;match=j+1;}}
   if(best>4e-12)match=0;
   run=(match&&match==previous)?run+1:(match?1:0);previous=match;
   if(run>=8){labels[id]=match;times[id]=offset+((real)k+u)*dt;s=v;break;}
  }
  s=v;
 }
 state[id]=s;lastmatch[id]=previous;streak[id]=run;
}
'''


def run(gpu,params,initial,refs,dt,duration=4000):
    n=len(initial);states=np.zeros((n,8));states[:,:3]=initial;params=np.ascontiguousarray(params);reference=np.zeros((4,4));reference[:,:2]=refs
    labels=np.zeros(n,np.int32);last=np.zeros(n,np.int32);streak=np.zeros(n,np.int32);times=np.full(n,np.nan)
    arrays=[states,params,reference,labels,last,streak,times];buffers=[gpu.buffer(x,i in [1,2]) for i,x in enumerate(arrays)]
    program=cl.Program(gpu.context,(ROOT/'kernels.cl').read_text()+KERNEL).build(options=['-DUSE_DOUBLE']);kernel=cl.Kernel(program,'capture')
    steps=round(duration/dt);tick=-1
    for offset in range(0,steps,200):
        kernel(gpu.queue,(n,),None,*buffers,np.int32(n),np.int32(min(200,steps-offset)),gpu.scalar(dt),gpu.scalar(offset*dt)).wait()
        t=int((offset+200)*dt/500)
        if t>tick:
            cl.enqueue_copy(gpu.queue,labels,buffers[3]).wait();print('capture',round((offset+200)*dt,1),'remaining',int((labels==0).sum()),flush=True);tick=t
            if np.all(labels!=0):break
    for a,b in zip(arrays,buffers):cl.enqueue_copy(gpu.queue,a,b).wait()
    return labels,times


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--perturb',action='store_true');p.add_argument('--transients',action='store_true');p.add_argument('--full',action='store_true')
    args=p.parse_args();perturb=args.perturb;out=args.output;gpu=GPU();r=json.loads((out/'four-cycles.json').read_text())
    params=np.array([r[k] for k in ['alpha','beta','m0','m1']]);refs=np.array([x['cycle']['q'][:2] for x in r['records']])
    raw=np.load(out/'overlap-basins.npz')
    if args.transients:
        indices=np.load(out/'basin-perturbations.npz')['indices'];initial=raw['initial'][indices]
        result=gpu.scan(np.tile(params,(len(initial),1)),initial_states=initial[:,None,:],dt=.000625,transient=0,duration=80)
        np.savez_compressed(out/'basin-transients.npz',indices=indices,initial=initial,**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        report=dict(settings=result['settings'],samples=len(initial),positive_four_window_screens=int((result['classes']==4).sum()),
            lle_quantiles={str(q):float(np.quantile(result['lle'],q)) for q in [0,.1,.25,.5,.75,.9,1.]},
            note='Finite initial-time instability, not an asymptotic chaotic attractor: the same parameter point has attracting cycles.')
        (out/'basin-transients.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(report,flush=True);return
    if perturb:
        rng=np.random.default_rng(20261030);indices=rng.choice(len(raw['initial']),1024,replace=False);base=raw['initial'][indices]
        eps=np.array([0.,1e-2,1e-3,1e-4,1e-5,1e-6,1e-8,1e-10]);initial=np.tile(base,(len(eps),1));initial[:,1]+=np.repeat(eps,len(base))
        labels,times=run(gpu,np.tile(params,(len(initial),1)),initial,refs,.000625)
        matrix=labels.reshape(len(eps),-1);fractions=np.mean(matrix!=matrix[0],axis=1)
        report=dict(epsilon=eps.tolist(),changed_fraction=fractions.tolist(),samples=len(base),dt=.000625,
                    note='Sensitivity at fixed numerical precision; the smallest perturbations need not resolve exact basin boundaries.')
        np.savez_compressed(out/'basin-perturbations.npz',indices=indices,epsilon=eps,labels=matrix,times=times.reshape(len(eps),-1))
        name='basin-perturbations.json'
    else:
        stride=1 if args.full else 16;dt=.000625 if args.full else .0003125
        initial=raw['initial'][::stride];labels,times=run(gpu,np.tile(params,(len(initial),1)),initial,refs,dt)
        tight=np.load(out/'overlap-basins-tight.npz')['labels']
        comparison=labels[::4] if args.full else labels
        target=tight if args.full else tight[::4]
        report=dict(samples=len(labels),dt=dt,counts={str(k):int((labels==k).sum()) for k in range(5)},
            match_baseline=float(np.mean(labels==raw['labels'][::stride])),match_half_step=float(np.mean(comparison==target)),
            capture_time_quantiles={str(q):float(np.nanquantile(times,q)) for q in [0,.1,.25,.5,.75,.9,.99,1.]})
        name='basin-capture-full' if args.full else 'basin-capture'
        np.savez_compressed(out/(name+'.npz'),initial=initial,labels=labels,times=times);name+='.json'
    (out/name).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(report,flush=True)


if __name__=='__main__':main()
