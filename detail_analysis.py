"""Dense beta transects, recurrence diagnostics, example orbits, and a basin slice."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pyopencl as cl
from chua_gpu import GPU, DEFAULT_SEEDS
from run import save_json, progress_printer

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"/"detailed"
EVENT_KERNEL=r"""
__kernel void peaks(__global real8 *state, __global const real4 *params,
                    __global real *previous_delta, __global real *values,
                    __global int *counts, __global int *status,
                    int n, int steps, int capacity, real dt) {
    int id=get_global_id(0);
    if(id>=n || status[id]) return;
    real8 s=state[id]; real4 p=params[id];
    real prev=previous_delta[id]; int count=counts[id];
    for(int k=0;k<steps;++k) {
        real8 next=rk4(s,p,dt);
        if(!isfinite(next.s0) || !isfinite(next.s1) || !isfinite(next.s2) ||
           fmax(fabs(next.s0),fmax(fabs(next.s1),fabs(next.s2)))>100) {status[id]=1;break;}
        real delta=next.s0-s.s0;
        if(prev>0 && delta<=0) {
            real curvature=delta-prev;
            real height=s.s0-(delta+prev)*(delta+prev)/(8*curvature);
            values[id*capacity+count%capacity]=height;
            ++count;
        }
        prev=delta; s=next;
    }
    state[id]=s; previous_delta[id]=prev; counts[id]=count;
}
"""

def collect_peaks(gpu,parameters,initial,duration=400,dt=.005,capacity=128):
    p=np.asarray(parameters,dtype=gpu.dtype)
    states=np.zeros((len(p),8),dtype=gpu.dtype)
    states[:,:3]=initial
    previous=np.zeros(len(p),dtype=gpu.dtype)
    values=np.full((len(p),capacity),np.nan,dtype=gpu.dtype)
    counts=np.zeros(len(p),dtype=np.int32)
    status=np.zeros(len(p),dtype=np.int32)
    buffers=[gpu.buffer(a,i==1) for i,a in enumerate([states,p,previous,values,counts,status])]
    program=cl.Program(gpu.context,(ROOT/"kernels.cl").read_text()+EVENT_KERNEL).build(options=["-DUSE_DOUBLE"])
    kernel=cl.Kernel(program,"peaks")
    steps=round(duration/dt)
    for offset in range(0,steps,200):
        kernel(gpu.queue,(len(p),),None,*buffers,np.int32(len(p)),np.int32(min(200,steps-offset)),
               np.int32(capacity),gpu.scalar(dt)).wait()
    for a,b in zip([states,p,previous,values,counts,status],buffers):
        cl.enqueue_copy(gpu.queue,a,b).wait()
    for j,count in enumerate(counts):
        if count>capacity:
            values[j]=np.roll(values[j],-int(count%capacity))
    values[status!=0]=np.nan
    return values,counts,status

def peak_period(values,tolerance=5e-4):
    """Smallest repeated peak pattern; zero means not resolved, not proof of chaos."""
    values=np.asarray(values)
    values=values[np.isfinite(values)][-96:]
    if len(values)<16 or np.ptp(values)<1e-10:
        # A limit cycle with constant maxima is period one; equilibrium is excluded upstream.
        return 1 if len(values)>=16 else 0
    for period in range(1,min(24,len(values)//4)+1):
        scale=max(1,float(np.max(np.abs(values))))
        if np.max(np.abs(values[period:]-values[:-period]))<tolerance*scale:
            return period
    return 0

def intervals(axis,mask,min_points=1):
    changes=np.diff(np.r_[False,mask,False].astype(int))
    return [{"start":float(axis[a]),"end":float(axis[b-1]),"nodes":int(b-a)}
            for a,b in zip(np.flatnonzero(changes==1),np.flatnonzero(changes==-1)) if b-a>=min_points]

def line(gpu,name,cfg,n=1601):
    beta=np.linspace(11,19,n) if name=="classic" else np.linspace(23,38,n)
    p=np.tile(cfg["parameters"],(n,1));p[:,1]=beta
    result=gpu.scan(p,dt=.005,transient=600,duration=1600,progress=progress_printer(name+" beta line"))
    peaks=np.full((n,128),np.nan)
    period=np.zeros(n,dtype=np.int16)
    valid=result["status"][:,0]==0
    if valid.any():
        vals,counts,status=collect_peaks(gpu,p[valid],result["final_states"][valid,0])
        peaks[valid]=vals
        for idx,row in zip(np.flatnonzero(valid),vals):
            if result["classes"][idx,0] in [2,3]:
                period[idx]=peak_period(row)
    np.savez_compressed(OUT/f"{name}_beta.npz",beta=beta,parameters=p,peaks=peaks,peak_period=period,
                        **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    # Require every tested seed to be near-zero/nonpositive; the peak recurrence refers to seed 0.
    regular=np.isin(result["classes"],[1,2]).all(axis=1)&(period>0)
    both=(result["classes"]==4)&result["visits_both_outer_regions"]
    one=(result["classes"]==4)&~result["visits_both_outer_regions"]
    meta={"configuration":cfg,"settings":result["settings"],"device":gpu.info,"timing":result["timing"],
          "beta_step":float(beta[1]-beta[0]),"sampled_parameters":n,
          "recurrent_regular_intervals":intervals(beta,regular,3),
          "double_scroll_screen_intervals":intervals(beta,both.any(axis=1),3),
          "single_scroll_screen_intervals":intervals(beta,one.any(axis=1),3),
          "peak_observation_duration":400,"peak_capacity":128,"peak_period_tolerance_relative":.0005,
          "scope":"Independent parameter points and three starts; no parameter continuation. "
                  "Repeated local maxima support periodicity but do not establish an exact mathematical period."}
    save_json(OUT/f"{name}_beta.json",meta)
    print(name,"largest recurrent intervals",sorted(meta["recurrent_regular_intervals"],key=lambda r:-r["nodes"])[:8],flush=True)

def trace_segment(gpu,p,initial,duration=150,dt=.005):
    state=np.asarray(initial)
    parts=[]
    steps=round(duration/dt)
    for offset in range(0,steps,1000):
        data=gpu.trace([p],[state],dt,min(1000,steps-offset))[:,0]
        parts.append(data);state=data[-1]
    return np.concatenate(parts)

def examples(gpu,cfg):
    summaries=[]
    for label,beta in [("double_scroll",28),("periodic_window",29),("single_scroll_pair",32),("regular_cycle",36)]:
        p=np.array(cfg["parameters"]);p[1]=beta
        seeds=[[.7,0,0],[-.7,0,0]]
        result=gpu.scan([p],seeds=seeds,dt=.005,transient=1000,duration=2400)
        if (result["status"]!=0).any():
            raise RuntimeError(f"Example {label} escaped")
        xyz=np.stack([trace_segment(gpu,p,state) for state in result["final_states"][0]])
        vals,counts,status=collect_peaks(gpu,np.tile(p,(2,1)),result["final_states"][0],duration=600)
        periods=[peak_period(v) if cls in [2,3] else 0 for v,cls in zip(vals,result["classes"][0])]
        np.savez_compressed(OUT/f"example_{label}.npz",xyz=xyz,parameters=p,peaks=vals,
                            lle=result["lle"],classes=result["classes"],lle_windows=result["lle_windows"])
        item={"name":label,"parameters":p.tolist(),"initial_states":seeds,
              "lle":result["lle"][0].tolist(),"classes":result["classes"][0].tolist(),
              "x_range_per_start":[[float(t[:,0].min()),float(t[:,0].max())] for t in xyz],
              "peak_period":periods,"settings":result["settings"],"trajectory_duration":150,
              "peak_duration":600}
        summaries.append(item)
        print("example",label,item["lle"],"peak periods",periods,"x range",item["x_range_per_start"],flush=True)
    save_json(OUT/"examples.json",summaries)

def basin(gpu,cfg,n=151):
    p=np.array(cfg["parameters"]);p[1]=32
    x=np.linspace(-2,2,n); y=np.linspace(-.5,.5,n)
    xy=np.stack(np.meshgrid(x,y,indexing="ij"),axis=-1).reshape(-1,2)
    seeds=np.column_stack([xy,np.zeros(len(xy))])
    result=gpu.scan([p],seeds=seeds,dt=.005,transient=600,duration=1600,progress=progress_printer("basin slice"))
    c=result["classes"][0];st=result["stats"][0]
    means=st[:,2]/(1600/.005)
    both=result["visits_both_outer_regions"][0]
    labels=np.full(len(seeds),3,dtype=np.uint8)
    labels[c==0]=0
    labels[(c==4)&~both&(means<0)]=1
    labels[(c==4)&~both&(means>0)]=2
    labels[(c==4)&both]=4
    np.savez_compressed(OUT/"basin.npz",x=x,y=y,initial_states=seeds,parameters=p,
                        basin=labels.reshape(n,n),classes=c.reshape(n,n),mean_x=means.reshape(n,n),
                        lle=result["lle"][0].reshape(n,n))
    values,counts=np.unique(labels,return_counts=True)
    save_json(OUT/"basin.json",{"parameters":p.tolist(),"grid":n,"fixed_initial_z":0,
        "settings":result["settings"],"counts":dict(zip(map(str,values.tolist()),counts.tolist())),
        "labels":{"0":"escaped","1":"left single-scroll chaos","2":"right single-scroll chaos",
                  "3":"other or unresolved","4":"visits both outer regions"},
        "scope":"A 2-D slice of initial-condition space at one fixed parameter tuple, not the 3-D parameter region.",
        "device":gpu.info,"timing":result["timing"]})
    print("basin",dict(zip(values,counts)),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=["all","lines","examples","basin"],default="all",nargs="?")
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((ROOT/"centers.json").read_text())
    gpu=GPU()
    if args.stage in ["all","lines"]:
        for name,conf in cfg.items():line(gpu,name,conf)
    if args.stage in ["all","examples"]:examples(gpu,cfg["reference"])
    if args.stage in ["all","basin"]:basin(gpu,cfg["reference"])

if __name__=="__main__":main()
