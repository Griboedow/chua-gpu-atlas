"""Command line: GPU trajectories, 3-D parameter surveys, and offline reports."""
from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from chua_gpu import GPU, DEFAULT_SEEDS, CLASS_NAMES

ROOT=Path(__file__).resolve().parent


def save_json(path,data):
    Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False),encoding="utf-8")


def progress_printer(label):
    last=[-100.0]
    def report(fraction,elapsed):
        if elapsed-last[0]>=10 or fraction==1:
            eta=elapsed*(1/fraction-1)
            print(f"{label}: {fraction:6.1%} | {elapsed:.1f}s elapsed | {eta:.1f}s remaining",flush=True)
            last[0]=elapsed
    return report


def make_grid(center, half_widths, n):
    if n<3 or n%2!=1:
        raise ValueError("Use an odd grid size >=3 so the reference point is sampled exactly")
    p=np.asarray(center,dtype=float)
    widths=np.asarray(half_widths,dtype=float)
    if p.shape!=(4,) or widths.shape!=(3,) or not np.isfinite(p).all() or not np.isfinite(widths).all():
        raise ValueError("center must have four finite parameters; half_widths three")
    if (widths<=0).any() or (p[:2]-widths[:2]<=0).any():
        raise ValueError("half_widths must be positive and alpha,beta ranges must stay positive")
    axes=[np.linspace(c-w,c+w,n) for c,w in zip(p[:3],widths)]
    xyz=np.stack(np.meshgrid(*axes,indexing="ij"),axis=-1).reshape(-1,3)
    return axes,np.column_stack([xyz,np.full(len(xyz),p[3])])


def save_scan(gpu, name, cfg, args, directory):
    axes,parameters=make_grid(cfg["parameters"],cfg["half_widths"],args.grid)
    print(f"Scanning {name}: {len(parameters):,} parameters x {len(DEFAULT_SEEDS)} seeds",flush=True)
    result=gpu.scan(parameters,dt=args.dt,transient=args.transient,duration=args.duration,
                    threshold=args.threshold,chunk_steps=args.chunk_steps,
                    progress=progress_printer(name))
    arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray)}
    candidate=arrays["classes"]==4
    any_seed=candidate.any(axis=1)
    all_seeds=candidate.all(axis=1)
    mixed=(arrays["classes"]!=arrays["classes"][:,:1]).any(axis=1)
    mask=any_seed.reshape((args.grid,)*3)
    touches=[bool(np.take(mask,k,axis=ax).any()) for ax in range(3) for k in (0,-1)]
    # Grid-node fractions are reported directly, not mislabeled as exact volumes.
    meta={"name":name,"configuration":cfg,"grid":args.grid,
          "created_utc":datetime.now(timezone.utc).isoformat(),
          "device":gpu.info,"settings":result["settings"],"timing":result["timing"],
          "parameter_order":["alpha","beta","m0","m1"],
          "varying_parameters":["alpha","beta","m0"],
          "classes":CLASS_NAMES,"grid_nodes":len(parameters),
          "candidate_nodes_any_seed":int(any_seed.sum()),
          "candidate_nodes_all_seeds":int(all_seeds.sum()),
          "mixed_class_nodes":int(mixed.sum()),
          "candidate_fraction_of_nodes":float(any_seed.mean()),
          "boundary_faces_with_candidates":dict(zip(
              ["alpha_min","alpha_max","beta_min","beta_max","m0_min","m0_max"],touches)),
          "scope":"Finite-time, finite-grid evidence from the listed initial conditions. "
                  "Not a rigorous existence region, basin census, or proof of absence elsewhere."}
    path=directory/f"{name}_scan.npz"
    np.savez_compressed(path,parameters=parameters,alpha=axes[0],beta=axes[1],m0=axes[2],
                        candidate_any=any_seed,candidate_all=all_seeds,**arrays)
    save_json(directory/f"{name}_scan.json",meta)
    print(f"Saved {path.name}: {any_seed.sum():,} candidate nodes; "
          f"{all_seeds.sum():,} for all seeds",flush=True)


def save_trajectory(gpu,name,cfg,args,directory):
    p=np.array(cfg["parameters"])
    result=gpu.scan([p],seeds=[DEFAULT_SEEDS[0]],dt=args.dt,transient=args.transient,
                    duration=args.duration,threshold=args.threshold,chunk_steps=args.chunk_steps)
    if result["status"][0,0]!=0:
        raise RuntimeError(f"Reference {name} escaped or failed numerically")
    state=result["final_states"][0,0]
    parts=[]
    steps=int(round(args.plot_duration/args.dt))
    for offset in range(0,steps,1000):
        part=gpu.trace([p],[state],args.dt,min(1000,steps-offset))[:,0,:]
        if not np.isfinite(part).all() or np.abs(part).max()>100:
            raise RuntimeError(f"Reference {name} escaped during plotted segment")
        parts.append(part)
        state=part[-1]
    trajectory=np.concatenate(parts)
    np.savez_compressed(directory/f"{name}_trajectory.npz",xyz=trajectory,
                        t=np.arange(1,steps+1)*args.dt,parameters=p)
    save_json(directory/f"{name}_trajectory.json",{
        "parameters":p.tolist(),"initial_condition":DEFAULT_SEEDS[0],
        "lle":float(result["lle"][0,0]),"lle_windows":result["lle_windows"][0,0].tolist(),
        "class":CLASS_NAMES[int(result["classes"][0,0])],"device":gpu.info,
        "settings":result["settings"],"plotted_duration":steps*args.dt,
        "plot_time_origin":"Relative to the end of transient + LLE observation",
        "configuration":cfg})
    print(f"{name} reference LLE = {result['lle'][0,0]:.6f}",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=["all","scan","trajectory","report","devices"],nargs="?",default="all")
    parser.add_argument("--config",type=Path,default=ROOT/"centers.json")
    parser.add_argument("--center",default="all",help="configuration key or all")
    parser.add_argument("--output",type=Path,default=ROOT/"outputs")
    parser.add_argument("--grid",type=int,default=41)
    parser.add_argument("--dt",type=float,default=0.01)
    parser.add_argument("--transient",type=float,default=300.0)
    parser.add_argument("--duration",type=float,default=800.0)
    parser.add_argument("--plot-duration",type=float,default=200.0)
    parser.add_argument("--threshold",type=float,default=0.01)
    parser.add_argument("--chunk-steps",type=int,default=200)
    parser.add_argument("--precision",choices=["float64","float32"],default="float64")
    parser.add_argument("--device",default="",help="substring of GPU name")
    args=parser.parse_args()
    if args.command=="report":
        from report import build_report
        build_report(args.output)
        return
    gpu=GPU(args.device,args.precision)
    print(json.dumps(gpu.info,indent=2),flush=True)
    if args.command=="devices":
        return
    if not np.isfinite(args.plot_duration) or args.plot_duration<=0:
        parser.error("--plot-duration must be positive and finite")
    configs=json.loads(args.config.read_text(encoding="utf-8"))
    if args.center!="all":
        configs={args.center:configs[args.center]}
    args.output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    for name,cfg in configs.items():
        if not name.replace("_","").replace("-","").isalnum():
            raise ValueError("configuration keys must be simple filenames")
        if args.command in ("all","trajectory"):
            save_trajectory(gpu,name,cfg,args,args.output)
        if args.command in ("all","scan"):
            save_scan(gpu,name,cfg,args,args.output)
    if args.command=="all":
        from report import build_report
        build_report(args.output)
    print(f"Completed in {time.perf_counter()-started:.1f}s",flush=True)


if __name__=="__main__":
    main()
