"""Independent numerical validation; writes a machine-readable evidence report."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from chua_gpu import GPU, rhs, DEFAULT_SEEDS
from run import save_json, progress_printer

ROOT=Path(__file__).resolve().parent


def cpu_lle(parameters, transient=100, duration=200):
    """Independent adaptive DOP853 integration, renormalized every 1 time unit."""
    state=np.r_[DEFAULT_SEEDS[0],np.array([1,2,3])/np.sqrt(14)]
    total=0.0
    for k in range(int(transient+duration)):
        solution=solve_ivp(rhs,(0,1),state,args=(parameters,),method="DOP853",
                           rtol=2e-9,atol=2e-11,max_step=.05)
        if not solution.success:
            raise RuntimeError(solution.message)
        state=solution.y[:,-1]
        norm=np.linalg.norm(state[3:])
        if k>=transient:
            total+=np.log(norm)
        state[3:]/=norm
    return float(total/duration)


def run_checks(directory):
    start=time.perf_counter()
    gpu=GPU()
    configs=json.loads((ROOT/"centers.json").read_text())
    checks=[]
    details=[]
    for name,cfg in configs.items():
        p=cfg["parameters"]
        initial=DEFAULT_SEEDS[0]
        errors=[]
        for dt in [.01,.005]:
            times=np.arange(1,round(2/dt)+1)*dt
            reference=solve_ivp(rhs,(0,2),initial,args=(p,),method="DOP853",t_eval=times,
                                rtol=1e-12,atol=1e-14,max_step=.002)
            if not reference.success:
                raise RuntimeError(reference.message)
            calculated=gpu.trace([p],[initial],dt,len(times))[:,0,:]
            errors.append(float(np.max(np.abs(reference.y.T-calculated))))
        checks.append({"name":f"{name}: short trajectory versus adaptive DOP853",
                       "passed":errors[0]<.01 and errors[1]<errors[0],
                       "max_abs_error_dt_0_01":errors[0],"max_abs_error_dt_0_005":errors[1]})
        estimate=gpu.scan([p],seeds=[initial],transient=100,duration=200)["lle"][0,0]
        reference_lle=cpu_lle(p)
        checks.append({"name":f"{name}: positive LLE using independent solvers",
                       "passed":estimate>.1 and reference_lle>.1 and abs(estimate-reference_lle)<.12,
                       "gpu_rk4":float(estimate),"cpu_dop853":reference_lle,
                       "note":"Long chaotic trajectories decorrelate; compare finite-time exponents, not endpoint positions."})
        print(f"{name}: short-path errors={errors}; LLE GPU={estimate:.4f}, DOP853={reference_lle:.4f}",flush=True)
    # Smooth linear control with an exact matrix-exponential solution.
    p=[9,100/7,.2,.2]
    jac=np.array([[-p[0]*(1+p[2]),p[0],0],[1,-1,1],[0,-p[1],0]])
    seed=np.array([.7,.1,-.2])
    numerical=gpu.trace([p],[seed],.005,400)[-1,0]
    error=float(np.max(np.abs(numerical-expm(2*jac)@seed)))
    checks.append({"name":"Linear control against exact matrix exponential","passed":error<1e-7,"max_abs_error":error})
    result=gpu.scan([p],seeds=[seed],transient=300,duration=800)
    eigenvalue=float(np.linalg.eigvals(jac).real.max())
    checks.append({"name":"Stable equilibrium: negative LLE matches eigenvalue",
                   "passed":int(result['classes'][0,0])==1 and abs(float(result['lle'][0,0])-eigenvalue)<.01,
                   "gpu_lle":float(result['lle'][0,0]),"expected_lle":eigenvalue})
    # An exact unstable equilibrium has positive tangent growth but is no strange attractor.
    eq=gpu.scan([configs["classic"]["parameters"]],seeds=[[0,0,0]],transient=100,duration=200)
    checks.append({"name":"Unstable stationary solution is not counted as chaos",
                   "passed":int(eq["classes"][0,0])!=4,"classification":int(eq["classes"][0,0])})
    # An escaped trajectory must have no finite reported LLE, even after partial windows.
    escaped=gpu.scan([[9,100/7,-.8,-5/7]],seeds=[DEFAULT_SEEDS[0]],transient=100,duration=200)
    checks.append({"name":"Escape detection masks truncated LLE","passed":int(escaped['classes'][0,0])==0 and bool(np.isnan(escaped['lle'][0,0]))})
    # Stratify by base candidate membership: reproducible and explicit, not an unbiased population estimate.
    rng=np.random.default_rng(20261029)
    for file in sorted(directory.glob("*_scan.npz")):
        name=file.stem[:-5]
        data=np.load(file)
        meta=json.loads(file.with_suffix(".json").read_text())
        groups=[]
        indices=[]
        for value in [False,True]:
            pool=np.flatnonzero(data["candidate_any"]==value)
            selected=rng.choice(pool,min(256,len(pool)),replace=False)
            indices.extend(selected.tolist())
            groups.extend([value]*len(selected))
        indices=np.asarray(indices)
        parameters=data["parameters"][indices]
        settings=meta["settings"]
        fine=gpu.scan(parameters,seeds=settings["seeds"],dt=settings["dt"]/2,
                      transient=settings["transient"]*2,duration=settings["duration"]*2,
                      threshold=settings["threshold"],escape=settings["escape"],
                      progress=progress_printer(f"{name} refinement"))
        fine_any=(fine["classes"]==4).any(axis=1)
        coarse=data["candidate_any"][indices]
        summary={"name":name,"sample_size":len(indices),
                 "sampling":"256 random base candidates and 256 random other grid nodes; fixed RNG seed 20261029",
                 "fine_settings":fine["settings"],"timing":fine["timing"],
                 "candidate_membership_agreement":float((fine_any==coarse).mean()),
                 "base_candidates_retained":int((coarse&fine_any).sum()),
                 "base_candidate_sample_count":int(coarse.sum()),
                 "other_nodes_becoming_candidates":int((~coarse&fine_any).sum()),
                 "other_node_sample_count":int((~coarse).sum()),
                 "interpretation":"This stratified comparison is a sensitivity diagnostic, not a full-grid certificate or confidence interval."}
        details.append(summary)
        np.savez_compressed(directory/f"{name}_refinement.npz",indices=indices,parameters=parameters,
                            base_candidate=coarse,fine_candidate=fine_any,
                            fine_lle=fine["lle"],fine_classes=fine["classes"],fine_windows=fine["lle_windows"])
        print(f"{name}: refinement membership agreement {summary['candidate_membership_agreement']:.1%}",flush=True)
    # Convert NumPy booleans before strict JSON serialization.
    for check in checks:
        check["passed"]=bool(check["passed"])
    evidence={"device":gpu.info,"checks":checks,"all_checks_passed":all(c["passed"] for c in checks),
              "refinement":details,"wall_seconds":time.perf_counter()-start}
    save_json(directory/"validation.json",evidence)
    print(json.dumps(evidence,indent=2),flush=True)
    if not evidence["all_checks_passed"]:
        raise SystemExit("A numerical validation check failed; inspect validation.json")


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"outputs")
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    run_checks(args.output)
