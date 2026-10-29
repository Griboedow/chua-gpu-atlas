"""Validate the new peak sampler and the interpretation of narrow regular windows."""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from chua_gpu import GPU,rhs
from detail_analysis import collect_peaks,peak_period,OUT
from run import save_json

gpu=GPU()
checks=[]
for name,p,expected in [("classic_five_peak_window",[9,14.32,-8/7,-5/7],5),
                        ("reference_six_peak_window",[15.6,29,-1.143,-.714],6)]:
    result=gpu.scan([p],dt=.0025,transient=1200,duration=3200)
    vals,counts,status=collect_peaks(gpu,[p],[result['final_states'][0,0]],dt=.0025,duration=400)
    gp=peak_period(vals[0])
    t=np.arange(0,300.0001,.0025)
    sol=solve_ivp(rhs,(0,300),result['final_states'][0,0],args=(p,),method="DOP853",
                  t_eval=t,rtol=1e-10,atol=1e-12,max_step=.03)
    if not sol.success:raise RuntimeError(sol.message)
    idx,_=find_peaks(sol.y[0]);cpu_period=peak_period(sol.y[0,idx])
    passed=bool((result['classes']==2).all() and gp==expected and cpu_period==expected)
    checks.append({"name":name,"passed":passed,"expected_peak_pattern":expected,
        "gpu_peak_pattern":gp,"dop853_peak_pattern":cpu_period,"refined_lle":result['lle'][0].tolist(),
        "refined_settings":result['settings'],"dop853_time":300})
    print(name,"GPU",gp,"DOP853",cpu_period,"LLE",result['lle'][0],flush=True)

# Compare the event buffer against peaks extracted from a retained trajectory,
# including chronological order after a ring-buffer wrap.
p=[15.6,36,-1.143,-.714]
data=np.load(OUT/'example_regular_cycle.npz')
initial=data['xyz'][0,-1]
values,counts,status=collect_peaks(gpu,[p],[initial],dt=.005,duration=120,capacity=16)
state=initial;parts=[]
for _ in range(24):
    part=gpu.trace([p],[state],.005,1000)[:,0]
    parts.append(part);state=part[-1]
x=np.concatenate(parts)[:,0]
idx,_=find_peaks(x)
previous=x[idx-1];current=x[idx];following=x[idx+1]
heights=current-(following-previous)**2/(8*(previous-2*current+following))
error=float(np.max(np.abs(heights[-16:]-values[0])))
checks.append({"name":"GPU peak buffer versus retained trajectory, with wrap",
               "passed":bool(counts[0]>16 and error<1e-7),"max_abs_error":error,
               "total_peaks":int(counts[0]),"buffer_capacity":16})

# Symmetry is an exact property of the odd vector field and the symmetric start grid.
basin=np.load(OUT/'basin.npz')['basin']
mirrored=basin[::-1,::-1]
swapped=np.where(basin==1,2,np.where(basin==2,1,basin))
checks.append({"name":"Basin slice respects central symmetry","passed":bool(np.array_equal(swapped,mirrored)),
               "unresolved_nodes":int((basin==3).sum())})

# Local-sheet sensitivity check: equal samples of candidates and regular nodes.
d=np.load(OUT/'window_zoom_scan.npz')
regular=np.isin(d['classes'],[1,2]).all(axis=1)&(d['classes']==2).any(axis=1)
rng=np.random.default_rng(20261029)
indices=np.r_[rng.choice(np.flatnonzero(regular),128,replace=False),
              rng.choice(np.flatnonzero(d['candidate_all']),128,replace=False)]
fine=gpu.scan(d['parameters'][indices],dt=.0025,transient=1200,duration=3200)
fine_regular=np.isin(fine['classes'],[1,2]).all(axis=1)&(fine['classes']==2).any(axis=1)
fine_chaos=(fine['classes']==4).all(axis=1)
diagnostic={"sampled_regular":128,"regular_retained":int(fine_regular[:128].sum()),
            "sampled_all_seed_chaos":128,"all_seed_chaos_retained":int(fine_chaos[128:].sum()),
            "settings":fine['settings'],
            "scope":"A stratified 256-node sensitivity check, not a certificate for the full local surface."}
np.savez_compressed(OUT/'zoom_refinement.npz',indices=indices,parameters=d['parameters'][indices],
                    fine_classes=fine['classes'],fine_lle=fine['lle'])
result={"checks":checks,"all_checks_passed":all(c['passed'] for c in checks),
        "zoom_refinement":diagnostic,"device":gpu.info}
save_json(OUT/'validation_details.json',result)
print(result,flush=True)
if not result['all_checks_passed']:raise SystemExit("Detailed validation failed")
