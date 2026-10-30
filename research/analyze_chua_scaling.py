"""Recorded finite-sample summaries and an exact affine-transient counterexample."""
import json
import hashlib
import numpy as np
from scipy.linalg import eig, expm
from chua_scaling_gpu import OUT, field_jac


def verify_resume_source():
    if not (OUT/'refined-clock-resume.npz').exists():return
    meta=json.loads((OUT/'refined-clock-resume.json').read_text())
    actual=hashlib.sha256((OUT/'refined.npz').read_bytes()).hexdigest()
    if meta.get('source_refined_sha256')!=actual:
        raise RuntimeError('Stale resumed dataset: its original refined NPZ hash differs. Do not merge results from different runs.')


def affine_control():
    source=np.load(OUT/'scout.npz');i=963;p=source['parameters'][i];alpha=p[0]
    xeq=(p[3]-p[2])/(1+p[3]);eq=np.array([xeq,0,-xeq/alpha])
    initial=source['boundary_states'][i,0];J=alpha*field_jac(eq,p)[1]
    roots,V=eig(J);coeff=np.linalg.solve(V,initial-eq)
    duration=200.
    # Triangle bound valid at every t in [0,duration], not just sampled states.
    uniform_x_bound=float(np.sum(np.abs(V[0]*coeff)*np.exp(np.maximum(roots.real,0)*duration)))
    endpoint=eq+expm(J*duration)@(initial-eq)
    data=dict(source_index=i,parameters=p.tolist(),outer_equilibrium=eq.tolist(),
              observation_initial_state=initial.tolist(),duration=duration,
              original_time_eigenvalues_real=roots.real.tolist(),original_time_eigenvalues_imag=roots.imag.tolist(),
              uniform_affine_x_deviation_bound=uniform_x_bound,
              entire_affine_solution_stays_above_one=bool(xeq-uniform_x_bound>1),
              exact_affine_endpoint=endpoint.tolist(),gpu_endpoint_error=float(np.max(abs(endpoint-source['final_states'][i]))),
              gpu_window_lle=source['lle_windows'][i].tolist(),
              interpretation='This entire observation is self-consistently confined to one affine branch. Its positive exponent describes an unstable oscillatory transient, not a strange attractor. Future switching remains possible.')
    if not data['entire_affine_solution_stays_above_one']:raise RuntimeError('Affine bound fails')
    (OUT/'affine-transient-963.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    return data


def summarize():
    verify_resume_source()
    data={}
    for name,maskkey in [('scout','preliminary_candidate'),('refined','retained_candidate'),('refined-clock-resume','retained_candidate')]:
        path=OUT/f'{name}.npz'
        if not path.exists():continue
        x=np.load(path);p=x['parameters'];ok=x['status']==0;positive=x[maskkey]
        rows=[]
        for m in np.unique(p[:,2]):
            for a in np.unique(p[:,0]):
                take=(p[:,2]==m)&(p[:,0]==a)
                if not take.any():continue
                rows.append(dict(alpha=float(a),m0=float(m),lanes=int(take.sum()),completed=int((take&ok).sum()),
                                 positive=int((take&positive).sum()),
                                 max_completed_lle=float(x['lle'][take&ok].max()) if np.any(take&ok) else None))
        data[name]=dict(total=len(p),status_counts={str(k):int(v) for k,v in zip(*np.unique(x['status'],return_counts=True))},
                        positive=int(positive.sum()),by_slope_alpha=rows)
    if (OUT/'refined-clock-resume.npz').exists():
        base=np.load(OUT/'refined.npz');resume=np.load(OUT/'refined-clock-resume.npz')
        status=base['status'].copy();positive=base['retained_candidate'].copy()
        status[resume['refined_indices']]=resume['status'];positive[resume['refined_indices']]=resume['retained_candidate']
        data['merged_refinement']=dict(total=len(status),completed=int((status==0).sum()),retained_positive=int(positive.sum()),
                                       status_counts={str(k):int(v) for k,v in zip(*np.unique(status,return_counts=True))})
    data['scope']='Finite sampled trajectories and finite observation intervals; a positive window screen does not certify an attractor. No inference of an infinite-alpha chaotic family.'
    (OUT/'summary.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(affine_control())
    print({k:{kk:vv for kk,vv in v.items() if kk!='by_slope_alpha'} for k,v in data.items() if isinstance(v,dict)})


if __name__=='__main__':summarize()
