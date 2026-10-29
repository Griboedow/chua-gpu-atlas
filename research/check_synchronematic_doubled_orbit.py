"""Independent shooting check that a regular window contains a doubled orbit."""
import json
from pathlib import Path
import numpy as np
from continue_synchronematic_orbit import Orbit,record


def main():
    path=Path('results/synchronematic/n4-periodic-branch')
    branch=json.loads((path/'branch.json').read_text());old=np.load('results/synchronematic/n4-transition/scan.npz')
    model=Orbit(branch['positions']);results=[]
    for eta in [1.078,1.076,1.082]:
        row=min(branch['records'],key=lambda r:abs(r['eta']-eta));w=np.array(row['winding'])
        ids=np.flatnonzero(abs(old['experiment_labels'][:,0]-eta)<1e-8)
        multiplier=2 if eta==1.078 else 1
        candidates=[]
        for i in ids:
            q=np.angle(np.exp(1j*(model.D@old['final_states'][i])))
            end,_,_=model.integrate(q,multiplier*row['period'],eta)
            candidates.append((np.linalg.norm(end-q-multiplier*2*np.pi*w),int(i),q))
        candidates.sort(key=lambda x:x[0]);distance,index,start=candidates[0]
        u,M,error,it=model.shoot(np.r_[start,multiplier*row['period']],eta,start,multiplier*w,tol=1e-11)
        output=record(model,u,eta,M,error,it,multiplier*w)
        half,_,_=model.integrate(u[:7],u[7]/2,eta)
        output['half_period_circular_distance']=float(np.linalg.norm(np.angle(np.exp(1j*(half-u[:7])))))
        output['source_transition_index']=index;output['initial_period_multiple']=multiplier
        output['seed_return_residuals']=[dict(index=i,residual=float(d)) for d,i,q in candidates]
        independent=Orbit(branch['positions'],rtol=2e-12,max_step=.025,method='Radau')
        end2,M2,_=independent.integrate(u[:7],u[7],eta)
        output['radau_closure_error']=float(np.max(abs(end2-u[:7]-multiplier*2*np.pi*w)))
        output['radau_monodromy_difference']=float(np.max(abs(M2-M)))
        times=np.linspace(0,u[7],501)
        # Only a short orbit trace is needed for the scientific figure.
        from scipy.integrate import solve_ivp
        sol=solve_ivp(lambda t,q:model.field(q,eta)[0],(0,u[7]),u[:7],t_eval=times,
                      method='DOP853',rtol=2e-12,atol=2e-14,max_step=.025)
        if not sol.success:raise RuntimeError(sol.message)
        np.savez_compressed(path/f'orbit-{eta:.3f}.npz',time=times,state=sol.y.T,eta=eta)
        results.append(output)
        (path/'checked-orbits.json').write_text(json.dumps(dict(changes_equations=True,orbits=results,
            caveat='Shooting/Floquet evidence in the explicit homotopy; not a validated mathematical proof or device calibration.'),indent=2))
        print(eta,'period',u[7],'radius',output['spectral_radius'],'half-distance',output['half_period_circular_distance'],flush=True)


if __name__=='__main__':main()
