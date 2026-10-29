"""Follow the doubled orbit and distinguish it from a double-covered parent."""
import json
from pathlib import Path
import numpy as np
from continue_synchronematic_orbit import Orbit,record


def main():
    path=Path('results/synchronematic/n4-periodic-branch')
    data=json.loads((path/'branch.json').read_text());checks=json.loads((path/'checked-orbits.json').read_text())
    seed=next(r for r in checks['orbits'] if r['eta']==1.078);model=Orbit(data['positions']);w=np.array(seed['winding'])
    records=[seed];failures=[]
    def save():
        (path/'doubled-branch.json').write_text(json.dumps(dict(changes_equations=True,records=sorted(records,key=lambda r:r['eta']),
            failures=failures,caveat='Numerical continuation. Nonzero half-period distance rejects a double cover of the parent. No claim of generic bifurcation proof.'),indent=2))
    for sign,limit in [(-1,1.07694),(1,1.08018)]:
        eta=1.078;u=np.r_[seed['state'],seed['period']]
        while sign*(limit-eta)>1e-10:
            trial=eta+sign*min(.0001,abs(limit-eta))
            try:new,M,error,it=model.shoot(u,trial,u[:7],w,tol=2e-10)
            except (RuntimeError,ValueError,np.linalg.LinAlgError) as e:
                failures.append(dict(eta=trial,message=str(e)));save();break
            half,_,_=model.integrate(new[:7],new[7]/2,trial)
            distance=float(np.linalg.norm(np.angle(np.exp(1j*(half-new[:7])))))
            if distance<1e-5:
                failures.append(dict(eta=trial,message='Converged to double-covered parent; not retained as doubled branch'));save();break
            row=record(model,new,trial,M,error,it,w);row['half_period_circular_distance']=distance
            records.append(row);save();eta=trial;u=new
            print('Doubled',eta,u[7],row['spectral_radius'],distance,flush=True)
    save()


if __name__=='__main__':main()
