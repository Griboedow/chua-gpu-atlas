"""Numerical fold candidates: periodic shooting plus a transverse unit multiplier."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import root
from continue_synchronematic_orbit import Orbit,record


def main():
    path=Path('results/synchronematic/n4-periodic-branch')
    data=json.loads((path/'branch.json').read_text());analysis=json.loads((path/'bifurcations.json').read_text())
    model=Orbit(data['positions']);rows=data['records'];result=[]
    # Keep the diagnostic from the exploratory fixed-parameter run explicit.
    data['large_adjacent_period_jumps']=[dict(lower_eta=a['eta'],upper_eta=b['eta'],period_jump=abs(a['period']-b['period']))
        for a,b in zip(rows,rows[1:]) if abs(a['period']-b['period'])>.1]
    data['caveat']='Fixed-parameter Newton converged to a different orbit below the lower fold. Flagged jumps must not be plotted as one connected branch. Tiny monodromy multipliers are unresolved.'
    (path/'branch.json').write_text(json.dumps(data,indent=2))
    for name,seed in [('lower',analysis['lower_fold_sample']),('upper',max(rows,key=lambda r:r['eta']))]:
        ref=np.array(seed['state']);normal=model.field(ref,seed['eta'])[0];normal/=np.linalg.norm(normal)
        basis=null_space(normal[None,:]);w=np.array(seed['winding']);calls=0
        def fun(z):
            nonlocal calls
            calls+=1;q=z[:7];T=z[7];eta=z[8]/100
            end,M,_=model.integrate(q,T,eta)
            f=model.field(end,eta)[0]
            P=basis.T@(np.eye(7)-np.outer(f,normal)/(normal@f))@M@basis
            eig=np.linalg.eigvals(P);critical=eig[np.argmin(abs(eig-1))]
            if abs(critical.imag)>1e-6:raise RuntimeError('Fold multiplier became complex')
            return np.r_[end-q-2*np.pi*w,normal@(q-ref),critical.real-1]
        z0=np.r_[ref,seed['period'],100*seed['eta']]
        solution=root(fun,z0,method='hybr',options={'xtol':1e-9,'eps':1e-12,'maxfev':160})
        z=solution.x;res=float(np.max(abs(fun(z))))
        if res>2e-8:raise RuntimeError(f'{name} fold solve failed: {solution.message}, residual {res}')
        end,M,_=model.integrate(z[:7],z[7],z[8]/100)
        row=record(model,z[:8],z[8]/100,M,res,calls,w);row['boundary']=name
        row['root_message']=str(solution.message);row['condition']='periodic closure, fixed phase, transverse Floquet multiplier +1'
        independent=Orbit(data['positions'],rtol=1e-12,max_step=.025,method='Radau')
        e2,M2,_=independent.integrate(z[:7],z[7],z[8]/100)
        row['radau_closure_error']=float(np.max(abs(e2-z[:7]-2*np.pi*w)))
        row['radau_monodromy_difference']=float(np.max(abs(M2-M)))
        f2=independent.field(e2,z[8]/100)[0];P2=basis.T@(np.eye(7)-np.outer(f2,normal)/(normal@f2))@M2@basis
        eigen2=np.linalg.eigvals(P2);row['radau_transverse_unit_error']=float(min(abs(eigen2-1)))
        result.append(row)
        (path/'folds.json').write_text(json.dumps(dict(changes_equations=True,folds=result,
            caveat='Numerical fold candidates, not interval-validated bifurcations or proof of generic nondegeneracy.'),indent=2))
        print(name,'fold',row['eta'],'period',row['period'],'residual',res,flush=True)


if __name__=='__main__':main()
