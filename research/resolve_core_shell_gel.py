"""Stiff CPU continuation of GPU step-budget cases; no chaos inference."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from core_shell_gel import rhs_jac,ROOT,save


def resolve(file,duration=100.):
    old=np.load(file);indices=np.flatnonzero(old['status']==3);rows=[];finals=[];params=[]
    for i in indices:
        p=old['parameters'][i];initial=old['final_states'][i];wall=time.perf_counter()
        def lower(t,y):return y[3]
        def upper(t,y):return 1-y[3]
        lower.terminal=upper.terminal=True;lower.direction=upper.direction=-1
        sol=solve_ivp(lambda t,y:rhs_jac(y,p)[0],(0,duration),initial,method='Radau',jac=lambda t,y:rhs_jac(y,p)[1],rtol=1e-9,atol=1e-11,max_step=.5,events=[lower,upper],dense_output=True)
        y=sol.y[:,-1];f,j=rhs_jac(y,p);eig=np.linalg.eigvals(j)
        late=sol.sol(np.linspace(max(0,sol.t[-1]-duration/4),sol.t[-1],2001)).T
        residual=float(np.max(abs(f)/(1+abs(y))));variation=float(np.max(np.ptp(late,axis=0)/(1+abs(y))))
        outcome='unresolved'
        if sol.success and sol.status==1:outcome='physical_domain_exit'
        elif sol.success and sol.t[-1]>=duration and residual<1e-8 and variation<1e-7 and np.max(eig.real)<-1e-5:outcome='locally_stable_equilibrium_observed'
        row=dict(source_index=int(i),outcome=outcome,success=bool(sol.success),elapsed=float(sol.t[-1]),scaled_residual=residual,scaled_late_variation=variation,eigenvalues=[[float(z.real),float(z.imag)] for z in eig],wall_seconds=time.perf_counter()-wall)
        rows.append(row);params.append(p);finals.append(y);print(row,flush=True)
        save(ROOT/(Path(file).stem+'-stiff-continuation'),dict(source=Path(file).name,source_indices=indices[:len(rows)],parameters=np.array(params),final_states=np.array(finals),rows=rows,method='Independent Radau, analytic Jacobian',rtol=1e-9,atol=1e-11,max_step=.5,additional_duration=duration,scope='Resolve stiff GPU budget cases by observed physical exit or stable equilibrium; oscillatory cases remain unresolved'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',default=str(ROOT/'scout.npz'));parser.add_argument('--duration',type=float,default=100.);a=parser.parse_args();resolve(a.source,a.duration)
