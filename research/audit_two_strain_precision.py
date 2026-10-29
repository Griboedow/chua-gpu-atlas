"""Probe very small populations and compare independent CPU integrators."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from two_strain import ROOT,decode,log_rhs,log_jac,rhs_jac


def run():
    source=np.load(ROOT/'slice-refined.npz');rows=[]
    for index in [17,32,39]:
        p=source['parameters'][index];start=source['final_log_states'][index]
        # This is a short solver comparison, not an independent chaos claim.
        first=solve_ivp(lambda t,q:log_rhs(q,p),(0,1000),start,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.5,dense_output=True)
        second=solve_ivp(lambda t,q:log_rhs(q,p),(0,1000),start,method='Radau',jac=lambda t,q:log_jac(q,p),rtol=1e-11,atol=1e-13,max_step=.5,dense_output=True)
        if not first.success or not second.success:raise RuntimeError('CPU solver failed')
        times=np.linspace(0,1000,2001)
        q1=first.sol(times).T;q2=second.sol(times).T
        y1=np.array([decode(q,p) for q in q1]);y2=np.array([decode(q,p) for q in q2])
        row=dict(refined_index=index,source_index=int(source['source_index'][index]),
            max_log_state_difference=float(np.max(abs(q1-q2))),max_physical_state_difference=float(np.max(abs(y1-y2))),
            minimum_log_prevalence=float(np.min(q1[:,1:3])),minimum_sampled_prevalence=float(y1[:,1:4].min()),
            duration=1000.,DOP853_rtol=1e-12,Radau_rtol=1e-11,max_step=.5)
        rows.append(row)
    # Deliberately tiny invader, within the literal deterministic simplex.
    p=source['parameters'][32];q=np.array([.48,-600.,-3.,.2,.03])
    sol=solve_ivp(lambda t,q:log_rhs(q,p),(0,200),q,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.25)
    direct=solve_ivp(lambda t,y:rhs_jac(y,p)[0],(0,200),decode(q,p),method='DOP853',rtol=1e-12,atol=1e-280,max_step=.25)
    if not sol.success or not direct.success:raise RuntimeError('Tiny-invader control failed')
    tiny=dict(initial_log_I1_over_mu=-600.,final_log_I1_over_mu=float(sol.y[1,-1]),
        direct_log_I1_over_mu=float(np.log(direct.y[1,-1]/p[0])),
        difference=float(abs(sol.y[1,-1]-np.log(direct.y[1,-1]/p[0]))),
        scope='Artificial solver control only; not a physically meaningful finite-population state')
    result=dict(comparisons=rows,tiny_invader=tiny,
        passed=bool(max(r['max_log_state_difference'] for r in rows)<1e-6 and tiny['difference']<1e-6),
        caveat='Pointwise agreement over this interval does not imply pathwise agreement for arbitrarily long chaotic trajectories')
    (ROOT/'precision-audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
    if not result['passed']:raise RuntimeError('Precision comparison needs investigation')


if __name__=='__main__':run()
