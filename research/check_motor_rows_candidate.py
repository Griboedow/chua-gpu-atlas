"""Independent full-spectrum CPU audit of one exact N-row GPU candidate.

Full original reduced state, no cluster reduction. Original normalized time.
Checkpointed QR, divergence and exact flow-vector controls. No novelty claim.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from motor_rows import rhs, field_jac, dimension, unpack, source_force_constraint_field

ROOT=Path(__file__).resolve().parents[1]


def run(n=3,index=783,warm=200.,observe=2000.,wall=600.):
    began=time.perf_counter();source=ROOT/f'results/motor-rows/gpu-scout-N{n}.npz';data=np.load(source)
    p=data['parameters'][index];d=dimension(n);x=data['final_states'][index].copy();initial=x.copy()
    out=ROOT/f'results/motor-rows/cpu-candidate-N{n}-{index}'
    def augmented(t,z):
        f,j=field_jac(z[:d],p,n)
        return np.r_[f,(j@z[d:d+d*d].reshape(d,d)).ravel(),np.trace(j)]
    # Independent algebra and a two-solver state/fundamental-matrix control.
    f,j=field_jac(x,p,n);eye=np.eye(d)
    cj=np.column_stack([np.imag(rhs(0,x.astype(complex)+1e-28j*eye[k],p,n))/1e-28 for k in range(d)])
    independent=source_force_constraint_field(x,p,n)
    c,s,X=unpack(x,n);closed_trace=-2*n-p[0]*(n-1)+p[1]*(1-1/n)*np.sum(c)
    check=dict(jacobian_complex_step_error=float(np.max(abs(j-cj))),
        independent_force_constraint_field_error=float(np.max(abs(f-independent))),
        analytic_trace_identity_error=float(abs(np.trace(j)-closed_trace)))
    ends={}
    for method in ('DOP853','Radau'):
        sol=solve_ivp(augmented,(0,2.),np.r_[x,eye.ravel(),0.],method=method,rtol=2e-12,atol=2e-14,max_step=.01)
        if not sol.success:raise RuntimeError(sol.message)
        ends[method]=sol.y[:,-1]
    original=solve_ivp(lambda t,z:source_force_constraint_field(z,p,n),(0,2.),x,
        method='DOP853',rtol=2e-12,atol=2e-14,max_step=.01)
    if not original.success:raise RuntimeError(original.message)
    check.update(duration=2.,DOP853_Radau_state_error=float(np.max(abs(ends['DOP853'][:d]-ends['Radau'][:d]))),
        DOP853_Radau_variational_error=float(np.max(abs(ends['DOP853'][d:d+d*d]-ends['Radau'][d:d+d*d]))),
        independent_force_constraint_state_error=float(np.max(abs(ends['DOP853'][:d]-original.y[:,-1]))))
    check['passed']=bool(max(check[k] for k in ('jacobian_complex_step_error','independent_force_constraint_field_error',
        'analytic_trace_identity_error','DOP853_Radau_state_error','DOP853_Radau_variational_error','independent_force_constraint_state_error'))<1e-8)
    report=dict(N=n,dimension=d,parameters=p.tolist(),source_row=index,initial_state=initial.tolist(),
        source_npz_exact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),short_control=check,status='running',
        settings=dict(warm=warm,observe=observe,time='Original normalized tau',rtol=1e-11,atol=1e-13,
            maxstep=.03,qr_max=.2,checkpoint_window=100.,wall_budget_seconds=wall),windows=[],
        max_local_trace_error=0.,max_local_flow_identity_error=0.,max_local_condition=0.,
        caveat='Existing exact 3N-1 dimensional model. Finite-time positive growth is not proof of a new attractor; 3D figures are projections.')
    basis=eye;elapsed=0.;logs=np.zeros(d);trace=0.;flowlog=0.;next_dt=.2;chunks=0;rejected=0
    times=[];states=[];sections=[];section_times=[];checkpoints=[];observed=0.;flow_initial=None
    def save():
        report.update(wall_seconds=time.perf_counter()-began,elapsed_original=elapsed,qr_chunks=chunks,rejected_qr_chunks=rejected,
            final_state=x.tolist(),final_basis=basis.tolist(),observation_log_sums=logs.tolist(),integrated_trace=float(trace),
            normalized_flow_log=float(flowlog),next_qr_interval=next_dt)
        out.with_suffix('.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        np.savez_compressed(out.with_suffix('.npz'),parameters=p,time=np.array(times),states=np.array(states),
            section_time=np.array(section_times),section_states=np.array(sections),checkpoint=np.array(checkpoints),
            final_state=x,final_basis=basis,observation_log_sums=logs,integrated_trace=trace,normalized_flow_log=flowlog,
            elapsed_original=elapsed,next_qr_interval=next_dt)
    save()
    if not check['passed']:report['status']='short controls failed';save();return
    phases=[('warm',min(100.,warm-i)) for i in np.arange(0,warm,100.)]+[('observation',min(100.,observe-i)) for i in np.arange(0,observe,100.)]
    for phase,length in phases:
        spent=0.;wlogs=np.zeros(d);wtrace=0.;wflow=0.;nf_start=np.linalg.norm(rhs(0,x,p,n));section_before=len(sections)
        if phase=='observation' and flow_initial is None:flow_initial=nf_start
        while spent<length:
            if time.perf_counter()-began>wall:
                report.update(status='wall budget incomplete',partial_window=dict(phase=phase,duration=spent,
                    log_sums=wlogs.tolist(),trace=float(wtrace),flow_log=float(wflow)));save();return
            dt=min(next_dt,length-spent);f0=rhs(0,x,p,n);n0=np.linalg.norm(f0)
            if n0<1e-14:report['status']='near-equilibrium flow; no neutral estimate';save();return
            def section(t,z):return z[2*n]
            section.direction=1
            sol=solve_ivp(augmented,(0,dt),np.r_[x,basis.ravel(),0.],method='DOP853',rtol=1e-11,atol=1e-13,
                max_step=min(.03,dt),events=section)
            if not sol.success:report.update(status='solver failure',solver_message=sol.message);save();return
            end=sol.y[:d,-1];matrix=sol.y[d:d+d*d,-1].reshape(d,d);new_basis,R=np.linalg.qr(matrix);diag=np.diag(R)
            error=np.inf if np.any(diag==0) else abs(np.log(abs(diag)).sum()-sol.y[-1,-1]);condition=np.linalg.cond(matrix)
            if error>1e-8 or condition>1e10:
                next_dt=dt/2;rejected+=1
                if next_dt<1e-6:report['status']='QR conditioning failure';save();return
                continue
            prop=matrix@(basis.T@(f0/n0));f1=rhs(0,end,p,n)
            identity=np.linalg.norm(prop*n0-f1)/(1+np.linalg.norm(f1));growth=np.log(abs(diag));flow_growth=np.log(np.linalg.norm(prop))
            report['max_local_trace_error']=max(report['max_local_trace_error'],float(error))
            report['max_local_flow_identity_error']=max(report['max_local_flow_identity_error'],float(identity))
            report['max_local_condition']=max(report['max_local_condition'],float(condition))
            x=end;basis=new_basis;elapsed+=dt;spent+=dt;chunks+=1
            if length-spent<1e-11:spent=length
            next_dt=min(.2,next_dt*1.05)
            if phase=='observation':
                logs+=growth;trace+=sol.y[-1,-1];flowlog+=flow_growth;wlogs+=growth;wtrace+=sol.y[-1,-1];wflow+=flow_growth
                for hit,hit_state in zip(sol.t_events[0],sol.y_events[0]):
                    if hit>1e-10:section_times.append(elapsed-dt+float(hit));sections.append(hit_state[:d])
            times.append(elapsed);states.append(x.copy())
        if phase=='observation':
            observed+=length
            w=dict(end_tau=elapsed,duration=length,spectrum=(wlogs/length).tolist(),mean_divergence=float(wtrace/length),
                sum_trace_error=float(abs(wlogs.sum()-wtrace)/length),flow_exponent=float(wflow/length),
                telescopic_flow_exponent=float(np.log(np.linalg.norm(rhs(0,x,p,n))/nf_start)/length),
                upward_X1_section_hits=len(sections)-section_before)
            report['windows'].append(w);report.update(completed_observation_tau=observed,cumulative_spectrum=(logs/observed).tolist(),
                mean_divergence=float(trace/observed),sum_trace_error=float(abs(logs.sum()-trace)/observed),
                flow_exponent=float(flowlog/observed),telescopic_flow_exponent=float(np.log(np.linalg.norm(rhs(0,x,p,n))/flow_initial)/observed),
                upward_X1_section_hits=len(sections))
            print(json.dumps(dict(elapsed=elapsed,spectrum=report['cumulative_spectrum'],sections=len(sections),
                trace_error=report['sum_trace_error'],wall=round(time.perf_counter()-began,1))),flush=True)
        checkpoints.append(np.r_[elapsed,x,basis.ravel(),logs,trace,flowlog,next_dt]);save()
    sample=np.array(states);mask=np.array(times)>warm;c=sample[mask,:n];s=sample[mask,n:2*n]
    xx=np.c_[sample[mask,2*n:],-sample[mask,2*n:].sum(axis=1)]
    pair=[]
    for i in range(n):
        for j in range(i+1,n):
            distance=np.sqrt((c[:,i]-c[:,j])**2+(s[:,i]-s[:,j])**2+(xx[:,i]-xx[:,j])**2)
            pair.append(dict(rows=[i+1,j+1],rms_distance=float(np.sqrt(np.mean(distance**2))),
                min_distance=float(distance.min()),max_distance=float(distance.max())))
    report.update(status='complete',sampled_domain=dict(max_density_disk_radius_squared=float(np.max(c*c+s*s)),
        invariant_disk_bound=float(1+p[2]**2),max_abs_sum_X=float(np.max(abs(xx.sum(axis=1))))),
        row_pair_distances=pair,clustering_caveat='Small instantaneous distance alone is not invariant row synchronization.')
    save()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,default=3);parser.add_argument('--index',type=int,default=783)
    parser.add_argument('--warm',type=float,default=200.);parser.add_argument('--observe',type=float,default=2000.)
    parser.add_argument('--wall',type=float,default=600.);a=parser.parse_args();run(a.n,a.index,a.warm,a.observe,a.wall)
