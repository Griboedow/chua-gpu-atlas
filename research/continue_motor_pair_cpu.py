"""Recover the previously unsaved QR frame, then extend one motor-pair audit.

Only the saved-endpoint 100+400 tau control is replayed. Its endpoint and QR
logs must match the immutable earlier result before any extension is accepted.
Every new 100-tau checkpoint retains a complete restartable QR state.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from check_motor_pair_gpu_candidate import field, physical

ROOT=Path(__file__).resolve().parents[1]


def run(index=444,duration=3000.,wall=1000.):
    started=time.perf_counter();source=ROOT/f'results/motor-pair/gpu-cpu-{index}.json'
    earlier=json.loads(source.read_text());reference=earlier['runs'][0]
    if reference['label']!='saved_endpoint' or reference['status']!='complete':raise RuntimeError('Completed original control required')
    p=np.array(earlier['parameters']);x=np.array(reference['start_balanced']);basis=np.eye(3);elapsed=0.;dt_next=.05
    logs=np.zeros(3);trace=0.;flowlog=0.;obs_initial_flow=None;chunks=0;rejected=0
    out=ROOT/f'results/motor-pair/gpu-cpu-{index}-extended';states=[];times=[];sections=[];section_times=[];checkpoints=[]
    report=dict(parameters=p.tolist(),source_id=index,source_control_exact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        frame_recovery='The first audit did not retain the final QR frame. Deterministically replay its saved-endpoint 100+400 tau path and verify state and QR logs before extension; no grid is rerun.',
        settings=dict(replayed_warm_tau=100.,replayed_observation_tau=400.,additional_observation_tau=duration,
            window_tau=100.,qr_max=.05,rtol=1e-11,atol=1e-13,maxstep=.01,wall_budget_seconds=wall),
        coordinate_norm='Euclidean Q=h/r-Y,Y,Z; original tau time',status='running',windows=[],
        max_local_trace_error=0.,max_local_flow_identity_error=0.,
        qualification='Finite-time evidence; convergence of the neutral exponent and adequate recurrence count are assessed, not assumed.')
    def save():
        report.update(wall_seconds=time.perf_counter()-started,elapsed_original=elapsed,qr_chunks=chunks,rejected_qr_chunks=rejected,
            final_balanced_state=x.tolist(),final_original_state=physical(x,p).tolist(),final_basis=basis.tolist(),
            observation_log_sums=logs.tolist(),integrated_trace=float(trace),integrated_normalized_flow_log=float(flowlog),
            next_qr_interval=dt_next)
        out.with_suffix('.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        np.savez_compressed(out.with_suffix('.npz'),parameters=p,time=np.array(times),balanced_states=np.array(states),
            section_time=np.array(section_times),section_states=np.array(sections),checkpoint=np.array(checkpoints),
            final_state=x,final_basis=basis,observation_log_sums=logs,integrated_trace=trace,
            normalized_flow_log=flowlog,next_qr_interval=dt_next,elapsed_original=elapsed)
    phases=[('warm',100.)]+[('replay',100.)]*4+[('extension',min(100.,duration-i)) for i in np.arange(0,duration,100.)]
    for phase,length in phases:
        spent=0.;wlogs=np.zeros(3);wtrace=0.;wflow=0.;flow_start=np.linalg.norm(field(x,p)[0]);start_section=len(sections)
        if phase!='warm' and obs_initial_flow is None:obs_initial_flow=flow_start
        while spent<length:
            if time.perf_counter()-started>wall:
                report.update(status='wall budget incomplete',partial_window=dict(phase=phase,duration=spent,
                    log_sums=wlogs.tolist(),trace=float(wtrace),flow_log=float(wflow)))
                save();return
            dt=min(dt_next,length-spent);f0=field(x,p)[0];n0=np.linalg.norm(f0)
            def aug(t,z):
                f,j=field(z[:3],p)
                return np.r_[f,(j@z[3:12].reshape(3,3)).ravel(),np.trace(j)]
            def section(t,z):return z[0]+z[1]
            section.direction=1
            sol=solve_ivp(aug,(0,dt),np.r_[x,basis.ravel(),0.],method='DOP853',rtol=1e-11,atol=1e-13,
                max_step=min(.01,dt),events=section)
            if not sol.success:
                report.update(status='solver failure',solver_message=sol.message);save();return
            end=sol.y[:3,-1];matrix=sol.y[3:12,-1].reshape(3,3);new_basis,R=np.linalg.qr(matrix);diag=np.diag(R)
            local_trace=np.inf if np.any(diag==0) else abs(np.log(abs(diag)).sum()-sol.y[-1,-1])
            if local_trace>1e-8 or np.linalg.cond(matrix)>1e10:
                dt_next=dt/2;rejected+=1
                if dt_next<1e-6:report['status']='QR conditioning failure';save();return
                continue
            prop=matrix@(basis.T@(f0/n0));f1=field(end,p)[0]
            identity=np.linalg.norm(prop*n0-f1)/(1+np.linalg.norm(f1))
            growth=np.log(abs(diag));flow_growth=np.log(np.linalg.norm(prop))
            report['max_local_trace_error']=max(report['max_local_trace_error'],float(local_trace))
            report['max_local_flow_identity_error']=max(report['max_local_flow_identity_error'],float(identity))
            x=end;basis=new_basis;elapsed+=dt;spent+=dt;chunks+=1
            if length-spent<1e-11:spent=length
            dt_next=min(.05,dt_next*1.05)
            if phase!='warm':
                logs+=growth;trace+=sol.y[-1,-1];flowlog+=flow_growth
                wlogs+=growth;wtrace+=sol.y[-1,-1];wflow+=flow_growth
                for hit,state in zip(sol.t_events[0],sol.y_events[0]):
                    if hit>1e-10:section_times.append(elapsed-dt+float(hit));sections.append(state[:3])
            times.append(elapsed);states.append(x.copy())
        if phase!='warm':
            w=dict(phase=phase,end_tau=elapsed,duration=length,spectrum=(wlogs/length).tolist(),
                mean_divergence=float(wtrace/length),sum_trace_error=float(abs(wlogs.sum()-wtrace)/length),
                flow_exponent=float(wflow/length),telescopic_flow_exponent=float(np.log(np.linalg.norm(field(x,p)[0])/flow_start)/length),
                upward_section_hits=len(sections)-start_section)
            report['windows'].append(w)
            observed=sum(v['duration'] for v in report['windows'])
            report.update(cumulative_observed_tau=observed,cumulative_spectrum=(logs/observed).tolist(),
                cumulative_sum_trace_error=float(abs(logs.sum()-trace)/observed),mean_divergence=float(trace/observed),
                cumulative_flow_exponent=float(flowlog/observed),
                cumulative_telescopic_flow_exponent=float(np.log(np.linalg.norm(field(x,p)[0])/obs_initial_flow)/observed),
                observed_section_hits=len(sections))
            print(json.dumps(dict(phase=phase,elapsed=elapsed,spectrum=report['cumulative_spectrum'],
                window=w['spectrum'],sections=len(sections),wall=round(time.perf_counter()-started,1))),flush=True)
        if phase=='replay' and len(report['windows'])==4:
            err=float(max(abs(x-np.array(reference['final_balanced_state']))))
            logs_err=float(max(abs(logs/400-np.array(reference['spectrum']))))
            report['replay_verification']=dict(endpoint_max_error=err,spectrum_max_error=logs_err,passed=err<1e-10 and logs_err<1e-10)
            if not report['replay_verification']['passed']:
                report['status']='replay mismatch; extension not started';save();return
        checkpoints.append(np.r_[elapsed,x,basis.ravel(),logs,trace,flowlog,dt_next])
        save()
    report['status']='complete'
    extra=[v for v in report['windows'] if v['phase']=='extension'];total=sum(v['duration'] for v in extra)
    report['extension_spectrum']=np.average([v['spectrum'] for v in extra],weights=[v['duration'] for v in extra],axis=0).tolist()
    report['extension_completed_tau']=total
    sample=np.array(states);report['sampled_domain']=dict(maximum_abs_Q=float(np.max(abs(sample[:,0]))),
        maximum_disk_radius_squared=float(np.max(sample[:,1]**2+(1-sample[:,2])**2)))
    save()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--index',type=int,default=444)
    parser.add_argument('--duration',type=float,default=3000.);parser.add_argument('--wall',type=float,default=1000.)
    a=parser.parse_args();run(a.index,a.duration,a.wall)
