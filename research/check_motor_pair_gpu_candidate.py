"""Independent bounded CPU full-spectrum audit of a saved exact-model scout.

Uses original transition-rate time tau, continuous QR with warmed bases,
trace quadrature and the exact flow-vector variational identity. Each complete
window is saved; a wall cutoff stays explicit. No positive FTLE is called a new
attractor. Original-start and resumed-endpoint controls are kept separate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from motor_pair import oriola_field, to_oriola, from_oriola

ROOT=Path(__file__).resolve().parents[1]


def field(x,p):
    q,y,z=x;rho,r=p;h=q+y
    return np.array([-rho*h,r*h*(1-z)-y,r*h*y-z]),np.array([
        [-rho,-rho,0.],[r*(1-z),r*(1-z)-1.,-r*h],[r*y,r*(q+2*y),-1.]])


def physical(x,p):
    return np.array([p[1]*(x[0]+x[1]),x[1],x[2]])


def balanced(y,p):
    return np.array([y[0]/p[1]-y[1],y[1],y[2]])


def short_control(start,p,duration=2.):
    def aug(t,z):
        f,j=field(z[:3],p)
        return np.r_[f,(j@z[3:].reshape(3,3)).ravel()]
    answers={}
    for method in ['DOP853','Radau']:
        sol=solve_ivp(aug,(0,duration),np.r_[start,np.eye(3).ravel()],method=method,
            rtol=2e-12,atol=2e-14,max_step=.002)
        if not sol.success:raise RuntimeError(sol.message)
        answers[method]=sol.y[:,-1]
    omega=.2;rho,r=p;mu=r/omega
    initial=to_oriola(physical(start,p),rho,r,omega)
    original=solve_ivp(lambda t,z:oriola_field(z,rho,mu,omega),(0,duration),initial,
        method='DOP853',rtol=2e-12,atol=2e-14,max_step=.002)
    if not original.success:raise RuntimeError(original.message)
    end=balanced(from_oriola(original.y[:,-1],rho,mu,omega),p)
    result=dict(duration_original=duration,DOP853_Radau_state_error=float(max(abs(answers['DOP853'][:3]-answers['Radau'][:3]))),
        DOP853_Radau_monodromy_error=float(max(abs(answers['DOP853'][3:]-answers['Radau'][3:]))),
        independent_Oriola_equation7_state_error=float(max(abs(end-answers['DOP853'][:3]))),
        method='DOP853/Radau in balanced memory coordinates, plus independent original equation (7), rtol2e-12,atol2e-14,maxstep.002')
    result['passed']=bool(result['DOP853_Radau_state_error']<1e-7 and result['DOP853_Radau_monodromy_error']<1e-6
        and result['independent_Oriola_equation7_state_error']<1e-7)
    return result


def run(index=444,warm=100.,observe=400.,wall_budget=300.,qr=.05):
    began=time.perf_counter();source=ROOT/'results/motor-pair/gpu-scout.npz';data=np.load(source)
    p=data['parameters'][index];out=ROOT/f'results/motor-pair/gpu-cpu-{index}'
    starts=[('saved_endpoint',data['balanced_states'][index]),('original_start',balanced(data['initial_states'][index],p))]
    results=[];all_times=[];all_states=[];all_labels=[];all_sections=[];all_section_times=[];all_section_labels=[]
    report=dict(source_npz_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_id=index,
        parameters=p.tolist(),time='Original dimensionless transition-rate time tau; no slowtime rescaling',
        coordinate_norm='Euclidean Q=h/r-Y,Y,Z',settings=dict(warm=warm,observe=observe,qr_max=qr,
            rtol=1e-11,atol=1e-13,max_step=.01,wall_budget=wall_budget),runs=results,
        caveat='Finite-time full-spectrum and independent solver evidence; not proof of an attractor or priority.')
    def save():
        report['wall_seconds']=time.perf_counter()-began
        out.with_suffix('.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        np.savez_compressed(out.with_suffix('.npz'),time=np.array(all_times),balanced_states=np.array(all_states),
            run_index=np.array(all_labels),section_time=np.array(all_section_times),section_states=np.array(all_sections),
            section_run_index=np.array(all_section_labels),parameters=p)
    for number,(label,start) in enumerate(starts):
        if time.perf_counter()-began>wall_budget-10:break
        check=short_control(start,p)
        item=dict(label=label,start_balanced=start.tolist(),short_control=check,status='running',windows=[],
            max_local_trace_error=0.,max_local_flow_identity_error=0.,rejected_qr_chunks=0)
        results.append(item);save()
        if not check['passed']:
            item['status']='short independent control failed';save();continue
        x=start.copy();basis=np.eye(3);elapsed=0.;logs=np.zeros(3);trace=0.;flowlog=0.;count=0
        initial_observed_flow=None;window_start_flow=None;window_logs=np.zeros(3);window_trace=0.;window_flowlog=0.
        candidate_dt=qr;phases=[warm]+[observe/4]*4
        for phase,length in enumerate(phases):
            phase_elapsed=0.;window_logs[:]=0.;window_trace=0.;window_flowlog=0.
            if phase>0:
                window_start_flow=np.linalg.norm(field(x,p)[0])
                if phase==1:initial_observed_flow=window_start_flow
            while phase_elapsed<length:
                if time.perf_counter()-began>wall_budget:
                    item['status']='wall budget incomplete';break
                dt=min(candidate_dt,length-phase_elapsed)
                fstart=field(x,p)[0];flow_start_norm=np.linalg.norm(fstart)
                def augmented(t,z):
                    f,j=field(z[:3],p)
                    return np.r_[f,(j@z[3:12].reshape(3,3)).ravel(),np.trace(j)]
                def section(t,z):return z[0]+z[1]
                section.direction=1
                sol=solve_ivp(augmented,(0,dt),np.r_[x,basis.ravel(),0.],method='DOP853',
                    rtol=1e-11,atol=1e-13,max_step=min(.01,dt),events=section)
                if not sol.success:raise RuntimeError(sol.message)
                end=sol.y[:3,-1];matrix=sol.y[3:12,-1].reshape(3,3)
                new_basis,R=np.linalg.qr(matrix);diag=np.diag(R)
                if np.any(diag==0):local_trace_error=np.inf
                else:local_trace_error=abs(np.log(abs(diag)).sum()-sol.y[-1,-1])
                if local_trace_error>1e-8 or np.linalg.cond(matrix)>1e10:
                    candidate_dt=dt/2;item['rejected_qr_chunks']+=1
                    if candidate_dt<1e-6:raise RuntimeError('Unable to control local QR conditioning')
                    continue
                propagated_flow=matrix@(basis.T@(fstart/flow_start_norm)) if flow_start_norm>0 else np.zeros(3)
                fend=field(end,p)[0]
                identity=np.linalg.norm(propagated_flow*flow_start_norm-fend)/(1+np.linalg.norm(fend))
                growth=np.log(abs(diag));flow_growth=np.log(np.linalg.norm(propagated_flow)) if flow_start_norm>0 else 0.
                item['max_local_trace_error']=max(item['max_local_trace_error'],float(local_trace_error))
                item['max_local_flow_identity_error']=max(item['max_local_flow_identity_error'],float(identity))
                x=end;basis=new_basis;elapsed+=dt;phase_elapsed+=dt;count+=1
                if length-phase_elapsed<1e-11:phase_elapsed=length
                candidate_dt=min(qr,candidate_dt*1.05)
                if phase>0:
                    logs+=growth;trace+=sol.y[-1,-1];flowlog+=flow_growth
                    window_logs+=growth;window_trace+=sol.y[-1,-1];window_flowlog+=flow_growth
                    for hit_time,hit_state in zip(sol.t_events[0],sol.y_events[0]):
                        if hit_time>1e-10:
                            all_section_times.append(elapsed-dt+float(hit_time));all_sections.append(hit_state[:3]);all_section_labels.append(number)
                all_times.append(elapsed);all_states.append(x.copy());all_labels.append(number)
            if item['status']=='wall budget incomplete':break
            if phase==0:
                item['warm_completed']=True;item['warm_final_state']=x.tolist();save()
            else:
                endpoint_flow=np.linalg.norm(field(x,p)[0])
                exact_flow=np.log(endpoint_flow/window_start_flow)/length
                value=dict(duration=length,spectrum=(window_logs/length).tolist(),mean_divergence=float(window_trace/length),
                    spectrum_sum_trace_error=float(abs(window_logs.sum()-window_trace)/length),
                    normalized_flow_exponent=float(window_flowlog/length),telescopic_flow_exponent=float(exact_flow),
                    flow_exponent_discrepancy=float(abs(window_flowlog/length-exact_flow)))
                item['windows'].append(value);save()
                print(json.dumps(dict(label=label,window=phase,**value,wall=round(time.perf_counter()-began,1))),flush=True)
        observation_done=sum(w['duration'] for w in item['windows'])
        item.update(elapsed_original=elapsed,qr_chunks=count,final_balanced_state=x.tolist(),
            final_original_state=physical(x,p).tolist())
        if len(item['windows'])==4:
            item['status']='complete';item['spectrum']=(logs/observe).tolist();item['mean_divergence']=float(trace/observe)
            item['spectrum_sum_trace_error']=float(abs(logs.sum()-trace)/observe)
            item['flow_exponent']=float(flowlog/observe)
            item['telescopic_flow_exponent']=float(np.log(np.linalg.norm(field(x,p)[0])/initial_observed_flow)/observe)
        else:
            item['completed_observation_windows_duration']=observation_done
            item['partial_phase_time']=phase_elapsed
        save()
    save()
    print(json.dumps(dict(parameters=p.tolist(),runs=[{'label':r['label'],'status':r['status'],'spectrum':r.get('spectrum')} for r in results],
        wall=report['wall_seconds']),indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--index',type=int,default=444)
    parser.add_argument('--warm',type=float,default=100.);parser.add_argument('--observe',type=float,default=400.)
    parser.add_argument('--wall',type=float,default=300.);parser.add_argument('--qr',type=float,default=.05)
    a=parser.parse_args();run(a.index,a.warm,a.observe,a.wall,a.qr)
