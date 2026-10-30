"""Finite survey along the published chemical-activity path, all constants fixed.

Alouges et al., arXiv:2409.03506v2, Eqs.(15),(16),(65), Table1 and Section3.5.
Uses the existing validated N-row GPU implementation without modifying it.
All requested durations are original normalized tau, not a slow-time clock.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from motor_rows import source_parameters, source_start

ROOT=Path(__file__).resolve().parents[1]
SEEDS=[20251030,20251031,20251032]


def mapping():
    reference=source_parameters(0.)
    low=source_parameters(-.1);high=source_parameters(.1)
    rho=np.sort(np.r_[np.linspace(1.,4.7,45),low['rho'],high['rho'],reference['rho']])
    assert len(rho)==48 and len(np.unique(rho))==48
    stiffness,eta,prefactor=9.5e-5,1e-7,1000.
    d=reference['d_source_convention'];activity=(stiffness/(eta*prefactor*rho))**2
    a0=stiffness/(eta*rho);A=d*np.sqrt(activity)/prefactor
    positivity=.5-np.sqrt(2)*abs(A)
    records=[]
    for index,rr in enumerate(rho):
        reconstructed=source_parameters(activity[index]/reference['activity_nN_nm']-1.)
        records.append(dict(tuple_id=index,rho=float(rr),r=float(reference['r']),b=1.,
            Omega_source_numeric_nN_nm=float(activity[index]),delta=float(activity[index]/reference['activity_nN_nm']-1.),
            a0_per_second=float(a0[index]),A=float(A[index]),B=float(A[index]),
            minimum_transition_rate_over_a0=float(positivity[index]),
            invariant_density_lower_bound=float(positivity[index]),invariant_density_upper_bound=float(1-positivity[index]),
            warm_physical_seconds=float(200/a0[index]),observed_physical_seconds=float(400/a0[index]),
            distance_to_linear_onset=float(rr-(reference['r']-1.)),
            reconstructed_rho_error=float(abs(reconstructed['rho']-rr)),
            reconstructed_r_error=float(abs(reconstructed['r']-reference['r']))))
    report=dict(source=dict(title='Some mathematical models for flagellar activation mechanisms',
        version='arXiv:2409.03506v2',url='https://arxiv.org/html/2409.03506v2',
        references=['Eq.(15): constrained N-row model','Eq.(16): first-harmonic transition rates',
            'Section3.1 and Eq.(65): a0=c sqrt(Omega), a1=b1=d Omega','Table1: constants','Section3.5/Fig9: N=8 starts and controls']),
        constants=dict(ell_m=1e-8,kBT_nN_nm=.0042668,eta_kg_per_second=eta,k_kg_per_second_squared=stiffness,
            U_joule=.042668e-18,c_source_prefactor=prefactor,Omega0_source_numeric_nN_nm=.064002,d_source_convention=float(d)),
        equations=dict(Omega='(k/(eta*c*rho))^2',a0='k/(eta*rho)',A_and_B='d*k/(eta*c^2*rho)',
            r='-2*pi^2*U*d/(eta*ell^2*c^2)',b='1',rate_nonnegativity='sqrt(A^2+B^2)<=1/2'),
        fixed_r=float(reference['r']),critical_rho=float(reference['r']-1.),
        rho_grid='45 equally spaced rho in [1,4.7] plus exact delta=-.1,+.1,0 controls:48 distinct values',
        source_controls=dict(below=low,above=high,critical=reference),
        positivity_all_passed=bool(np.all(positivity>0)),
        rho_rate_nonnegativity_lower_bound=float(2*np.sqrt(2)*abs(d)*stiffness/(eta*prefactor**2)),
        tuples=records,seeds=SEEDS,initial_prescription='Source equilibrium Fourier densities c=s=1; seeded absolute interior shifts uniform[-.01,.01]nm with fixed endpoints; take differences and normalize by 2pi/ell.',
        scope='Only chemical activity varies; all Table1 constants and Eq65 d stay fixed. The interval extends well beyond illustrated delta=+/-0.1 and is not an experimental calibration.',
        time='Original normalized tau=a0*t. For every lane s=tau; physical seconds differ with activity.',
        control_scope='Exact source control tuples are included, but source-style seeds are our explicit reproducibility choices.',
        script_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest())
    assert report['positivity_all_passed']
    return rho,report


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    m=json.loads((folder/'mapping.json').read_text());fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,n in enumerate((3,8)):
        path=folder/f'N{n}.npz'
        if not path.exists():continue
        a=np.load(path);valid=a['status']==0
        for seed in range(3):
            subset=(a['start_id']==seed)&valid
            axes[0,col].plot(a['parameters'][subset,0],a['lle'][subset],'o-',ms=3,lw=.7,label=f'Seed {SEEDS[seed]}')
            amp=(a['late_max'][:,2*n:]-a['late_min'][:,2*n:]).max(axis=1)
            axes[1,col].plot(a['parameters'][subset,0],amp[subset],'o-',ms=3,lw=.7,label=f'Seed {SEEDS[seed]}')
        for ax in axes[:,col]:
            ax.axvline(m['critical_rho'],c='black',ls='--',lw=.8,label='Linear onset' if ax is axes[0,col] else None)
            for key in ('below','above'):ax.axvline(m['source_controls'][key]['rho'],c='gray',ls=':',lw=.7)
            ax.set_xlabel('rho (activity increases toward the left)');ax.grid(alpha=.15)
        axes[0,col].axhline(0,c='gray',lw=.7);axes[0,col].set_title(f'N={n}: {valid.sum()}/144 completed');axes[0,col].legend(fontsize=7)
        axes[0,col].set_ylabel('Leading finite-time exponent per tau')
        axes[1,col].set_ylabel('Largest retained-X peak-to-peak range\n(last 100 tau; omitted last X not used)')
    fig.suptitle('Published activity path, fixed r=4.75514605 and b=1\n200 discarded +400 observed tau; screening evidence, not an asymptotic phase diagram',fontsize=13)
    fig.savefig(folder/'activity-path.png',dpi=170);plt.close(fig)


def run(folder,wall=600.,prepare_only=False):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);rho,m=mapping()
    (folder/'mapping.json').write_text(json.dumps(m,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if prepare_only:return
    validation=ROOT/'results/motor-rows/gpu-validation.json'
    if not json.loads(validation.read_text())['all_passed']:raise RuntimeError('Existing source and short GPU controls must pass')
    from motor_rows_gpu import MotorRowsGPU,save
    began=time.perf_counter();summary=dict(mapping='mapping.json',source_gpu_validation_exact_sha256=hashlib.sha256(validation.read_bytes()).hexdigest(),
        total_wall_budget=wall,settings=dict(transient=200.,duration=400.,time_mode='original',rtol=1e-9,
            step_factor=1.,max_steps_per_lane=8000000),runs=[],
        qualification='A finite source-activity-path screen. Subthreshold, near-onset and incomplete trajectories do not prove absence of chaos.')
    for n in (3,8):
        remaining=wall-(time.perf_counter()-began)
        if remaining<=1:
            summary['runs'].append(dict(N=n,status='not started: total wall budget'));break
        parameters=np.repeat(np.c_[rho,np.full(48,m['fixed_r']),np.ones(48)],3,axis=0)
        source_starts=[source_start(n,seed) for seed in SEEDS]
        starts=np.tile(np.array([v[0] for v in source_starts]),(48,1))
        result=MotorRowsGPU(n).scan(parameters,starts,transient=200.,duration=400.,rtol=1e-9,
            step_factor=1.,max_steps=8000000,max_wall=min(wall/2.,remaining),time_mode='original')
        result.update(tuple_id=np.repeat(np.arange(48),3),start_id=np.tile(np.arange(3),48),
            activity_path=dict(mapping='mapping.json',seeds=SEEDS,
                absolute_initial_positions_nm=[v[1].tolist() for v in source_starts],
                constants_fixed=True,fixed_r=m['fixed_r'],fixed_b=1.,
                scope=m['scope'],time=m['time']))
        save(folder/f'N{n}.npz',result)
        valid=result['status']==0;pre=valid&(result['lle']>.005)&np.all(result['lle_windows']>.001,axis=1)
        records=[]
        for index in np.flatnonzero(pre):
            records.append(dict(lane=int(index),tuple_id=int(result['tuple_id'][index]),seed=SEEDS[int(result['start_id'][index])],
                parameters=result['parameters'][index].tolist(),mean_lle=float(result['lle'][index]),
                windows=result['lle_windows'][index].tolist(),mean_divergence=float(result['mean_divergence_original'][index]),
                delta=m['tuples'][int(result['tuple_id'][index])]['delta']))
        entry=dict(N=n,lanes=len(valid),completed=int(valid.sum()),status_counts={str(j):int(np.sum(result['status']==j)) for j in range(5)},
            preliminary_count=int(pre.sum()),candidate_records=records,wall_seconds=result['wall_seconds'],
            maximum_completed_lle=float(np.max(result['lle'][valid])) if valid.any() else None,
            exact_source_control_rows={key:np.flatnonzero(result['parameters'][:,0]==m['source_controls'][key]['rho']).tolist() for key in ['below','above','critical']})
        summary['runs'].append(entry);summary['total_wall_seconds']=time.perf_counter()-began
        (folder/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print(json.dumps(entry),flush=True)
    render(folder)


def refine(folder,wall=600.):
    """Carry every completed preliminary state and tangent; no scout rerun."""
    from motor_rows_gpu import MotorRowsGPU,save
    folder=Path(folder);began=time.perf_counter();report=dict(runs=[],time_mode='original',
        selection='Every scout lane with mean>.005 and all four windows>.001',
        script_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest())
    for n in (3,8):
        source=folder/f'N{n}.npz';old=np.load(source)
        ids=np.flatnonzero((old['status']==0)&(old['lle']>.005)&np.all(old['lle_windows']>.001,axis=1))
        if not len(ids):continue
        remaining=wall-(time.perf_counter()-began)
        if remaining<=1:report['runs'].append(dict(N=n,status='not started: wall budget'));break
        result=MotorRowsGPU(n).scan(old['parameters'][ids],old['final_states'][ids],tangents=old['final_tangents'][ids],
            transient=300.,duration=1200.,rtol=1e-10,step_factor=.5,max_steps=8000000,
            max_wall=min(wall/2.,remaining),time_mode='original')
        result.update(source_lane_indices=ids,source_npz_exact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            source_npz=source.relative_to(ROOT).as_posix(),original_start_id=old['start_id'][ids],original_tuple_id=old['tuple_id'][ids],
            activity_path=dict(constants_fixed=True,fixed_r=float(old['parameters'][0,1]),fixed_b=1.,
                scope='Published activity path extended beyond illustrated delta=+/-0.1; not apparatus calibration'))
        save(folder/f'N{n}-refined.npz',result)
        valid=result['status']==0;kept=valid&(result['lle']>.005)&np.all(result['lle_windows']>.001,axis=1)
        entry=dict(N=n,lanes=len(ids),completed=int(valid.sum()),retained=int(kept.sum()),
            retained_source_lane_indices=ids[kept].tolist(),status_counts={str(j):int(np.sum(result['status']==j)) for j in range(5)},
            wall_seconds=result['wall_seconds'])
        report['runs'].append(entry);report['total_wall_seconds']=time.perf_counter()-began
        (folder/'refinement-summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(entry),flush=True)


def cpu_check(folder,index=69,warm=200.,observe=1000.,wall=600.,extend=False):
    """One independent original-time full23D QR check, with short Radau Jv."""
    from scipy.integrate import solve_ivp
    from motor_rows import rhs,field_jac,dimension,unpack,source_force_constraint_field
    folder=Path(folder);began=time.perf_counter();n=8;d=dimension(n);source=folder/'N8.npz';old=np.load(source)
    p=old['parameters'][index];x=old['final_states'][index].copy();out=folder/f'cpu-N8-row{index}'
    previous=None;previous_data=None
    if extend:
        previous_path=out.with_suffix('.json');previous=json.loads(previous_path.read_text());previous_data=np.load(out.with_suffix('.npz'))
        if previous['status']!='complete':raise RuntimeError('Only a completed CPU audit may be extended by this command')
        assert np.array_equal(p,previous_data['parameters'])
        assert np.array_equal(previous_data['final_state'],np.array(previous['final_state']))
        assert np.array_equal(previous_data['final_basis'],np.array(previous['final_basis']))
        assert np.array_equal(previous_data['log_sums'],np.array(previous['log_sums']))
        x=previous_data['final_state'].copy();out=folder/f'cpu-N8-row{index}-extended';warm=0.;observe=3000.
    f,j=field_jac(x,p,n);eye=np.eye(d);v=np.sin(np.arange(d)+1.);v/=np.linalg.norm(v)
    complex_j=np.column_stack([np.imag(rhs(0,x.astype(complex)+1e-28j*eye[k],p,n))/1e-28 for k in range(d)])
    c,s,xx=unpack(x,n);exact_trace=-2*n-p[0]*(n-1)+p[1]*(1-1/n)*np.sum(c)
    control=dict(jacobian_complex_step_error=float(np.max(abs(j-complex_j))),
        independent_force_field_error=float(np.max(abs(f-source_force_constraint_field(x,p,n)))),
        trace_formula_error=float(abs(np.trace(j)-exact_trace)))
    def short_rhs(t,z):
        ff,jj=field_jac(z[:d],p,n);return np.r_[ff,jj@z[d:]]
    if previous is None:
        answers={}
        for method in ('DOP853','Radau'):
            sol=solve_ivp(short_rhs,(0,2.),np.r_[x,v],method=method,rtol=2e-12,atol=2e-14,max_step=.01)
            if not sol.success:raise RuntimeError(sol.message)
            answers[method]=sol.y[:,-1]
        original=solve_ivp(lambda t,z:source_force_constraint_field(z,p,n),(0,2.),x,
            method='DOP853',rtol=2e-12,atol=2e-14,max_step=.01)
        if not original.success:raise RuntimeError(original.message)
        control.update(duration_tau=2.,DOP853_Radau_state_error=float(np.max(abs(answers['DOP853'][:d]-answers['Radau'][:d]))),
            DOP853_Radau_tangent_error=float(np.max(abs(answers['DOP853'][d:]-answers['Radau'][d:]))),
            independent_force_state_error=float(np.max(abs(answers['DOP853'][:d]-original.y[:,-1]))),
            scope='Full23D state plus one analytic tangent; full23x23 basis is used in the longer QR calculation')
        control['passed']=all(control[key]<1e-8 for key in ['jacobian_complex_step_error','independent_force_field_error',
            'trace_formula_error','DOP853_Radau_state_error','DOP853_Radau_tangent_error','independent_force_state_error'])
    else:
        control=dict(inherited_short_control=previous['short_control'],continuation_endpoint_algebra=control,
            passed=bool(previous['short_control']['passed'] and max(control.values())<1e-8),
            scope='No repeated integration of the old control; inherited solver controls and fresh endpoint algebra')
    report=dict(N=n,dimension=d,source_lane=index,parameters=p.tolist(),initial_state=x.tolist(),short_control=control,
        source_npz_exact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),status='running',windows=[],
        settings=dict(warm=warm,observe=observe,rtol=1e-11,atol=1e-13,maxstep=.03,qr_max=.2,checkpoint=100.,wall_budget=wall),
        time='Original normalized tau',max_local_trace_error=0.,max_local_flow_identity_error=0.,max_local_condition=0.,
        script_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        qualification='23D existing model on fixed-constant activity path. Finite-time positive exponents require convergence checks; no priority claim.')
    basis=eye;elapsed=0.;next_dt=.2;logs=np.zeros(d);trace=0.;flowlog=0.;count=0;rejected=0;observed=0.;flow_initial=None
    times=[];states=[];sections=[];section_times=[];checkpoints=[]
    if previous is not None:
        basis=previous_data['final_basis'].copy();elapsed=float(previous_data['elapsed_original']);next_dt=float(previous_data['next_qr_interval'])
        logs=previous_data['log_sums'].copy();trace=float(previous_data['integrated_trace']);flowlog=float(previous_data['normalized_flow_log'])
        observed=previous['observed_tau'];count=previous['qr_chunks'];rejected=previous['rejected_qr_chunks']
        times=previous_data['time'].tolist();states=list(previous_data['states'].copy());sections=list(previous_data['section_states'].copy())
        section_times=previous_data['section_time'].tolist();checkpoints=list(previous_data['checkpoint'].copy())
        warm_index=int(np.argmin(abs(previous_data['time']-previous['settings']['warm'])))
        flow_initial=np.linalg.norm(rhs(0,previous_data['states'][warm_index],p,n))
        flow_check=float(abs(np.log(np.linalg.norm(rhs(0,x,p,n))/flow_initial)/observed-previous['telescopic_flow_exponent']))
        if flow_check>1e-12:raise RuntimeError('Unable to recover the exact original observation flow norm')
        report['windows']=list(previous['windows'])
        report['continuation']=dict(input_metadata_exact_sha256=hashlib.sha256(previous_path.read_bytes()).hexdigest(),
            input_npz_exact_sha256=hashlib.sha256(previous_path.with_suffix('.npz').read_bytes()).hexdigest(),
            state_basis_log_sums_exactly_matched=True,original_observed_tau=observed,original_elapsed_tau=elapsed,
            recovered_initial_flow_rate_error=flow_check,existing_histories_carried=True,
            additional_requested_tau=observe,original_record_remains_immutable=True)
        report['max_local_trace_error']=previous['max_local_trace_error'];report['max_local_flow_identity_error']=previous['max_local_flow_identity_error']
        report['max_local_condition']=previous['max_local_condition']
    def save_cpu():
        report.update(wall_seconds=time.perf_counter()-began,elapsed_original=elapsed,qr_chunks=count,rejected_qr_chunks=rejected,
            final_state=x.tolist(),final_basis=basis.tolist(),log_sums=logs.tolist(),integrated_trace=float(trace),
            normalized_flow_log=float(flowlog),next_qr_interval=next_dt)
        out.with_suffix('.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        np.savez_compressed(out.with_suffix('.npz'),parameters=p,time=np.array(times),states=np.array(states),
            section_time=np.array(section_times),section_states=np.array(sections),checkpoint=np.array(checkpoints),
            final_state=x,final_basis=basis,log_sums=logs,integrated_trace=trace,normalized_flow_log=flowlog,
            elapsed_original=elapsed,next_qr_interval=next_dt)
    save_cpu()
    if not control['passed']:report['status']='short control failed';save_cpu();return
    def augmented(t,z):
        ff,jj=field_jac(z[:d],p,n)
        return np.r_[ff,(jj@z[d:d+d*d].reshape(d,d)).ravel(),np.trace(jj)]
    phases=[('warm',min(100.,warm-i)) for i in np.arange(0,warm,100.)]+[('observation',min(100.,observe-i)) for i in np.arange(0,observe,100.)]
    for phase,length in phases:
        spent=0.;wlogs=np.zeros(d);wtrace=0.;wflow=0.;flow_start=np.linalg.norm(rhs(0,x,p,n));old_sections=len(sections)
        if phase=='observation' and flow_initial is None:flow_initial=flow_start
        while spent<length:
            if time.perf_counter()-began>wall:
                report.update(status='wall budget incomplete',partial_window=dict(phase=phase,duration=spent,logs=wlogs.tolist(),trace=float(wtrace),flow_log=float(wflow)))
                save_cpu();return
            dt=min(next_dt,length-spent);f0=rhs(0,x,p,n);nf=np.linalg.norm(f0)
            if nf<1e-14:report['status']='near equilibrium; no neutral-flow estimate';save_cpu();return
            def section(t,z):return z[2*n]
            section.direction=1
            sol=solve_ivp(augmented,(0,dt),np.r_[x,basis.ravel(),0.],method='DOP853',rtol=1e-11,atol=1e-13,
                max_step=min(.03,dt),events=section)
            if not sol.success:report.update(status='solver failure',message=sol.message);save_cpu();return
            end=sol.y[:d,-1];matrix=sol.y[d:d+d*d,-1].reshape(d,d);new_basis,R=np.linalg.qr(matrix);diag=np.diag(R)
            error=np.inf if np.any(diag==0) else abs(np.log(abs(diag)).sum()-sol.y[-1,-1]);condition=np.linalg.cond(matrix)
            if error>1e-8 or condition>1e10:
                next_dt=dt/2;rejected+=1
                if next_dt<1e-6:report['status']='QR conditioning failure';save_cpu();return
                continue
            prop=matrix@(basis.T@(f0/nf));f1=rhs(0,end,p,n);identity=np.linalg.norm(prop*nf-f1)/(1+np.linalg.norm(f1))
            growth=np.log(abs(diag));fg=np.log(np.linalg.norm(prop))
            report['max_local_trace_error']=max(report['max_local_trace_error'],float(error))
            report['max_local_flow_identity_error']=max(report['max_local_flow_identity_error'],float(identity))
            report['max_local_condition']=max(report['max_local_condition'],float(condition))
            x=end;basis=new_basis;elapsed+=dt;spent+=dt;count+=1
            if length-spent<1e-11:spent=length
            next_dt=min(.2,next_dt*1.05)
            if phase=='observation':
                logs+=growth;trace+=sol.y[-1,-1];flowlog+=fg;wlogs+=growth;wtrace+=sol.y[-1,-1];wflow+=fg
                for hit,z in zip(sol.t_events[0],sol.y_events[0]):
                    if hit>1e-10:section_times.append(elapsed-dt+float(hit));sections.append(z[:d])
            times.append(elapsed);states.append(x.copy())
        if phase=='observation':
            observed+=length
            report['windows'].append(dict(end_tau=elapsed,duration=length,spectrum=(wlogs/length).tolist(),
                mean_divergence=float(wtrace/length),sum_trace_error=float(abs(wlogs.sum()-wtrace)/length),
                flow_exponent=float(wflow/length),telescopic_flow_exponent=float(np.log(np.linalg.norm(rhs(0,x,p,n))/flow_start)/length),
                upward_X1_returns=len(sections)-old_sections))
            report.update(observed_tau=observed,spectrum=(logs/observed).tolist(),mean_divergence=float(trace/observed),
                sum_trace_error=float(abs(logs.sum()-trace)/observed),flow_exponent=float(flowlog/observed),
                telescopic_flow_exponent=float(np.log(np.linalg.norm(rhs(0,x,p,n))/flow_initial)/observed),
                upward_X1_returns=len(sections))
            print(json.dumps(dict(CPU_N8_tau=elapsed,leading_spectrum=report['spectrum'][:10],
                trace_error=report['sum_trace_error'],returns=len(sections),wall=round(time.perf_counter()-began,1))),flush=True)
        checkpoints.append(np.r_[elapsed,x,basis.ravel(),logs,trace,flowlog,next_dt]);save_cpu()
    sample=np.array(states);active=np.array(times)>(previous['settings']['warm'] if previous is not None else warm);cs=sample[active,:n];ss=sample[active,n:2*n]
    positions=np.c_[sample[active,2*n:],-sample[active,2*n:].sum(axis=1)];pairs=[]
    for i in range(n):
        for j in range(i+1,n):
            distance=np.sqrt((cs[:,i]-cs[:,j])**2+(ss[:,i]-ss[:,j])**2+(positions[:,i]-positions[:,j])**2)
            pairs.append(dict(rows=[i+1,j+1],minimum=float(distance.min()),maximum=float(distance.max()),rms=float(np.sqrt(np.mean(distance**2)))))
    report.update(status='complete',row_pair_distances=pairs,sampled_domain=dict(max_density_disk_radius_squared=float(np.max(cs*cs+ss*ss)),
        disk_bound=float(1+p[2]**2),max_abs_sum_X=float(np.max(abs(positions.sum(axis=1))))))
    save_cpu()


def lifetime(folder,wall=600.):
    """Fixed entire 32-lane retained N8 cohort, ten blocks, no pruning/warmup."""
    from motor_rows_gpu import MotorRowsGPU,save
    folder=Path(folder);source=folder/'N8-refined.npz';original=np.load(source)
    selected=np.flatnonzero((original['status']==0)&(original['lle']>.005)&np.all(original['lle_windows']>.001,axis=1))
    parameters=original['parameters'][selected].copy();states=original['final_states'][selected].copy();tangents=original['final_tangents'][selected].copy()
    ids=original['source_lane_indices'][selected].copy();gpu=MotorRowsGPU(8);began=time.perf_counter()
    report=dict(source_npz_exact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cohort_source_lanes=ids.tolist(),
        cohort_size=len(ids),blocks=[],requested_additional_tau=5000.,block_tau=500.,observation_window_tau=125.,
        no_pruning=True,no_warmup=True,settings=dict(rtol=1e-10,step_factor=.5,max_steps_per_lane=8000000,wall_budget=wall),
        carry='Exact saved states and aligned tangent directions; scan normalizes the carried vectors at each new block',
        qualification='Finite persistence, not an asymptotic attractor proof. A lost threshold is not a verified periodic capture.',
        script_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest())
    window_history=[];mean_history=[];status_history=[]
    for block in range(1,11):
        remaining=wall-(time.perf_counter()-began)
        if remaining<=1:report['status']='wall budget before next block';break
        result=gpu.scan(parameters,states,tangents=tangents,transient=0.,duration=500.,rtol=1e-10,step_factor=.5,
            max_steps=8000000,max_wall=remaining,time_mode='original')
        result.update(source_lane_indices=ids,block_number=block,requested_absolute_tau_from_original_start=2100.+block*500.,
            baseline_source_npz_exact_sha256=report['source_npz_exact_sha256'],no_warmup=True,no_pruning=True)
        save(folder/f'N8-lifetime-block-{block:02d}.npz',result)
        valid=result['status']==0;positive=valid&(result['lle']>.005)&np.all(result['lle_windows']>.001,axis=1)
        entry=dict(block=block,completed=int(valid.sum()),positive_screen=int(positive.sum()),
            positive_source_lanes=ids[positive].tolist(),status_counts={str(j):int(np.sum(result['status']==j)) for j in range(5)},
            wall_seconds=result['wall_seconds'],additional_completed_common_tau=500.*(block if valid.all() else block-1))
        report['blocks'].append(entry);window_history.append(result['lle_windows'].copy());mean_history.append(result['lle'].copy());status_history.append(result['status'].copy())
        states=result['final_states'].copy();tangents=result['final_tangents'].copy()
        report.update(total_wall_seconds=time.perf_counter()-began,status='complete' if block==10 and valid.all() else 'checkpoint')
        (folder/'lifetime-summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        np.savez_compressed(folder/'N8-lifetime-history.npz',source_lane_indices=ids,parameters=parameters,
            lle_windows=np.array(window_history),lle=np.array(mean_history),status=np.array(status_history),
            final_states=states,final_tangents=tangents)
        print(json.dumps(entry),flush=True)
        if not valid.all():report['status']='incomplete block preserved; cohort not pruned';break
    report['total_wall_seconds']=time.perf_counter()-began
    (folder/'lifetime-summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def postprocess(folder,gif=False):
    """Saved-data-only spectra, persistence and unreduced row geometry."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation,PillowWriter
    folder=Path(folder);report=json.loads((folder/'cpu-N8-row69-extended.json').read_text())
    data=np.load(folder/'cpu-N8-row69-extended.npz');life=np.load(folder/'N8-lifetime-history.npz')
    if report['status']!='complete' or life['lle_windows'].shape!=(10,32,4):raise RuntimeError('Completed records required')
    t=data['time'];state=data['states'];mask=t>200.;t=t[mask];z=state[mask];c=z[:,:8];s=z[:,8:16];X=np.c_[z[:,16:],-z[:,16:].sum(axis=1)]
    pairs=[(i,j) for i in range(8) for j in range(i+1,8)]
    distances=np.stack([np.sqrt((c[:,i]-c[:,j])**2+(s[:,i]-s[:,j])**2+(X[:,i]-X[:,j])**2) for i,j in pairs],axis=1)
    nearest=np.argmin(distances,axis=1);cluster={}
    for threshold in [1e-3,1e-2,.05]:
        counts=[]
        for near in distances<threshold:
            parent=list(range(8))
            def find(i):
                while parent[i]!=i:i=parent[i]
                return i
            for yes,(i,j) in zip(near,pairs):
                if yes:parent[find(i)]=find(j)
            counts.append(len({find(i) for i in range(8)}))
        cluster[str(threshold)]={str(k):int(np.sum(np.array(counts)==k)) for k in range(1,9)}
    windows=report['windows'];rates=np.array([w['spectrum'] for w in windows]);weights=np.array([w['duration'] for w in windows]);tt=np.cumsum(weights)
    cumulative=np.cumsum(rates*weights[:,None],axis=0)/tt[:,None];blocks=[]
    for i in range(0,len(windows),10):
        blocks.append(dict(observed_tau_range=[i*100,(i+10)*100],spectrum=np.average(rates[i:i+10],weights=weights[i:i+10],axis=0).tolist()))
    retained=np.all(life['lle_windows']>.001,axis=(0,2))&np.all(life['lle']>.005,axis=0)
    selected=int(np.flatnonzero(life['source_lane_indices']==69)[0])
    result=dict(cpu_spectrum=report['spectrum'],cpu_observed_tau=report['observed_tau'],cpu_thousand_tau_blocks=blocks,
        cpu_upward_X1_returns=report['upward_X1_returns'],cpu_sum_trace_error=report['sum_trace_error'],
        cpu_flow_rate=report['flow_exponent'],cpu_exact_telescopic_flow_rate=report['telescopic_flow_exponent'],
        CPU_first_six_positive_in_each_thousand_tau_block=bool(all(np.all(np.array(v['spectrum'])[:6]>0) for v in blocks)),
        neutral_mode_qualification='Seventh QR estimate is near zero and the flow identity independently supplies a near-zero finite-time rate. Do not count a tiny finite positive seventh estimate as a further asymptotic expanding direction.',
        fixed_cohort=dict(lanes=32,blocks=10,extra_tau=5000.,total_tau_from_original_start=7100.,all_completed=bool(np.all(life['status']==0)),
            all40windows_and10means_retained=int(retained.sum()),retained_source_lanes=life['source_lane_indices'][retained].tolist(),
            lost_screen_source_lanes=life['source_lane_indices'][~retained].tolist(),
            minimum_window_among_retained=float(life['lle_windows'][:,retained].min()),
            representative69_mean_over_extra5000=float(life['lle'][:,selected].mean()),
            representative69_window_range=[float(life['lle_windows'][:,selected].min()),float(life['lle_windows'][:,selected].max())],
            lost_screen_mean_history=life['lle'][:,~retained].tolist(),
            caveat='No pruning. Lost positivity is not a verified periodic capture. Finite persistence is not a lifetime distribution or an asymptotic proof.'),
        row_structure=dict(pair_metric='Euclidean normalized(c_i,s_i,X_i)',pairs=[[i+1,j+1] for i,j in pairs],
            all_pair_rms=np.sqrt(np.mean(distances**2,axis=0)).tolist(),all_pair_max=distances.max(axis=0).tolist(),
            all_pair_min=distances.min(axis=0).tolist(),nearest_pair_fraction=[float(np.mean(nearest==j)) for j in range(28)],
            nearest_pair_switches=int(np.sum(np.diff(nearest)!=0)),sampled_component_count_histograms=cluster,
            sample_count=len(t),sample_spacing='About .2 original tau',
            qualification='Threshold graph connected components describe instantaneous proximity only, including transitive chains. No fixed synchronized cluster reduction is imposed or inferred.'),
        source_hashes={})
    for name in ['mapping.json','N3.npz','N8.npz','N3-refined.npz','N8-refined.npz','cpu-N8-row69.json','cpu-N8-row69.npz',
            'cpu-N8-row69-extended.json','cpu-N8-row69-extended.npz','N8-lifetime-history.npz']:
        b=(folder/name).read_bytes();result['source_hashes'][name]=hashlib.sha256(b if name.endswith('.npz') else b.replace(b'\r\n',b'\n')).hexdigest()
    result['source_hash_policy']='Exact NPZ bytes; text normalized CRLF to LF'
    (folder/'diagnostics-summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    fig=plt.figure(figsize=(14,9),layout='constrained');fig.suptitle('Eight-row motor model on the published activity path: numerical hyperchaotic candidate\nFixed source constants; all 3D views are projections of the full 23D state',fontsize=14)
    for panel,(values,labels) in enumerate([(np.c_[X[:,:2],s[:,0]],['X1','X2','s1']),
            (np.c_[c[:,0],s[:,0],X[:,0]],['c1','s1','X1']),(s[:,:3],['s1','s2','s3'])],1):
        ax=fig.add_subplot(2,3,panel,projection='3d');ax.scatter(*values[::2].T,s=.25,c='#296e87',alpha=.35,depthshade=False)
        ax.set(xlabel=labels[0],ylabel=labels[1],zlabel=labels[2],title='Physical normalized-coordinate projection');ax.view_init(24,35)
    ax=fig.add_subplot(234);sec=data['section_states'];st=data['section_time'];dots=ax.scatter(sec[:,17],sec[:,8],c=st,s=4,cmap='viridis');fig.colorbar(dots,ax=ax,label='tau')
    ax.set(xlabel='X2 on X1=0',ylabel='s1 on X1=0',title=f'Upward Poincare section: {len(sec)} returns');ax.grid(alpha=.15)
    ax=fig.add_subplot(235)
    for j in range(8):ax.plot(tt,cumulative[:,j],lw=1,label=f'{j+1}')
    ax.axhline(0,c='black',lw=.7);ax.set(xlabel='Observed original tau',ylabel='Finite-time QR exponent',title='Leading eight directions; full basis has 23');ax.legend(ncol=4,fontsize=7);ax.grid(alpha=.15)
    ax=fig.add_subplot(236);spectrum=np.array(report['spectrum']);colors=['#bc5429' if j<6 else '#34738c' if j==6 else '#8a959b' for j in range(23)]
    ax.bar(np.arange(1,24),spectrum,color=colors);ax.axhline(0,c='black',lw=.7);ax.set(xlabel='QR direction',ylabel='Finite-time exponent',title='Complete spectrum after 4,000 tau');ax.grid(alpha=.15)
    fig.savefig(folder/'N8-geometry.png',dpi=170);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained');fig.suptitle('Finite persistence and full-row diagnostics; no surviving lanes discarded',fontsize=14)
    ax=axs[0,0];im=ax.imshow(life['lle'].T,aspect='auto',origin='lower',extent=[0,5000,-.5,31.5],cmap='viridis');fig.colorbar(im,ax=ax,label='Leading FTLE per tau')
    ax.set(xlabel='Additional original tau',ylabel='Fixed cohort lane (source ID)',yticks=np.arange(32),yticklabels=life['source_lane_indices'],title='All 32 lanes, ten blocks of 500 tau');ax.tick_params(axis='y',labelsize=7)
    ax=axs[0,1];zoom=t>=t[-1]-50
    ax.plot(t[zoom],distances[zoom].min(axis=1),label='Closest pair',lw=.8);ax.plot(t[zoom],np.median(distances[zoom],axis=1),label='Median pair',lw=.8)
    ax.set(xlabel='tau',ylabel='Normalized full-row distance',title='Final 50 tau: rows remain distinct');ax.legend(fontsize=8);ax.grid(alpha=.15)
    ax=axs[1,0];rms=np.zeros((8,8))
    for value,(i,j) in zip(result['row_structure']['all_pair_rms'],pairs):rms[i,j]=rms[j,i]=value
    im=ax.imshow(rms,origin='lower',cmap='viridis');fig.colorbar(im,ax=ax,label='RMS full-row distance');ax.set(xticks=np.arange(8),yticks=np.arange(8),xticklabels=np.arange(1,9),yticklabels=np.arange(1,9),xlabel='Row',ylabel='Row',title='All 4,000 observed tau; no fixed cluster collapse')
    ax=axs[1,1];disk=c*c+s*s
    ax.fill_between(t[zoom],disk[zoom].min(axis=1),disk[zoom].max(axis=1),alpha=.4,color='#296e87',label='Range across 8 rows')
    ax.axhline(2,c='black',ls='--',lw=.8);ax.set(xlabel='tau',ylabel='c_i squared + s_i squared',title='Invariant density disks (bound 2)');ax.legend(fontsize=8);ax.grid(alpha=.15)
    fig.savefig(folder/'N8-persistence.png',dpi=170);plt.close(fig)
    if gif:
        fig=plt.figure(figsize=(5.5,5.5));ax=fig.add_subplot(111,projection='3d');ax.scatter(X[::2,0],X[::2,1],s[::2,0],s=.45,c='#296e87',alpha=.4,depthshade=False)
        ax.set(xlabel='X1',ylabel='X2',zlabel='s1');fig.suptitle('Published activity path: 23D motor candidate\nRotation of a physical-coordinate projection',fontsize=10)
        def update(frame):ax.view_init(22+16*np.sin(frame*2*np.pi/60),frame*6);return []
        FuncAnimation(fig,update,frames=60,interval=80).save(folder/'N8-rotation.gif',writer=PillowWriter(fps=12),dpi=95);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=ROOT/'results/motor-rows/activity-path')
    parser.add_argument('--wall',type=float,default=600.);parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--render-only',action='store_true');parser.add_argument('--refine',action='store_true')
    parser.add_argument('--cpu-check',type=int,metavar='N8_LANE');parser.add_argument('--cpu-extend',type=int,metavar='N8_LANE')
    parser.add_argument('--lifetime',action='store_true');parser.add_argument('--postprocess',action='store_true');parser.add_argument('--gif',action='store_true');args=parser.parse_args()
    if args.postprocess:postprocess(args.output_dir,args.gif)
    elif args.cpu_extend is not None:cpu_check(args.output_dir,args.cpu_extend,wall=args.wall,extend=True)
    elif args.lifetime:lifetime(args.output_dir,args.wall)
    elif args.cpu_check is not None:cpu_check(args.output_dir,args.cpu_check,wall=args.wall)
    elif args.refine:refine(args.output_dir,args.wall)
    elif args.render_only:render(args.output_dir)
    else:run(args.output_dir,args.wall,args.prepare_only)
