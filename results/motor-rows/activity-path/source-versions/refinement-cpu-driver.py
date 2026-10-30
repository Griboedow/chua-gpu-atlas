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


def cpu_check(folder,index=69,warm=200.,observe=1000.,wall=600.):
    """One independent original-time full23D QR check, with short Radau Jv."""
    from scipy.integrate import solve_ivp
    from motor_rows import rhs,field_jac,dimension,unpack,source_force_constraint_field
    folder=Path(folder);began=time.perf_counter();n=8;d=dimension(n);source=folder/'N8.npz';old=np.load(source)
    p=old['parameters'][index];x=old['final_states'][index].copy();out=folder/f'cpu-N8-row{index}'
    f,j=field_jac(x,p,n);eye=np.eye(d);v=np.sin(np.arange(d)+1.);v/=np.linalg.norm(v)
    complex_j=np.column_stack([np.imag(rhs(0,x.astype(complex)+1e-28j*eye[k],p,n))/1e-28 for k in range(d)])
    c,s,xx=unpack(x,n);exact_trace=-2*n-p[0]*(n-1)+p[1]*(1-1/n)*np.sum(c)
    control=dict(jacobian_complex_step_error=float(np.max(abs(j-complex_j))),
        independent_force_field_error=float(np.max(abs(f-source_force_constraint_field(x,p,n)))),
        trace_formula_error=float(abs(np.trace(j)-exact_trace)))
    def short_rhs(t,z):
        ff,jj=field_jac(z[:d],p,n);return np.r_[ff,jj@z[d:]]
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
    report=dict(N=n,dimension=d,source_lane=index,parameters=p.tolist(),initial_state=x.tolist(),short_control=control,
        source_npz_exact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),status='running',windows=[],
        settings=dict(warm=warm,observe=observe,rtol=1e-11,atol=1e-13,maxstep=.03,qr_max=.2,checkpoint=100.,wall_budget=wall),
        time='Original normalized tau',max_local_trace_error=0.,max_local_flow_identity_error=0.,max_local_condition=0.,
        script_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        qualification='23D existing model on fixed-constant activity path. Finite-time positive exponents require convergence checks; no priority claim.')
    basis=eye;elapsed=0.;next_dt=.2;logs=np.zeros(d);trace=0.;flowlog=0.;count=0;rejected=0;observed=0.;flow_initial=None
    times=[];states=[];sections=[];section_times=[];checkpoints=[]
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
    sample=np.array(states);active=np.array(times)>warm;cs=sample[active,:n];ss=sample[active,n:2*n]
    positions=np.c_[sample[active,2*n:],-sample[active,2*n:].sum(axis=1)];pairs=[]
    for i in range(n):
        for j in range(i+1,n):
            distance=np.sqrt((cs[:,i]-cs[:,j])**2+(ss[:,i]-ss[:,j])**2+(positions[:,i]-positions[:,j])**2)
            pairs.append(dict(rows=[i+1,j+1],minimum=float(distance.min()),maximum=float(distance.max()),rms=float(np.sqrt(np.mean(distance**2)))))
    report.update(status='complete',row_pair_distances=pairs,sampled_domain=dict(max_density_disk_radius_squared=float(np.max(cs*cs+ss*ss)),
        disk_bound=float(1+p[2]**2),max_abs_sum_X=float(np.max(abs(positions.sum(axis=1))))))
    save_cpu()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=ROOT/'results/motor-rows/activity-path')
    parser.add_argument('--wall',type=float,default=600.);parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--render-only',action='store_true');parser.add_argument('--refine',action='store_true')
    parser.add_argument('--cpu-check',type=int,metavar='N8_LANE');args=parser.parse_args()
    if args.cpu_check is not None:cpu_check(args.output_dir,args.cpu_check,wall=args.wall)
    elif args.refine:refine(args.output_dir,args.wall)
    elif args.render_only:render(args.output_dir)
    else:run(args.output_dir,args.wall,args.prepare_only)
