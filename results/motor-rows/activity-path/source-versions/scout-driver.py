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


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=ROOT/'results/motor-rows/activity-path')
    parser.add_argument('--wall',type=float,default=600.);parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--render-only',action='store_true');args=parser.parse_args()
    if args.render_only:render(args.output_dir)
    else:run(args.output_dir,args.wall,args.prepare_only)
