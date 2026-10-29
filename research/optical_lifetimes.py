"""Capture into independently checked optical periodic orbits.

Finite observation can support persistent chaos but cannot prove an attractor.
The candidate model uses an assumed logistic illumination profile.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from smooth_operators import SmoothOperatorsGPU, rhs_and_jac

ROOT=Path('results/light-operators/lifetimes')


def cycle_check(p,start,duration=400):
    def event(t,y):return y[0]-p[2]
    event.direction=1
    sol=solve_ivp(lambda t,y:rhs_and_jac(y,p)[0],(0,duration),start,
        method='DOP853',rtol=1e-12,atol=1e-14,max_step=.005,events=event)
    if not sol.success:raise RuntimeError(sol.message)
    section=sol.y_events[0];times=sol.t_events[0]
    result=dict(parameters=p.tolist(),initial_state=start.tolist(),duration=duration,
        method='DOP853',rtol=1e-12,atol=1e-14,max_step=.005,crossings=len(times))
    result['return_errors']={str(q):float(np.max(abs(section[-5*q:]-section[-6*q:-q])))
        for q in range(1,min(33,len(section)//6))}
    orders=[int(q) for q,error in result['return_errors'].items() if error<1e-8]
    if not orders:
        result['periodic']=False
        return result,section,times
    q=min(orders);period=times[-1]-times[-1-q];y=section[-1-q]
    def tangent(t,z):
        f,j=rhs_and_jac(z[:6],p);return np.r_[f,(j@z[6:].reshape(6,6)).ravel()]
    ref=solve_ivp(tangent,(0,period),np.r_[y,np.eye(6).ravel()],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.0025)
    if not ref.success:raise RuntimeError(ref.message)
    last=ref.y[:6,-1];m=ref.y[6:,-1].reshape(6,6);f=rhs_and_jac(last,p)[0]
    if abs(f[0])<1e-7:raise RuntimeError('Near-grazing periodic control')
    normal=np.eye(6)[0];derivative=((np.eye(6)-np.outer(f,normal)/f[0])@m)[1:,1:]
    eig=np.linalg.eigvals(derivative);closure=float(np.max(abs(last-y)))
    result.update(periodic=True,section_return_order=q,period=float(period),section_state=y.tolist(),
        closure_error=closure,poincare_eigenvalues=[[float(z.real),float(z.imag)] for z in eig],
        spectral_radius=float(max(abs(eig))),stable_numerically=bool(closure<1e-8 and max(abs(eig))<1),
        caveat='Numerical periodic-orbit and Floquet check, not a rigorous existence or basin proof')
    return result,section,times


def controls():
    source=np.load('results/light-operators/smooth-width-refined.npz');rows=[]
    for width in [.003,.005]:
        candidates=np.flatnonzero((source['reference_index']==1026)&(source['parameters'][:,4]==width))
        index=candidates[np.argmin(source['lle_windows'][candidates,-1])]
        result,section,times=cycle_check(source['parameters'][index],source['final_states'][index])
        result.update(source_file='results/light-operators/smooth-width-refined.npz',source_index=int(index))
        rows.append(result);np.savez_compressed(ROOT/f'cycle-w{width:g}.npz',sections=section,crossing_times=times)
        print(json.dumps(result,indent=2),flush=True)
    data=json.loads(Path('results/light-operators/smooth-coexisting-periodic-control.json').read_text())
    cpu=json.loads(Path('results/light-operators/smooth-cpu-1026-w004-regular.json').read_text())
    section=np.load('results/light-operators/smooth-section-w004-regular.npz')['sections'][-1]
    rows.append(dict(parameters=cpu['parameters']+[.04],section_state=section.tolist(),
        periodic=True,stable_numerically=True,**data['periodic_control'],
        source_file='results/light-operators/smooth-coexisting-periodic-control.json'))
    (ROOT/'periodic-controls.json').write_text(json.dumps(rows,indent=2))


def validate_chunks():
    gpu=SmoothOperatorsGPU();p=np.array([[74.1567393414624,1.784180421282182,.4139517339970917,.25683679452165964,.005]])
    start=np.array([[.3,.4,.2,.4,-.2,.5]])
    common=dict(transient=0,rtol=1e-10,step_factor=.5)
    whole=gpu.scan(p,start,duration=2,**common)
    first=gpu.scan(p,start,duration=1,**common)
    second=gpu.scan(p,first['final_states'],duration=1,balanced_tangents=first['balanced_tangents'],**common)
    result=dict(state_error=float(np.max(abs(whole['final_states']-second['final_states']))),
        tangent_error=float(np.max(abs(whole['balanced_tangents']-second['balanced_tangents']))),
        lle_error=float(abs(whole['lle'][0]-(first['lle'][0]+second['lle'][0])/2)))
    result['passed']=bool(max(result.values())<1e-7 and all(d['status'][0]==0 for d in [whole,first,second]))
    (ROOT/'chunk-validation.json').write_text(json.dumps(result,indent=2));print(result)
    if not result['passed']:raise RuntimeError('Chunk continuation validation failed')


def run(starts_per_width=32,blocks=40,chunk=500.,rtol=1e-9,step_factor=.5):
    if starts_per_width<1 or blocks<1 or chunk<=0:raise ValueError('Positive cohort and observation settings required')
    controls=json.loads((ROOT/'periodic-controls.json').read_text())
    if not all(c.get('stable_numerically') for c in controls):raise RuntimeError('Need checked stable periodic controls first')
    seed=2026102969;rng=np.random.default_rng(seed);gpu=SmoothOperatorsGPU()
    params=np.repeat(np.array([c['parameters'] for c in controls]),starts_per_width,axis=0)
    n=len(params);initial=rng.uniform(.05,.95,(n,6));initial[:,[1,4]]=0
    states=initial.copy();tangents=None;history=[];status=np.zeros(n,np.int32)
    capture_time=np.full(n,np.nan);near_streak=np.zeros(n,np.int32);previous_near=np.full(n,np.nan)
    wall=time.perf_counter();snapshots=[]
    for block in range(blocks):
        active=np.flatnonzero(status==0)
        if not len(active):break
        result=gpu.scan(params[active],states[active],transient=0,duration=chunk,rtol=rtol,
            step_factor=step_factor,max_wall_seconds=180,
            balanced_tangents=None if tangents is None else tangents[active])
        if tangents is None:tangents=np.zeros((n,6))
        states[active]=result['final_states'];tangents[active]=result['balanced_tangents']
        failed=result['status']!=0;status[active[failed]]=10+result['status'][failed]
        exponents=np.full(n,np.nan);windows=np.full((n,4),np.nan);residuals=np.full(n,np.nan)
        exponents[active]=result['lle'];windows[active]=result['lle_windows']
        for c in controls:
            group=np.flatnonzero((status==0)&(params[:,4]==c['parameters'][4]))
            if not len(group):continue
            probe=gpu.scan(params[group],states[group],transient=0,duration=c['period'],rtol=1e-10,step_factor=.5,max_wall_seconds=30)
            scale=np.array([1,np.sqrt(c['parameters'][0]),1,1,np.sqrt(c['parameters'][0]),1])
            residuals[group]=np.max(abs((probe['final_states']-states[group])/scale),axis=1)
            near=(probe['status']==0)&(residuals[group]<1e-7)&(np.max(abs(windows[group]),axis=1)<.01)
            previous_near[group[near&(near_streak[group]==0)]]=(block+1)*chunk
            near_streak[group]=np.where(near,near_streak[group]+1,0)
            captured=group[near_streak[group]>=2];status[captured]=1;capture_time[captured]=previous_near[captured]
        history.append(dict(exponents=exponents,windows=windows,residuals=residuals,status=status.copy()))
        snapshots.append(states.copy())
        settings=dict(seed=seed,starts_per_width=starts_per_width,requested_blocks=blocks,completed_blocks=block+1,
            chunk_duration=chunk,rtol=rtol,atol=rtol*.001,step_factor=step_factor,
            capture_rule='Known-period flow residual <1e-7 in velocity-balanced max norm and all four chunk exponent magnitudes <.01 at two consecutive checkpoints. Capture time records first qualifying checkpoint, confirmation one chunk later.',
            caveat='Numerical capture diagnostic; unclassified survivors need not be chaotic. Checkpoint timing is interval resolved; finite persistence cannot establish an attractor.',
            status_codes={'0':'still observed','1':'captured near checked stable periodic orbit','11':'thermal failure','12':'numerical failure','13':'step budget unresolved','14':'wall budget unresolved'},
            wall_seconds=time.perf_counter()-wall,device=gpu.device.name)
        np.savez_compressed(ROOT/'ensemble.npz',parameters=params,initial_states=initial,final_states=states,balanced_tangents=tangents,
            status=status,capture_time=capture_time,checkpoint_times=np.arange(1,block+2)*chunk,
            exponents=np.array([h['exponents'] for h in history]),lle_windows=np.array([h['windows'] for h in history]),
            period_residuals=np.array([h['residuals'] for h in history]),status_history=np.array([h['status'] for h in history]),
            state_history=np.array(snapshots))
        (ROOT/'ensemble.json').write_text(json.dumps(settings,indent=2))
        counts=[dict(width=float(c['parameters'][4]),captured=int(np.sum((params[:,4]==c['parameters'][4])&(status==1))),
            active=int(np.sum((params[:,4]==c['parameters'][4])&(status==0)))) for c in controls]
        print(json.dumps(dict(time=(block+1)*chunk,counts=counts,wall_seconds=round(time.perf_counter()-wall,1))),flush=True)


if __name__=='__main__':
    ROOT.mkdir(parents=True,exist_ok=True)
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true');parser.add_argument('--validate',action='store_true');parser.add_argument('--run',action='store_true')
    parser.add_argument('--starts',type=int,default=32);parser.add_argument('--blocks',type=int,default=40);args=parser.parse_args()
    if args.controls:controls()
    elif args.validate:validate_chunks()
    elif args.run:run(args.starts,args.blocks)
