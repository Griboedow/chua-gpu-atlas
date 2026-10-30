"""FP64 GPU scouting of the exact collective-motor three-state model.

No model extension: internal Q=h/r-Y improves conditioning. Tangents use the
Euclidean (Q,Y,Z) norm. Scaled computational time s=min(1,rho)*tau gives longer
physical observation to slow cases; every reported exponent is converted back
to original tau. Neither finite positivity nor budget exhaustion classifies an
asymptotic attractor. Raw partial state, clocks, logs and statuses are retained.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp

from motor_pair import field_jac, original_field, SOURCE_CONTROLS, equilibrium_density_start

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_TANGENT=np.array([1., .37, -.21]);DEFAULT_TANGENT/=np.linalg.norm(DEFAULT_TANGENT)


def to_balanced(y,p):
    y=np.asarray(y);return np.array([y[0]/p[1]-y[1],y[1],y[2]])


def from_balanced(x,p):
    return np.array([p[1]*(x[0]+x[1]),x[1],x[2]])


def balanced_field_jac(x,p):
    q,y,z=x;rho,r=p;H=q+y
    return np.array([-rho*H,r*H*(1-z)-y,r*H*y-z]),np.array([
        [-rho,-rho,0],[r*(1-z),r*(1-z)-1,-r*H],[r*y,r*(q+2*y),-1]])


class MotorPairGPU:
    def __init__(self):
        devices=[d for platform in cl.get_platforms() for d in platform.get_devices()
            if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('An FP64 GPU is required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        self.kernel=cl.Kernel(cl.Program(self.context,Path(__file__).with_suffix('.cl').read_text()).build(),'advance')

    def scan(self,parameters,starts,transient=200.,duration=400.,rtol=1e-8,step_factor=1.,
             max_steps=8000000,max_wall=600.,balanced_starts=None,balanced_tangents=None):
        params=np.ascontiguousarray(parameters,dtype=float);starts=np.ascontiguousarray(starts,dtype=float)
        n=len(params)
        if params.shape!=(n,2) or starts.shape!=(n,3) or n<1:raise ValueError('Invalid input shapes')
        if not np.isfinite(params).all() or np.any(params<=0) or not np.isfinite(starts).all():raise ValueError('Nonfinite or nonpositive parameters')
        if transient<0 or min(duration,rtol,step_factor,max_wall,max_steps)<=0:raise ValueError('Invalid numerical settings')
        p=np.zeros((n,4));p[:,:2]=params;p[:,2]=np.minimum(1.,params[:,0])
        state=np.zeros((n,8));state[:,:3]=np.array([to_balanced(y,pa) for y,pa in zip(starts,params)])
        if balanced_starts is not None:state[:,:3]=np.asarray(balanced_starts)
        if not np.isfinite(state).all() or np.any(abs(state[:,0])>1.+1e-10) or np.any(state[:,1]**2+(1-state[:,2])**2>1.+1e-10):raise ValueError('Initial states must lie in the declared invariant cylinder')
        state[:,3:6]=DEFAULT_TANGENT
        if balanced_tangents is not None:
            tangent=np.asarray(balanced_tangents,dtype=float)
            if tangent.shape!=(n,3) or not np.isfinite(tangent).all() or np.any(np.linalg.norm(tangent,axis=1)==0):raise ValueError('Invalid carried tangent')
            state[:,3:6]=tangent/np.linalg.norm(tangent,axis=1)[:,None]
        clock=np.zeros((n,4));clock[:,1]=np.minimum(.001,.1/(1+params.sum(axis=1)))*p[:,2]
        logs=np.zeros((n,4));stats=np.zeros((n,16));stats[:,[0,2,4,9,11,13]]=np.inf;stats[:,[1,3,5,10,12,14]]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[state,p,clock,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        wall=time.perf_counter();last=wall;calls=0
        while True:
            self.kernel(self.queue,(n,),None,*buffers,np.int32(n),np.int32(128),np.float64(transient),np.float64(duration),
                np.float64(rtol),np.float64(rtol*.001),np.float64(step_factor),np.int32(max_steps)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clock,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clock[:,0]<transient+duration)
                if time.perf_counter()-last>25:
                    print(json.dumps(dict(mean_scaled_progress=float(np.mean(clock[:,0])/(transient+duration)),
                        complete=int(np.sum((status==0)&~active)), failed=int(np.sum(status!=0)), wall=round(time.perf_counter()-wall,1))),flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-wall>max_wall:
                    status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for array,buffer in zip(arrays,buffers):cl.enqueue_copy(self.queue,array,buffer).wait()
        rates=p[:,2]
        windows=logs/(duration/4.)*rates[:,None];windows[status!=0]=np.nan
        final=np.array([from_balanced(x,pa) for x,pa in zip(state[:,:3],params)])
        return dict(parameters=params,initial_states=starts,final_states=final,balanced_states=state[:,:3].copy(),
            balanced_tangents=state[:,3:6].copy(),lle=windows.mean(axis=1),lle_windows=windows,
            status=status,elapsed_scaled=clock[:,0].copy(),elapsed_original=clock[:,0]/rates,
            raw_state=state,raw_clock=clock,raw_log_sums=logs,stats=stats,step_counts=counts,
            mean_divergence_original=stats[:,8]/duration*rates,
            settings=dict(transient_scaled=transient,observation_scaled=duration,rtol=rtol,atol=rtol*.001,
                step_factor=step_factor,max_steps_per_lane=max_steps,max_wall_seconds=max_wall,
                time_mapping='s=min(1,rho)*tau; exponent_tau=min(1,rho)*exponent_s',
                state_norm='Euclidean (Q=h/r-Y,Y,Z); final_states retain original (h,Y,Z)',
                stability_cap='h_tau<=step_factor*min(.05,.5/local_Jacobian_row_bound)',
                invariant_check='No clipping; reject accepted step when |Q| or Y^2+(1-Z)^2 exceeds its unit bound by more than 1e-6'),
            stats_columns=['min_Q','max_Q','min_Y','max_Y','min_Z','max_Z','maximum_disk_radius_squared','maximum_abs_Q',
                'observed_integrated_divergence','late_min_H','late_max_H','late_min_Y','late_max_Y','late_min_Z','late_max_Z','unused'],
            status_codes={'0':'completed','1':'invariant-cylinder violation','2':'numerical failure','3':'step budget incomplete','4':'wall budget incomplete'},
            device=self.device.name,wall_seconds=time.perf_counter()-wall,
            caveat='Exact previously published model. Finite-time positivity is preliminary; incomplete lanes are unresolved, not negative results.')


def save(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,**{k:v for k,v in data.items() if isinstance(v,np.ndarray)})
    metadata={k:v for k,v in data.items() if not isinstance(v,np.ndarray)}
    valid=data['status']==0
    preliminary=valid&(data['lle']>.01)&np.all(data['lle_windows']>.005,axis=1)
    metadata.update(lanes=len(valid),completed=int(valid.sum()),preliminary_count=int(preliminary.sum()),
        preliminary_ids=np.flatnonzero(preliminary).tolist(),status_counts={str(i):int(np.sum(data['status']==i)) for i in range(5)},
        source_hash_policy='Canonical LF for Python/OpenCL source; exact bytes for NPZ data',
        source_canonical_lf_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in [ROOT/'research/motor_pair.py',Path(__file__),Path(__file__).with_suffix('.cl')]})
    path.with_suffix('.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:metadata[k] for k in ['lanes','completed','preliminary_count','status_counts','wall_seconds']},indent=2),flush=True)


def validate():
    gpu=MotorPairGPU();records=[];rng=np.random.default_rng(2026103082)
    algebra=[];jacobian=[]
    for _ in range(128):
        p=np.array([10**rng.uniform(-2,2),10**rng.uniform(-.3,4.5)])
        x=rng.uniform(-.5,.5,3);x[2]+=.9
        y=from_balanced(x,p);original=field_jac(y,p)[0];f,j=balanced_field_jac(x,p)
        algebra.append(float(max(abs(f-np.array([original[0]/p[1]-original[1],original[1],original[2]])))/(1+max(abs(f)))))
        transform=np.array([[1/p[1],-1,0],[0,1,0],[0,0,1.]])
        inverse=np.array([[p[1],p[1],0],[0,1,0],[0,0,1.]])
        jacobian.append(float(np.max(abs(j-transform@field_jac(y,p)[1]@inverse))/(1+np.max(abs(j)))))
    cases=[(.05,40.,2.),(.7,400.,.15),(4.5,400.,.15),(6.,7600.,.015),(6.,26815.8,.005),(.01,.5,2.),(100.,.5,.2),(1.,30000.,.005)]
    for rho,r,duration in cases:
        p=np.array([rho,r]);rate=min(1.,rho)
        y=equilibrium_density_start(rho) if (rho,r)==(.05,40.) else from_balanced(np.array([.2,.2,.5]),p)
        x=to_balanced(y,p)
        def fun(t,z):
            f,j=balanced_field_jac(z[:3],p);return np.r_[f,j@z[3:]]
        answers={}
        for method in ['DOP853','Radau']:
            sol=solve_ivp(fun,(0,duration),np.r_[x,DEFAULT_TANGENT],method=method,
                rtol=2e-12,atol=2e-14,max_step=min(.005,duration/100))
            if not sol.success:raise RuntimeError(sol.message)
            end=sol.y[:,-1];norm=np.linalg.norm(end[3:]);answers[method]=(end[:3],end[3:]/norm,np.log(norm)/duration)
        result=gpu.scan(p[None],y[None],transient=0,duration=duration*rate,rtol=1e-10,step_factor=.5,max_wall=30)
        row=dict(parameters=p.tolist(),duration_original=duration,status=int(result['status'][0]),
            solver_state_max_difference=float(np.max(abs(answers['DOP853'][0]-answers['Radau'][0]))),
            solver_lle_difference=float(abs(answers['DOP853'][2]-answers['Radau'][2])),
            gpu_cpu_state_max_difference=float(np.max(abs(result['balanced_states'][0]-answers['DOP853'][0]))),
            gpu_cpu_tangent_max_difference=float(np.max(abs(result['balanced_tangents'][0]-answers['DOP853'][1]))),
            gpu_cpu_lle_difference=float(abs(result['lle'][0]-answers['DOP853'][2])),
            cpu_lle=float(answers['DOP853'][2]),gpu_lle=float(result['lle'][0]))
        row['passed']=bool(row['status']==0 and row['solver_state_max_difference']<2e-8 and row['gpu_cpu_state_max_difference']<2e-8
            and row['gpu_cpu_tangent_max_difference']<2e-8 and row['solver_lle_difference']<1e-6*(1+abs(row['cpu_lle']))
            and row['gpu_cpu_lle_difference']<2e-6*(1+abs(row['cpu_lle'])))
        records.append(row);print(json.dumps(row),flush=True)
    report=dict(algebra_max=max(algebra),jacobian_transform_max=max(jacobian),rows=records,
        all_passed=all(r['passed'] for r in records) and max(algebra)<1e-12 and max(jacobian)<1e-12,
        note='Short controls validate original-time conversion, balanced field/tangent, DOP853 and Radau agreement; not a chaos or source-periodic classification.')
    path=ROOT/'results/motor-pair/gpu-validation.json';path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    if not report['all_passed']:raise RuntimeError('GPU validation failed')


def survey(output,wall=600.,rho_points=24,r_points=40):
    validation=json.loads((ROOT/'results/motor-pair/gpu-validation.json').read_text())
    if not validation['all_passed']:raise RuntimeError('Short validation must pass first')
    regular=json.loads((ROOT/'results/motor-pair/gpu-regular-control-gate.json').read_text())
    if not regular['all_passed']:raise RuntimeError('Periodic and theorem-backed negative controls must pass first')
    starts_balanced=np.array([[.01,0.,0.],[.2,.2,.5],[-.5,.3,1.2]])
    rho_values=np.geomspace(.01,100.,rho_points);r_values=np.geomspace(.5,30000.,r_points)
    tuples=np.array([[rho,r] for rho in rho_values for r in r_values])
    parameters=np.repeat(tuples,3,axis=0);balanced=np.tile(starts_balanced,(len(tuples),1))
    initial=np.array([from_balanced(x,p) for x,p in zip(balanced,parameters)])
    result=MotorPairGPU().scan(parameters,initial,max_wall=wall,balanced_starts=balanced)
    result.update(parameter_tuple_id=np.repeat(np.arange(len(tuples)),3),start_id=np.tile(np.arange(3),len(tuples)),
        survey=dict(rho_range=[.01,100.],r_range=[.5,30000.],rho_points=rho_points,r_points=r_points,
            spacing='logarithmic rectangle, three declared invariant-cylinder starts per tuple',
            balanced_starts=starts_balanced.tolist(),source='Exact Oriola 2014 model; broad exploratory dimensionless ratios, not apparatus calibration'))
    save(output,result)


def regular_controls():
    source=json.loads((ROOT/'results/motor-pair/source-audit.json').read_text())['control']
    cycle=source['shooting'];p=np.array([.05,40.]);state=np.array(cycle['section_state'])
    if not cycle['optimizer_success'] or cycle['radau_closure_max_abs']>1e-4:
        raise RuntimeError('Source periodic control is not sufficiently resolved even for the declared bounded-accuracy gate')
    q=to_balanced(state,p);flow=balanced_field_jac(q,p)[0];flow/=np.linalg.norm(flow)
    gpu=MotorPairGPU()
    periodic=gpu.scan(p[None],state[None],transient=0,duration=cycle['period']*.05,
        rtol=1e-11,step_factor=.5,max_wall=60,balanced_tangents=flow[None])
    save(ROOT/'results/motor-pair/gpu-periodic-control.npz',periodic)
    closure=float(max(abs(periodic['balanced_states'][0]-q)))
    periodic_passed=bool(periodic['status'][0]==0 and closure<5e-5 and abs(periodic['lle'][0])<5e-5)
    parameters=np.array([[rho,.5] for rho in [.01,.05,1.,100.]])
    initial=np.array([from_balanced(np.array([.2,.2,.5]),pa) for pa in parameters])
    negative=gpu.scan(parameters,initial,max_wall=90)
    save(ROOT/'results/motor-pair/gpu-negative-controls.npz',negative)
    spectra=np.array([np.linalg.eigvals(balanced_field_jac(np.zeros(3),pa)[1]).real.max() for pa in parameters])
    discrepancy=np.abs(negative['lle']-spectra)
    negative_passed=bool(np.all(negative['status']==0) and np.all(negative['lle_windows']<0) and
        np.all(discrepancy<.02*(1+np.abs(spectra))))
    report=dict(all_passed=periodic_passed and negative_passed,
        source_control_reference='results/motor-pair/source-audit.json: control.shooting',
        source_shooting_snapshot=cycle,
        source_shooting_content_sha256=hashlib.sha256(json.dumps(cycle,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        periodic=dict(passed=periodic_passed,period_original=cycle['period'],gpu_balanced_closure=closure,
            gpu_original_time_flow_tangent_mean=float(periodic['lle'][0]),
            qualification='One-period GPU consistency check against a bounded-accuracy published-parameter periodic control, not machine-precision proof of orbit closure'),
        negative=dict(passed=negative_passed,theorem='The exact model has global origin attraction for r<=1, rho>0',
            parameters=parameters.tolist(),equilibrium_spectral_abscissae=spectra.tolist(),
            gpu_original_time_leading_estimates=negative['lle'].tolist(),spectral_discrepancy=discrepancy.tolist()))
    (ROOT/'results/motor-pair/gpu-regular-control-gate.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2),flush=True)
    if not report['all_passed']:raise RuntimeError('Regular control gate failed')


def refine(source,output,wall=600.):
    data=np.load(source);ids=np.flatnonzero((data['status']==0)&(data['lle']>.01)&np.all(data['lle_windows']>.005,axis=1))
    if not len(ids):raise RuntimeError('No completed preliminary lanes to refine')
    result=MotorPairGPU().scan(data['parameters'][ids],data['final_states'][ids],transient=400.,duration=800.,
        rtol=1e-10,step_factor=.5,max_wall=wall,balanced_starts=data['balanced_states'][ids],balanced_tangents=data['balanced_tangents'][ids])
    result.update(source_ids=ids,source_npz=str(source),source_npz_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest())
    save(output,result)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--validate',action='store_true');parser.add_argument('--survey',action='store_true')
    parser.add_argument('--regular-controls',action='store_true')
    parser.add_argument('--refine',type=Path);parser.add_argument('--output',type=Path)
    parser.add_argument('--wall',type=float,default=600.);parser.add_argument('--rho-points',type=int,default=24);parser.add_argument('--r-points',type=int,default=40)
    args=parser.parse_args()
    if args.validate:validate()
    elif args.regular_controls:regular_controls()
    elif args.survey:survey(args.output or ROOT/'results/motor-pair/gpu-scout.npz',args.wall,args.rho_points,args.r_points)
    elif args.refine:refine(args.refine,args.output or ROOT/'results/motor-pair/gpu-refined.npz',args.wall)
    else:parser.print_help()
