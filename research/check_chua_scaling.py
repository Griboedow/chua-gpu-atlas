"""Independent full-field DOP853 spectrum/sections in exact scaled coordinates."""
import argparse
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from chua_scaling_gpu import OUT, field_jac, save


def check(index,duration=40.,transient=10.,max_fast_step=.15,rtol=2e-11,wall=360.,label=None,align_tangent=True):
    source=np.load(OUT/'scout.npz');p=source['parameters'][index];y=source['final_states'][index].copy()
    alpha=p[0];start=time.perf_counter();deadline=start+wall
    eps=1/alpha;kappa=p[1];m0=p[2];m1=p[3]
    def rhs(s,z):
        if time.perf_counter()>deadline:raise TimeoutError('CPU wall budget; completed QR intervals only are reported')
        x,v,Z=z[:3];slope=m0 if abs(x)<1 else m1
        out=np.empty(13)
        out[:3]=[v-(1+m1)*x-(m0-m1)*max(-1,min(1,x)),Z+eps*(x-v),-kappa*v]
        out[3:6]=-(1+slope)*z[3:6]+z[6:9]
        out[6:9]=eps*(z[3:6]-z[6:9])+z[9:12]
        out[9:12]=-kappa*z[6:9];out[12]=-1-slope-eps
        return out
    def state_rhs(s,z):
        if time.perf_counter()>deadline:raise TimeoutError('CPU transient wall budget')
        x,v,Z=z
        return [v-(1+m1)*x-(m0-m1)*max(-1,min(1,x)),Z+eps*(x-v),-kappa*v]
    records=[];st=[];sy=[];tt=[];yy=[];sums=np.zeros(3);divergence=0.
    Q=np.linalg.qr(np.array([[1.,3.,-2.],[2.,-2.,1.],[3.,1.,4.]]))[0] if align_tangent else np.eye(3)
    s=0.;status='completed'
    def section(s,z):return z[1]
    section.direction=1
    initial=y.copy();observation_initial=None;flownorm0=0.
    windowlogs=np.zeros((4,3));windowtime=np.zeros(4)
    rng=np.random.default_rng(20261030);field_error=0.
    for _ in range(64):
        trial=rng.normal(size=13)*3
        f,j=field_jac(trial[:3],p)
        expected=np.r_[f,(j@trial[3:12].reshape(3,3)).ravel(),np.trace(j)]
        field_error=max(field_error,float(np.max(abs(rhs(0,trial)-expected)/(1+abs(expected)))))
    if field_error>1e-12:raise RuntimeError('Optimized augmented field identity failed')
    try:
        if transient:
            remaining=transient*alpha
            while remaining>1e-12:
                length=min(5.,remaining) if align_tangent else remaining
                warm=solve_ivp(rhs if align_tangent else state_rhs,(0,length),np.r_[y,Q.ravel(),0.] if align_tangent else y,
                               method='DOP853',rtol=rtol,atol=rtol*.001,max_step=max_fast_step)
                if not warm.success:raise RuntimeError(warm.message)
                y=warm.y[:3,-1]
                if align_tangent:Q=np.linalg.qr(warm.y[3:12,-1].reshape(3,3))[0]
                remaining-=length
        observation_initial=y.copy();flow0=field_jac(y,p)[0];flownorm0=np.linalg.norm(flow0)
        # Fast QR spacing prevents loss of strongly contracting outer direction.
        qr=5.;target=duration*alpha
        quarter=duration*alpha/4;nextboundary=quarter;window=0
        while s<target*(1-1e-14):
            end=min(s+qr,target,nextboundary)
            # Autonomous field: resetting the local clock removes accumulated
            # large-s spacing constraints near a nonsmooth Jacobian crossing.
            sol=solve_ivp(rhs,(0,end-s),np.r_[y,Q.ravel(),0.],method='DOP853',rtol=rtol,atol=rtol*.001,
                          max_step=max_fast_step,events=section,dense_output=True)
            if not sol.success:raise RuntimeError(sol.message)
            endpoint=sol.y[:,-1];Q,R=np.linalg.qr(endpoint[3:12].reshape(3,3));lg=np.log(np.abs(np.diag(R)))
            sums+=lg;divergence+=endpoint[-1];windowlogs[window]+=lg;windowtime[window]+=(end-s)/alpha
            # Fixed fast-time samples retain physical coordinates, not warped projection.
            grid=np.linspace(0,end-s,max(2,int(np.ceil((end-s)/.2))+1))[:-1]
            tt.extend(((s+grid)/alpha).tolist());yy.extend(sol.sol(grid)[:3].T.tolist())
            if len(sol.t_events[0]):st.extend(((s+sol.t_events[0])/alpha).tolist());sy.extend(sol.y_events[0][:,:3].tolist())
            y=endpoint[:3];s=end
            if abs(s-nextboundary)<1e-8:
                records.append(dict(original_elapsed=s/alpha,spectrum=(sums/(s/alpha)).tolist()))
                window+=1;nextboundary=quarter*(window+1)
    except (TimeoutError,RuntimeError) as exc:status=str(exc)
    elapsed=s/alpha
    spectrum=sums/elapsed if elapsed else np.full(3,np.nan)
    # The neutral-flow direction can have a nonzero finite-time endpoint rate.
    flowrate=(np.log(np.linalg.norm(field_jac(y,p)[0]))-np.log(flownorm0))/elapsed if elapsed and flownorm0>0 else None
    out=dict(source_index=index,parameters=p,initial_state=initial,observation_initial_state=observation_initial,
             final_state=y,spectrum=spectrum,qr_history=records,
             window_spectra=np.divide(windowlogs,windowtime[:,None],out=np.full((4,3),np.nan),where=windowtime[:,None]>0),
             window_elapsed=windowtime,elapsed_original_time=elapsed,requested_duration=duration,requested_transient=transient,
             integrated_divergence_original_rate=float(divergence/elapsed) if elapsed else None,
             sum_spectrum_minus_trace=float(spectrum.sum()-divergence/elapsed) if elapsed else None,
             finite_time_neutral_flow_rate=flowrate,status=status,wall_seconds=time.perf_counter()-start,
             time=np.asarray(tt),states=np.asarray(yy),section_time=np.asarray(st),section_states=np.asarray(sy),
             settings=dict(method='DOP853 full piecewise field/Jacobian, fast QR interval5 with local clock reset, y=0 upward sections',
                           rtol=rtol,atol=rtol*.001,max_fast_step=max_fast_step,wall_budget=wall,
                           tangent_initialization='QR basis whose first vector is proportional to (1,2,3); QR evolves during transient' if align_tangent else 'Identity Q after state-only warm-up; historical finite-alignment diagnostic'),
             scope='Independent CPU integrator; same exact scaled equations. Finite-time spectrum, not a proof of a strange attractor.',
             optimized_augmented_field_identity_relative_error=field_error,
             time_unit='original dimensionless time; state=(x,y,Z=z/alpha)')
    name=label or f'cpu-{index}'
    save(OUT/name,out)
    print({k:out[k] for k in ['source_index','spectrum','elapsed_original_time','sum_spectrum_minus_trace','finite_time_neutral_flow_rate','status','wall_seconds']},flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('index',type=int);ap.add_argument('--duration',type=float,default=40.)
    ap.add_argument('--transient',type=float,default=10.);ap.add_argument('--step',type=float,default=.15)
    ap.add_argument('--rtol',type=float,default=2e-11);ap.add_argument('--wall',type=float,default=360.)
    ap.add_argument('--label');ap.add_argument('--state-only-warm',action='store_true');a=ap.parse_args()
    check(a.index,a.duration,a.transient,a.step,a.rtol,a.wall,a.label,not a.state_only_warm)
