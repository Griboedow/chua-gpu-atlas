"""Continue a known optical cycle using full-state shooting, not new scouts.

The illumination law is the explicitly assumed smooth logistic profile.
Fixed-width continuation may stop at a fold; no missing branch is filled in.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from smooth_operators import rhs_and_jac

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/light-operators/cycle-continuation'


def shoot(p,y,period,method='DOP853',rtol=2e-11,max_step=.01):
    def rhs(t,z):
        f,j=rhs_and_jac(z[:6],p)
        return np.r_[f,(j@z[6:].reshape(6,6)).ravel()]
    sol=solve_ivp(rhs,(0,float(period)),np.r_[y,np.eye(6).ravel()],
                  method=method,rtol=rtol,atol=rtol*.01,max_step=max_step)
    if not sol.success:raise RuntimeError(sol.message)
    return sol.y[:6,-1],sol.y[6:,-1].reshape(6,6)


def close(p,seed,period,method='DOP853',rtol=2e-11,max_step=.01):
    scale=np.array([1,np.sqrt(p[0]),1,1,np.sqrt(p[0]),1])
    y=np.array(seed,copy=True);y[0]=p[2];T=float(period);history=[]
    for it in range(15):
        end,M=shoot(p,y,T,method,rtol,max_step)
        residual=(end-y)/scale;norm=float(max(abs(residual)));history.append(norm)
        if norm<2e-10:break
        f=rhs_and_jac(end,p)[0]
        jac=np.column_stack([((M-np.eye(6))/scale[:,None]*scale[None,:])[:,1:],f/scale])
        delta=np.linalg.solve(jac,-residual)
        accepted=False
        for fac in [1.,.5,.25,.125,.0625]:
            trial=y.copy();trial[1:]+=fac*delta[:5]*scale[1:];tt=T+fac*delta[-1]
            if not .3<tt<6:continue
            e,_=shoot(p,trial,tt,method,rtol,max_step)
            if max(abs((e-trial)/scale))<norm:
                y,T=trial,tt;accepted=True;break
        if not accepted:raise RuntimeError(f'Newton line search stopped: {norm:g}')
    else:raise RuntimeError(f'Newton iteration budget: {norm:g}')
    f=rhs_and_jac(end,p)[0]
    if y[1]<=0 or abs(f[0])<1e-7:raise RuntimeError('Wrong section orientation or near grazing')
    section=((np.eye(6)-np.outer(f,np.eye(6)[0])/f[0])@M)[1:,1:]
    eig=np.linalg.eigvals(section)
    return dict(parameters=p.tolist(),section_state=y.tolist(),period=T,
        closure_balanced_max=norm,newton_history=history,
        multipliers=[[float(v.real),float(v.imag)] for v in eig],
        spectral_radius=float(max(abs(eig))),flow_multiplier_error=float(np.linalg.norm(M@rhs_and_jac(y,p)[0]-f)/np.linalg.norm(f)),
        liouville_logdet_error=float(abs(np.linalg.slogdet(M)[1]+2*(p[1]+1)*T)),
        method=method,rtol=rtol,atol=rtol*.01,max_step=max_step,
        stability='numerically stable' if max(abs(eig))<1 else 'unstable',
        caveat='Numerical closure and local Floquet stability; no rigorous existence, basin, global branch, or chaos-crisis proof.')


def run(branch='from-narrow',step=.0005,wall_budget=300.):
    OUT.mkdir(parents=True,exist_ok=True)
    controls=json.loads((ROOT/'results/light-operators/lifetimes/periodic-controls.json').read_text())
    ref=controls[0 if branch=='from-narrow' else 2]
    p=np.array(ref['parameters']);seed=np.array(ref['section_state']);T=ref['period']
    widths=(np.arange(.003,.0500001,step) if branch=='from-narrow' else np.arange(.04,.0029999,-step))
    records=[];failure=None;wall=time.perf_counter()
    for width in widths:
        if time.perf_counter()-wall>wall_budget:
            failure=dict(width=float(width),reason='wall budget before next width');break
        p[4]=width
        try:r=close(p,seed,T)
        except Exception as exc:
            failure=dict(width=float(width),reason=f'{type(exc).__name__}: {exc}');break
        seed=np.array(r['section_state']);T=r['period'];records.append(r)
        data=dict(branch=branch,requested_widths=widths.tolist(),step=step,records=records,failure=failure,
                  wall_seconds=time.perf_counter()-wall,wall_budget=wall_budget,
                  source_control_width=ref['parameters'][4],model='Assumed logistic optical-pair model')
        (OUT/f'{branch}.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf8')
        print(f'{branch}: w={width:.6f}, T={T:.9f}, rho={r["spectral_radius"]:.7f}, {time.perf_counter()-wall:.1f}s',flush=True)
    data=dict(branch=branch,requested_widths=widths.tolist(),step=step,records=records,failure=failure,
              wall_seconds=time.perf_counter()-wall,wall_budget=wall_budget,
              source_control_width=ref['parameters'][4],model='Assumed logistic optical-pair model')
    (OUT/f'{branch}.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print(json.dumps(dict(branch=branch,completed=len(records),failure=failure,wall=data['wall_seconds'])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--branch',choices=['from-narrow','from-broad'],default='from-narrow')
    parser.add_argument('--step',type=float,default=.0005);parser.add_argument('--wall',type=float,default=300.)
    a=parser.parse_args();run(a.branch,a.step,a.wall)
