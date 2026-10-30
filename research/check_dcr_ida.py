"""Short, bounded IDA controls of the literal strict-contact author-code variant."""
import argparse
from dataclasses import asdict
import time
import numpy as np
import sksundae
from sksundae.ida import IDA
from dcr import Parameters,ROOT,field,residual,residual_jacobian,geometry,save


def run(alpha=1.,beta=1.,duration=.2,step=1e-4,tol=1e-7,analytic=False,wall=120.,label='ida-control'):
    p=Parameters(alpha=alpha,beta=beta);rng=np.random.default_rng(2026103021)
    y0=np.r_[0.,rng.uniform(0,2*np.pi,p.n)];yp0=field(y0,p)
    start=time.perf_counter();calls=0;last=[0.,y0.copy(),yp0.copy()]
    def res(t,y,yp,out):
        nonlocal calls
        calls+=1
        if time.perf_counter()-start>wall:raise TimeoutError('CPU callback wall budget')
        out[:]=residual(y,yp,p)
        if calls%100==0:last[:]=[float(t),y.copy(),yp.copy()]
    def jac(t,y,yp,r,cj,out):out[:]=residual_jacobian(y,yp,cj,p)
    options=dict(rtol=tol,atol=tol,max_step=step,max_num_steps=200000,linsolver='dense')
    if analytic:options['jacfn']=jac
    grid=np.linspace(0,duration,max(3,int(np.ceil(duration/.001))+1))
    solver=IDA(res,**options)
    try:
        solution=solver.solve(grid,y0,yp0)
        data=dict(time=solution.t,states=solution.y,derivatives=solution.yp,
                  success=bool(solution.success),status=int(solution.status),message=solution.message,
                  nfev=int(solution.nfev),njev=int(solution.njev))
        counts=np.array([geometry(y,p)[0].sum() for y in solution.y])
        data['contacts']=counts
    except Exception as exc:
        data=dict(success=False,status='exception',message=f'{type(exc).__name__}: {exc}',
                  last_callback_time=last[0],last_callback_state=last[1],last_callback_derivative=last[2])
    data.update(parameters=asdict(p),initial_state=y0,initial_derivative=yp0,seed=2026103021,
                settings=dict(requested_duration=duration,max_step=step,rtol=tol,atol=tol,analytic_branch_jacobian=analytic,
                              wall_budget=wall,contact_policy='Literal strict author-code predicates; no capture/Filippov convention added',solver='SUNDIALS IDA via scikit-sundae '+sksundae.__version__),
                callback_calls=calls,wall_seconds=time.perf_counter()-start,
                scope='Short source-solver control, not reproduction of a long Fig5 regime or a chaos test.')
    save(ROOT/label,data)
    print({k:v for k,v in data.items() if not isinstance(v,np.ndarray) and k not in ['settings','parameters']},flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--alpha',type=float,default=1);ap.add_argument('--beta',type=float,default=1)
    ap.add_argument('--duration',type=float,default=.2);ap.add_argument('--step',type=float,default=1e-4)
    ap.add_argument('--tol',type=float,default=1e-7);ap.add_argument('--analytic',action='store_true')
    ap.add_argument('--wall',type=float,default=120);ap.add_argument('--label',default='ida-control');a=ap.parse_args()
    run(a.alpha,a.beta,a.duration,a.step,a.tol,a.analytic,a.wall,a.label)
