"""Literal four-state Shinbo et al. model, DOI 10.1002/syst.70050.

The paper analyzes its two-state quasi-steady reduction. Here substrate
relaxations remain dynamic; lambda_g is an unspecified positive source rate.
No finite-bath depletion, new feedback, or concentration clipping is added.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1] / 'results' / 'core-shell-gel'
COLUMNS = ['A', 'C', 'gamma', 'rho', 'B', 'epsilon', 'lambda_g', 'Hill_n']
HC, K1, K2, KA = 5.62e-5, 5e-6, 2e-9, 1.404e-5
PMIN, PMAX, ELASTIC, PHIC, MH, X0 = .1992, .8544, 25., .3, .4, 1.


def parameters(lam=25.):
    return np.array([3., 110., 25., 5.68, 7.2, .5, lam, 4.568])


def response(x, n):
    # Algebraically regular at x=0. Negative trial stages are rejected by GPU.
    if x < 0:
        raise ValueError('Negative proton concentration')
    aa, bb, kk = K2/HC, HC/K1, KA/HC
    den = x + aa + bb*x*x
    fu = x/den
    dfu = (aa-bb*x*x)/(den*den)
    xn, kn = x**n, kk**n
    gate = PMIN+(PMAX-PMIN)*xn/(xn+kn)
    dg = (PMAX-PMIN)*n*kn*x**(n-1)/(xn+kn)**2
    return fu*gate, dfu*gate+fu*dg


def rhs_jac(y, p):
    # y=(G,U,x,phi), G=k_Rg*g/(k_Dh*h_c), U=2*k_Ru*u/(k_Dh*h_c).
    G,U,x,phi = y
    A,C,gamma,rho,B,eps,lam,n = p
    f,df = response(x,n)
    target=A*np.exp(-B*phi)
    field=np.array([lam*(target-G), C-(gamma+rho*f)*U,
                    X0-x+G-f*U, (phi-PHIC-ELASTIC*(phi-PHIC)**3+MH*(x-1))/eps])
    jac=np.array([[-lam,0,0,-lam*B*target], [0,-gamma-rho*f,-rho*df*U,0],
                  [1,-f,-1-df*U,0], [0,0,MH/eps,(1-3*ELASTIC*(phi-PHIC)**2)/eps]])
    return field,jac


def start(p,x=1.,phi=.2):
    f,_=response(x,p[7])
    return np.array([p[0]*np.exp(-p[4]*phi),p[1]/(p[2]+p[3]*f),x,phi])


def reduced_rhs(t,y,p):
    x,phi=y;f,_=response(x,p[7])
    return [X0-x+p[0]*np.exp(-p[4]*phi)-p[1]*f/(p[2]+p[3]*f),
            (phi-PHIC-ELASTIC*(phi-PHIC)**3+MH*(x-1))/p[5]]


def save(name,data):
    name=Path(name);name.parent.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(str(name)+'.npz',**arrays)
    meta={k:v for k,v in data.items() if k not in arrays}
    meta['array_shapes']={k:list(v.shape) for k,v in arrays.items()}
    Path(str(name)+'.json').write_text(json.dumps(meta,indent=2,allow_nan=False),encoding='utf-8')


def validate():
    from core_shell_gel_gpu import GelGPU
    rng=np.random.default_rng(2026103001);errors={'jacobian':0.,'quasi_steady_identity':0.}
    p=parameters()
    for _ in range(64):
        y=np.array([rng.uniform(.02,4),rng.uniform(.1,8),rng.uniform(.05,3),rng.uniform(.05,.7)])
        f,j=rhs_jac(y,p);h=1e-6
        fd=np.column_stack([(rhs_jac(y+np.eye(4)[i]*h,p)[0]-rhs_jac(y-np.eye(4)[i]*h,p)[0])/(2*h) for i in range(4)])
        errors['jacobian']=max(errors['jacobian'],float(np.max(abs(j-fd)/(1+abs(j)))))
        q=start(p,y[2],y[3]);ff=rhs_jac(q,p)[0]
        errors['quasi_steady_identity']=max(errors['quasi_steady_identity'],float(max(abs(ff[:2]))),float(max(abs(ff[2:]-reduced_rhs(0,q[2:],p)))))
    y=start(p);v=np.array([1.,.3,-.2,.4]);v/=np.linalg.norm(v)
    def aug(t,z):
        f,j=rhs_jac(z[:4],p);return np.r_[f,j@z[4:]]
    cpu=solve_ivp(aug,(0,3),np.r_[y,v],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.02)
    rad=solve_ivp(lambda t,y:rhs_jac(y,p)[0],(0,3),y,method='Radau',jac=lambda t,y:rhs_jac(y,p)[1],rtol=1e-11,atol=1e-13)
    gpu=GelGPU();r=gpu.scan(p[None,:],y[None,:],transient=0,duration=3,rtol=1e-10,max_step=.02)
    errors.update(gpu_state=float(max(abs(r['final_states'][0]-cpu.y[:4,-1]))),gpu_lle=float(abs(r['lle'][0]-np.log(np.linalg.norm(cpu.y[4:,-1]))/3)),dop853_vs_radau=float(max(abs(cpu.y[:4,-1]-rad.y[:,-1]))))
    out=dict(errors=errors,passed=bool(max(errors.values())<1e-6 and r['status'][0]==0 and cpu.success and rad.success),device=gpu.device.name,random_cases=64)
    ROOT.mkdir(parents=True,exist_ok=True);(ROOT/'validation.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
    if not out['passed']:raise RuntimeError('Gel validation failed')


def source_control():
    p=parameters();wall=time.perf_counter()
    sol=solve_ivp(lambda t,y:reduced_rhs(t,y,p),(0,200),[1.,.2],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.03,dense_output=True)
    t=np.arange(100,200,.01);states=sol.sol(t).T
    save(ROOT/'source-reduced-control',dict(parameters=p,times=t,states=states,success=bool(sol.success),method='DOP853',rtol=1e-12,atol=1e-14,max_step=.03,initial_state=[1.,.2],time_range=[0,200],wall_seconds=time.perf_counter()-wall,scope='Literal Table 1 and Eqs. 10-11, two-dimensional control; initial condition chosen by us'))
    print({'success':sol.success,'range':np.ptp(states,axis=0).tolist()},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source-control',action='store_true');a=parser.parse_args()
    if a.source_control:source_control()
    else:validate()
