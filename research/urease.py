"""Exact four-species Straube et al. urease network (arXiv:2212.14503v1).

Concentrations are divided by a fixed 1e-4 M; time remains in seconds.
No quasi-steady-state approximation, water-ionization extension or clipping.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1] / 'results' / 'urease'
SCALE = 1e-4
COLUMNS = ['k2', 'k2r', 'KM', 'KE1', 'KE2', 'kS', 'kH', 'k', 'vmax', 'Sext', 'Hext']
TRANSFORM = np.array([[1,0,0,0],[0,-1,1,0],[0,1,0,1],[0,1,0,0]],dtype=float)
INVERSE = np.linalg.inv(TRANSFORM)


def parameters():
    return np.array([4.3e10, 24., .003, 5e-6, 2e-9, .0014, .009, .0014, 1.85e-4, 3.8e-4, 1.3e-4])


def rate(s, h, p):
    km, e1, e2 = p[2:5] / SCALE
    den = h + h*h/e1 + e2
    activity = h/den
    vmax = p[8]/SCALE
    v = vmax*s/(km+s)*activity
    vs = vmax*km/(km+s)**2*activity
    vh = vmax*s/(km+s)*(e2-h*h/e1)/den**2
    return v, vs, vh


def rhs_jac(y, p):
    s,h,a,q = y
    k2,k2r,_,_,_,ks,kh,k,_,se,he = p
    v,vs,vh = rate(s,h,p)
    fast = k2*SCALE
    j = fast*a*h-k2r*q
    f = np.array([-v+ks*(se/SCALE-s), -j+kh*(he/SCALE-h), 2*v-j-k*a, j-k*q])
    jac = np.array([[-vs-ks,-vh,0,0], [0,-fast*a-kh,-fast*h,k2r],
                    [2*vs,2*vh-fast*a,-fast*h-k,k2r], [0,fast*a,fast*h,-k2r-k]])
    return f,jac


def manifold_rhs_jac(z,p):
    """Exact (S,H,P) restriction on N=2*Sext, only when k=kS."""
    if p[7] != p[5]:
        raise ValueError('Exact invariant-plane restriction requires k=kS')
    s,h,a = z
    y = np.r_[z,2*p[9]/SCALE-2*s-a]
    f,j = rhs_jac(y,p)
    lift = np.array([[1,0,0],[0,1,0],[0,0,1],[-2,0,-1]])
    return f[:3],j[:3]@lift


def transport_rhs_jac(z,p):
    """Exact linear coordinates (S, D=P-H, Z=H+Q, H).

    Only H has the large protonation term. The three other equations retain
    their slow balances without cancellation between two large reaction rates.
    """
    s,d,total,h=z
    v,vs,vh=rate(s,h,p)
    k2,k2r,_,_,_,ks,kh,k,_,se,he=p
    fast=k2*SCALE
    f=np.array([-v+ks*(se/SCALE-s),2*v-k*d+(kh-k)*h-kh*he/SCALE,
                kh*he/SCALE-k*total+(k-kh)*h,
                -fast*(d+h)*h+k2r*(total-h)+kh*(he/SCALE-h)])
    j=np.array([[-vs-ks,0,0,-vh],[2*vs,-k,0,2*vh+kh-k],
                [0,0,-k,k-kh],[0,-fast*h,k2r,-fast*(d+2*h)-k2r-kh]])
    return f,j


def save(name, data):
    name=Path(name);name.parent.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(str(name)+'.npz',**arrays)
    meta={k:v for k,v in data.items() if k not in arrays}
    meta['array_shapes']={k:list(v.shape) for k,v in arrays.items()}
    Path(str(name)+'.json').write_text(json.dumps(meta,indent=2,allow_nan=False),encoding='utf-8')


def solve(y,p,duration,method='Radau',rtol=1e-10,atol=1e-13,max_step=2.,reduced=False):
    fun=manifold_rhs_jac if reduced else rhs_jac
    out=solve_ivp(lambda t,z:fun(z,p)[0],(0,duration),y,method=method,
                  jac=lambda t,z:fun(z,p)[1],rtol=rtol,atol=atol,max_step=max_step,dense_output=True)
    if not out.success:
        raise RuntimeError(out.message)
    return out


def validate():
    p=parameters();rng=np.random.default_rng(2026103005)
    err={'jacobian':0.,'mass_balance_scaled':0.,'divergence_relative':0.,'linear_transform_field':0.,'linear_transform_jacobian':0.}
    for _ in range(64):
        y=10**rng.uniform(-4,1,4)
        f,j=rhs_jac(y,p);vs=rate(y[0],y[1],p)[1]
        # Complex steps avoid subtracting huge protonation terms to recover
        # much smaller substrate derivatives in this strongly stiff field.
        delta=1e-20
        fd=np.column_stack([rhs_jac(y.astype(complex)+1j*np.eye(4)[i]*delta,p)[0].imag/delta for i in range(4)])
        err['jacobian']=max(err['jacobian'],float(np.max(abs(j-fd)/(1+abs(j)))))
        expected=p[5]*(2*p[9]/SCALE-2*y[0]-y[2]-y[3])
        err['mass_balance_scaled']=max(err['mass_balance_scaled'],float(abs(2*f[0]+f[2]+f[3]-expected)/(1+max(abs(f)))))
        div=-vs-p[5]-p[6]-p[1]-2*p[7]-p[0]*SCALE*(y[1]+y[2])
        err['divergence_relative']=max(err['divergence_relative'],float(abs(np.trace(j)-div)/(1+abs(div))))
        ft,jt=transport_rhs_jac(TRANSFORM@y,p)
        err['linear_transform_field']=max(err['linear_transform_field'],float(np.max(abs(ft-TRANSFORM@f))/(1+max(abs(f)))))
        err['linear_transform_jacobian']=max(err['linear_transform_jacobian'],float(np.max(abs(jt-TRANSFORM@j@INVERSE))/(1+np.max(abs(j)))))
    y=np.array([.5,.1,0.,0.])
    rad=solve(y,p,1200);bdf=solve(y,p,1200,method='BDF',rtol=1e-11,atol=1e-14,max_step=.5)
    t=np.linspace(0,1200,2401);yr=rad.sol(t).T;yb=bdf.sol(t).T
    err['radau_bdf_state']=float(np.max(abs(yr-yb)/(1+abs(yr))))
    expected=2*p[9]/SCALE+(2*y[0]+y[2]+y[3]-2*p[9]/SCALE)*np.exp(-p[5]*t)
    err['exact_mass_trajectory']=float(np.max(abs(2*yr[:,0]+yr[:,2]+yr[:,3]-expected)))
    ym=np.array([.5,.1,.01,2*p[9]/SCALE-1.-.01])
    full=solve(ym,p,1000);red=solve(ym[:3],p,1000,reduced=True)
    t3=np.linspace(0,1000,2001);yf=full.sol(t3).T;zr=red.sol(t3).T
    lift=np.column_stack([zr,2*p[9]/SCALE-2*zr[:,0]-zr[:,2]])
    err['four_vs_exact_three_state']=float(np.max(abs(yf-lift)/(1+abs(yf))))
    passed=bool(err['jacobian']<1e-7 and max(v for k,v in err.items() if k!='jacobian')<1e-5 and yr.min()>-1e-12 and yf.min()>0)
    out=dict(parameters=p,errors=err,passed=passed,random_cases=64,state_scale_molar=SCALE,time_unit='s',minimum_scaled_state=float(yr.min()),
             methods='Complex-step Jacobian; Radau full four-state vs BDF; separate exact three-state manifold Radau',
             initial_product_policy='P=Q=0 chosen for off-manifold control; source Section2 does not specify them')
    save(ROOT/'validation',out);print({k:v for k,v in out.items() if not isinstance(v,np.ndarray)},flush=True)
    if not passed:raise RuntimeError('Urease validation failed')


def reference():
    p=parameters();y=np.array([.5,.1,0.,0.]);wall=time.perf_counter()
    sol=solve(y,p,12000,rtol=1e-10,atol=1e-13)
    t=np.arange(0,12000,.5);states=sol.sol(t).T
    save(ROOT/'reference',dict(parameters=p,initial_state=y,times=t,states=states,state_scale_molar=SCALE,time_unit='s',
        method='Radau, analytic Jacobian',rtol=1e-10,atol=1e-13,max_step=2.,duration=12000,
        parameter_columns=COLUMNS,state_columns=['S','H','P','Q'],
        initial_condition_scope='Source S(0),H(0); product P(0)=Q(0)=0 chosen by us, not attributed to source',
        wall_seconds=time.perf_counter()-wall))
    last=states[t>=6000]
    print({'wall_seconds':time.perf_counter()-wall,'late_pH_range':[-np.log10(last[:,1].max()*SCALE),-np.log10(last[:,1].min()*SCALE)],'late_scaled_ranges':np.ptp(last,axis=0).tolist()},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--reference',action='store_true');args=parser.parse_args()
    if args.reference:reference()
    else:validate()
