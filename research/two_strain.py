"""Independent implementation of Gavish arXiv:2608.13995v1, Eq. 2a-e.

State integration uses log prevalence to preserve very small invader pools.
Variational equations use the original compartment coordinates, scaled only
by a CONSTANT diagonal matrix. No log-coordinate exponent is reported.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path('results/two-strain')
COLUMNS = ['mu', 'beta1', 'beta2', 'gamma1', 'gamma2', 'gamma12', 'sigma12', 'eta2']


def parameters(r2=1.58, advantage=1.33, r1=2.5, mu=.00025):
    g1, g2, sigma, eta = 1., 1.1, 1.1, 1.1
    return np.array([mu, r1*(g1+mu), r2*(g2+mu), g1, g2,
                     sigma*eta*g2/advantage, sigma, eta])


def rhs_jac(y, p):
    """Original physical state order S,I1,I2,I12,R1 (fractions)."""
    mu, b1, b2, g1, g2, g12, sigma, eta = p
    s, i1, i2, i12, r1 = y
    force = b2*(i2+eta*i12)
    f = np.array([mu*(1-s)-b1*i1*s-force*s,
        (b1*s-mu-g1)*i1, force*s-(mu+g2)*i2,
        sigma*force*r1-(mu+g12)*i12, g1*i1-sigma*force*r1-mu*r1])
    j = np.array([[-mu-b1*i1-force, -b1*s, -b2*s, -b2*eta*s, 0],
        [b1*i1, b1*s-mu-g1, 0, 0, 0],
        [force, 0, b2*s-mu-g2, b2*eta*s, 0],
        [0, 0, sigma*b2*r1, sigma*b2*eta*r1-mu-g12, sigma*force],
        [0, g1, -sigma*b2*r1, -sigma*b2*eta*r1, -sigma*force-mu]])
    return f, j


def coexistence(p):
    """Closed form Eq. 12, checked against the literal field."""
    mu, b1, b2, g1, g2, g12, sigma, eta = p
    a, b = b1/(g1+mu), b2/(g2+mu)
    r = (g12+mu)*(a-b)/((g2+mu)*eta*sigma*a*b)
    invasion = b/a*(1+g1*(g2+mu)*eta*sigma*(a-1)/((g1+mu)*(g12+mu)))
    c = (g1+mu)*(g12+mu)/((g1+mu)*(g12+mu)*(a-b)+g1*eta*(g2+mu)*b)
    i2 = mu*c/((g2+mu)*sigma)*(invasion-1)
    i1 = (1+sigma*(a-1))*mu*r/(g1+(g1+mu)*sigma*a*r)
    result = np.array([1/a, i1, i2, (a-b)/(eta*b)*i2, r])
    if min(result) <= 0 or sum(result) >= 1:
        raise ValueError('No interior coexistence equilibrium at these parameters')
    return result


def initial(p, kind):
    if kind == 'coexistence':
        y = coexistence(p)
        old = y[1]
        y[1] *= .9
        return y/(1-.1*old)
    if kind == 'invasion':
        mu, _, b2, _, g2, *_ = p
        r2 = b2/(g2+mu)
        y = np.array([1/r2, 0, mu*(r2-1)/((g2+mu)*r2), 0, 0])
        y /= 1+1e-5
        y[1] = 1e-5/(1+1e-5)
        return y
    raise ValueError(kind)


def encode(y, p):
    s, i1, i2, i12, r1 = y
    if min(y) < 0 or sum(y) > 1+1e-12 or i1 <= 0 or i2+i12 <= 0:
        raise ValueError('Need positive total prevalence in the population simplex')
    return np.array([s, np.log(i1/p[0]), np.log((i2+i12)/p[0]), i12/(i2+i12), r1])


def decode(q, p):
    s, u, v, ratio, r1 = q
    i1, total = p[0]*np.exp(u), p[0]*np.exp(v)
    return np.array([s, i1, total*(1-ratio), total*ratio, r1])


def log_rhs(q, p):
    mu, b1, b2, g1, g2, g12, sigma, eta = p
    s, u, v, ratio, r1 = q
    i1, total = mu*np.exp(u), mu*np.exp(v)
    a = 1+(eta-1)*ratio
    return np.array([mu*(1-s)-b1*i1*s-b2*total*a*s,
        b1*s-mu-g1,
        b2*a*(s+sigma*r1)-mu-g2-(g12-g2)*ratio,
        b2*a*(sigma*r1-ratio*(s+sigma*r1))+(g2-g12)*ratio*(1-ratio),
        g1*i1-sigma*b2*total*a*r1-mu*r1])


def scale(p):
    return np.array([1., p[0], p[0], p[0], 1.])


def balanced_jac(q, p):
    j = rhs_jac(decode(q, p), p)[1]
    d = scale(p)
    return j*d[None, :]/d[:, None]


def log_jac(q,p):
    """Analytic Jacobian of log_rhs, useful for well-conditioned full QR.

    Its exponents must be labelled by coordinates. Equivalence to physical
    asymptotic exponents requires control of the coordinate transformation.
    """
    mu,b1,b2,g1,g2,g12,sigma,eta=p
    s,u,v,r,r1=q;i1=mu*np.exp(u);total=mu*np.exp(v)
    e=eta-1;a=1+e*r;force=b2*total*a;h=s+sigma*r1;d=g12-g2
    return np.array([[-mu-b1*i1-force,-b1*i1*s,-force*s,-b2*total*e*s,0],
        [b1,0,0,0,0],
        [b2*a,0,0,b2*e*h-d,b2*a*sigma],
        [-b2*a*r,0,0,b2*(e*(sigma*r1-r*h)-a*h)-d*(1-2*r),b2*a*sigma*(1-r)],
        [0,g1*i1,-sigma*force*r1,-sigma*b2*total*e*r1,-sigma*force-mu]])


def save(path, data):
    path = Path(path)
    if path.suffix == '.npz':path=path.with_suffix('')
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(path, **arrays)
    info = {k:v for k,v in data.items() if not isinstance(v,np.ndarray)}
    Path(str(path)+'.json').write_text(json.dumps(info,indent=2,allow_nan=False))


def validate():
    rng = np.random.default_rng(2026102970)
    errors = dict(jacobian=0., transformed_field=0., equilibrium=0., population_balance=0.,log_jacobian=0.,trace_coordinate_identity=0.)
    for _ in range(32):
        p = parameters(rng.uniform(1.25,2.1), rng.uniform(1.2,1.8))
        y = rng.dirichlet(np.ones(6))[:5]
        f,j = rhs_jac(y,p)
        eps=1e-6
        fd=np.column_stack([(rhs_jac(y+np.eye(5)[i]*eps,p)[0]-rhs_jac(y-np.eye(5)[i]*eps,p)[0])/(2*eps) for i in range(5)])
        errors['jacobian']=max(errors['jacobian'],float(np.max(abs(fd-j))))
        q=encode(y,p);dq=log_rhs(q,p)
        logfd=np.column_stack([(log_rhs(q+np.eye(5)[i]*eps,p)-log_rhs(q-np.eye(5)[i]*eps,p))/(2*eps) for i in range(5)])
        errors['log_jacobian']=max(errors['log_jacobian'],float(np.max(abs(logfd-log_jac(q,p)))))
        errors['trace_coordinate_identity']=max(errors['trace_coordinate_identity'],float(abs(np.trace(j)-np.trace(log_jac(q,p))-dq[1]-2*dq[2])))
        total=y[2]+y[3]
        reconstructed=np.array([dq[0],y[1]*dq[1],y[2]*dq[2]-total*dq[3],y[3]*dq[2]+total*dq[3],dq[4]])
        errors['transformed_field']=max(errors['transformed_field'],float(np.max(abs(f-reconstructed))))
        ce=coexistence(p)
        errors['equilibrium']=max(errors['equilibrium'],float(np.max(abs(rhs_jac(ce,p)[0]))))
        recovered=1-y.sum()
        balance=f.sum()+p[4]*y[2]+p[5]*y[3]-p[0]*recovered
        errors['population_balance']=max(errors['population_balance'],abs(float(balance)))
    p=parameters(1.9);y=initial(p,'invasion');q=encode(y,p)
    direct=solve_ivp(lambda t,y:rhs_jac(y,p)[0],(0,100),y,method='DOP853',rtol=1e-12,atol=1e-15,max_step=.1)
    log=solve_ivp(lambda t,q:log_rhs(q,p),(0,100),q,method='Radau',rtol=1e-11,atol=1e-13,max_step=.5)
    errors['independent_physical_vs_log']=float(np.max(abs(direct.y[:,-1]-decode(log.y[:,-1],p))))
    from two_strain_gpu import TwoStrainGPU
    gpu=TwoStrainGPU()
    result=gpu.scan(np.array([p]),np.array([q]),transient=0,duration=100,rtol=1e-10,max_step=.25)
    z0=np.r_[q,np.array([1.,.3,.4,-.2,.5])/np.sqrt(1.54)]
    def tangent(t,z):return np.r_[log_rhs(z[:5],p),balanced_jac(z[:5],p)@z[5:]]
    cpu=solve_ivp(tangent,(0,100),z0,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.1)
    errors['gpu_log_state']=float(np.max(abs(result['final_log_states'][0]-cpu.y[:5,-1])))
    errors['gpu_physical_lle']=float(abs(result['lle'][0]-np.log(np.linalg.norm(cpu.y[5:,-1]))/100))
    passed=max(errors.values())<1e-7 and result['status'][0]==0 and direct.success and log.success and cpu.success
    out=dict(errors=errors,passed=bool(passed),random_cases=32,
        tangent_coordinates='Original five compartments with fixed infection scale 1/mu; not log variables',device=gpu.device.name)
    ROOT.mkdir(parents=True,exist_ok=True)
    (ROOT/'validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
    if not passed:raise RuntimeError('Model/solver validation failed')


def source_controls():
    rows=[]
    for r2,kind in [(1.58,'coexistence'),(1.9,'coexistence'),(1.9,'invasion')]:
        p=parameters(r2);q0=encode(initial(p,kind),p)
        wall=time.perf_counter()
        sol=solve_ivp(lambda t,q:log_rhs(q,p),(0,50000),q0,method='DOP853',rtol=1e-10,atol=1e-12,max_step=1.,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        times=np.arange(30000,50000.1,.5)
        q=sol.sol(times).T;y=np.array([decode(v,p) for v in q])
        out=dict(parameters=p,initial_log_state=q0,times=times,log_states=q,states=y,final_log_state=sol.y[:,-1],
            r2=r2,start_kind=kind,method='DOP853',rtol=1e-10,atol=1e-12,max_step=1.,
            min_log_prevalence=float(np.min(q[:,1:3])),min_compartment=float(y.min()),
            max_population_sum=float(y.sum(axis=1).max()),wall_seconds=time.perf_counter()-wall,
            scope='Source-parameter trajectory control; not a chaos diagnosis; q_advantage=1.33 interpreted literally')
        name=f'source-r2-{r2:g}-{kind}'
        save(ROOT/name,out)
        rows.append({k:v for k,v in out.items() if not isinstance(v,np.ndarray)})
        print(json.dumps(rows[-1]),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true');args=parser.parse_args()
    if args.controls:source_controls()
    else:validate()
