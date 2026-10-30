"""Literal four-state S5-S13 closure of Case, Angilella & Motter (2020).

State (Q1,Q2,Q3,r); parameters are the two dimensionless, constant inlet
pressures. All geometric/fitted coefficients retain their published values.
Source time tau corresponds to 0.25 physical seconds. No flow clipping.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1] / 'results/microfluidic'
DEFAULT = np.array([62500., 62500.])
L = np.array([24., 144., 24., 240., 120.])
A = np.array([[1/L[0]+1.08/L[2]+1/L[4], -1/L[2]],
              [-1.08/L[2], 1/L[1]+1/L[3]+1/L[2]]])
AI = np.linalg.inv(A)


def components(s, p):
    x, y, z, r = s
    eta = .12 + 8e-6*x
    xi = (5.5**3*eta-25)/5.5
    resistance = 42+7.35*(r-4.1)**2
    mixing = z-z*z/x
    a3 = -resistance*z/24+15.75*mixing
    a4 = -(31+.2565*(y+z))*(y+z)/12
    a5 = -(x-z)+12*mixing
    pressures = AI @ np.array([p[0]/24-x-a3-a5, p[1]/144-y-a4+a3])
    return pressures, np.array([a3,a4,a5]), eta, xi, resistance


def field_jac(s, p):
    x,y,z,r = s
    (u,v), (a3,a4,a5), eta,xi,resistance = components(s,p)
    f = np.array([12*((p[0]-u)/24-x), 12*((p[1]-v)/144-y),
                  12*((1.08*u-v)/24+a3),
                  10.9*(z-25-eta*(r-5.5)**3+xi*(r-5.5))])
    ex,ey,ez,er = np.eye(4)
    dm = np.array([z*z/(x*x),0,1-2*z/x,0])
    da3 = -resistance/24*ez-14.7*(r-4.1)*z/24*er+15.75*dm
    da4 = -(31+2*.2565*(y+z))/12*(ey+ez)
    da5 = -ex+ez+12*dm
    du,dv = AI @ np.vstack([-ex-da3-da5, -ey-da4+da3])
    j = np.vstack([12*(-du/24-ex), 12*(-dv/144-ey),
                   12*((1.08*du-dv)/24+da3),
                   10.9*(8e-6*(-(r-5.5)**3+5.5**2*(r-5.5))*ex
                         +ez+(-3*eta*(r-5.5)**2+xi)*er)])
    return f,j


def original_field(s,p):
    """Independent five-flow pressure-constraint reconstruction (Q4,Q5 explicit)."""
    x,y,z,r = s
    q4,q5 = y+z,x-z
    rr = 42+7.35*(r-4.1)**2
    # Unknowns P3,P4 are solved anew; no shared inverse or reduced field.
    b3 = -rr*z/24+15.75*z*q5/x
    b4 = -(31+.2565*q4)*q4/12
    b5 = -q5+12*z*q5/x
    mat = np.array([[1/24+1.08/24+1/120,-1/24],
                    [-1.08/24,1/144+1/240+1/24]])
    u,v = np.linalg.solve(mat, [p[0]/24-x-b3-b5,p[1]/144-y-b4+b3])
    five = 12*np.array([(p[0]-u)/24-x,(p[1]-v)/144-y,
                        (1.08*u-v)/24+b3,v/240+b4,u/120+b5])
    vortex = 10.9*(z-r*(25/5.5+(.12+8e-6*x)*(r-5.5)*(r-11)))
    return np.r_[five[:3],vortex],five


def save(name,data):
    name=Path(name);name.parent.mkdir(parents=True,exist_ok=True)
    arrays={k:v for k,v in data.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(str(name)+'.npz',**arrays)
    meta={k:v for k,v in data.items() if k not in arrays}
    meta['array_shapes']={k:list(v.shape) for k,v in arrays.items()}
    Path(str(name)+'.json').write_text(json.dumps(meta,indent=2,allow_nan=False),encoding='utf-8')


def validate():
    rng=np.random.default_rng(2026103016)
    errors=dict(jacobian=0.,five_flow_reconstruction=0.,flow_constraints=0.)
    for _ in range(128):
        x=rng.uniform(10,3000);s=np.array([x,rng.uniform(1,300),rng.uniform(.01,.9)*x,rng.uniform(.1,10)])
        p=rng.uniform(1e3,2e5,2);f,j=field_jac(s,p)
        ref=np.column_stack([field_jac(s.astype(complex)+1e-30j*np.eye(4)[i],p)[0].imag/1e-30 for i in range(4)])
        errors['jacobian']=max(errors['jacobian'],float(np.max(abs(j-ref)/(1+abs(j)))))
        other,five=original_field(s,p)
        errors['five_flow_reconstruction']=max(errors['five_flow_reconstruction'],float(np.max(abs(f-other))/(1+np.max(abs(f)))))
        errors['flow_constraints']=max(errors['flow_constraints'],float(max(abs(five[0]-five[2]-five[4]),abs(five[1]-five[3]+five[2]))/(1+np.max(abs(five)))))
    passed=max(errors.values())<1e-10
    save(ROOT/'validation',dict(passed=bool(passed),errors=errors,random_states=128,seed=2026103016))
    print(dict(passed=passed,errors=errors),flush=True)
    if not passed:raise RuntimeError('Microfluidic algebra/Jacobian validation failed')


def source_control():
    p=DEFAULT.copy();starts=np.array([[650,50,25,5.],[650,50,25,3.],[550,50,20,2.]])
    records=[]
    for k,initial in enumerate(starts):
        sol=solve_ivp(lambda t,s:field_jac(s,p)[0],(0,100),initial,
                      method='DOP853',rtol=1e-10,atol=1e-12,max_step=.01,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        times=np.linspace(90,100,10001);states=sol.sol(times).T
        record=dict(start=k,success=bool(sol.success),state_minimum=states.min(axis=0).tolist(),
                    state_maximum=states.max(axis=0).tolist(),endpoint=sol.y[:,-1].tolist())
        records.append(record);print(record,flush=True)
        save(ROOT/f'source-start{k}',dict(parameters=p,initial_state=initial,times=times,states=states,record=record,
             solver=dict(method='DOP853',rtol=1e-10,atol=1e-12,max_step=.01,duration=100),
             scope='Our three initial conditions, not author-supplied; source Fig.S5 parameter setting 250 Pa each.'))


def independent_control():
    initial=np.array([650.,50.,25.,5.]);p=DEFAULT.copy()
    a=solve_ivp(lambda t,s:field_jac(s,p)[0],(0,10),initial,method='DOP853',
                rtol=1e-11,atol=1e-13,max_step=.005,dense_output=True)
    b=solve_ivp(lambda t,s:original_field(s,p)[0],(0,10),initial,method='Radau',
                rtol=1e-11,atol=1e-13,max_step=.005,dense_output=True)
    if not (a.success and b.success):raise RuntimeError('Independent integration failed')
    times=np.linspace(0,10,10001);x=a.sol(times).T;y=b.sol(times).T
    error=float(np.max(abs(x-y)/(1+abs(y))))
    eig=np.linalg.eigvals(field_jac(x[-1],p)[1])
    passed=bool(a.success and b.success and error<1e-8)
    save(ROOT/'independent-source-control',dict(passed=passed,relative_sup_error=error,
         parameters=p,initial_state=initial,times=times,dop853_states=x,
         independent_radau_states=y,equilibrium_eigenvalues=eig,
         settings=dict(duration=10,rtol=1e-11,atol=1e-13,max_step=.005),
         scope='Independent five-flow pressure reconstruction Radau vs reduced-field DOP853, at source 250Pa. No cycle reproduced.'))
    print(dict(independent_control_passed=passed,relative_sup_error=error),flush=True)
    if not passed:raise RuntimeError('Independent microfluidic control failed')


if __name__=='__main__':
    validate();source_control();independent_control()
