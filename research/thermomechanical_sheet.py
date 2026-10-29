"""Independent algebra audit of Shen's 2025 sheet model; no chaos claim."""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


def field_jacobian(state, parameters):
    """Time s=t/tau, state (X,V,Y), parameters (K,G,B,S,D)."""
    X,V,Y=state;K,G,B,S,D=parameters
    f=np.array([V,G*Y-K*(1-S*Y)*X-B*V-D*V**3,1-Y-X])
    jac=np.array([[0.,1.,0.],[-K*(1-S*Y),-B-3*D*V**2,G+K*S*X],[-1.,0.,-1.]])
    return f,jac


def equilibrium(parameters):
    K,G,B,S,D=parameters
    Y=brentq(lambda y:K*(1-S*y)*(1-y)-G*y,0.,min(1.,1/S) if S>0 else 1.)
    return np.array([1-Y,0.,Y])


def dimensionless(m,c,k0,alpha,gamma,mu,tau,beta,delta):
    return np.array([k0*tau*tau/m,gamma*mu*tau*tau/m,c*tau/m,-alpha*delta/k0,
                     beta*delta*delta/(m*tau*mu*mu)])


def jerk_field(state, parameters):
    """Same flow in (Y,Y',Y''); exact affine change, not a new model."""
    y,dy,ddy=state;K,G,B,S,D=parameters
    dddy=-(1+B)*ddy-(B+K*(1-S*y))*dy-G*y+K*(1-S*y)*(1-y)-D*(dy+ddy)**3
    return np.array([dy,ddy,dddy])


def main():
    rng=np.random.default_rng(2026102942);coordinate_errors=[];jacobian_errors=[];jerk_errors=[]
    for _ in range(32):
        m,c,k0,gamma,mu,tau,beta,delta=10**rng.uniform(-.5,.5,8)
        alpha=-10**rng.uniform(-1.,0.)
        p=dimensionless(m,c,k0,alpha,gamma,mu,tau,beta,delta)
        A,velocity,theta=rng.uniform(.1,.5,3)
        scale=np.array([mu/delta,tau*mu/delta,1/delta])
        raw=np.array([velocity,(gamma*theta-(k0+alpha*theta)*A-c*velocity-beta*velocity**3)/m,
                      (delta-theta-mu*A)/tau])
        state=np.array([A,velocity,theta])*scale
        f,j=field_jacobian(state,p)
        coordinate_errors.append(float(np.max(abs(f-raw*scale*tau))))
        eps=1e-6
        finite=np.column_stack([(field_jacobian(state+eps*e,p)[0]-field_jacobian(state-eps*e,p)[0])/(2*eps) for e in np.eye(3)])
        jacobian_errors.append(float(np.max(abs(j-finite))))
        # Y'=1-Y-X, Y''=-Y'-V, Y'''=-Y''-V'.
        dy=1-state[2]-state[0];ddy=-dy-state[1]
        jerk_errors.append(float(np.max(abs(jerk_field([state[2],dy,ddy],p)-[dy,ddy,-ddy-f[1]]))))
    # Dimensionless illustrative values; not claimed to reproduce a paper figure.
    m=c=k0=gamma=mu=tau=delta=1.;alpha=-.1;beta=.1
    p=dimensionless(m,c,k0,alpha,gamma,mu,tau,beta,delta)
    eq=equilibrium(p);A=eq[0]*delta/mu;theta=eq[2]*delta
    printed=gamma*delta/(k0+alpha*delta+gamma*mu)
    printed_residual=(k0+alpha*(delta-mu*printed))*printed-gamma*(delta-mu*printed)
    f,j=field_jacobian(eq,p);K,G,B,S,D=p
    effective=K*(1-S*eq[2]);coupling=G+K*S*eq[0]
    polynomial=np.array([1.,1+B,effective+B,effective+coupling])
    error=max(coordinate_errors);jac_error=max(jacobian_errors)
    if error>1e-9 or jac_error>1e-5 or max(jerk_errors)>1e-9 or np.max(abs(np.poly(j)-polynomial))>1e-10:
        raise RuntimeError('Coordinate, Jacobian, or characteristic polynomial audit failed')
    result=dict(source_doi='10.1088/0256-307X/42/9/090801',
                model='Printed Eqs 1-2 with the cubic velocity damping explicitly added in the numerical section',
                coordinate_validation_max_error=error,jacobian_fd_max_error=jac_error,validation_samples=32,
                jerk_coordinate_validation_max_error=max(jerk_errors),
                example_parameters=p.tolist(),example_is_paper_calibration=False,
                exact_equilibrium=eq.tolist(),exact_equilibrium_residual=float(np.linalg.norm(f)),
                printed_equilibrium_displacement=printed,printed_force_balance_residual=printed_residual,
                corrected_displacement=A,corrected_temperature_rise=theta,
                characteristic_polynomial=polynomial.tolist(),
                stability_margin=float(B*(effective+B+1)-coupling),
                caveat='Algebra audit only. No trajectory survey, physical calibration, or chaos conclusion.')
    target=Path('results/thermomechanical-sheet');target.mkdir(exist_ok=True)
    (target/'algebra-audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
