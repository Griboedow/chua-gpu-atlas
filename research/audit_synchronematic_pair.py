"""A first integral for two fixed-center phase/orientation oscillators.

The integral is local to invariant angular chambers and is not a novelty claim.
Both source-coefficient signs are checked independently against the full field.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from synchronematic import COEFFICIENTS,REPOSITORY_NUMPY_COEFFICIENTS,geometry,flow_and_jac


def primitive(delta,c):
    a=c[1,1]-c[1,3];b=4*c[1,3];cc=c[0,1]-c[0,3];d=4*c[0,3]
    if cc<=0 or cc+d<=0 or d==0:raise ValueError('Closed form requires c>0,c+d>0,d nonzero')
    turns=np.floor(delta/np.pi);u=delta-turns*np.pi
    angle=np.arctan2(np.sqrt(cc)*np.sin(u),np.sqrt(cc+d)*np.cos(u))+turns*np.pi
    return b/d*delta+(a-b*cc/d)*angle/np.sqrt(cc*(cc+d))


def invariant(state,psi,c):
    theta=state[:2];chi=state[2:]-psi
    return np.log(abs(np.sin(chi[0])/np.sin(chi[1])))+primitive(theta[1]-theta[0],c)


def gradient(state,psi,c):
    delta=state[1]-state[0];a=state[2:]-psi
    ratio=(c[1,1]+c[1,3]*(4*np.cos(delta)**2-1))/(c[0,1]+c[0,3]*(4*np.cos(delta)**2-1))
    return np.array([-ratio,ratio,1/np.tan(a[0]),-1/np.tan(a[1])])


def main():
    rng=np.random.default_rng(2026102950);records=[]
    for name,c in [('printed',COEFFICIENTS),('repository_numpy',REPOSITORY_NUMPY_COEFFICIENTS)]:
        errors=[];fd=[]
        for _ in range(256):
            psi=rng.uniform(-np.pi,np.pi);distance=rng.uniform(1,3)
            positions=np.array([[0.,0.],distance*np.array([np.cos(psi),np.sin(psi)])]);pairs=geometry(positions)
            state=np.r_[0,rng.uniform(.1,np.pi-.1),psi+rng.uniform(.15,np.pi-.15,2)]
            f,j=flow_and_jac(state,pairs,c);g=gradient(state,psi,c)
            errors.append(float(abs(g@f)));h=2e-6
            approx=np.array([(invariant(state+h*e,psi,c)-invariant(state-h*e,psi,c))/(2*h) for e in np.eye(4)])
            fd.append(float(max(abs(approx-g))))
        # Short, regular interior trajectory: avoid log singularities at chamber boundaries.
        positions=np.array([[0.,0.],[2.,0.]]);pairs=geometry(positions);start=np.array([0.,1.1,.7,1.9])
        sol=solve_ivp(lambda t,y:flow_and_jac(y,pairs,c)[0],(0,12),start,method='DOP853',
                      rtol=1e-12,atol=1e-14,max_step=.02,t_eval=np.linspace(0,12,601))
        if not sol.success:raise RuntimeError(sol.message)
        values=np.array([invariant(y,0,c) for y in sol.y.T]);drift=float(max(abs(values-values[0])))
        row=dict(variant=name,coefficients=c.tolist(),algebra_samples=256,
                 maximum_directional_identity_error=max(errors),maximum_gradient_fd_error=max(fd),
                 trajectory_duration=12,trajectory_invariant_drift=drift,
                 minimum_abs_orientation_sine=float(abs(np.sin(sol.y[2:])).min()),
                 minimum_phase_distance_from_boundary=float(np.minimum(sol.y[1]-sol.y[0],np.pi-(sol.y[1]-sol.y[0])).min()))
        if max(errors)>1e-10 or max(fd)>1e-7 or drift>1e-8:raise RuntimeError('First integral audit failed')
        records.append(row)
    result=dict(identity='I=log|sin(a)/sin(b)|+integral S_chi(delta)/S_theta(delta) ddelta',
                assumptions='Two fixed centers; symmetric pair weight; no noise; Fourier coupling as implemented; regular invariant angular chamber',
                scope='First integral reduces regular recurrent dynamics to two-dimensional level surfaces. No claim that three or four oscillators share this integral.',
                novelty='Independent derivation for screening; not claimed as a new theorem or an exhaustive prior-art finding.',records=records)
    out=Path('results/synchronematic/n2-first-integral');out.mkdir(parents=True,exist_ok=True)
    (out/'audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
