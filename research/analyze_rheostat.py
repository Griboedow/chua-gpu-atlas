"""Independent periodic controls and exact linear-stability audit."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from rheostat import field_jacobian, critical_gain


def crossing(section, p, method='DOP853'):
    basis=np.array([[0.,0.],[1.,0.],[0.,1.]])
    initial=np.r_[0.,section,basis.ravel()]
    def rhs(t,z):
        f,j=field_jacobian(z[:3],p)
        return np.r_[f,(j@z[3:].reshape(3,2)).ravel()]
    def event(t,z):return z[0]
    event.direction=1;event.terminal=True
    kwargs=dict(method=method,rtol=2e-11,atol=2e-13,max_step=.03)
    kick=solve_ivp(rhs,(0,1e-5),initial,**kwargs)
    sol=solve_ivp(rhs,(1e-5,200),kick.y[:,-1],events=event,**kwargs)
    if not sol.success or not len(sol.t_events[0]):raise RuntimeError('No upward return')
    t=sol.t_events[0][0];z=sol.y_events[0][0];f,_=field_jacobian(z[:3],p)
    if f[0]<=1e-9:raise RuntimeError('Degenerate return')
    matrix=z[3:].reshape(3,2)
    matrix-=np.outer(f,matrix[0]/f[0])
    return z[1:3],t,matrix[1:]


def close_orbit(p, guess, label):
    def residual(s):
        end,_,jac=crossing(s,p)
        return end-s,jac-np.eye(2)
    solved=root(residual,guess,jac=True,tol=1e-10)
    s=solved.x;end,T,jac=crossing(s,p)
    rad,Tr,jr=crossing(s,p,'Radau')
    values=np.linalg.eigvals(jac)
    sol=solve_ivp(lambda t,y:field_jacobian(y,p)[0],(0,T),np.r_[0,s],
                  method='DOP853',rtol=1e-12,atol=1e-14,max_step=T/1000,dense_output=True)
    times=np.linspace(0,T,2001);trace=sol.sol(times).T
    assert max(abs(end-s))<1e-8
    assert max(abs(rad-s))<1e-7
    assert s[0]>1e-6 and max(abs(values))<1
    return dict(label=label,parameters=list(p),section=s.tolist(),period=float(T),
                transverse_multipliers=[[float(x.real),float(x.imag)] for x in values],
                closure_error=float(max(abs(end-s))),radau_closure_error=float(max(abs(rad-s))),
                radau_period_error=float(abs(Tr-T)),radau_jacobian_error=float(np.max(abs(jr-jac))),
                amplitude_X=float(np.ptp(trace[:,0])/2),amplitude_V=float(np.ptp(trace[:,1])/2),
                trace_file=label+'.npz'), times,trace


def main():
    out=Path('results/rheostat/controls');out.mkdir(parents=True,exist_ok=True)
    results=[];c,r=.2,1.
    for excess in [.1,.03,.01]:
        g=critical_gain(c,r)+excess
        effective=(-(r+1/r)+np.sqrt((r+1/r)**2+4*g))/2
        velocity=3*np.pi/8*(effective-c)
        result,t,y=close_orbit([c,r,g],[velocity,effective*velocity],f'near-threshold-{excess:g}')
        result['excess_gain']=excess
        result['leading_amplitude_X']=3*np.pi/8*excess/(np.sqrt(1+c*r)*(1/r+2*c+r))
        result['amplitude_to_leading_ratio']=result['amplitude_X']/result['leading_amplitude_X']
        np.savez_compressed(out/result['trace_file'],time=t,state=y)
        results.append(result);print(result,flush=True)
    data=np.load('results/rheostat/scout/scan.npz');index=int(np.argmax(data['lle']))
    p=data['parameters'][index];start=data['final_states'][index]
    def event(t,y):return y[0]
    event.direction=1
    sol=solve_ivp(lambda t,y:field_jacobian(y,p)[0],(0,100),start,method='DOP853',rtol=1e-11,atol=1e-13,events=event,max_step=.02)
    assert sol.success and len(sol.t_events[0])>2
    result,t,y=close_orbit(p,sol.y_events[0][-1,1:],'strongest-scout-control')
    result['scout_index']=index;result['scout_lle']=float(data['lle'][index])
    np.savez_compressed(out/result['trace_file'],time=t,state=y)
    results.append(result);print(result,flush=True)
    # Exact spectral crossing, distinct from a nonlinear Hopf proof.
    rr=np.geomspace(.01,100,400);g=1.
    eigenvalues=np.array([np.linalg.eigvals(field_jacobian([0,0,0],[c,x,g])[1]) for x in rr])
    roots=np.sort(np.roots([c,c*c-g,c]))
    algebra=[]
    for rv in [roots[0],1.,roots[1]]:
        j=field_jacobian([0,0,0],[c,rv,g])[1]
        algebra.append(float(np.max(abs(np.poly(j)-[1,c+rv,1+c*rv,rv*(1+g)]))))
    np.savez_compressed(out/'linear-boundary.npz',rate=rr,max_real=eigenvalues.real.max(axis=1),gain_critical=critical_gain(c,rr))
    summary=dict(orbits=results,linear_example=dict(c=c,g=g,rate_boundaries=roots.tolist(),polynomial_error=max(algebra)),
        amplitude_derivation='Leading harmonic balance: A_X ~ (3*pi/8)*(g-g_H)/(sqrt(1+c*r)*(1/r+2*c+r)); not an interval-validated bifurcation theorem.',
        interpretation='Stable numerically closed periodic controls; tiny contracting multipliers may have limited relative precision.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2))


if __name__=='__main__':main()
