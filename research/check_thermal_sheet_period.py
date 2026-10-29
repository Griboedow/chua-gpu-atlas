"""Independent shooting/Floquet control for a small-positive scout example."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from thermomechanical_sheet import equilibrium,field_jacobian


def main():
    base=Path('results/thermomechanical-sheet');meta=json.loads((base/'refined/summary.json').read_text())
    data=np.load(base/'refined/scan.npz');index=meta['parent_indices'].index(5906)
    p=data['parameters'][index];start=data['final_states'][index];eq=equilibrium(p)
    def field(t,y):return field_jacobian(y,p)[0]
    def section(t,y):return y[0]-eq[0]
    section.direction=1
    sol=solve_ivp(field,(0,5),start,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.005,events=section)
    if not sol.success or len(sol.t_events[0])<3:raise RuntimeError('No repeated section crossings')
    seed=sol.y_events[0][-1];period=float(np.median(np.diff(sol.t_events[0][-6:])))
    scale=np.array([.01,.2,.001])
    def endpoint(z,variational=False):
        s=np.array([eq[0],z[0],z[1]])
        if variational:
            def rhs(t,y):
                f,j=field_jacobian(y[:3],p)
                return np.r_[f,(j@y[3:].reshape(3,3)).ravel()]
            s=np.r_[s,np.eye(3).ravel()]
        else:rhs=field
        out=solve_ivp(rhs,(0,z[2]),s,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.005,dense_output=True)
        if not out.success:raise RuntimeError(out.message)
        return out
    def residual(z):return (endpoint(z).y[:3,-1]-[eq[0],z[0],z[1]])/scale
    fit=least_squares(residual,[seed[1],seed[2],period],bounds=([0,0,.8*period],[np.inf,1,1.2*period]),
                      xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60)
    orbit=endpoint(fit.x,True);floquet=np.linalg.eigvals(orbit.y[3:,-1].reshape(3,3))
    times=np.linspace(0,fit.x[2],1001);states=orbit.sol(times)[:3].T
    report=dict(source_scout_index=5906,parameters=p.tolist(),period=float(fit.x[2]),
                section_state=[float(eq[0]),float(fit.x[0]),float(fit.x[1])],
                scaled_closure_error=float(np.linalg.norm(fit.fun)),shooting_success=bool(fit.success),
                floquet_multipliers=[dict(real=float(v.real),imag=float(v.imag),absolute=float(abs(v))) for v in floquet],
                neutral_multiplier_error=float(np.min(abs(floquet-1))),
                sampled_min_X=float(states[:,0].min()),sampled_min_Y=float(states[:,2].min()),
                sampled_max_Y=float(states[:,2].max()),sampled_min_stiffness=float((1-p[3]*states[:,2]).min()),
                settings=dict(method='DOP853',rtol=1e-12,atol=1e-14,max_step=.005),
                interpretation='Numerical periodic-orbit control, not a validated existence or stability theorem.')
    if report['scaled_closure_error']>1e-7 or not fit.success:raise RuntimeError('Shooting failed')
    (base/'periodic-control.json').write_text(json.dumps(report,indent=2))
    np.savez_compressed(base/'periodic-control.npz',time=times,state=states,parameters=p)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
