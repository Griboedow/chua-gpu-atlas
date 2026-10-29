"""Locate flip crossings and follow the lower fold with pseudo-arclength.

The family changes the published interaction law. Bifurcations reported here
are numerical candidates; this is not an interval-validated existence proof.
"""
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import brentq

from continue_synchronematic_orbit import Orbit, record


def transverse(model,q,T,eta,mono,normal):
    z=null_space(normal[None,:]);end,_,_=model.integrate(q,T,eta)
    f=model.field(end,eta)[0]
    section=(np.eye(7)-np.outer(f,normal)/(normal@f))@mono
    return np.linalg.eigvals(z.T@section@z)


def main():
    path=Path('results/synchronematic/n4-periodic-branch')
    data=json.loads((path/'branch.json').read_text());model=Orbit(data['positions'])
    rows=data['records'];seed=min(rows,key=lambda r:abs(r['eta']-1.074))
    reference=np.array(seed['state']);w=np.array(seed['winding']);normal=model.field(reference,1.074)[0];normal/=np.linalg.norm(normal)
    result=dict(changes_equations=True,eta_scaling=100,flip_crossings=[],lower_fold_arc=[],validation={},
                warning='Numerical bifurcation evidence in a modified interaction law, not a proof of new physical chaos.')
    def save():(path/'bifurcations.json').write_text(json.dumps(result,indent=2))
    for bracket in [(1.07675,1.077),(1.08,1.0805)]:
        cache={}
        def target(eta):
            row=min(rows,key=lambda r:abs(r['eta']-eta));u=np.r_[row['state'],row['period']]
            u,mono,error,it=model.shoot(u,eta,reference,w,tol=1e-11)
            eig=transverse(model,u[:7],u[7],eta,mono,normal)
            closest=eig[np.argmin(abs(eig+1))]
            if abs(closest.imag)>1e-7:raise RuntimeError('Flip multiplier is not real')
            cache[eta]=(u,mono,error,it,closest)
            return closest.real+1
        eta=brentq(target,*bracket,xtol=2e-12)
        target(eta);u,mono,error,it,mult=cache[eta]
        row=record(model,u,eta,mono,error,it,w)
        h=1e-6;row['multiplier_eta_derivative']=(target(eta+h)-target(eta-h))/(2*h)
        row['crossing_bracket']=list(bracket)
        tight=Orbit(data['positions'],rtol=2e-12,max_step=.025)
        end2,M2,_=tight.integrate(u[:7],u[7],eta)
        row['tighter_closure_error']=float(np.max(abs(end2-u[:7]-2*np.pi*w)))
        row['tighter_monodromy_difference']=float(np.max(abs(M2-mono)))
        result['flip_crossings'].append(row);save()
        print('Flip crossing',eta,'T',u[7],'slope',row['multiplier_eta_derivative'],flush=True)
    # Derivatives and coordinate controls at an existing regular orbit.
    u=np.r_[reference,seed['period']];end,mono,sensitivity=model.integrate(u[:7],u[7],1.074,True)
    direction=np.arange(1,8,dtype=float);direction/=np.linalg.norm(direction);h=2e-6
    plus=model.integrate(u[:7]+h*direction,u[7],1.074)[0];minus=model.integrate(u[:7]-h*direction,u[7],1.074)[0]
    ep=model.integrate(u[:7],u[7],1.074+h)[0];em=model.integrate(u[:7],u[7],1.074-h)[0]
    independent=Orbit(data['positions'],rtol=2e-11,max_step=.025,method='Radau')
    er,mr,_=independent.integrate(u[:7],u[7],1.074)
    result['validation']=dict(coordinate_inverse_error=float(np.max(abs(model.D@model.E-np.eye(7)))),
        monodromy_directional_fd_error=float(np.max(abs((plus-minus)/(2*h)-mono@direction))),
        eta_sensitivity_fd_error=float(np.max(abs((ep-em)/(2*h)-sensitivity))),
        dop853_radau_state_difference=float(np.max(abs(er-end))),
        dop853_radau_monodromy_difference=float(np.max(abs(mr-mono))))
    save();print('Validation',result['validation'],flush=True)
    # Fixed phase section and weighted eta: continuation can turn at a fold.
    row=min(rows,key=lambda r:abs(r['eta']-1.0735));eta=row['eta']
    u,_,_,_=model.shoot(np.r_[row['state'],row['period']],eta,reference,w,tol=1e-11)
    z=np.r_[u,100*eta]
    def equations(z):
        q=z[:7];T=z[7];eta=z[8]/100
        end,M,s=model.integrate(q,T,eta,True)
        r=np.r_[end-q-2*np.pi*w,normal@(q-reference)]
        J=np.zeros((8,9));J[:7,:7]=M-np.eye(7);J[:7,7]=model.field(end,eta)[0];J[:7,8]=s/100;J[7,:7]=normal
        return r,J,M
    _,J,_=equations(z);tangent=null_space(J)[:,0]
    if tangent[-1]>0:tangent=-tangent
    ds=.008;wall=time.perf_counter()
    for index in range(80):
        pred=z+ds*tangent;trial=pred.copy()
        for iteration in range(12):
            r,J,M=equations(trial);r=np.r_[r,tangent@(trial-pred)]
            if np.linalg.norm(r,np.inf)<2e-10:break
            step=np.linalg.solve(np.vstack([J,tangent]),-r)
            trial+=step
            if np.linalg.norm(step)>.2:raise RuntimeError('Arclength corrector left local neighborhood')
        else:raise RuntimeError('Arclength corrector failed')
        z=trial;new=null_space(J)[:,0]
        if new@tangent<0:new=-new
        tangent=new
        rr=record(model,z[:8],z[8]/100,M,float(np.linalg.norm(r,np.inf)),iteration+1,w)
        rr['eta_arclength_derivative']=float(tangent[-1]/100);result['lower_fold_arc'].append(rr);save()
        if index%10==0:print('Arc',index,rr['eta'],rr['period'],rr['spectral_radius'],flush=True)
        if index>10 and z[8]/100>1.0741:break
    result['arc_wall_seconds']=time.perf_counter()-wall
    result['lower_fold_sample']=min(result['lower_fold_arc'],key=lambda r:r['eta'])
    save();print('Saved bifurcation analysis',flush=True)


if __name__=='__main__':main()
