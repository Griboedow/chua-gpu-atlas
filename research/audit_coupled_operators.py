"""Structural identities and reproducible local parameter/start perturbations."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit
from coupled_operators import OperatorsGPU, flow, initial_modes, save_result, next_event, segment_matrix


def symmetry(y):
    return np.array([y[3],y[4],y[5],1-y[0],-y[1],1-y[2]])


def symmetry_derivative(f):
    return np.r_[f[3:],-f[:3]]


def smooth_equilibrium(p,width):
    k,b,dm,dp=p
    xm=brentq(lambda x:expit((expit((dm-x)/width)-dp)/width)-x,0.,1.,xtol=5e-15)
    xp=expit((dm-xm)/width)
    am=xm*(1-xm)/width;ap=xp*(1-xp)/width
    block=np.array([[0,1,0],[-k,-b,k],[0,0,-1.]])
    jac=np.zeros((6,6));jac[:3,:3]=block;jac[3:,3:]=block;jac[2,3]=am;jac[5,0]=-ap
    eig=np.linalg.eigvals(jac)
    polynomial=np.array([1,b+1,k+b,k]);char=np.polymul(polynomial,polynomial);char[-1]+=k*k*am*ap
    residual=max(abs(np.polyval(char,eig)))/(np.linalg.norm(char)*max(1.,max(abs(eig))**6))
    return dict(equilibrium=[xm,0.,xm,xp,0.,xp],gain_product=am*ap,
                eigenvalues=[[float(z.real),float(z.imag)] for z in eig],characteristic_scaled_residual=float(residual))


def reference_period():
    source=np.load('results/light-operators/coupled-scout.npz');p=source['parameters'][0];y=source['final_states'][0].copy();modes=initial_modes(y,p)
    for _ in range(100):
        y,modes,dt,salt,which=next_event(y,p,modes,256)
        if which==0 and y[1]>0:break
    initial=y.copy();matrix=np.eye(6);elapsed=0.;events=0
    for _ in range(100):
        y,modes,dt,salt,which=next_event(y,p,modes,256)
        matrix=salt@segment_matrix(p,dt)@matrix;elapsed+=dt;events+=1
        if which==0 and y[1]>0:break
    f=flow(y,p,modes);normal=np.array([1.,0,0,0,0,0])
    derivative=((np.eye(6)-np.outer(f,normal)/f[0])@matrix)[1:,1:]
    eig=np.linalg.eigvals(derivative)
    return dict(parameters=p.tolist(),section_state=initial.tolist(),return_state=y.tolist(),
                period=elapsed,event_count=events,closure_error=float(max(abs(y-initial))),
                poincare_eigenvalues=[[float(z.real),float(z.imag)] for z in eig],
                spectral_radius=float(max(abs(eig))),caveat='Numerical stable periodic control of the continuous-heat interpretation; very small multipliers have limited relative precision')


def audit():
    rng=np.random.default_rng(2026102965);errors=[];four=[];jerk=[]
    for _ in range(256):
        p=np.array([10**rng.uniform(1,4),rng.uniform(.2,4),*rng.uniform(.01,.99,2)])
        y=rng.uniform(-.2,1.2,6);p2=p.copy();p2[2:]=[p[3],1-p[2]]
        f=flow(y,p,initial_modes(y,p));g=flow(symmetry(y),p2,initial_modes(symmetry(y),p2))
        errors.append(float(np.max(abs(g-symmetry_derivative(f)))))
        z=y.copy()
        for _ in range(4):z=symmetry(z)
        four.append(float(np.max(abs(z-y))))
        modes=initial_modes(y,p);k,b=p[:2]
        for j in (0,3):
            third=k*(f[j+2]-f[j])-b*f[j+1]
            jerk.append(float(abs(third+(b+1)*f[j+1]+(k+b)*f[j]+k*y[j]-k*modes[j//3])))
    source=np.load('results/light-operators/coupled-scout.npz');p=source['parameters'][1026];k,b=p[:2]
    critical_gain=b*b*(k+b+1)**2/(k*(b+1)**3)
    width=brentq(lambda w:smooth_equilibrium(p,w)['gain_product']-critical_gain,.005,.5,xtol=1e-14)
    smooth=dict(reference_index=1026,parameters=p.tolist(),width_005=smooth_equilibrium(p,.005),
                critical_gain_product=critical_gain,linear_Hopf_candidate_width=width,
                angular_frequency=float(np.sqrt(k/(b+1))),
                characteristic='[(lambda+1)(lambda^2+b*lambda+k)]^2+k^2*a_minus*a_plus',
                caveat='A simple imaginary pair is a linear Hopf candidate, not a nonlinear bifurcation proof or chaos boundary; logistic optical law is an explicit model modification')
    result=dict(symmetry='T(xm,vm,ym,xp,vp,yp)=(xp,vp,yp,1-xm,-vm,-ym+1); (dm,dp)->(dp,1-dm)',
                symmetry_max_error=max(errors),fourth_power_error=max(four),jerk_identity_error=max(jerk),
                exact_symmetry_at_fixed_parameters='C4 when dm=dp=1/2; otherwise maps parameter-related systems',
                jerk_pair='P(D)xm=k H(xp-dp); P(D)xp=k H(dm-xm); P(D)=(D+1)(D^2+bD+k)',
                divergence='-2*(b+1); each transverse saltation determinant is 1',
                boundedness='Thermal fractions in [0,1]; stable forced mechanical filters for k,b>0',
                no_branch_equilibrium='For 0<dm,dp<1, branch equilibria require xm=xp and xp=1-xm with xm,xp in {0,1}, a contradiction. This does not exclude a Filippov equilibrium on the intersection of switch surfaces.',
                smooth_equilibrium_analysis=smooth,
                source_parameter_periodic_control=reference_period(),
                passed=max(errors)<1e-10 and max(four)<1e-14 and max(jerk)<1e-9)
    Path('results/light-operators/coupled-structure.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
    if not result['passed']:raise RuntimeError('Structural identity failed')


def robustness():
    source=np.load('results/light-operators/coupled-scout.npz');rng=np.random.default_rng(2026102966)
    params=[];starts=[];reference=[];jitters=[]
    for index in [1026,3054]:
        for amount in [0.,.001,.01,.03]:
            p=np.tile(source['parameters'][index],(32,1))*np.exp(rng.uniform(-amount,amount,(32,4)))
            y=rng.uniform(.05,.95,(32,6));y[:,[1,4]]=0.
            params.extend(p);starts.extend(y);reference.extend([index]*32);jitters.extend([amount]*32)
    data=OperatorsGPU().scan(np.asarray(params),np.asarray(starts),transient=4000,events=16000,subdivisions=256)
    data['reference_index']=np.array(reference);data['log_parameter_perturbation']=np.array(jitters)
    save_result('results/light-operators/coupled-robustness.npz',data)
    summarize_robustness()


def summarize_robustness():
    data=np.load('results/light-operators/coupled-robustness.npz')
    good=(data['status']==0)&(data['lle']>.005)&np.all(data['lle_windows']>.002,axis=1)
    rows=[]
    for i in [1026,3054]:
        for amount in [0.,.001,.01,.03]:
            mask=(data['reference_index']==i)&(data['log_parameter_perturbation']==amount)
            rows.append(dict(reference=i,log_parameter_jitter=amount,retained=int(good[mask].sum()),runs=int(mask.sum())))
    out=dict(method='Independent random starts; all four parameters multiplied by exp(uniform(-j,j)); fractions are not basin probabilities',groups=rows)
    Path('results/light-operators/coupled-robustness-summary.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--robustness',action='store_true');p.add_argument('--summarize',action='store_true');args=p.parse_args()
    if args.robustness:robustness()
    elif args.summarize:summarize_robustness()
    else:audit()
