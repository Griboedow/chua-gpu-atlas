"""Independent checks of the adaptive GPU implementation before wide surveying."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from adaptive_gpu import AdaptiveGPU
from chua_gpu import rhs
from explore_extent import starting_states
from explore_expanded import ROOT,OUT
from run import save_json
from independent_reference import reference_lle

def main():
    OUT.mkdir(parents=True,exist_ok=True);g=AdaptiveGPU();checks=[]
    p=np.array([[9,100/7,-8/7,-5/7],[255.693129732,7216.58622261,-3.477304427,-5/7],[1024,150000,-29,-5/7]])
    starts=starting_states(p)[:,:1,:]
    result=g.scan(p,initial_states=starts,transient=0,duration=.2,rtol=2e-8,atol=2e-11)
    for i,point in enumerate(p):
        reference=solve_ivp(lambda t,s:rhs(t,s,point),[0,.2],starts[i,0],method='DOP853',
            rtol=1e-12,atol=1e-13,max_step=.01/np.sqrt(point[1])).y[:,-1]
        error=float(np.max(abs(reference-result['final_states'][i,0])/np.maximum(1,abs(reference))))
        checks.append({'name':'Short state against DOP853','parameters':point.tolist(),'scaled_error':error,'passed':error<2e-5})
    linear=np.array([[2,3,-.5,-.5]])
    initial=np.array([[[.7,.1,-.3]]]);a,b,m0,m1=linear[0]
    jac=np.array([[-a*(1+m0),a,0],[1,-1,1],[0,-b,0]])
    r=g.scan(linear,initial_states=initial,transient=0,duration=4,rtol=1e-9,atol=1e-12)
    error=float(np.max(abs(r['final_states'][0,0]-expm(jac*4)@initial[0,0])))
    checks.append({'name':'Smooth linear exact matrix exponential','max_error':error,'passed':error<1e-7})
    r=g.scan(p[:1],initial_states=np.zeros((1,1,3)),transient=10,duration=40)
    checks.append({'name':'Stationary unstable origin is not classified as chaos','class':int(r['classes'][0,0]),'passed':bool(r['classes'][0,0]!=4)})
    # Revisit the old extreme with the corrected full-field CPU reference.
    old=np.load(ROOT/'outputs/extent/verified.npz');idx=int(np.argmin(abs(old['parameters'][:,0]-p[1,0])))
    point=old['parameters'][idx];seed=int(np.argmax(np.where(old['classes'][idx]==4,old['lle'][idx],-np.inf)))
    initial=old['final_states'][idx,seed]
    adaptive=g.scan(point[None,:],initial_states=initial[None,None,:],transient=40,duration=160,rtol=1e-8,atol=1e-11)
    reference=reference_lle(point,initial,transient=40,duration=160)
    value=float(adaptive['lle'][0,0]);sign=bool((adaptive['lle_windows'][0,0]>0).all() and min(reference['windows'])>0)
    checks.append({'name':'Positive LLE at the previous difficult far-field example','parameters':point.tolist(),
        'adaptive_lle':value,'adaptive_windows':adaptive['lle_windows'][0,0].tolist(),'full_field_dop853':reference,
        'passed':sign,'scope':'Agreement of positive sign at matched settings; not a precision certificate for arbitrary parameters.'})
    evidence={'device':g.info,'checks':checks,'all_checks_passed':all(x['passed'] for x in checks)}
    save_json(OUT/'adaptive_validation.json',evidence)
    print(json.dumps(evidence,indent=2),flush=True)
    if not evidence['all_checks_passed']:raise SystemExit('Adaptive validation failed')

if __name__=='__main__':main()
