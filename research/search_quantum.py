"""Reproducible GPU exploration; run from the repository root."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from scipy.integrate import solve_ivp,quad
from quantum_neuron import QuantumGPU,conductance,rhs,PARAMETERS

def save(path,result):
    path.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path/'scan.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    summary={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    valid=result['status']==0
    candidate=valid&(result['lle_windows'].min(axis=1)>.002)&(result['lle']>.005)&((result['stats'][:,7]-result['stats'][:,6])>1e-5)
    summary.update(parameter_columns=PARAMETERS,trajectories=len(valid),valid=int(valid.sum()),candidate_indices=np.flatnonzero(candidate).tolist(),
                   max_lle=float(np.nanmax(result['lle'])),status_counts={str(k):int((result['status']==k).sum()) for k in np.unique(result['status'])},
                   criterion='mean LLE > 0.005; all four windows > 0.002; late V range > 1e-5; physical-domain status 0',
                   interpretation='Finite-time screen, not proof of chaos, absence, or publication novelty.')
    (path/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)

def controls(gpu):
    p=np.array([2.,5.,.1,7.,1.,1.,.8,0.]);s=np.array([.05,0.,.95,.1])
    # Derivative cross-check against finite differences of the state-only RHS.
    tangent=np.array([1.,2.,3.,4.]);eps=1e-6
    analytic=rhs(0,np.r_[s,tangent],p)[4:]
    numerical=(rhs(0,s+eps*tangent,p)-rhs(0,s-eps*tangent,p))/(2*eps)
    jac_error=float(np.max(abs(analytic-numerical)))
    grid=np.linspace(-3,3,61);f,_=conductance(grid,512)
    reference=[]
    for x in grid:
        reference.append([quad(lambda u: np.exp(-u*u)/np.sqrt(np.pi)*factor(u)*2*np.exp(-abs((.5*u+x)/.13))/(1+np.exp(-2*abs((.5*u+x)/.13))),-12,12,epsabs=2e-13,epsrel=2e-13)[0]
                          for factor in [lambda u:1,lambda u:np.sqrt(2)*u,lambda u:2*u*u]])
    quadrature_error=float(np.max(abs(f-reference)))
    testx=np.random.default_rng(45).uniform(-5,5,1000)
    ff,df=conductance(testx,1024)
    table_error=float(np.max(abs(gpu.spline(testx)-ff)))
    derivative_error=float(np.max(abs(gpu.spline(testx,1)-df)))
    pp=np.tile(p,(3,1));pp[1,:4]=[5,50,.03,25];pp[2,5]=0
    starts=np.tile(s,(3,1))
    short=gpu.scan(pp,starts,transient=0,duration=5,rtol=1e-9,atol=1e-12)
    state_errors=[]
    for i in range(3):
        sol=solve_ivp(lambda t,y:rhs(t,y,pp[i],1024),(0,5),starts[i],method='DOP853',rtol=2e-12,atol=2e-14)
        if not sol.success:raise RuntimeError(sol.message)
        state_errors.append(float(np.max(abs(short['final_states'][i]-sol.y[:,-1]))))
    result=dict(jacobian_error=jac_error,quadrature_vs_quad=quadrature_error,table_value_error=table_error,
                table_derivative_error=derivative_error,gpu_vs_direct_cpu_state_errors=state_errors,
                short_status=short['status'].tolist(),backend=gpu.info)
    result['passed']=bool(jac_error<1e-7 and quadrature_error<1e-9 and table_error<1e-9 and derivative_error<1e-6 and max(state_errors)<1e-6 and (short['status']==0).all())
    path=Path('results/candidate-research');path.mkdir(parents=True,exist_ok=True)
    (path/'validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2),flush=True)
    if not result['passed']:raise RuntimeError('Numerical controls failed')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--power',type=int,default=12)
    parser.add_argument('--validate',action='store_true');parser.add_argument('--refine',type=Path)
    parser.add_argument('--output',type=Path,default=Path('results/candidate-research/scout'))
    parser.add_argument('--transient',type=float,default=200);parser.add_argument('--duration',type=float,default=400)
    parser.add_argument('--rtol',type=float,default=1e-7);parser.add_argument('--x0-variable',action='store_true')
    args=parser.parse_args();gpu=QuantumGPU()
    if args.validate:controls(gpu);raise SystemExit()
    if args.refine:
        old=np.load(args.refine/'scan.npz');meta=json.loads((args.refine/'summary.json').read_text())
        indices=np.array(meta['candidate_indices'],int);p=old['parameters'][indices];starts=old['starts'][indices]
        if not len(p):raise SystemExit('No candidates to refine')
        args.output.mkdir(parents=True,exist_ok=True);np.save(args.output/'parent_indices.npy',indices)
    else:
        u=qmc.Sobol(7,scramble=True,seed=20261029).random_base2(args.power)
        p=np.zeros((len(u),8))
        for i,(lo,hi) in enumerate([(.05,12),(1,200),(.01,1),(.3,30),(.1,2)]):p[:,i]=lo*(hi/lo)**u[:,i]
        p[:,5]=.3+.7*u[:,5];p[:,6]=-.9+1.8*u[:,6] if args.x0_variable else .8
        # Three physical states: thermal/low voltage, coherent, thermal/high voltage.
        p=np.repeat(p,3,axis=0);starts=np.zeros((len(p),4))
        starts[0::3,2]=p[0::3,5];starts[0::3,3]=.01*p[0::3,0]
        starts[1::3,:3]=[.6,0,.6];starts[1::3,3]=.35*p[1::3,0]
        starts[2::3,2]=p[2::3,5];starts[2::3,3]=.95*p[2::3,0]
    print(f'Scanning {len(p)} trajectories on {gpu.info["device"]}',flush=True)
    result=gpu.scan(p,starts,transient=args.transient,duration=args.duration,rtol=args.rtol,atol=args.rtol/1000,
                    progress=lambda fraction,seconds:print(f'{fraction:.1%}, {seconds:.1f}s',flush=True))
    save(args.output,result)
