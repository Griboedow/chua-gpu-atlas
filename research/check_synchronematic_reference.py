"""Independent adaptive DOP853/QR spectrum on the global-phase quotient."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space,block_diag
from synchronematic import geometry,flow_and_jac,COEFFICIENTS

def spectrum(positions,start,transient=300,duration=1000,rtol=1e-10,coefficients=None):
    n=len(positions);d=2*n;pairs=geometry(positions);basis=block_diag(null_space(np.ones((1,n))),np.eye(n));dim=d-1
    state=np.r_[start,basis.ravel()];sums=np.zeros((4,dim));divergence=np.zeros(4);qerrors=[];trace_t=[];trace_s=[];history=[]
    def fun(t,y):
        f,j=flow_and_jac(y[:d],pairs) if coefficients is None else flow_and_jac(y[:d],pairs,coefficients)
        return np.r_[f,(j@y[d:d+d*dim].reshape(d,dim)).ravel(),np.trace(j)]
    state=np.r_[state,0.];wall=time.perf_counter();last=wall
    for k in range(int(transient+duration)):
        # Short QR intervals prevent loss of the strongly contracting directions.
        for sub in range(4):
            t=k+sub*.25;state[-1]=0
            sol=solve_ivp(fun,(t,t+.25),state,method='DOP853',rtol=rtol,atol=rtol*.01,max_step=.1)
            if not sol.success:raise RuntimeError(sol.message)
            state=sol.y[:,-1];matrix=state[d:d+d*dim].reshape(d,dim)
            matrix[:n]-=matrix[:n].mean(axis=0)
            q,r=np.linalg.qr(matrix);logs=np.log(np.abs(np.diag(r)))
            state[d:d+d*dim]=q.ravel();state[:d]=np.angle(np.exp(1j*state[:d]))
            if k>=transient:
                w=min(3,int(4*(k-transient)/duration));sums[w]+=logs;divergence[w]+=state[-1]
                trace_t.append(t+.25-transient);trace_s.append(state[:d].copy())
        if k>=transient and (k-transient)%10==9:history.append([k+1-transient,*list(sums.sum(axis=0)/(k+1-transient))])
        if time.perf_counter()-last>20:
            print(f'CPU spectrum {(k+1)/(transient+duration):.1%}, {time.perf_counter()-wall:.1f}s',flush=True);last=time.perf_counter()
    windows=sums/(duration/4);spectrum=windows.mean(axis=0)
    result={'positions':positions.tolist(),'start':start.tolist(),'transient':transient,'duration':duration,
            'rtol':rtol,'atol':rtol*.01,'method':'DOP853, analytic full Jacobian, thin QR every 0.25, common phase removed',
            'spectrum':spectrum.tolist(),'spectrum_windows':windows.tolist(),
            'mean_divergence':float(divergence.sum()/duration),'spectrum_sum':float(spectrum.sum()),
            'divergence_identity_error':float(abs(spectrum.sum()-divergence.sum()/duration)),
            'coefficients':(COEFFICIENTS if coefficients is None else np.asarray(coefficients)).tolist(),
            'wall_seconds':time.perf_counter()-wall,'model_limit':'Specified Fourier approximation; no noise or position dynamics; not a hardware prediction.'}
    return result,np.array(trace_t),np.array(trace_s),np.array(history)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=Path('results/synchronematic/n4-scout'))
    parser.add_argument('--index',type=int,default=43);parser.add_argument('--transient',type=int,default=300)
    parser.add_argument('--duration',type=int,default=1000);parser.add_argument('--rtol',type=float,default=1e-10)
    parser.add_argument('--output',type=Path,default=Path('results/synchronematic/n4-cpu'))
    args=parser.parse_args();a=np.load(args.source/'scan.npz');i=args.index
    coefficients=a['coefficients'][i] if 'coefficients' in a else None
    result,t,s,h=spectrum(a['positions'][i],a['starts'][i],args.transient,args.duration,args.rtol,coefficients)
    args.output.mkdir(exist_ok=True,parents=True)
    result['parent']={'source':str(args.source),'index':i}
    (args.output/'spectrum.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(args.output/'trace.npz',time=t,state=s,history=h,positions=a['positions'][i])
    print(json.dumps(result,indent=2),flush=True)
