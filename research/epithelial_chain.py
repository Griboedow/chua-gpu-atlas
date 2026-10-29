"""Discrete fixed-end cell train of Parmar et al. (PRX Life 2025).

Exact normalization of discrete Eqs. 1--4, checked against arXiv:2408.02806v3.
State: N-1 internal vertex displacements u, N rest-length deficits a, N
relative myosin fluctuations c. Time s=t/sqrt(tau_a*tau_c), where the two
times are already in mechanical relaxation units. Parameters rho,ratio,alpha,beta.
No author code is used. No cell replacement, clipping or collision law is added.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.stats import qmc
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from chua_gpu import GPU


def field_jacobian(y,p):
    n=(len(y)+1)//3; m=n-1
    rho,r,alpha,beta=p
    u=np.r_[0.,y[:m],0.]; a=y[m:m+n]; c=y[m+n:]
    strain=np.diff(u);tc=np.tanh(c);te=np.tanh(strain)
    f=np.r_[rho*(np.diff(u,2)+np.diff(a)),(-a+alpha*tc)/r,r*(-c+beta*te)]
    jac=np.zeros((3*n-1,3*n-1))
    jac[:m,:m]=rho*(np.diag(np.full(m,-2.))+np.diag(np.ones(m-1),1)+np.diag(np.ones(m-1),-1))
    for j in range(m):jac[j,m+j]=-rho;jac[j,m+j+1]=rho
    for i in range(n):
        jac[m+i,m+i]=-1/r;jac[m+i,m+n+i]=alpha*(1-tc[i]**2)/r
        jac[m+n+i,m+n+i]=-r
        if i<n-1:jac[m+n+i,i]=r*beta*(1-te[i]**2)
        if i>0:jac[m+n+i,i-1]=-r*beta*(1-te[i]**2)
    return f,jac


def threshold(n,rho,ratio):
    q2=4*np.sin(np.arange(1,n)*np.pi/(2*n))**2
    b=ratio+1/ratio; z=rho[...,None]*q2
    return np.min(b[...,None]*(z+b[...,None]+1/z),axis=-1)


class ChainGPU(GPU):
    def __init__(self,n):
        super().__init__('','float64'); self.n=n;self.d=3*n-1
        source=f'#define N {n}\n#define D {self.d}\n'+Path(__file__).with_suffix('.cl').read_text()
        self.kernel=cl.Kernel(cl.Program(self.context,source).build(),'advance')
        self.info.update(model='Parmar et al. 2025 discrete fixed-end chain',cells=n,integrator='FP64 RK4, analytic tangent')

    def scan(self,p,starts,transient=250.,duration=500.,dt=.005):
        p=np.ascontiguousarray(p,float);starts=np.ascontiguousarray(starts,float);count=len(p)
        if p.shape!=(count,4) or starts.shape!=(count,self.d) or (p<=0).any() or not np.isfinite(p).all() or not np.isfinite(starts).all():raise ValueError('Invalid ensemble')
        if transient<0 or min(duration,dt)<=0:raise ValueError('Invalid times')
        state=np.zeros((count,2*self.d));state[:,:self.d]=starts
        direction=np.sin(np.arange(self.d)+1);state[:,self.d:]=direction/np.linalg.norm(direction)
        logs=np.zeros((count,4));domain=np.tile([np.inf,np.inf,np.inf,0.],(count,1));status=np.zeros(count,np.int32)
        arrays=[state,p,logs,domain,status];buffers=[self.buffer(a) for a in arrays]
        begin=time.perf_counter();gpu_seconds=0.
        for window,target in [(-1,transient)]+[(i,duration/4) for i in range(4)]:
            steps=int(np.ceil(target/dt));h=target/steps if steps else dt
            for first in range(0,steps,64):
                event=self.kernel(self.queue,(count,),None,*buffers,np.int32(count),np.int32(min(64,steps-first)),np.float64(h),np.int32(window))
                event.wait();gpu_seconds+=(event.profile.end-event.profile.start)*1e-9
            print('window',window,'elapsed',round(time.perf_counter()-begin,2),flush=True)
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,starts=starts,final_states=state[:,:self.d],lle=windows.mean(axis=1),lle_windows=windows,status=status,domain=domain,
                    settings=dict(cells=self.n,transient=transient,duration=duration,max_dt=dt,time='s=t/sqrt(tau_a*tau_c)'),
                    backend=self.info,timing=dict(wall_seconds=time.perf_counter()-begin,gpu_seconds=gpu_seconds),
                    parameter_columns=['rho=sqrt(tau_a*tau_c)','ratio=sqrt(tau_a/tau_c)','alpha','beta'],
                    domain_columns=['min_cell_length','min_rest_length','min_concentration_ratio','max_abs_state'],
                    domain_note='Minima over observation only; positive lengths/concentration necessary, not sufficient, for physical validity')


def save(result,out,extra):
    out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'scan.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    candidate=(result['status']==0)&(result['lle']>.005)&(result['lle_windows'].min(axis=1)>.002)
    clear=(result['domain'][:,:3]>0).all(axis=1)
    meta={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    meta.update(extra,trajectories=len(candidate),candidate_indices=np.flatnonzero(candidate).tolist(),domain_clear_candidates=np.flatnonzero(candidate&clear).tolist(),
                domain_clear_count=int(clear.sum()),status_counts={str(i):int((result['status']==i).sum()) for i in np.unique(result['status'])},
                criterion='Mean > .005 and all four windows > .002, in normalized time; finite-time screen only')
    (out/'summary.json').write_text(json.dumps(meta,indent=2))
    print({k:meta[k] for k in ['trajectories','candidate_indices','domain_clear_candidates','domain_clear_count','timing']},flush=True)


def validate(gpu,out):
    rng=np.random.default_rng(2026102961);d=gpu.d;n=gpu.n
    # Avoid the defective uniform chemical block at ratio=1 when comparing
    # unstructured floating-point eigenvalues; its splitting is O(sqrt(eps)).
    p=np.array([[1.,1.2,4.,3.],[4.,.25,33.4,1.],[.3,3.,5.,5.]])
    starts=rng.normal(0,.02,(3,d));errors=[];lle_errors=[];jac_errors=[];eig_errors=[]
    result=gpu.scan(p,starts,0,2.,.001)
    direction=np.sin(np.arange(d)+1);direction/=np.linalg.norm(direction)
    for pp,s,end,lle in zip(p,starts,result['final_states'],result['lle']):
        f,j=field_jacobian(s,pp);eps=1e-6
        fd=np.column_stack([(field_jacobian(s+eps*e,pp)[0]-field_jacobian(s-eps*e,pp)[0])/(2*eps) for e in np.eye(d)])
        jac_errors.append(float(np.max(abs(j-fd))))
        def rhs(t,z):
            f,j=field_jacobian(z[:d],pp);return np.r_[f,j@z[d:]]
        sol=solve_ivp(rhs,(0,2),np.r_[s,direction],method='DOP853',rtol=1e-12,atol=1e-14)
        if not sol.success:raise RuntimeError(sol.message)
        errors.append(float(max(abs(sol.y[:d,-1]-end))));lle_errors.append(float(abs(np.log(np.linalg.norm(sol.y[d:,-1]))/2-lle)))
        rho,r,alpha,beta=pp;predicted=[-1/r,-r]
        for q2 in 4*np.sin(np.arange(1,n)*np.pi/(2*n))**2:
            z=rho*q2;predicted.extend(np.roots([1,z+r+1/r,1+z*(r+1/r),z*(1+alpha*beta)]))
        observed=np.linalg.eigvals(field_jacobian(np.zeros(d),pp)[1])
        eig_errors.append(float(max(min(abs(observed-z)) for z in predicted)))
    audit=dict(jacobian_fd_errors=jac_errors,gpu_dop853_state_errors=errors,gpu_dop853_lle_errors=lle_errors,linear_modal_errors=eig_errors,
               passed=bool(max(errors)<1e-7 and max(lle_errors)<1e-7 and max(jac_errors)<1e-7 and max(eig_errors)<1e-8 and (result['status']==0).all()))
    out.mkdir(parents=True,exist_ok=True);(out/'validation.json').write_text(json.dumps(audit,indent=2));print(audit,flush=True)
    if not audit['passed']:raise RuntimeError('Validation failed')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cells',type=int,default=16);ap.add_argument('--power',type=int,default=9)
    ap.add_argument('--validate',action='store_true');ap.add_argument('--dt',type=float,default=.005)
    ap.add_argument('--transient',type=float,default=250);ap.add_argument('--duration',type=float,default=500)
    ap.add_argument('--refine',type=Path);ap.add_argument('--top',type=int,default=16);ap.add_argument('--indices',nargs='+',type=int)
    ap.add_argument('--output',type=Path,default=Path('results/epithelial-chain/n16-scout'));args=ap.parse_args()
    gpu=ChainGPU(args.cells)
    if args.validate:return validate(gpu,Path('results/epithelial-chain'))
    if args.refine:
        old=np.load(args.refine/'scan.npz');ids=np.array(args.indices,int) if args.indices else np.argsort(-old['lle'])[:args.top]
        p=old['parameters'][ids];starts=old['starts'][ids];extra=dict(parent=str(args.refine),parent_indices=ids.tolist())
    else:
        rows=qmc.Sobol(4,scramble=True,seed=2026102961).random_base2(args.power)
        rho=10**(-.5+1.5*rows[:,0]);r=10**(-.6+1.2*rows[:,1]);beta=10**(-.7+1.*rows[:,2])
        gain=threshold(args.cells,rho,r)*10**(-.05+1.05*rows[:,3]);p=np.repeat(np.column_stack([rho,r,gain/beta,beta]),2,axis=0)
        rng=np.random.default_rng(2026102961);common=rng.normal(0,.03,(2,gpu.d));common[:,:args.cells-1]=0
        starts=np.tile(common,(len(rows),1));extra=dict(seed=2026102961,sampling='Sobol rho [.316,10], ratio [.251,3.981], beta [.20,2], gain/finite-chain-threshold [.891,10], logarithmic; two common starts',source='10.1103/PRXLife.3.013002; arXiv:2408.02806v3 discrete Eqs.1--4')
    save(gpu.scan(p,starts,args.transient,args.duration,args.dt),args.output,extra)


if __name__=='__main__':main()
