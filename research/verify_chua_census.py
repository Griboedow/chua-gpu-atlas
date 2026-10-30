"""Close and check selected cycles from the layer-count transect independently."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from trace_chua_window import Chua
from chua_layer_bifurcations import from_orbit


def one(model,beta,start):
    warm=model.integrate(beta,start,1600.,dense=False)
    late=model.integrate(beta,warm.y[:,-1],500.)
    times,states=late.t_events[2],late.y_events[2]
    rec=dict(beta=float(beta),status='no_recurrence',positive_returns=len(times))
    for lag in range(1,min(60,len(times)//3)):
        errors=np.max(abs(states[-lag:]-states[-2*lag:-lag]))
        if errors>2e-6:continue
        guess=np.r_[states[-1,1:],times[-1]-times[-1-lag]]
        family,q=from_orbit(model,beta,guess,False)
        r=family.record(beta,q,True)
        T=r['period'];u=np.r_[0.,q[:2]]
        orbit=model.integrate(beta,u,T)
        peak_count=len(orbit.t_events[3])
        sym_error=float(np.max(abs(orbit.sol(np.linspace(0,T/2,1001))+orbit.sol(np.linspace(T/2,T,1001)))))
        other=Chua(rtol=2e-12,step=.02,method='Radau').integrate(beta,u,T)
        def rhs(t,v):
            J=model.jacobian(beta,abs(v[0])<1)
            return np.r_[model.field(beta,v[:3]),(J@v[3:].reshape(3,3)).ravel()]
        var=solve_ivp(rhs,[0,T],np.r_[u,np.eye(3).ravel()],method='DOP853',rtol=2e-12,atol=2e-13,max_step=.02)
        f=model.field(beta,var.y[:3,-1]);M=var.y[3:,-1].reshape(3,3)
        P=((np.eye(3)-np.outer(f,[1.,0.,0.])/f[0])@M)[1:,1:]
        multipliers=np.linalg.eigvals(P);rho=float(max(abs(multipliers)))
        closure=float(np.max(abs(other.y[:,-1]-u)))
        rec.update(status='verified_attracting_cycle' if var.success and closure<2e-7 and rho<.9999 and r['region_violation']<1e-7 else 'unresolved',
                   cycle=r,regions=family.regions,targets=family.targets,peak_count=peak_count,
                   section_period=lag,half_symmetry_error=sym_error,radau_closure=closure,
                   independent_variational_rho=rho,independent_variational_success=var.success,
                   independent_variational_closure=float(np.max(abs(var.y[:3,-1]-u))),
                   variation_vs_affine_rho=abs(rho-r['spectral_radius']))
        return rec,orbit
    return rec,None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results/chua-layer-census')
    p.add_argument('--betas',type=float,nargs='*',default=[14.259875,14.2908125,14.3036875,14.30525,14.30725,14.32,14.34990625,14.354])
    p.add_argument('--label',default='cycles')
    p.add_argument('--source',default='line')
    p.add_argument('--refine-existing',action='store_true')
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    if args.refine_existing:
        path=out/(args.label+'.json');report=json.loads(path.read_text())
        for r in report['records']:
            if r['status']!='unresolved' or 'cycle' not in r:continue
            model=Chua(rtol=2e-13,step=.01,method='Radau')
            u=np.r_[0.,r['cycle']['q'][:2]];T=r['cycle']['period']
            other=model.integrate(r['beta'],u,T)
            closure=float(np.max(abs(other.y[:,-1]-u)))
            r['radau_refinement']=dict(rtol=2e-13,atol=2e-14,max_step=.01,closure=closure,
                                      original_closure=r['radau_closure'])
            if r['independent_variational_success'] and closure<2e-7 and r['independent_variational_rho']<.9999 and r['cycle']['region_violation']<1e-7:
                r['status']='verified_attracting_cycle'
            print('refinement',r['beta'],closure,r['status'],flush=True)
        report['refinement_code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
        return
    scout=np.load(out/(args.source+'.npz'));model=Chua(rtol=2e-12,step=.025)
    records=[];arrays={}
    for j,beta in enumerate(args.betas):
        i=int(np.argmin(abs(scout['beta']-beta)));start=scout['final_states'][i,0]
        try:r,orbit=one(model,beta,start)
        except (RuntimeError,ValueError,OverflowError) as e:r=dict(beta=beta,status='failed',error=str(e));orbit=None
        records.append(r)
        if orbit is not None:
            t=np.linspace(0,r['cycle']['period'],4001);arrays[f'case{j}_t']=t;arrays[f'case{j}_xyz']=orbit.sol(t).T
        print(json.dumps({k:v for k,v in r.items() if k not in ['cycle','regions','targets']}),flush=True)
        (out/(args.label+'.json')).write_text(json.dumps(dict(records=records,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n',encoding='utf-8',newline='\n')
        np.savez_compressed(out/(args.label+'.npz'),**arrays)


if __name__=='__main__':main()
