"""Independently check two attracting cycles where the continued sheets overlap."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from trace_chua_window import Chua
from chua_layer_bifurcations import Flights

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/chua-layer-census')
    args=p.parse_args();out=args.output
    alpha=9.04;m0=-8/7+.004;beta=14.352
    model=Chua(alpha=alpha,m0=m0,rtol=2e-12,step=.02)
    first=json.loads((ROOT/'results/chua-layers/extended.json').read_text())['surfaces']
    second=json.loads((out/'second-surfaces.json').read_text())
    report=dict(alpha=alpha,beta=beta,m0=m0,m1=-5/7,cycles=[]);arrays={}
    for j,rows in enumerate([first,second]):
        row=min(rows,key=lambda r:abs(r['alpha']-alpha)+abs(r['m0']-m0))
        edge=row['symmetry'];q=np.array(edge['q'])
        regions=[r['region'] for r in edge['arcs']]
        targets=[int(round(r['end'][0])) for r in edge['arcs'][:-1]]
        family=Flights(model,regions,targets,True)
        for b in np.linspace(edge['beta'],beta,40)[1:]:q=family.solve(b,q)
        rec=family.record(beta,q,True);T=rec['period'];u=np.r_[0.,q[:2]]
        radau=Chua(alpha=alpha,m0=m0,rtol=2e-12,step=.02,method='Radau').integrate(beta,u,T)
        orbit=model.integrate(beta,u,T)
        def rhs(t,v):
            J=model.jacobian(beta,abs(v[0])<1)
            return np.r_[model.field(beta,v[:3]),(J@v[3:].reshape(3,3)).ravel()]
        var=solve_ivp(rhs,[0,T],np.r_[u,np.eye(3).ravel()],method='DOP853',rtol=2e-12,atol=2e-13,max_step=.02)
        f=model.field(beta,var.y[:3,-1]);M=var.y[3:,-1].reshape(3,3)
        P=((np.eye(3)-np.outer(f,[1.,0.,0.])/f[0])@M)[1:,1:]
        rho=float(max(abs(np.linalg.eigvals(P))))
        trials=[]
        for delta in [[.001,-.0005,.0007],[-.001,.0005,-.0007]]:
            start=u+delta;sol=model.integrate(beta,start,500.,dense=False)
            error=float(np.max(abs(sol.y_events[2][-1]-u)))
            trials.append(dict(initial=start.tolist(),duration=500.,last_section_error=error))
        item=dict(family=j+1,cycle=rec,maxima=len(orbit.t_events[3]),radau_closure=float(max(abs(radau.y[:,-1]-u))),
                  independent_rho=rho,variation_success=var.success,variation_closure=float(max(abs(var.y[:3,-1]-u))),
                  perturbed_starts=trials)
        report['cycles'].append(item)
        t=np.linspace(0,T,4001);arrays[f'family{j+1}_t']=t;arrays[f'family{j+1}_xyz']=orbit.sol(t).T
        print('coexistence',j+1,T,rho,item['radau_closure'],trials,flush=True)
    report['verified']=all(r['cycle']['region_violation']<1e-7 and r['radau_closure']<2e-7 and
        r['variation_success'] and r['independent_rho']<.9999 and
        all(t['last_section_error']<2e-6 for t in r['perturbed_starts']) for r in report['cycles'])
    (out/'coexistence.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    np.savez_compressed(out/'coexistence.npz',**arrays)


if __name__=='__main__':main()
