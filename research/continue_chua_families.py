"""Find local symmetric stability windows and map the second broad sheet."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import root

from trace_chua_window import Chua
from chua_layer_bifurcations import Flights,from_orbit


def boundary(f,beta,q,target):
    def fun(w):
        r,J,P,aux=f.evaluate(w[-1],w[:-1])
        return np.r_[r,target*target-target*np.trace(P)+f.sign*np.exp(aux['logdet'])]
    w=root(fun,np.r_[q,beta],tol=1e-10).x
    # The characteristic equation is more sensitive than state closure.
    # Polish using a larger central-difference stencil, retaining the analytic
    # closing Jacobian. Do not relax the acceptance threshold after a failure.
    for iteration in range(8):
        residual=fun(w)
        if max(abs(residual))<3e-10:break
        J=np.empty((len(w),len(w)))
        J[:-1,:-1]=f.evaluate(w[-1],w[:-1])[1]
        for k in range(len(w)):
            h=1e-6*max(1,abs(w[k]));e=np.zeros(len(w));e[k]=h
            derivative=(fun(w+e)-fun(w-e))/(2*h)
            J[-1,k]=derivative[-1]
            if k==len(w)-1:J[:-1,k]=derivative[:-1]
        change=np.linalg.solve(J,-residual)
        for scale in [1.,.5,.25,.125]:
            trial=w+scale*change
            if max(abs(fun(trial)))<max(abs(residual)):
                w=trial;break
        else:break
    error=float(max(abs(fun(w))))
    if error>1e-8:raise RuntimeError(f'Boundary residual {error}')
    return w[-1],w[:-1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results/chua-layer-census')
    p.add_argument('--surfaces-only',action='store_true')
    p.add_argument('--boundaries-only',action='store_true')
    p.add_argument('--families-only',action='store_true')
    args=p.parse_args();out=args.output
    model=Chua(rtol=2e-12,step=.025)
    report=dict(families=[],second_surfaces=[])
    if not args.surfaces_only and not args.boundaries_only:
        records=[]
        for name in ['cycles','cycles-extra','adding-cycles']:
            path=out/(name+'.json')
            if path.exists():records+=json.loads(path.read_text())['records']
        for r in records:
            if r['status']!='verified_attracting_cycle' or r['half_symmetry_error']>2e-5:continue
            b=r['beta'];v=[*r['cycle']['q'][:2],r['cycle']['period']]
            rec=dict(beta=b,peak_count=r['peak_count'],period=r['cycle']['period'],status='unresolved')
            try:
                family,q=from_orbit(model,b,v,True)
                rec.update(regions=family.regions,targets=family.targets)
                for key,target in [('fold',1),('symmetry',-1)]:
                    beta,qq=boundary(family,b,q,target)
                    rr=family.record(beta,qq,True)
                    if rr['region_violation']>1e-7 or min(qq[2:])<=0:raise RuntimeError('Invalid affine itinerary')
                    rec[key]=rr
                if abs(rec['fold']['beta']-rec['symmetry']['beta'])>.01 and r['peak_count']>5:
                    raise RuntimeError('Unexpected large jump to another closing root')
                rec['status']='two_symmetric_stability_edges'
            except (RuntimeError,ValueError,OverflowError) as e:rec['error']=str(e)
            report['families'].append(rec)
            print('family',b,r['peak_count'],rec['status'],[(k,rec[k]['beta']) for k in ['fold','symmetry'] if k in rec],flush=True)
        (out/'family-edges.json').write_text(json.dumps(report['families'],indent=2)+'\n',encoding='utf-8',newline='\n')
    if args.families_only:return
    # Recover the two templates from the verified beta=14.354 asymmetric orbit.
    r=next(r for r in json.loads((out/'cycles.json').read_text())['records'] if abs(r['beta']-14.354)<1e-10)
    full_regions,full_targets=r['regions'],r['targets']
    half,qq=from_orbit(model,r['beta'],[*r['cycle']['q'][:2],r['cycle']['period']],True)
    full=Flights(model,full_regions,full_targets,False)
    base={}
    for key,target in [('fold',1),('symmetry',-1),('period_doubling',-1)]:
        f=full if key=='period_doubling' else half
        seed=np.array(r['cycle']['q']) if key=='period_doubling' else qq
        b,q=boundary(f,r['beta'],seed,target)
        base[key]=f.record(b,q,True)
        if base[key]['region_violation']>1e-7 or min(q[2:])<=0:
            raise RuntimeError('Invalid central second-family itinerary')
    (out/'second-boundaries.json').write_text(json.dumps(base,indent=2)+'\n',encoding='utf-8',newline='\n')
    if args.boundaries_only:return
    for alpha in np.linspace(8.94,9.06,7):
        for m0 in np.linspace(-8/7-.006,-8/7+.006,7):
            row=dict(alpha=float(alpha),m0=float(m0))
            for key,target in [('fold',1),('symmetry',-1),('period_doubling',-1)]:
                b,q=base[key]['beta'],np.array(base[key]['q'])
                for weight in np.linspace(0,1,25)[1:]:
                    local=Chua(alpha=float(9+weight*(alpha-9)),m0=float(-8/7+weight*(m0+8/7)))
                    symmetric=key!='period_doubling'
                    f=Flights(local,half.regions if symmetric else full_regions,half.targets if symmetric else full_targets,symmetric)
                    b,q=boundary(f,b+(2.14*(alpha-9)-15*(m0+8/7))/24,q,target)
                rr=f.record(b,q,True)
                if rr['region_violation']>1e-7 or min(q[2:])<=0:raise RuntimeError('Invalid surface itinerary')
                row[key]=rr
            report['second_surfaces'].append(row)
        print('second surface',alpha,flush=True)
        (out/'second-surfaces.json').write_text(json.dumps(report['second_surfaces'],indent=2)+'\n',encoding='utf-8',newline='\n')
    (out/'second-surfaces.json').write_text(json.dumps(report['second_surfaces'],indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
