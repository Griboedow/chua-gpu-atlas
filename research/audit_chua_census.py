"""Audit census evidence, boundary degeneracies, and independent multipliers."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from trace_chua_window import Chua
from chua_layer_bifurcations import Flights

ROOT=Path(__file__).resolve().parents[1]


def independent(r):
    model=Chua(alpha=r['alpha'],m0=r['m0'],m1=r['m1'])
    beta=r['beta'];u=np.r_[0.,r['q'][:2]]
    def rhs(t,v):
        return np.r_[model.field(beta,v[:3]),(model.jacobian(beta,abs(v[0])<1)@v[3:].reshape(3,3)).ravel()]
    sol=solve_ivp(rhs,[0,r['period']],np.r_[u,np.eye(3).ravel()],method='DOP853',rtol=2e-12,atol=2e-13,max_step=.02)
    f=model.field(beta,sol.y[:3,-1]);M=sol.y[3:,-1].reshape(3,3)
    P=((np.eye(3)-np.outer(f,[1.,0.,0.])/f[0])@M)[1:,1:]
    eig=np.linalg.eigvals(P);large=eig[np.argmax(abs(eig))]
    return dict(success=sol.success,closure=float(max(abs(sol.y[:3,-1]-u))),
                dominant=float(large.real),expected=r['full_multipliers'][0][0],
                multiplier_error=float(abs(large-r['full_multipliers'][0][0])))


def fold_controls(r):
    f=Flights(Chua(alpha=r['alpha'],m0=r['m0'],m1=r['m1']),
              [a['region'] for a in r['arcs']],[round(a['end'][0]) for a in r['arcs'][:-1]],True)
    b=r['beta'];q=np.array(r['q']);U,S,V=np.linalg.svd(f.evaluate(b,q)[1]);left=U[:,-1];right=V[-1]
    checks=[]
    for h in [1e-4,3e-5,1e-5]:
        residual=f.evaluate(b,q)[0]
        rb=(f.evaluate(b+h,q)[0]-f.evaluate(b-h,q)[0])/(2*h)
        rqq=(f.evaluate(b,q+h*right)[0]-2*residual+f.evaluate(b,q-h*right)[0])/h**2
        checks.append(dict(step=h,parameter=float(left@rb),quadratic=float(left@rqq)))
    return dict(singular_values=S.tolist(),coefficients=checks)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/chua-layer-census')
    out=p.parse_args().output
    analysis=json.loads((out/'census-analysis.json').read_text())
    records=[r for name in ['cycles','cycles-extra','adding-cycles'] for r in json.loads((out/(name+'.json')).read_text())['records']]
    verified=[r for r in records if r['status']=='verified_attracting_cycle']
    surfaces=[r[k] for r in json.loads((out/'second-surfaces.json').read_text()) for k in ['fold','symmetry','period_doubling']]
    central=json.loads((out/'second-boundaries.json').read_text())
    picks=[('second_'+k,v) for k,v in central.items()]
    picks += [(str(r['maxima'])+'_'+k,r[k]) for r in analysis['adding_sequence'] for k in ['fold','symmetry']]
    checks=[];folds={}
    for name,r in picks:
        check=dict(name=name,**independent(r));checks.append(check)
        if name.endswith('fold'):folds[name]=fold_controls(r)
        print(name,check,flush=True)
    fold=central['fold'];f=Flights(Chua(),[a['region'] for a in fold['arcs']],
                                 [round(a['end'][0]) for a in fold['arcs'][:-1]],True)
    q=np.array(fold['q']);direction=np.linalg.svd(f.evaluate(fold['beta'],q)[1])[2][-1]
    partners=[f.record(fold['beta']-1e-4,f.solve(fold['beta']-1e-4,q+s*.01*direction),True) for s in [-1,1]]
    coexist=json.loads((out/'coexistence.json').read_text())
    metrics=dict(verified_cycles=len(verified),distinct_transect_anchors=len(analysis['anchors']),
        max_radau_closure=max(r.get('radau_refinement',{}).get('closure',r['radau_closure']) for r in verified),
        max_variational_closure=max(r['independent_variational_closure'] for r in verified),
        max_variation_affine_rho_difference=max(r['variation_vs_affine_rho'] for r in verified),
        max_cycle_rho=max(r['independent_variational_rho'] for r in verified),
        max_cycle_region_violation=max(r['cycle']['region_violation'] for r in verified),
        max_surface_residual=max(r['residual'] for r in surfaces),
        max_surface_region_violation=max(r['region_violation'] for r in surfaces),
        min_surface_crossing_speed=min(r['minimum_breakpoint_speed'] for r in surfaces),
        max_boundary_closure=max(r['closure'] for r in checks),
        max_boundary_multiplier_error=max(r['multiplier_error'] for r in checks))
    gates=dict(all_variational_integrations_succeeded=all(r['independent_variational_success'] for r in verified) and all(r['success'] for r in checks),
        radau_closure=metrics['max_radau_closure']<2e-7,
        variational_closure=metrics['max_variational_closure']<2e-7,
        independent_multiplier_agreement=metrics['max_variation_affine_rho_difference']<2e-5,
        attracting_cycles=metrics['max_cycle_rho']<.9999,
        cycle_itineraries=metrics['max_cycle_region_violation']<1e-7,
        surface_closure=metrics['max_surface_residual']<1e-8,
        surface_itineraries=metrics['max_surface_region_violation']<1e-7,
        surface_crossings=metrics['min_surface_crossing_speed']>.1,
        boundary_closure=metrics['max_boundary_closure']<2e-7,
        boundary_multipliers=metrics['max_boundary_multiplier_error']<2e-5,
        fold_nondegeneracy=all(min(abs(x['parameter']) for x in c['coefficients'])>1e-4 and min(abs(x['quadratic']) for x in c['coefficients'])>1e-4 for c in folds.values()),
        distinct_second_fold_partners=abs(partners[0]['period']-partners[1]['period'])>1e-3,
        stable_and_saddle_partner=min(r['spectral_radius'] for r in partners)<1<max(r['spectral_radius'] for r in partners),
        gap_screens=all('witness_beta' in r and np.min(r['four_window_lle'])>.01 for r in analysis['gap_witnesses']),
        coexistence=coexist['verified'])
    report=dict(status='passed' if all(gates.values()) else 'failed',gates=gates,metrics=metrics,
                independent_boundaries=checks,fold_controls=folds,second_fold_partners=partners,
                scope='Numerical checks, not interval proofs. Gap screens do not exclude hidden periodic basins.')
    (out/'audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(report['status'],metrics,gates,flush=True)
    if not all(gates.values()):raise SystemExit(2)


if __name__=='__main__':main()
