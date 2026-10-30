"""Check every continued point against matrix exponentials and selected ODE flows."""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_layer_bifurcations import Flights
from audit_chua_census import independent


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output;rows=[];features=[];arrays={}
    registry=json.loads((out/'family-registry.json').read_text())['families'];lookup={r['id']:r for r in registry}
    for path in sorted(out.glob('*-sheet.json')):
        sheet=json.loads(path.read_text());original=[]
        for node in sheet['nodes']:
            model=Chua(alpha=node['alpha'],m0=node['m0']);f=Flights(model,node.get('regions',sheet['regions']),node.get('targets',sheet['targets']),sheet.get('symmetric',True))
            for key in ['core','fold','symmetry']:
                if key not in node:continue
                r=node[key];check=f.record(r['beta'],np.array(r['q']))
                original.append(dict(key=key,residual=check['residual'],rho_difference=abs(check['spectral_radius']-r['spectral_radius'])))
        inside=[n for n in sheet['nodes'] if 14.24<=n['core']['beta']<=14.4]
        assert inside,'No representative inside original atlas'
        selected=[min(inside,key=lambda n:n['m0']),max(inside,key=lambda n:n['m0'])]
        controls=[]
        for node in selected:
            r=node['core'];check=independent(r);controls.append(dict(alpha=node['alpha'],m0=node['m0'],beta=r['beta'],**check))
        root=lookup[sheet['id']]['core'];model=Chua(rtol=2e-12,step=.02)
        u=np.r_[0.,root['q'][:2]];T=root['period'];orbit=model.integrate(root['beta'],u,T)
        t=np.linspace(0,T,8001);xyz=orbit.sol(t).T
        feature=dict(id=sheet['id'],reference_beta=root['beta'],period=T,maxima=len(orbit.t_events[3]),
            symmetric=lookup[sheet['id']]['symmetric'],section_returns=lookup[sheet['id']]['section_returns'],
            x_min=float(xyz[:,0].min()),x_max=float(xyz[:,0].max()),x_mean=float(np.trapezoid(xyz[:,0],t)/T),
            x_rms=float(np.sqrt(np.trapezoid(xyz[:,0]**2,t)/T)),rho=root['spectral_radius'],
            relaxation_time=float(-T/np.log(root['spectral_radius'])),
            both_outer_regions=bool(xyz[:,0].min()<-1 and xyz[:,0].max()>1))
        features.append(feature);arrays[sheet['id']+'_t']=t[::2];arrays[sheet['id']+'_xyz']=xyz[::2]
        row=dict(id=sheet['id'],nodes=len(sheet['nodes']),inside_nodes=len(inside),complete=sheet['complete'],
            max_expm_closure=max(r['residual'] for r in original),max_expm_rho_difference=max(r['rho_difference'] for r in original),
            max_recorded_region_violation=max(n['core']['region_violation'] for n in sheet['nodes']),independent_controls=controls)
        rows.append(row);(out/'atlas-audit-partial.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(row['id'],row['nodes'],'expm closure',row['max_expm_closure'],'rho difference',row['max_expm_rho_difference'],
              'ODE error',max(c['multiplier_error'] for c in controls),flush=True)
    controls=[c for r in rows for c in r['independent_controls']]
    gates=dict(all_families_present=len(rows)==len(registry),each_family_intersects_original_box=all(r['inside_nodes']>0 for r in rows),
               original_expm_closure=max(r['max_expm_closure'] for r in rows)<2e-8,
               original_expm_multipliers=max(r['max_expm_rho_difference'] for r in rows)<2e-5,
               sampled_region_validity=max(r['max_recorded_region_violation'] for r in rows)<1e-7,
               independent_solvers=all(c['success'] for c in controls),
               independent_closure=max(c['closure'] for c in controls)<2e-7,
               independent_multipliers=max(c['multiplier_error'] for c in controls)<2e-5,
               independent_attraction=all(abs(c['dominant'])<.9999 for c in controls))
    report=dict(status='passed' if all(gates.values()) else 'failed',gates=gates,families=rows,total_nodes=sum(r['nodes'] for r in rows),
        incomplete_parameter_grids=[r['id'] for r in rows if not r['complete']],
        scope='Numerical continuation checks, not interval proofs or exhaustive family enumeration.')
    (out/'atlas-audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    (out/'family-behavior.json').write_text(json.dumps(features,indent=2)+'\n',encoding='utf-8',newline='\n')
    np.savez_compressed(out/'family-orbits.npz',**arrays)
    print(report['status'],gates,flush=True)
    if not all(gates.values()):raise SystemExit(2)


if __name__=='__main__':main()
