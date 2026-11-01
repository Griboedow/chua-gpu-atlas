"""Close denser family samples around the visible atlas box for README figures.

Interpolation selects where to spend solves; every retained point is closed
and checked independently with the original affine matrix-exponential flow.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from trace_chua_window import Chua
from chua_fast_flights import FastFlights
from chua_layer_bifurcations import Flights
from bridge_chua_itineraries import local_step
from audit_chua_census import independent

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/chua-layer-atlas'


def slim(r):
    return {k:v for k,v in r.items() if k not in ['arcs','section_matrix']}


def normalize(sheet,node,target):
    symmetric=sheet.get('symmetric',True)
    regions=node.get('regions',sheet['regions']);targets=node.get('targets',sheet['targets'])
    f=FastFlights(Chua(alpha=node['alpha'],m0=node['m0']),regions,targets,symmetric)
    old=node['representative']
    if sheet.get('target_multiplier')==target:beta,q=old['beta'],np.array(old['q'])
    else:beta,q=f.boundary(old['beta'],np.array(old['q']),target)
    r=f.record(beta,q,True)
    if r['region_violation']>1e-7 or abs(r['period']/old['period']-1)>.01 or r['spectral_radius']>=.999:
        raise RuntimeError('Normalization did not retain a nearby attracting orbit')
    return dict(alpha=node['alpha'],m0=node['m0'],representative=r,regions=regions,targets=targets)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--families',nargs='*');p.add_argument('--resume',action='store_true')
    p.add_argument('--p65-half-multiplier',type=float,help='Override the existing P65 display target, e.g. 0.5 to replay the retained preliminary attempt.')
    args=p.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    aa=np.linspace(8.94,9.06,17);mm=np.linspace(-8/7-.006,-8/7+.006,17)
    registry=json.loads((SOURCE/'family-registry.json').read_text())['families'];summaries=[]
    for family in sorted(registry,key=lambda r:r['id']):
        id=family['id']
        if args.families and id not in args.families:continue
        path=out/(id+'.json')
        if args.resume and path.exists():
            summaries.append(json.loads(path.read_text())['summary']);continue
        start=time.monotonic();sheet=json.loads((SOURCE/(id+'-sheet.json')).read_text())
        old=sheet['nodes'];coords=np.array([[(n['alpha']-9)/.06,(n['m0']+8/7)/.006] for n in old])
        interp=LinearNDInterpolator(coords,[n['representative']['beta'] for n in old])
        x,y=np.meshgrid((aa-9)/.06,(mm+8/7)/.006,indexing='ij');pred=interp(x,y);selected=np.zeros((17,17),bool)
        # Include a full cell around each visible interval, for subsequent
        # geometric clipping. Missing source coverage is never extrapolated.
        for i in range(16):
            for j in range(16):
                cell=pred[i:i+2,j:j+2]
                if np.isfinite(cell).all() and cell.min()<=14.402 and cell.max()>=14.238:selected[i:i+2,j:j+2]=True
        # P65 already has a well-conditioned interior representative at a
        # negative half-return multiplier. Moving it closer to the fold solely
        # to match the catalogue's identity convention is unnecessary for a view.
        cache={};nodes=[];failures=[];target=sheet['target_multiplier'] if id=='P65' else .5 if family['symmetric'] else -.5
        if id=='P65' and args.p65_half_multiplier is not None:target=args.p65_half_multiplier
        for i,j in zip(*np.where(selected)):
            point=np.array([(aa[i]-9)/.06,(mm[j]+8/7)/.006]);nearest=int(np.argmin(np.sum((coords-point)**2,axis=1)))
            try:
                if nearest not in cache:cache[nearest]=normalize(sheet,old[nearest],target)
                seed=cache[nearest];template=dict(symmetric=family['symmetric'],target_multiplier=target,regions=seed['regions'],targets=seed['targets'])
                if abs(seed['alpha']-aa[i])<1e-12 and abs(seed['m0']-mm[j])<1e-12:node=seed
                else:node=local_step(template,seed,float(aa[i]),float(mm[j]))
                r=node['representative'];f=Flights(Chua(alpha=r['alpha'],m0=r['m0']),node['regions'],node['targets'],family['symmetric'])
                check=f.record(r['beta'],np.array(r['q']));error=abs(check['spectral_radius']-r['spectral_radius'])
                if r['region_violation']>1e-7 or check['residual']>2e-8 or error>2e-5 or check['spectral_radius']>=.999:
                    raise RuntimeError('Independent affine check failed')
                nodes.append(dict(index=[int(i),int(j)],representative=slim(r),regions=node['regions'],targets=node['targets'],
                    expm_closure=check['residual'],expm_rho_error=error,new_parameter_location=bool(i%2 or j%2)))
            except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:
                failures.append(dict(index=[int(i),int(j)],error=str(e)))
        assert nodes,id+' has no checked samples'
        controls=[]
        if id in ['E09','E13','P65','A05','S08']:
            inside=[n for n in nodes if 14.24<=n['representative']['beta']<=14.4 and n['new_parameter_location']]
            for n in [min(inside,key=lambda n:n['representative']['m0']),max(inside,key=lambda n:n['representative']['m0'])]:
                c=independent(n['representative']);controls.append(dict(index=n['index'],**c))
        summary=dict(id=id,nodes=len(nodes),selected_nodes=int(selected.sum()),new_parameter_locations=sum(n['new_parameter_location'] for n in nodes),
            inside_nodes=sum(14.24<=n['representative']['beta']<=14.4 for n in nodes),failures=len(failures),
            max_expm_closure=max(n['expm_closure'] for n in nodes),max_expm_rho_error=max(n['expm_rho_error'] for n in nodes),
            independent_controls=controls,elapsed_seconds=time.monotonic()-start)
        report=dict(id=id,alpha=aa.tolist(),m0=mm.tolist(),symmetric=family['symmetric'],target_multiplier=target,
            nodes=nodes,failures=failures,summary=summary,
            scope='Closed attracting representatives at constant transverse multiplier; interpolation only selects sample locations. Rendering facets do not certify attraction between vertices.')
        path.write_text(json.dumps(report,separators=(',',':'))+'\n',encoding='utf-8',newline='\n');summaries.append(summary)
        print(id,'closed',len(nodes),'new locations',summary['new_parameter_locations'],'failures',len(failures),'seconds',round(summary['elapsed_seconds'],1),flush=True)
    if args.families:return
    controls=[c for r in summaries for c in r['independent_controls']]
    gates=dict(all_families=len(summaries)==32,affine_closure=all(r['max_expm_closure']<2e-8 for r in summaries),
        selected_coverage=all(r['nodes']==r['selected_nodes'] and not r['failures'] for r in summaries),
        affine_multipliers=all(r['max_expm_rho_error']<2e-5 for r in summaries),
        full_field_controls=len(controls)==10 and all(c['success'] and c['closure']<2e-7 and c['multiplier_error']<2e-5 and abs(c['dominant'])<.999 for c in controls))
    report=dict(status='passed' if all(gates.values()) else 'failed',gates=gates,families=summaries,
        total_nodes=sum(r['nodes'] for r in summaries),new_parameter_locations=sum(r['new_parameter_locations'] for r in summaries),
        inside_nodes=sum(r['inside_nodes'] for r in summaries),failures=sum(r['failures'] for r in summaries),
        scope='Targeted display refinement of the existing 32-family catalogue, not a new exhaustive search.')
    (out/'audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    if not all(gates.values()):raise SystemExit('Refinement audit failed')
    print('PASS',report['total_nodes'],'samples;',report['new_parameter_locations'],'new parameter locations',flush=True)


if __name__=='__main__':main()
