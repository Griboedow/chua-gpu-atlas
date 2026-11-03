"""Precision repeats and an independent eigendecomposition of new orbit flights."""
import argparse
import json
import mpmath as mp
import numpy as np
from scipy.signal import find_peaks
from chua_precision_flights import PrecisionFlights
from analyze_chua_return_geometry import sections
from study_chua_return_geometry import OUT, read, digest


def edge_checks(id):
    source=OUT/f'{id}-edges.json';row=read(OUT/f'{id}-reference.json')
    previous=read(source);f=PrecisionFlights(row['regions'],row['targets'],True,95)
    controls={}
    for name,target in [('fold',1),('symmetry',-1)]:
        old=previous[name];new=f.boundary(old['beta'],old['q'],target)
        controls[name]=dict(root=new,beta_change=mp.nstr(abs(mp.mpf(new['beta'])-mp.mpf(old['beta'])),20))
    width=mp.mpf(controls['fold']['root']['beta'])-mp.mpf(controls['symmetry']['root']['beta'])
    relative=abs(width/mp.mpf(previous['width'])-1)
    result=dict(id=id,working_digits=95,source=source.name,source_sha256=digest(source),
        controls=controls,width=mp.nstr(width,85),relative_width_change=mp.nstr(relative,20),
        passed=bool(relative<mp.mpf('1e-40') and all(mp.mpf(c['beta_change'])<mp.mpf('1e-50') for c in controls.values())),
        scope='Numerical precision repeat of both window endpoints; not interval arithmetic')
    (OUT/f'{id}-edge-precision-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(id,'edge relative width change',result['relative_width_change'],'passed',result['passed'],flush=True)
    if not result['passed']:raise SystemExit('Edge precision check failed')


def full_field_checks(id):
    """Replay the separate whole-period and section-by-section controls."""
    from audit_chua_layer_predictions import full_field, local_return_variation
    path=OUT/f'{id}-independent-orbit.json'
    if path.exists():
        print('Retained',path.name,flush=True);return
    row=read(OUT/f'{id}-reference.json')
    record=dict(representative=row,regions=row['regions'],targets=row['targets'])
    result=dict(whole_period=full_field(record,'DOP853',.00125),
                section_variation=local_return_variation(record,.005))
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(id,'whole-period closure',result['whole_period']['closure'],flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edges',action='store_true')
    parser.add_argument('--full-field',action='store_true')
    args=parser.parse_args()
    for id in ['P95','P105']:
        if args.edges:
            edge_checks(id);continue
        if args.full_field:
            full_field_checks(id);continue
        row=read(OUT/f'{id}-reference.json')
        f=PrecisionFlights(row['regions'],row['targets'],True,95)
        refined=f.boundary(row['beta'],row['q'],'.5')
        beta=mp.mpf(row['beta']);q=[mp.mpf(x) for x in row['q']]
        mats={}
        for r in set(row['regions']):
            slope=-mp.mpf(8)/7 if r==0 else -mp.mpf(5)/7
            A=mp.matrix([[-9*(1+slope),9,0],[1,-1,1],[0,-beta,0]])
            c=mp.matrix([r*mp.mpf('1.5'),0,-r*mp.mpf('1.5')])
            # m1-m0 = 3/7 and 1+m1 = 2/7, hence outer equilibria ±3/2.
            eig,V=mp.eig(A);W=V**-1
            mats[r]=(A,c,eig,V,W)
        initial=mp.matrix([0,q[0],q[1]]);u=initial.copy();M=mp.eye(3)
        largest_target=mp.mpf(0);violation=mp.mpf(0);values=[]
        for k,(r,dt) in enumerate(zip(row['regions'],q[2:])):
            A,c,eig,V,W=mats[r]
            E=V*mp.diag([mp.exp(l*dt) for l in eig])*W
            step=V*mp.diag([mp.exp(l*dt/128) for l in eig])*W
            sample=u.copy()
            for j in range(128):
                x=mp.re(sample[0]);values.append(float(x))
                violation=max(violation,abs(x)-1 if r==0 else 1-r*x)
                sample=c+step*(sample-c)
            u=c+E*(u-c);M=E*M
            if k<len(row['targets']):largest_target=max(largest_target,abs(u[0]-row['targets'][k]))
        velocity=mats[row['regions'][-1]][0]*(u-mats[row['regions'][-1]][1])
        project=mp.eye(3)
        for i in range(3):project[i,0]-=velocity[i]/velocity[0]
        P=-project*M;trace=mp.re(P[1,1]+P[2,2])
        closed=max(abs(v) for v in u+initial)
        xyz=np.r_[values,-np.array(values),values[0]]
        peaks=len(find_peaks(xyz)[0])
        points,t=sections(row,row)
        distances=np.linalg.norm(points[:,None,:]-points[None,:,:],axis=2)
        distances+=np.eye(len(points))*1e10
        result=dict(id=id,working_digits=95,beta_change=float(abs(mp.mpf(refined['beta'])-beta)),
            maximum_shooting_coordinate_change=float(max(abs(mp.mpf(a)-mp.mpf(b)) for a,b in zip(refined['q'],row['q']))),
            eigen_flow_half_closure=mp.nstr(closed,14),maximum_target_error=mp.nstr(largest_target,14),
            sampled_region_violation=mp.nstr(violation,14),signed_half_trace=mp.nstr(trace,24),
            full_dominant_multiplier_error=float(abs(trace**2-mp.mpf('.25'))),
            sampled_full_period_peak_count=peaks,half_period_normalized_section_count=len(points),
            minimum_distinct_normalized_section_distance=float(distances.min()),
            scope='Numerical high-precision eigen-flow reconstruction independent of mp.expm; path checks are sampled, not interval enclosures')
        result['passed']=bool(closed<mp.mpf('1e-40') and largest_target<mp.mpf('1e-40') and violation<mp.mpf('1e-40') and peaks==int(id[1:]) and distances.min()>1e-8 and result['full_dominant_multiplier_error']<1e-30)
        (OUT/f'{id}-precision-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(result,flush=True)
        if not result['passed']:raise SystemExit('New orbit check failed')


if __name__=='__main__':main()
