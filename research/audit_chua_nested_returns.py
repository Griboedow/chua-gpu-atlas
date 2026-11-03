"""Independent reconstruction and full-field controls of nested Chua cycles."""
import argparse
import time
import mpmath as mp
import numpy as np
from scipy.signal import find_peaks
from test_chua_nested_returns import OUT,read,save,digest
from chua_block_return import BlockShooting
from audit_chua_layer_predictions import full_field
from audit_chua_return_families import section_variation


def exponential_reconstruction(row,digits=90):
    mp.mp.dps=digits;beta=mp.mpf(row['beta']);q=[mp.mpf(v) for v in row['q']]
    initial=mp.matrix([0,*q[:2]]);state=initial.copy();product=mp.eye(2);mats={}
    for reg in set(row['regions']):
        slope=-mp.mpf(8)/7 if reg==0 else -mp.mpf(5)/7
        A=mp.matrix([[-9*(1+slope),9,0],[1,-1,1],[0,-beta,0]])
        c=mp.matrix([reg*mp.mpf('1.5'),0,-reg*mp.mpf('1.5')]);mats[reg]=A,c
    violation=mp.mpf(0);target_error=mp.mpf(0);speed=mp.inf;values=[];xyz=[];times=[];elapsed=mp.mpf(0);logdiv=mp.mpf(0)
    for i,(reg,dt) in enumerate(zip(row['regions'],q[2:])):
        A,c=mats[reg];E=mp.expm(A*dt);step=mp.expm(A*dt/128);sample=state.copy()
        for j in range(128):
            x=sample[0];values.append(float(x));violation=max(violation,abs(x)-1 if reg==0 else 1-reg*x)
            if j%2==0:xyz.append([float(v) for v in sample]);times.append(float(elapsed+j*dt/128))
            sample=c+step*(sample-c)
        state=c+E*(state-c);v=A*(state-c);speed=min(speed,abs(v[0]));logdiv+=sum(A[j,j] for j in range(3))*dt
        H=mp.matrix([[E[a,b]-v[a]*E[0,b]/v[0] for b in [1,2]] for a in [1,2]])
        product=H*product;target=row['targets'][i] if i<len(row['targets']) else 0
        target_error=max(target_error,abs(state[0]-target));elapsed+=dt
    trace=-(product[0,0]+product[1,1]);det=-mp.exp(logdiv)
    mu=(trace+mp.sqrt(trace*trace-4*det))/2;closure=max(abs(v) for v in state+initial)
    points=np.array(xyz);times=np.array(times);half=float(elapsed)
    plot=dict(xyz=np.vstack([points,-points,points[:1]]),t=np.r_[times,times+half,2*half])
    result=dict(method='mpmath matrix-exponential flights, independent of eigen-flow event shooting',working_digits=digits,
        half_closure=mp.nstr(closure,25),maximum_target_error=mp.nstr(target_error,25),
        sampled_region_violation=mp.nstr(violation,25),minimum_crossing_speed=float(speed),
        signed_half_multiplier=mp.nstr(mu,30),full_multiplier_error=float(abs(mu**2-mp.mpf('.25'))),
        full_peak_count=len(find_peaks(np.r_[values,-np.array(values),values[0]])[0]),
        path_samples_per_flight=128,scope='Numerical fixed-duration reconstruction and sampled region checks; not interval arithmetic')
    return result,plot


def audit(row,calibration):
    id=row['id'];path=OUT/f'{id}-audit.json'
    if path.exists():print('Retained',path.name,flush=True);return
    started=time.monotonic();k=row['block_count'];edge=read(OUT/f'{id}-edges.json');shoot=BlockShooting(calibration,k,85);controls={}
    for name,old,mu in [('reference',row,'.5'),('fold',edge['fold'],1),('symmetry',edge['symmetry'],-1)]:
        w=[mp.mpf(v) for p in old['block_coordinates'] for v in p]+[mp.mpf(old['normalized_coordinate'])]
        new=shoot.solve(w,mu)
        controls[name]=dict(root=new,beta_change=mp.nstr(abs(mp.mpf(new['beta'])-mp.mpf(old['beta'])),25))
        print(id,name,'precision change',controls[name]['beta_change'],flush=True)
    width=abs(mp.mpf(controls['fold']['root']['beta'])-mp.mpf(controls['symmetry']['root']['beta']))
    relative=float(abs(width/mp.mpf(edge['width'])-1))
    precise_path=OUT/f'{id}-precision.json'
    save(precise_path,dict(id=id,working_digits=85,controls=controls,relative_width_change=relative,
        source_sha256=digest(OUT/f'{id}-reference.json')))
    independent,plot=exponential_reconstruction(row)
    np.savez_compressed(OUT/f'{id}-orbit.npz',**plot)
    save(OUT/f'{id}-exponential.json',independent)
    whole=[full_field(dict(representative=row),'DOP853',.005)]
    if k==7:whole.append(full_field(dict(representative=row),'DOP853',.0025))
    print(id,'whole-field',whole,flush=True)
    local=section_variation(row,.005);points=np.array(local['points']);dist=np.linalg.norm(points[:,None,:]-points[None,:,:],axis=2);np.fill_diagonal(dist,np.inf)
    local_controls=[local]
    if local['multiplier_error']>2e-5:local_controls.append(section_variation(row,.0025));local=local_controls[-1]
    gates=dict(precision=max(float(v['beta_change']) for v in controls.values())<1e-38 and relative<1e-25,
        independent_closure=float(independent['half_closure'])<1e-35 and float(independent['maximum_target_error'])<1e-35,
        path=float(independent['sampled_region_violation'])<1e-35 and independent['minimum_crossing_speed']>.01,
        independent_stability=independent['full_multiplier_error']<1e-28,
        peak_counts=independent['full_peak_count']==row['forecast_peaks'] and all(w['peaks']==row['forecast_peaks'] for w in whole),
        section_counts=len(points)==7*k and all(w['section_returns']==7*k for w in whole),
        primitive=bool(dist.min()>1e-12),
        section_field_control=local['maximum_next_point_error']<1e-10 and local['multiplier_error']<2e-4)
    result=dict(id=id,source_sha256=digest(OUT/f'{id}-reference.json'),precision=precise_path.name,
        exponential_reconstruction=independent,whole_period=whole,section_controls=local_controls,
        minimum_distinct_section_distance=float(dist.min()),gates=gates,passed=all(gates.values()),
        elapsed_seconds=time.monotonic()-started,
        scope='Whole-period drift is reported separately. Independent high-precision closure, path checks and section-field derivatives determine the numerical orbit acceptance.')
    save(path,result);print(id,'AUDIT',gates,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ids',nargs='*');args=p.parse_args()
    calibration=read(OUT/'calibration.json')
    for path in sorted(OUT.glob('*-reference.json')):
        row=read(path)
        if args.ids and row['id'] not in args.ids:continue
        parent=next(c for c in calibration['parents'] if c['id']==row['parent']);audit(row,parent)


if __name__=='__main__':main()
