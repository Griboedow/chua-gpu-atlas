"""Independent orbit, derivative and precision controls for the new sequences."""
import argparse
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from test_chua_return_families import OUT, read, save, digest, fixed_root, parent_at
from chua_precision_flights import PrecisionFlights
from analyze_chua_return_geometry import sections
from audit_chua_layer_predictions import full_field
from check_chua_return_geometry import variational
from trace_chua_window import Chua
from study_chua_companion_saddles import path_checks
from sample_chua_return_normal_form import event_return


def section_variation(row,step=.005):
    points,times=sections(row,row);model=Chua();beta=float(row['beta'])
    mp.mp.dps=60;product=mp.eye(2);det=mp.mpf(1);checks=[]
    for i,p in enumerate(points):
        def rhs(t,v):
            return np.r_[model.field(beta,v[:3]),(model.jacobian(beta,abs(v[0])<1)@v[3:].reshape(3,3)).ravel()]
        def event(t,v):return v[0]
        event.direction=-1;event.terminal=True
        sol=solve_ivp(rhs,[0,30],np.r_[0.,p,np.eye(3).ravel()],method='DOP853',
                      rtol=2.5e-14,atol=2.5e-15,max_step=step,events=event)
        assert sol.success and len(sol.t_events[0])==1
        u=sol.y[:3,-1];M=sol.y[3:,-1].reshape(3,3);v=model.field(beta,u)
        H=-((np.eye(3)-np.outer(v,[1.,0.,0.])/v[0])@M)[1:,1:]
        mat=mp.matrix(H.tolist());product=mat*product;det*=mp.det(mat)
        checks.append(dict(index=i,duration=float(sol.t[-1]),matrix=H.tolist(),
            next_point_error=float(max(abs(-u[1:]-points[(i+1)%len(points)])))))
    tr=product[0,0]+product[1,1];mu=(tr+mp.sqrt(tr*tr-4*det))/2
    # A single Floquet direction is transported through every section. This
    # trajectory is not a curve of independently optimized operator norms.
    P=np.array(product.tolist(),float);eig,V=np.linalg.eig(P)
    direction=np.real(V[:,np.argmax(abs(eig))]);direction/=np.linalg.norm(direction)
    amplitudes=[1.];factors=[]
    for c in checks:
        direction=np.array(c['matrix'])@direction
        amplitude=float(np.linalg.norm(direction));factors.append(amplitude/amplitudes[-1]);amplitudes.append(amplitude)
    return dict(method='DOP853 full-field section-by-section variational integration',max_step=step,
        points=points.tolist(),times=times.tolist(),checks=checks,
        signed_half_multiplier=float(mp.re(mu)),full_multiplier=float(mp.re(mu**2)),
        multiplier_error=float(abs(mu**2-mp.mpf('.25'))),
        maximum_next_point_error=max(c['next_point_error'] for c in checks),
        period_difference=abs(2*sum(c['duration'] for c in checks)-float(row['period'])),
        floquet_direction_amplitudes=amplitudes,floquet_direction_factors=factors,
        contracting_returns=sum(v<1 for v in factors),
        scope='Euclidean section (y,z) norm; one transported Floquet direction; next-section errors are not whole-period drift')


def audit(row):
    id=row['id'];path=OUT/f'{id}-audit.json'
    if path.exists():return
    edge=read(OUT/f'{id}-edges.json');f=PrecisionFlights(row['regions'],row['targets'],True,90)
    precision={}
    for name,old,mu in [('reference',row,.5),('fold',edge['fold'],1),('symmetry',edge['symmetry'],-1)]:
        new=f.boundary(old['beta'],old['q'],mu)
        precision[name]=dict(root=new,beta_change=mp.nstr(abs(mp.mpf(new['beta'])-mp.mpf(old['beta'])),20))
    width=abs(mp.mpf(precision['fold']['root']['beta'])-mp.mpf(precision['symmetry']['root']['beta']))
    relative_width=float(abs(width/mp.mpf(edge['width'])-1))
    # Independent eigendecomposition, rather than mp.expm, for each affine flow.
    beta=mp.mpf(row['beta']);q=[mp.mpf(v) for v in row['q']];mats={}
    for reg in set(row['regions']):
        slope=-mp.mpf(8)/7 if reg==0 else -mp.mpf(5)/7
        A=mp.matrix([[-9*(1+slope),9,0],[1,-1,1],[0,-beta,0]])
        c=mp.matrix([reg*mp.mpf('1.5'),0,-reg*mp.mpf('1.5')]);eig,V=mp.eig(A)
        mats[reg]=(A,c,eig,V,V**-1)
    initial=mp.matrix([0,*q[:2]]);u=initial.copy();M=mp.eye(3)
    violation=mp.mpf(0);target_error=mp.mpf(0);values=[];speeds=[]
    for k,(reg,dt) in enumerate(zip(row['regions'],q[2:])):
        A,c,eig,V,W=mats[reg];E=V*mp.diag([mp.exp(l*dt) for l in eig])*W
        step=V*mp.diag([mp.exp(l*dt/128) for l in eig])*W;sample=u.copy()
        for j in range(128):
            x=mp.re(sample[0]);values.append(float(x))
            violation=max(violation,abs(x)-1 if reg==0 else 1-reg*x);sample=c+step*(sample-c)
        u=c+E*(u-c);M=E*M;speeds.append(abs((A*(u-c))[0]))
        if k<len(row['targets']):target_error=max(target_error,abs(u[0]-row['targets'][k]))
    closure=max(abs(v) for v in u+initial);v=mats[0][0]*u;project=mp.eye(3)
    for i in range(3):project[i,0]-=v[i]/v[0]
    P=-project*M;trace=mp.re(P[1,1]+P[2,2])
    peaks=len(find_peaks(np.r_[values,-np.array(values),values[0]])[0])
    eigen=dict(closure=mp.nstr(closure,20),target_error=mp.nstr(target_error,20),
        sampled_region_violation=mp.nstr(violation,20),minimum_crossing_speed=float(min(speeds)),
        half_trace=mp.nstr(trace,25),full_multiplier_error=float(abs(trace**2-mp.mpf('.25'))),peaks=peaks)
    record=dict(representative=row,regions=row['regions'],targets=row['targets'])
    print(id,'90-digit and eigen-flow controls complete',flush=True)
    whole=[full_field(record,'DOP853',.005)]
    if row['n']==5:whole.append(full_field(record,'DOP853',.0025))
    local=section_variation(row)
    points=np.array(local['points']);dist=np.linalg.norm(points[:,None,:]-points[None,:,:],axis=2)
    np.fill_diagonal(dist,np.inf)
    gates=dict(precision=max(float(v['beta_change']) for v in precision.values())<1e-45 and relative_width<1e-35,
        eigen_closure=float(closure)<1e-40 and float(target_error)<1e-40,
        path=float(violation)<1e-40 and min(speeds)>.01,
        eigen_stability=eigen['full_multiplier_error']<1e-30,
        counts=peaks==row['forecast_peaks'] and all(w['peaks']==peaks and w['section_returns']==2*row['n']+1 for w in whole),
        primitive=bool(dist.min()>1e-8),
        full_field=all(w['closure']<2e-6 and w['half_inversion_error']<2e-6 for w in whole),
        section_control=local['maximum_next_point_error']<1e-10 and local['multiplier_error']<2e-5)
    save(path,dict(id=id,source_sha256=digest(OUT/f'{id}-reference.json'),working_digits=90,
        precision=precision,relative_width_change=relative_width,eigen_reconstruction=eigen,
        whole_period=whole,section_variation=local,minimum_distinct_section_distance=float(dist.min()),gates=gates,passed=all(gates.values())))
    print(id,'audit',gates,'contractions',local['contracting_returns'],'factors',local['floquet_direction_factors'],flush=True)


def companion(row):
    id=row['id'];path=OUT/f'{id}-companion.json'
    if path.exists():return
    mp.mp.dps=85;edge=read(OUT/f'{id}-edges.json');q=mp.matrix([mp.mpf(v) for v in row['q']])
    fold=mp.matrix([mp.mpf(v) for v in edge['fold']['q']])
    saddle,ev=fixed_root(row,row['beta'],2*fold-q,85)
    f=PrecisionFlights(row['regions'],row['targets'],True,85)
    checks=path_checks(f,mp.mpf(row['beta']),ev)
    d=mp.matrix([mp.mpf(v) for v in saddle['q'][:2]])-q[:2,:]
    out=dict(id=id,companion=saddle,path_checks=checks,section_distance=mp.nstr(mp.norm(d),30),
        half_multiplier_error=abs(float(saddle['signed_map_multiplier'])-1.5),
        passed=bool(float(saddle['signed_map_multiplier'])>1 and mp.norm(d)>mp.mpf('1e-14') and checks['sampled_region_violation']<1e-9))
    save(path,out);print(id,'companion',saddle['signed_map_multiplier'],flush=True)


def normal_form(row,digits):
    id=row['id'];path=OUT/f'{id}-normal-form-d{digits}.json'
    if path.exists():return
    saddle=read(OUT/f'{id}-companion.json')['companion'];mp.mp.dps=digits
    q=mp.matrix([mp.mpf(v) for v in row['q'][:2]]);d=mp.matrix([mp.mpf(v) for v in saddle['q'][:2]])-q
    dd=mp.fdot(d,d);samples=[]
    for ss in ['-.5','-.25','0','.25','.5','.75','1','1.25','1.5']:
        s=mp.mpf(ss);value,checks=event_return(row,q+s*d,digits)
        coordinate=mp.fdot(value-q,d)/dd;pred=(s+s*s)/2
        samples.append(dict(s=float(s),coordinate=mp.nstr(coordinate,30),predicted=float(pred),
            error=float(abs(coordinate-pred)),off_secant=float(mp.norm(value-q-coordinate*d)/mp.sqrt(dd)),checks=checks))
    save(path,dict(id=id,digits=digits,samples=samples,maximum_error=max(v['error'] for v in samples),
        maximum_off_secant=max(v['off_secant'] for v in samples)))
    print(id,'normal form',max(v['error'] for v in samples),flush=True)


def gain_control(row):
    id=row['id'];path=OUT/f'{id}-gain-control.json'
    if path.exists():return
    gain=read(OUT/f'{id}-gain-s8.json');gain['beta']=row['beta']
    with np.load(OUT/f'{id}-gain-s8.npz') as z:data={k:z[k] for k in z.files}
    checks=[variational(gain,data,s) for s in [.025,.0125]]
    save(path,dict(id=id,checks=checks,passed=all(c['relative_gain_error']<2e-5 and c['maximum_local_state_error']<1e-10 for c in checks)))
    print(id,'gain errors',[c['relative_gain_error'] for c in checks],flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['orbits','companions','normal','gain'])
    p.add_argument('--ids',nargs='*');p.add_argument('--digits',type=int,default=75);args=p.parse_args()
    for path in sorted(OUT.glob('*-reference.json')):
        row=read(path)
        if args.ids and row['id'] not in args.ids:continue
        if args.stage=='orbits':audit(row)
        elif args.stage=='companions':companion(row)
        elif args.stage=='normal':normal_form(row,args.digits)
        else:gain_control(row)


if __name__=='__main__':main()
