"""Calibrate, freeze, then test two additional Chua period-adding sequences."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
import mpmath as mp
import numpy as np
from scipy.optimize import minimize
from chua_precision_flights import PrecisionFlights
from chua_fast_flights import FastFlights
from trace_chua_window import Chua
from study_chua_return_geometry import ROOT, read, digest, mesh, amplification
from refine_chua_transient_maxima import arcs_for, transport

OUT=ROOT/'results/chua-return-families'
SOURCES=ROOT/'results/chua-knot-atlas'
SEQUENCES=[dict(id='S11',references=['S09','S10'],base_peaks=5),
           dict(id='E05',references=['E03','E04'],base_peaks=7)]


def save(path,row):
    path.write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8',newline='\n')


def fixed_root(row,beta,guess,digits=70):
    f=PrecisionFlights(row['regions'],row['targets'],row.get('symmetric',True),digits)
    b=mp.mpf(str(beta));q=mp.matrix([mp.mpf(str(v)) for v in guess]);history=[]
    for iteration in range(30):
        ev=f.evaluate(b,list(q),False);error=max(abs(v) for v in ev['residual'])
        history.append(mp.nstr(error,15))
        if error<mp.mpf(10)**(-(digits-20)):break
        change=mp.lu_solve(ev['jacobian'],-ev['residual'])
        for scale in [1,mp.mpf('.5'),mp.mpf('.25'),mp.mpf('.125'),mp.mpf('.0625')]:
            trial=q+scale*change
            if min(trial[2:])>0 and max(abs(v) for v in f.evaluate(b,list(trial),False)['residual'])<error:
                q=trial;break
        else:raise RuntimeError('Fixed-beta line search stalled: '+str(history))
    else:raise RuntimeError('Fixed-beta Newton exhausted: '+str(history))
    tr=ev['trace'];det=ev['determinant']
    large=(tr+mp.sign(tr)*mp.sqrt(tr*tr-4*det))/2
    result=dict(beta=mp.nstr(b,digits-5),q=[mp.nstr(v,digits-5) for v in q],
        period=mp.nstr((2 if f.sign==-1 else 1)*sum(q[2:]),digits-5),
        signed_map_multiplier=mp.nstr(large,30),full_multiplier=mp.nstr(large**(2 if f.sign==-1 else 1),30),
        residual=mp.nstr(error,15),newton_residuals=history,working_digits=digits,
        regions=row['regions'],targets=row['targets'],symmetric=row.get('symmetric',True))
    return result,ev


def parent_at(id,beta):
    row=read(SOURCES/f'{id}-d70-n192.json')
    f=FastFlights(Chua(),row['regions'],row['targets'],True)
    q=np.array(row['q'],float)
    for b in np.linspace(float(row['beta']),float(beta),81)[1:]:q=f.solve(float(b),q)
    check=f.record(float(beta),q,True)
    assert check['region_violation']<1e-8 and check['full_multipliers'][0][0]>1
    result,ev=fixed_root(row,beta,q)
    result['sampled_path_check']={k:check[k] for k in ['region_violation','minimum_breakpoint_speed']}
    return result


def calibrate():
    result=[]
    for spec in SEQUENCES:
        refs=[]
        for id in spec['references']:
            source=SOURCES/f'{id}-d70-n192.json';row=read(source)
            f=PrecisionFlights(row['regions'],row['targets'],True,70)
            fold=f.boundary(row['beta'],row['q'],1)
            symmetry=f.boundary(row['beta'],row['q'],-1)
            width=abs(mp.mpf(fold['beta'])-mp.mpf(symmetry['beta']))
            gain_source=ROOT/f'results/chua-return-geometry/{id}-optimized.json'
            g=read(gain_source)['maximum']['gain']
            refs.append(dict(id=id,source=source.relative_to(ROOT).as_posix(),source_sha256=digest(source),
                gain_source=gain_source.relative_to(ROOT).as_posix(),gain_source_sha256=digest(gain_source),
                representative=row,fold=fold,symmetry=symmetry,width=mp.nstr(width,55),gain=g))
        parent=parent_at(spec['id'],refs[-1]['representative']['beta'])
        result.append(dict(**spec,calibration=refs,parent=parent))
        print(spec['id'],'widths',[r['width'] for r in refs],'gains',[r['gain'] for r in refs],
              'expansion',parent['full_multiplier'],'period',parent['period'],flush=True)
    save(OUT/'calibration.json',dict(sequences=result,scope='Previously known long roots and proposed short saddles; no new long roots solved'))


def freeze():
    path=OUT/'forecasts-frozen.json'
    if path.exists():raise SystemExit('Frozen forecasts already exist')
    calibration=read(OUT/'calibration.json');targets=[]
    for seq in calibration['sequences']:
        a,b=seq['calibration'];lam=float(seq['parent']['full_multiplier']);tp=float(seq['parent']['period'])
        w1,w2=float(a['width']),float(b['width']);g1,g2=a['gain'],b['gain']
        bw=(w1/lam**2-w2)/(lam-1);aw=w2-bw
        ag=(g2-g1)/(1-1/lam);bg=g2-ag
        for n in [3,4,5]:
            k=n-2;t=lam**-k
            targets.append(dict(id=f"{seq['id']}-n{n}",organizer=seq['id'],n=n,
                peaks=seq['base_peaks']*(2*n+1),section_returns=2*n+1,
                expansion=lam,parent_period=tp,
                leading=dict(width=w2*t*t,gain=g2/t),
                corrected=dict(width=t*t*(aw+bw*t),gain=ag/t+bg),
                period=float(b['representative']['period'])+2*k*tp,
                seed_beta=float(b['fold']['beta'])+(float(b['fold']['beta'])-float(a['fold']['beta']))*sum(lam**-j for j in range(1,k+1)),
                coefficients=dict(A_width=aw,B_width=bw,A_gain=ag,B_gain=bg),
                gates=dict(width_relative=.25,gain_relative=.20,period_absolute=.5),
                companion_half_multiplier=1.5))
    save(path,dict(created_utc=datetime.now(timezone.utc).isoformat(),targets=targets,
        source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in [OUT/'protocol.md',OUT/'calibration.json']},
        scope='Forecasts for new long orbits, using two existing members per sequence. Locations are corrected by shooting, not scored as predictions.'))
    for t in targets:print(t['id'],t['peaks'],t['corrected'],t['period'],flush=True)


def insert_parent(root,parent):
    """Insert a full short-saddle orbit at the nearest positive x=0 crossing."""
    beta=float(root['beta']);q=np.array(root['q'],float);regions=root['regions']
    hreg=parent['regions'];hd=np.array(parent['q'][2:],float)
    preg=hreg[:-1]+[0]+[-r for r in hreg[1:]]
    pd=np.r_[hd[:-1],hd[-1]+hd[0],hd[1:]]
    ptarg=[preg[j] or preg[j+1] for j in range(len(preg)-1)]
    pf=FastFlights(Chua(),preg,ptarg,False)
    pq=pf.solve(beta,np.r_[np.array(parent['q'][:2],float),pd])
    model=Chua(rtol=2e-12,step=.02)
    sol=model.integrate(beta,[0.,*q[:2]],float(root['period'])/2)
    ts=sol.t_events[2];us=sol.y_events[2]
    inds=np.flatnonzero((ts>1e-6)&(ts<float(root['period'])/2-1e-6))
    i=inds[np.argmin(np.linalg.norm(us[inds,1:]-pq[:2],axis=1))]
    t=ts[i];ends=np.cumsum(q[2:]);starts=np.r_[0.,ends[:-1]];k=int(np.searchsorted(ends,t))
    assert regions[k]==0 and preg[0]==preg[-1]==0
    durations=np.r_[q[2:2+k],t-starts[k]+pq[2],pq[3:-1],pq[-1]+ends[k]-t,q[3+k:]]
    newregions=regions[:k]+preg+regions[k+1:]
    targets=[newregions[j] or newregions[j+1] for j in range(len(newregions)-1)]
    return newregions,targets,np.r_[q[:2],durations],float(t)


def search():
    forecast=read(OUT/'forecasts-frozen.json');calibration=read(OUT/'calibration.json')
    for seq in calibration['sequences']:
        prior=dict(seq['calibration'][-1]['fold'],regions=seq['calibration'][-1]['representative']['regions'])
        for target in [t for t in forecast['targets'] if t['organizer']==seq['id']]:
            id=target['id'];path=OUT/f'{id}-reference.json'
            if path.exists():
                prior=dict(read(OUT/f'{id}-edges.json')['fold'],regions=read(path)['regions']);continue
            started=time.monotonic()
            try:
                parent=parent_at(seq['id'],prior['beta'])
                regions,targets,guess,t=insert_parent(prior,parent)
                f=PrecisionFlights(regions,targets,True,70)
                # Predictor retries affect the search only; the published
                # width, gain and period forecasts remain untouched.
                attempts=[];fold=None
                fast=FastFlights(Chua(),regions,targets,True)
                for scale in [1.,.5,1.5,0.,2.]:
                    seed=float(prior['beta'])+scale*(target['seed_beta']-float(prior['beta']))
                    try:
                        b,q=fast.boundary(seed,guess,1)
                        trial=f.boundary(str(b),q,1,max_iterations=40)
                        attempts.append(dict(seed_scale=scale,status='converged',beta=trial['beta']))
                        fold=trial;break
                    except Exception as exc:
                        attempts.append(dict(seed_scale=scale,status='failed',error=str(exc)))
                save(OUT/f'{id}-search-retries.json',dict(attempts=attempts))
                if fold is None:raise RuntimeError('All predictor retries failed')
                symmetry=f.boundary(fold['beta'],fold['q'],-1,max_iterations=40)
                rep=f.boundary(fold['beta'],fold['q'],'.5',max_iterations=40)
                rep.update(id=id,symmetric=True,regions=regions,targets=targets,organizer=seq['id'],
                           n=target['n'],forecast_peaks=target['peaks'])
                width=abs(mp.mpf(fold['beta'])-mp.mpf(symmetry['beta']))
                save(path,rep)
                save(OUT/f'{id}-edges.json',dict(id=id,fold=fold,symmetry=symmetry,width=mp.nstr(width,60),
                    insertion_time=t,forecast_sha256=digest(OUT/'forecasts-frozen.json'),elapsed_seconds=time.monotonic()-started))
                prior=dict(fold,regions=regions)
                print(id,'width',mp.nstr(width,12),'T',rep['period'],'seconds',time.monotonic()-started,flush=True)
            except Exception as exc:
                failure=OUT/f'{id}-search-failure.json'
                if not failure.exists():save(failure,dict(id=id,error=str(exc),elapsed_seconds=time.monotonic()-started))
                print('FAILED',id,str(exc),flush=True);break


def refined_gain(row,data,rough):
    arcs,period=arcs_for(row);t=data['t'];i=rough['maximizing_start_index'];k=rough['maximizing_steps'];n=len(t)
    tt=np.r_[t,t+period,2*period];start=tt[i];end=tt[i+k]
    ds=2*max(tt[i+1]-tt[i],tt[i]-tt[i-1] if i else tt[n]-tt[n-1])
    de=2*max(tt[i+k+1]-tt[i+k],tt[i+k]-tt[i+k-1]) if i+k+1<len(tt) else ds
    bounds=[(max(0,start-ds),min(period,start+ds)),(end-de,end+de)]
    def objective(x):
        return -np.log(transport(arcs,period,*x)) if 0<x[1]-x[0]<=period else 1e3
    sol=minimize(objective,[start,end],method='Powell',bounds=bounds,options={'xtol':1e-10,'ftol':1e-12,'maxiter':80})
    return dict(success=bool(sol.success),message=str(sol.message),maximum=transport(arcs,period,*sol.x,True),
                bounds=bounds,scope='Local refinement of a sampled gain maximum, not a certified global supremum')


def gains(subdivisions,ids):
    for path in sorted(OUT.glob('*-reference.json')):
        row=read(path);id=row['id']
        if ids and id not in ids:continue
        target=OUT/f'{id}-gain-s{subdivisions}.json'
        if target.exists():continue
        data,diag=mesh(row,subdivisions,85);curves,rough=amplification(data)
        record=dict(id=id,source=path.name,source_sha256=digest(path),diagnostics=diag,**rough)
        record['refined']=refined_gain(row,data,rough) if subdivisions==8 else None
        record['passed']=diag['closure']<1e-30 and diag['multiplier_error']<1e-20
        np.savez_compressed(target.with_suffix('.npz'),**data,**curves)
        save(target,record)
        print(id,'mesh',subdivisions,'gain',rough['maximum_gain'],'refined',record['refined']['maximum']['gain'] if record['refined'] else None,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['calibrate','freeze','search','gains'])
    p.add_argument('--subdivisions',type=int,default=8)
    p.add_argument('--ids',nargs='*')
    args=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if args.stage=='calibrate':calibrate()
    elif args.stage=='freeze':freeze()
    elif args.stage=='search':search()
    else:gains(args.subdivisions,args.ids)


if __name__=='__main__':main()
