"""Freeze quadratic-map forecasts, then search nested long-block Chua cycles."""
import argparse
from datetime import datetime,timezone
import time
import mpmath as mp
import numpy as np
from scipy.optimize import brentq
from study_chua_return_geometry import ROOT,read,digest
from test_chua_return_families import save,fixed_root
from chua_precision_flights import PrecisionFlights
from chua_fast_flights import FastFlights
from trace_chua_window import Chua
from chua_block_return import BlockShooting

OUT=ROOT/'results/chua-nested-returns'
PARENTS=[dict(id='S06',source='results/chua-knot-atlas/S06-d70-n192.json',
              companion_source='results/chua-return-geometry/S06-companion-d70.json'),
         dict(id='S11-n3',source='results/chua-return-families/S11-n3-reference.json',
              companion_source='results/chua-return-families/S11-n3-companion.json')]


def quadratic(k):
    def critical(M):
        x=0.
        for _ in range(k):x=M-x*x
        return x
    grid=np.linspace(1.5,2,20001);v=critical(grid);brackets=[(a,b) for a,b,u,vv in zip(grid[:-1],grid[1:],v[:-1],v[1:]) if u*vv<0]
    root=brentq(critical,*brackets[-1],xtol=1e-15)
    mp.mp.dps=70;root=mp.findroot(lambda M:iterate(M,mp.mpf(0),k)[0],mp.mpf(str(root)))
    def target(mu,seed):return mp.findroot(lambda x,M:(iterate(M,x,k)[0]-x,iterate(M,x,k)[1]-mu),list(seed),tol=mp.mpf('1e-60'))
    rep=target(mp.mpf('.5'),(mp.mpf(0),root));fold=target(mp.mpf(1),rep);flip=target(mp.mpf(-1),rep)
    orbit=[];x=rep[0]
    for _ in range(k):orbit.append(mp.nstr(x,60));x=rep[1]-x*x
    return dict(k=k,superstable=mp.nstr(root,60),representative=mp.nstr(rep[1],60),
        fold=mp.nstr(fold[1],60),flip=mp.nstr(flip[1],60),orbit=orbit,
        width=mp.nstr(abs(flip[1]-fold[1]),60))


def iterate(M,x,k):
    derivative=1
    for _ in range(k):derivative*=-2*x;x=M-x*x
    return x,derivative


def continued_parent(row,old,beta):
    fast=FastFlights(Chua(),row['regions'],row['targets'],True);q=np.array(old['q'],float)
    for b in np.linspace(float(old['beta']),float(beta),61)[1:]:q=fast.solve(float(b),q)
    r,e=fixed_root(row,beta,q,70);check=fast.record(float(beta),np.array(r['q'],float),True)
    r['path_check']={k:check[k] for k in ['region_violation','minimum_breakpoint_speed']}
    assert check['region_violation']<1e-8
    return r


def freeze():
    frozen=OUT/'forecasts-frozen.json'
    if frozen.exists():raise SystemExit('Refusing to overwrite frozen forecasts')
    qs=[quadratic(k) for k in [3,5,7]];parents=[];targets=[]
    for spec in PARENTS:
        row=read(ROOT/spec['source']);comp=read(ROOT/spec['companion_source'])['companion']
        f=PrecisionFlights(row['regions'],row['targets'],True,70)
        fold=f.boundary(row['beta'],row['q'],1);sym=f.boundary(row['beta'],row['q'],-1)
        width=mp.mpf(fold['beta'])-mp.mpf(sym['beta']);assert width>0
        cal=dict(**spec,representative=row,companion=comp,fold=fold,symmetry=sym,width=mp.nstr(width,60),
            source_sha256=digest(ROOT/spec['source']),companion_source_sha256=digest(ROOT/spec['companion_source']))
        parents.append(cal)
        for model in qs:
            k=model['k'];beta=mp.mpf(fold['beta'])-width*(mp.mpf(model['representative'])+mp.mpf('.25'))
            saddle=continued_parent(row,comp,mp.nstr(beta,65))
            t=dict(id=f"{spec['id']}-k{k}",parent=spec['id'],block_count=k,peaks=35*k,section_returns=7*k,
                quadratic=model,beta=mp.nstr(beta,60),width=mp.nstr(width*mp.mpf(model['width']),60),
                period=mp.nstr(k*mp.mpf(saddle['period']),60),parent_at_predicted_beta=saddle,
                gates=dict(location_in_widths=1.,width_relative=.3,period_absolute=.5))
            targets.append(t);print(t['id'],'beta',t['beta'],'width',t['width'],'period',t['period'],flush=True)
    save(OUT/'calibration.json',dict(parents=parents))
    save(frozen,dict(created_utc=datetime.now(timezone.utc).isoformat(),targets=targets,
        source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in [OUT/'protocol.md',OUT/'calibration.json']},
        scope='Quadratic-map forecasts and continued long parents only; no child Chua orbit searched'))


def search(digits,ids):
    forecasts=read(OUT/'forecasts-frozen.json');cal=read(OUT/'calibration.json')
    for t in forecasts['targets']:
        id=t['id'];path=OUT/f'{id}-reference.json'
        if ids and id not in ids:continue
        if path.exists() and (OUT/f'{id}-edges.json').exists():continue
        parent=next(c for c in cal['parents'] if c['id']==t['parent']);shoot=BlockShooting(parent,t['block_count'],digits)
        guess=[]
        for X in t['quadratic']['orbit']:guess.extend([-.5-2*mp.mpf(X),0])
        guess.append(mp.mpf(t['quadratic']['representative']));start=time.monotonic()
        try:
            if path.exists():rep=read(path)
            else:
                rep=shoot.solve(guess,'.5');rep.update(id=id,parent=t['parent'],forecast_peaks=t['peaks'],forecast_sha256=digest(OUT/'forecasts-frozen.json'))
                save(path,rep)
            print(id,'representative',rep['beta'],'T',rep['period'],'seconds',time.monotonic()-start,flush=True)
            w=[mp.mpf(v) for p in rep['block_coordinates'] for v in p]+[mp.mpf(rep['normalized_coordinate'])]
            fold=shoot.solve(w,1);flip=shoot.solve(w,-1)
            width=abs(mp.mpf(fold['beta'])-mp.mpf(flip['beta']))
            save(OUT/f'{id}-edges.json',dict(id=id,fold=fold,symmetry=flip,width=mp.nstr(width,digits-5),
                forecast_sha256=digest(OUT/'forecasts-frozen.json'),elapsed_seconds=time.monotonic()-start))
            print(id,'width',mp.nstr(width,12),'seconds',time.monotonic()-start,flush=True)
        except Exception as exc:
            fail=OUT/f'{id}-search-failure-d{digits}.json'
            if not fail.exists():save(fail,dict(id=id,error=str(exc),elapsed_seconds=time.monotonic()-start))
            print('FAILED',id,str(exc),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['freeze','search']);p.add_argument('--digits',type=int,default=60);p.add_argument('--ids',nargs='*');args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.stage=='freeze':freeze()
    else:search(args.digits,args.ids)


if __name__=='__main__':main()
