"""Interval certificates for pulse-line intersections with the local basin boundary."""
import argparse
import json
import time
from decimal import Decimal
from flint import arb,arb_mat,ctx
from chua_research_followup import ROOT,OUT,dump
from certify_chua_basin_separator import (refine_root,local_block,DynamicReturn,
    col,eye,matrix,midpoint,origin_ellipsoid,capture_cycle,captured,propagate,
    continuous_origin_capture)
from chua_interval_map import bounds

DIR=OUT/'trigger'


def context(digits):
    ctx.dps=digits
    seeds=json.loads((ROOT/'results/chua-hidden-boundary/cycle-seeds.json').read_text())['full_records']
    f,p,pbox,J=refine_root(seeds[0],digits)
    V,W,d,local,umax,smax=local_block(f,p,pbox,J)
    return f,p,pbox,V,W,d,local


def iterate_ray(amplitude,direction,count,digits,pbox,W):
    state=amplitude*direction;derivative=direction;cells=0;time=arb(0)
    for _ in range(count):
        f=DynamicReturn(state.entries(),digits)
        state,J,arcs=f.evaluate(state.entries(),True)
        derivative=J*derivative
        cells+=sum(a['path']['accepted_cells'] for a in arcs)
        time+=sum(a['duration'] for a in arcs)
    return W*(state-pbox),W*derivative,cells,time


def local(digits):
    f,p,pbox,V,W,d,cert=context(digits);L=arb('.02');rows=[]
    numeric=json.loads((DIR/'thresholds.json').read_text())['records']
    for row in numeric:
        direction=col([arb(str(v)) for v in row['direction']])
        center=arb(str(row['approximants'][-1]['amplitude']))
        # The wider interval proves a unique local crossing; its two inner
        # endpoints locate the crossing much more tightly.
        broad=arb(center,arb('2e-10'));left=center-arb('1e-11');right=center+arb('1e-11')
        entire,D,cells,T=iterate_ray(broad,direction,6,digits,pbox,W)
        ql,_,lc,_=iterate_ray(left,direction,6,digits,pbox,W)
        qr,_,rc,_=iterate_ray(right,direction,6,digits,pbox,W)
        contained=all(abs(x)<d for x in entire.entries())
        transverse=D[0,0]-L*abs(D[1,0])>0
        negative=ql[0,0]+L*abs(ql[1,0])<0
        positive=qr[0,0]-L*abs(qr[1,0])>0
        passed=bool(contained and transverse and negative and positive)
        rec=dict(angle=row['angle'],direction=[str(v) for v in row['direction']],
            half_returns=6,bracket=[left.mid().str(digits,radius=False),right.mid().str(digits,radius=False)],
            broad_interval=bounds(broad),image=[bounds(v) for v in entire.entries()],
            derivative=[bounds(v) for v in D.entries()],left_image=[bounds(v) for v in ql.entries()],
            right_image=[bounds(v) for v in qr.entries()],path_cells=cells+lc+rc,
            gates=dict(local_domain=bool(contained),transverse=bool(transverse),negative_side=bool(negative),positive_side=bool(positive)),passed=passed)
        rows.append(rec);print('pulse intersection',row['angle'],passed,rec['gates'],flush=True)
    dump(DIR/f'local-certificate-{digits}.json',dict(passed=all(r['passed'] for r in rows),working_digits=digits,
        local_block=cert,records=rows,
        scope='Unique intersection on each displayed local amplitude interval, with opposite basin destinations. Does not exclude earlier or later crossings along the ray.'))
    assert all(r['passed'] for r in rows)


def central_segment(digits,amplitude='.299',duration='30'):
    f,p,pbox,V,W,d,cert=context(digits)
    A,_,_=f.matrices(None,0);P,pc=origin_ellipsoid(A);level=arb('.99')/P.inv()[0,0]
    initial=col([0,arb(amplitude),0]);T=arb(duration);maximum=arb(0)
    pending=[(T*k/1200,T*(k+1)/1200,0) for k in range(1200)];cells=0
    while pending:
        left,right,depth=pending.pop();t=left.union(right)
        x=((A*t).exp()*initial)[0,0]
        if abs(x)<1:
            maximum=maximum.max(abs(x));cells+=1
        else:
            assert depth<18,'Central-flight enclosure unresolved'
            mid=(left+right)/2;pending.extend([(left,mid,depth+1),(mid,right,depth+1)])
    end=(A*T).exp()*initial;energy=(end.transpose()*P*end)[0,0]
    assert energy<level
    return dict(passed=True,upper_amplitude=amplitude,duration=duration,time_cells=cells,
        maximum_absolute_x=bounds(maximum),final_energy=bounds(energy),capture_level=bounds(level),
        argument='All smaller nonnegative amplitudes scale this central trajectory. At T they enter the invariant origin ellipsoid.')


def below(digits,start_band=0,end_band=None):
    f,p,pbox,V,W,d,cert=context(digits)
    seeds=json.loads((ROOT/'results/chua-hidden-boundary/cycle-seeds.json').read_text())['full_records']
    fc,pc,cbox,Jc=refine_root(seeds[1],digits);target,_=capture_cycle(fc,pc,Jc)
    A,_,_=f.matrices(None,0);P,pcert=origin_ellipsoid(A);level=arb('.99')/P.inv()[0,0]
    localdata=json.loads((DIR/f'local-certificate-{digits}.json').read_text())
    row=next(r for r in localdata['records'] if r['angle']==0)
    stop=arb(row['bracket'][0]);start=arb('.299');central=central_segment(digits)
    # Geometric intervals keep their relative distance to the boundary bounded.
    center=arb(str(json.loads((DIR/'thresholds.json').read_text())['records'][3]['approximants'][-1]['amplitude']))
    edges=[start]
    while center-edges[-1]>arb('2e-11'):
        edges.append(center-(center-edges[-1])*arb('.8'))
    edges.append(stop);edges=sorted(set(x.mid().str(80,radius=False) for x in edges),key=lambda x:float(x))
    # Remove any tiny geometric overshoot past the prescribed stopping point.
    edges=[arb(x) for x in edges if arb(x)<=stop]
    if edges[-1]!=stop:edges.append(stop)
    rows=[];failures=[]
    chunk=end_band is not None or start_band!=0
    output=DIR/(f'below-part-{digits}-{start_band}-{end_band}.json' if chunk else f'below-certificate-{digits}.json')
    def capture(a,b,depth=0):
        interval=a.union(b);centre=col([interval.mid(),0]);basis=eye(2);radii=[interval.rad(),arb(0)]
        count=0;cells=0;elapsed=0.;label=None
        try:
            for k in range(100):
                box=centre+basis*col([arb(0,r) for r in radii])
                localcoords=W*(box-pbox)
                if all(abs(v)<d for v in localcoords.entries()) and localcoords[0,0]+arb('.02')*abs(localcoords[1,0])<0:
                    label='origin_by_separator';break
                label=captured(centre,basis,radii,target,P,level)
                if label:break
                try:centre,basis,radii,c,dt=propagate(centre,basis,radii,digits)
                except (ValueError,AssertionError):
                    flow=continuous_origin_capture(centre,basis,radii,P,level)
                    label='origin_by_continuous_flow';cells+=flow.get('steps',0);break
                cells+=c;elapsed+=dt;count+=1
            if label not in ['origin','origin_by_separator','origin_by_continuous_flow']:
                raise ValueError('Not captured at origin: '+str(label))
            rows.append(dict(lower=a.mid().str(80,radius=False),upper=b.mid().str(80,radius=False),
                label=label,returns=count,path_cells=cells,time=elapsed,depth=depth))
        except (AssertionError,ValueError,ZeroDivisionError) as exc:
            failures.append(dict(lower=a.mid().str(25,radius=False),upper=b.mid().str(25,radius=False),depth=depth,error=str(exc)))
            if depth>=12:raise
            mid=(a+b)/2;capture(a,mid,depth+1);capture(mid,b,depth+1)
    for k,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        if k<start_band or (end_band is not None and k>=end_band):continue
        capture(a,b)
        dump(output,dict(passed=False,partial=True,central=central,records=rows,rejected_enclosures=failures,
            first_band=start_band,last_completed_band=k,total_bands=len(edges)-1,working_digits=digits))
        print('below',k,'/',len(edges)-1,'leaves',len(rows),'A',float(b.mid()),flush=True)
    if chunk:
        data=json.loads(output.read_text());data['chunk_complete']=True;dump(output,data);return
    ordered=sorted(rows,key=lambda r:float(r['lower']))
    # These are serialized labels of shared endpoints. Arb equality would test
    # equality of uncertain balls, not continuity of the subdivision tree.
    assert Decimal(ordered[0]['lower'])==Decimal(start.mid().str(80,radius=False))
    assert Decimal(ordered[-1]['upper'])==Decimal(stop.mid().str(80,radius=False))
    assert all(Decimal(a['upper'])==Decimal(b['lower']) for a,b in zip(ordered[:-1],ordered[1:]))
    dump(DIR/f'below-certificate-{digits}.json',dict(passed=True,partial=False,working_digits=digits,
        central=central,records=ordered,rejected_enclosures=failures,
        coverage=dict(adjoining_intervals=len(ordered),first=ordered[0]['lower'],last=ordered[-1]['upper']),
        threshold_bracket=row['bracket'],scope='Every amplitude from zero to the lower bracket endpoint converges to the origin; the upper endpoint converges to the hidden cycle by the local certificate.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['local','below']);ap.add_argument('--digits',type=int,default=80);ap.add_argument('--start-band',type=int,default=0);ap.add_argument('--end-band',type=int);a=ap.parse_args()
    if a.phase=='local':local(a.digits)
    else:below(a.digits,a.start_band,a.end_band)
