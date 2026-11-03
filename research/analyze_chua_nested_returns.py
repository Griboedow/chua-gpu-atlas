"""Test long-block ancestry and score the original nested-cycle forecasts."""
import argparse
import mpmath as mp
import numpy as np
from test_chua_nested_returns import OUT,ROOT,read,save,digest,continued_parent
from test_chua_return_families import parent_at
from analyze_chua_return_geometry import sections
from chua_block_return import BlockReturn,BlockShooting
from chua_precision_flights import PrecisionFlights


def side_coordinate(short,points):
    f=BlockReturn(short,65);ev=f.evaluate(mp.mpf(short['beta']),[mp.mpf(v) for v in short['q'][:2]])
    eig,V=np.linalg.eig(np.array(ev['matrix'].tolist(),float).T)
    left=V[:,np.argmax(abs(eig))].real;left*=np.sign(left[0]);left/=np.linalg.norm(left)
    u=(points-np.array(short['q'][:2],float))@left
    return u,left


def word(u):return ''.join('+' if v>0 else '-' for v in u)


def primitive_word_length(s):
    return next(k for k in range(1,len(s)+1) if len(s)%k==0 and s==s[:k]*(len(s)//k))


def geometry(row,cal):
    id=row['id'];path=OUT/f'{id}-geometry.json'
    if path.exists():return read(path)
    companion=continued_parent(cal['representative'],cal['companion'],row['beta'])
    other=continued_parent(cal['representative'],cal['representative'],row['beta'])
    short=parent_at('S08' if row['parent']=='S06' else 'S11',row['beta'])
    p,t=sections(row,row);pc,tc=sections(companion,companion);po,to=sections(other,other)
    assert len(p)==7*row['block_count'] and len(pc)==len(po)==7
    scale=np.sqrt(np.mean(np.sum((pc-po)**2,axis=1)))
    distance=[float(np.sqrt(np.mean(np.sum((p[j:j+7]-pc)**2,axis=1)))/scale) for j in range(0,len(p),7)]
    other_distance=[float(np.sqrt(np.mean(np.sum((p[j:j+7]-po)**2,axis=1)))/scale) for j in range(0,len(p),7)]
    delta=po[0]-pc[0];block_coordinate=(p[::7]-pc[0])@delta/np.dot(delta,delta)
    u,left=side_coordinate(short,p);parent_u,_=side_coordinate(short,pc)
    symbols=word(u);parent_symbols=word(parent_u)
    signs=np.sign(u);same_side=int(np.sum(signs==np.roll(signs,1)))
    result=dict(id=id,parent=companion,continued_attracting_branch=other,short_parent=short,
        points=p.tolist(),times=t.tolist(),parent_points=pc.tolist(),other_parent_points=po.tolist(),
        block_rms_distance_to_parent=distance,block_rms_distance_to_other_parent=other_distance,
        normalization_parent_rms_separation=float(scale),block_start_coordinate=block_coordinate.tolist(),
        blocks_with_distance_below_one_tenth=sum(d<.1 for d in distance),
        closeness_threshold_scope='Descriptive display threshold, not a preregistered acceptance gate',
        short_saddle_coordinate=u.tolist(),parent_short_saddle_coordinate=parent_u.tolist(),left_covector=left.tolist(),
        side_word=symbols,parent_side_word=parent_symbols,side_word_minimal_period=primitive_word_length(symbols),
        exact_parent_word_repetition=symbols==parent_symbols*row['block_count'],
        same_side_transitions=same_side,minimum_side_coordinate=float(min(abs(u))),
        period_minus_parent_multiple=float(row['period'])-row['block_count']*float(companion['period']),
        source_sha256=digest(OUT/f'{id}-reference.json'))
    save(path,result);print(id,'parent word',parent_symbols,'child',symbols,'close blocks',result['blocks_with_distance_below_one_tenth'],'min distance',min(distance),flush=True)
    return result


def controls(calibration):
    path=OUT/'negative-controls.json'
    if path.exists():return
    source=ROOT/'results/chua-return-geometry/P105-reference.json';ordinary=read(source)
    short=parent_at('S08',ordinary['beta']);p,t=sections(ordinary,ordinary);u,left=side_coordinate(short,p);s=word(u)
    result=dict(ordinary_105=dict(source=source.relative_to(ROOT).as_posix(),source_sha256=digest(source),beta=ordinary['beta'],
        points=p.tolist(),times=t.tolist(),short_saddle_coordinate=u.tolist(),side_word=s,
        side_word_minimal_period=primitive_word_length(s),same_side_transitions=int(np.sum(np.sign(u)==np.roll(np.sign(u),1))),
        scope='Each orbit is compared with its own same-parameter short saddle; the binary word alone is not an ancestry proof'),repeated_parent=[])
    for cal in calibration['parents']:
        row=cal['representative'];f=PrecisionFlights(row['regions'],row['targets'],True,70);mu=mp.mpf('.5')**(mp.mpf(1)/3)
        parent=f.boundary(row['beta'],row['q'],mu);parent.update(regions=row['regions'],targets=row['targets'],symmetric=True)
        points,times=sections(parent,parent)
        shoot=BlockShooting(cal,3,65);q=mp.matrix([mp.mpf(v) for v in parent['q'][:2]])
        a=shoot.W*(q-shoot.origin);M=(shoot.beta_fold-mp.mpf(parent['beta']))/shoot.width-mp.mpf('.25')
        w=mp.matrix([*a,*a,*a,M]);res,P,det=shoot.residual(w,mp.mpf('.5'))
        tiled=np.tile(points,(3,1));d=np.linalg.norm(tiled[:,None,:]-tiled[None,:,:],axis=2);np.fill_diagonal(d,np.inf)
        record=dict(parent=cal['id'],parent_root=parent,apparent_peak_count=105,apparent_period=3*float(parent['period']),
            block_closure_residual=mp.nstr(max(abs(v) for v in res),20),
            apparent_signed_half_multiplier=mp.nstr(P[0,0]+P[1,1],20),minimum_distinct_section_distance=float(d.min()),
            rejected_as_nonprimitive=bool(d.min()<=1e-12),
            scope='Three traversals of a 35-peak parent adjusted so even the apparent half multiplier equals 1/2')
        assert record['rejected_as_nonprimitive'];result['repeated_parent'].append(record)
    save(path,result);print('Negative control: ordinary 105 word',s,'period',primitive_word_length(s),flush=True)


def branch_continuation(calibration):
    path=OUT/'parent-branch-continuation.json'
    if path.exists():return
    rows=[]
    for cal in calibration['parents']:
        source=cal['representative'];f=PrecisionFlights(source['regions'],source['targets'],True,60);old=source;roots=[]
        for target in ['.5','.75','1','1.25','1.5','2','3','4']:
            r=f.boundary(old['beta'],old['q'],target);roots.append(r);old=r
        rows.append(dict(parent=cal['id'],roots=roots))
        print(cal['id'],'branch continued through the fold to half multiplier 4',flush=True)
    save(path,dict(rows=rows,scope='Numerical continuation of one 35-peak periodic-orbit family through its saddle-node'))


def scores():
    forecast=read(OUT/'forecasts-frozen.json');rows=[]
    for t in forecast['targets']:
        id=t['id'];r=read(OUT/f'{id}-reference.json');e=read(OUT/f'{id}-edges.json');mp.mp.dps=70
        width=mp.mpf(e['width']);delta=mp.mpf(r['beta'])-mp.mpf(t['beta'])
        errors=dict(location_absolute=mp.nstr(delta,25),location_in_widths=float(delta/width),
            width_relative=float(width/mp.mpf(t['width'])-1),period_absolute=float(mp.mpf(r['period'])-mp.mpf(t['period'])))
        passed=dict(location=abs(errors['location_in_widths'])<=t['gates']['location_in_widths'],
            width=abs(errors['width_relative'])<=t['gates']['width_relative'],period=abs(errors['period_absolute'])<=t['gates']['period_absolute'])
        rows.append(dict(id=id,parent=t['parent'],k=t['block_count'],peaks=t['peaks'],beta=r['beta'],width=e['width'],period=r['period'],errors=errors,passed=passed))
    save(OUT/'scores.json',dict(forecast_sha256=digest(OUT/'forecasts-frozen.json'),rows=rows,
        passed_by_quantity={key:sum(r['passed'][key] for r in rows) for key in ['location','width','period']}))
    print('Forecast passes',{key:sum(r['passed'][key] for r in rows) for key in ['location','width','period']},flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--branches',action='store_true');args=parser.parse_args()
    cal=read(OUT/'calibration.json')
    if args.branches:branch_continuation(cal);return
    scores();rows=[]
    for path in sorted(OUT.glob('*-reference.json')):
        r=read(path);parent=next(c for c in cal['parents'] if c['id']==r['parent']);rows.append(geometry(r,parent))
    save(OUT/'geometry-summary.json',dict(rows=rows));controls(cal)


if __name__=='__main__':main()
