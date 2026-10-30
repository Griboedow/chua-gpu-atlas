"""Continue identified symmetric window cores through the original atlas box."""
import argparse,json,time
from collections import deque
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_fast_flights import FastFlights

ROOT=Path(__file__).resolve().parents[1]


def registry():
    old=ROOT/'results/chua-layer-census'
    raw=json.loads((old/'family-edges.json').read_text())
    unique=[r for r in raw if abs(r['beta']-14.259953125)>1e-9]
    second=json.loads((old/'second-boundaries.json').read_text())
    unique.append(dict(beta=14.354,peak_count=5,fold=second['fold'],symmetry=second['symmetry'],
        regions=[a['region'] for a in second['fold']['arcs']],targets=[round(a['end'][0]) for a in second['fold']['arcs'][:-1]]))
    unique.sort(key=lambda r:r['beta'])
    for k,r in enumerate(unique):r['id']=f'S{k+1:02d}'
    return unique


def advance(template,source,alpha,m0,depth=0):
    row=dict(alpha=float(alpha),m0=float(m0));last=[]
    try:
        f=FastFlights(Chua(alpha=alpha,m0=m0),template['regions'],template['targets'],True)
        # Shift with measured family-specific slopes when supplied; a close
        # predictor is particularly important for the high-period windows.
        for key,target in [('fold',1),('symmetry',-1)]:
            previous=source[key];q=np.array(previous['q'])
            old=FastFlights(Chua(alpha=source['alpha'],m0=source['m0']),template['regions'],template['targets'],True)
            w=np.r_[q,previous['beta']];forcing=np.zeros(len(w))
            for param,delta,h in [('alpha',alpha-source['alpha'],1e-6),('m0',m0-source['m0'],1e-7)]:
                if delta==0:continue
                values=[]
                for sign in [-1,1]:
                    kwargs=dict(alpha=source['alpha'],m0=source['m0']);kwargs[param]+=sign*h
                    trial=FastFlights(Chua(**kwargs),template['regions'],template['targets'],True)
                    values.append(trial.boundary_system(w,target))
                forcing+=delta*(values[1]-values[0])/(2*h)
            prediction=w-np.linalg.solve(old.boundary_system(w,target,True),forcing)
            b,q=f.boundary(prediction[-1],prediction[:-1],target);r=f.record(b,q,True)
            if r['region_violation']>1e-7:raise RuntimeError('Itinerary violation')
            if abs(r['period']/previous['period']-1)>.08:raise RuntimeError('Period jump')
            row[key]=r
        width=row['fold']['beta']-row['symmetry']['beta']
        if width<=1e-10:raise RuntimeError('Nonpositive or unresolved core width')
        if abs(width/(source['fold']['beta']-source['symmetry']['beta'])-1)>.7:
            raise RuntimeError('Core width jumped')
        # Solve inside the core from the symmetry boundary, away from the
        # fold's double root, and explicitly check attraction.
        b=row['symmetry']['beta']+.45*width
        q=f.solve(b,np.array(row['symmetry']['q']));r=f.record(b,q,True)
        if r['spectral_radius']>=.9999 or r['region_violation']>1e-7:raise RuntimeError('Core not attracting')
        row['core']=r;row['width']=width;row['status']='continued';return row
    except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:
        if depth>=7:raise RuntimeError(str(e))
        midalpha=(alpha+source['alpha'])/2;midm0=(m0+source['m0'])/2
        mid=advance(template,source,midalpha,midm0,depth+1)
        return advance(template,mid,alpha,m0,depth+1)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--nodes',type=int,default=9)
    p.add_argument('--families',nargs='*');p.add_argument('--resume',action='store_true')
    p.add_argument('--adding',action='store_true')
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    aa=np.linspace(8.94,9.06,args.nodes);mm=np.linspace(-8/7-.006,-8/7+.006,args.nodes)
    templates=registry()
    if args.adding:
        templates=[dict(id='P'+str(r['peak_count']),beta=r['core']['beta'],peak_count=r['peak_count'],fold=r['fold'],symmetry=r['symmetry'],
                        regions=r['regions'],targets=r['targets']) for r in json.loads((out/'extended-adding.json').read_text()) if r['status']=='verified_attracting_cycle']
    (out/('adding-registry.json' if args.adding else 'symmetric-registry.json')).write_text(json.dumps(templates,indent=2)+'\n',encoding='utf-8',newline='\n')
    for template in templates:
        if args.families and template['id'] not in args.families:continue
        path=out/(template['id']+'-sheet.json')
        if args.resume and path.exists() and json.loads(path.read_text()).get('complete'):continue
        start=time.monotonic();template['gradient']={}
        for key,target in [('fold',1),('symmetry',-1)]:
            edge=template[key];q=np.array(edge['q']);grad=[]
            base=FastFlights(Chua(),template['regions'],template['targets'],True)
            w=np.r_[q,edge['beta']];J=base.boundary_system(w,target,True)
            for param,h in [('alpha',1e-6),('m0',1e-7)]:
                vals=[]
                for sign in [-1,1]:
                    model=Chua(**{param:(9 if param=='alpha' else -8/7)+sign*h})
                    f=FastFlights(model,template['regions'],template['targets'],True)
                    vals.append(f.boundary_system(w,target))
                grad.append(float(np.linalg.solve(J,-(vals[1]-vals[0])/(2*h))[-1]))
            template['gradient'][key]=grad
        center=args.nodes//2;initial=dict(alpha=9.,m0=-8/7,fold=template['fold'],symmetry=template['symmetry'])
        centerrow=advance(template,initial,9.,-8/7);centerrow['index']=[center,center]
        found={(center,center):centerrow};queue=deque([(center,center)]);failures=[];attempted={(center,center)}
        while queue:
            i,j=queue.popleft()
            for ii,jj in [(i-1,j),(i+1,j),(i,j-1),(i,j+1)]:
                if not (0<=ii<len(aa) and 0<=jj<len(mm)) or (ii,jj) in attempted:continue
                attempted.add((ii,jj))
                try:
                    row=advance(template,found[i,j],aa[ii],mm[jj]);row['index']=[ii,jj]
                    found[ii,jj]=row;queue.append((ii,jj))
                except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:
                    failures.append(dict(index=[ii,jj],alpha=float(aa[ii]),m0=float(mm[jj]),error=str(e)))
            report=dict(id=template['id'],anchor_beta=template['beta'],peak_count=template['peak_count'],
                alpha=aa.tolist(),m0=mm.tolist(),gradient=template['gradient'],regions=template['regions'],targets=template['targets'],
                nodes=list(found.values()),failures=failures,complete=len(found)==len(aa)*len(mm),elapsed_seconds=time.monotonic()-start)
            path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(template['id'],template['peak_count'],'nodes',len(found),'failures',len(failures),'seconds',round(time.monotonic()-start,1),flush=True)


if __name__=='__main__':main()
