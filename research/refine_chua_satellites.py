"""Continue attracting representatives at fixed transverse multiplier."""
import argparse,json,time
from collections import deque
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_fast_flights import FastFlights

ROOT=Path(__file__).resolve().parents[1]


def advance(template,source,alpha,m0,depth=0):
    try:
        old=FastFlights(Chua(alpha=source['alpha'],m0=source['m0']),template['regions'],template['targets'],template['symmetric'])
        f=FastFlights(Chua(alpha=alpha,m0=m0),template['regions'],template['targets'],template['symmetric'])
        r=source['representative'];w=np.r_[r['q'],r['beta']];forcing=np.zeros(len(w));target=template['target']
        for param,delta,h in [('alpha',alpha-source['alpha'],1e-6),('m0',m0-source['m0'],1e-7)]:
            if delta==0:continue
            values=[]
            for sign in [-1,1]:
                kwargs=dict(alpha=source['alpha'],m0=source['m0']);kwargs[param]+=sign*h
                trial=FastFlights(Chua(**kwargs),template['regions'],template['targets'],template['symmetric'])
                values.append(trial.boundary_system(w,target))
            forcing+=delta*(values[1]-values[0])/(2*h)
        prediction=w-np.linalg.solve(old.boundary_system(w,target,True),forcing)
        b,q=f.boundary(prediction[-1],prediction[:-1],target);r=f.record(b,q,True)
        if r['region_violation']>1e-7 or r['spectral_radius']>=.9999:raise RuntimeError('Invalid or non-attracting orbit')
        if abs(r['period']/source['representative']['period']-1)>.08:raise RuntimeError('Period jump')
        return dict(alpha=float(alpha),m0=float(m0),representative=r,status='continued')
    except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:
        if depth>=6:raise RuntimeError(str(e))
        mid=advance(template,source,(alpha+source['alpha'])/2,(m0+source['m0'])/2,depth+1)
        return advance(template,mid,alpha,m0,depth+1)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--resume',action='store_true');p.add_argument('--adding',action='store_true')
    p.add_argument('--families',nargs='*')
    p.add_argument('--registry',type=Path)
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True);templates=[]
    if args.registry:
        templates=json.loads(args.registry.read_text())
    elif args.adding:
        for item in json.loads((out/'extended-adding.json').read_text()):
            if item['status']!='verified_attracting_cycle':continue
            r=item['representative'];templates.append(dict(id='P'+str(item['peak_count']),peak_count=item['peak_count'],
               regions=item['regions'],targets=item['targets'],symmetric=True,target=r['map_multipliers'][0][0],representative=r))
    else:
        for name in ['cycles','cycles-extra','adding-cycles']:
            for item in json.loads((ROOT/'results/chua-layer-census'/(name+'.json')).read_text())['records']:
                if item['status']!='verified_attracting_cycle' or item['half_symmetry_error']<2e-5 or abs(item['beta']-14.354)<1e-9:continue
                r=item['cycle'];templates.append(dict(peak_count=item['peak_count'],regions=item['regions'],targets=item['targets'],
                   symmetric=False,target=r['map_multipliers'][0][0],representative=r))
        templates.sort(key=lambda r:r['representative']['beta'])
        for k,t in enumerate(templates):t['id']=f'A{k+1:02d}'
    aa=np.linspace(8.94,9.06,9);mm=np.linspace(-8/7-.006,-8/7+.006,9)
    for template in templates:
        if args.families and template['id'] not in args.families:continue
        path=out/(template['id']+'-sheet.json')
        if args.resume and path.exists() and json.loads(path.read_text())['complete']:continue
        started=time.monotonic();initial=dict(alpha=9.,m0=-8/7,representative=template['representative']);center=advance(template,initial,9.,-8/7);center['index']=[4,4]
        found={(4,4):center};queue=deque([(4,4)]);tried={(4,4)};failures=[]
        while queue:
            i,j=queue.popleft()
            for ii,jj in [(i-1,j),(i+1,j),(i,j-1),(i,j+1)]:
                if not (0<=ii<9 and 0<=jj<9) or (ii,jj) in tried:continue
                tried.add((ii,jj))
                try:
                    row=advance(template,found[i,j],aa[ii],mm[jj]);row['index']=[ii,jj];found[ii,jj]=row;queue.append((ii,jj))
                except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:
                    failures.append(dict(index=[ii,jj],error=str(e)))
            report=dict(id=template['id'],peak_count=template['peak_count'],anchor_beta=template['representative']['beta'],symmetric=template['symmetric'],
                regions=template['regions'],targets=template['targets'],target_multiplier=template['target'],alpha=aa.tolist(),m0=mm.tolist(),
                nodes=list(found.values()),failures=failures,complete=len(found)==81,elapsed_seconds=time.monotonic()-started,
                scope='Stable representative surface at fixed transverse multiplier; not a complete stability-volume boundary.')
            path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(template['id'],template['peak_count'],'nodes',len(found),'failures',len(failures),'seconds',round(time.monotonic()-started,1),flush=True)


if __name__=='__main__':main()
