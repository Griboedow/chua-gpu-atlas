"""Continue cycles through changes of diode-breakpoint itinerary."""
import argparse,json,time
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_event_map import EventMap
from chua_layer_bifurcations import from_orbit
from chua_fast_flights import FastFlights


class ContinuationBudget(Exception):
    pass


deadline = float('inf')


def step(source,alpha,m0,count,target,depth=0):
    if time.monotonic() > deadline:
        raise ContinuationBudget('Per-node continuation time limit; point remains unresolved')
    model=Chua(alpha=source['alpha'],m0=source['m0']);old=EventMap(model,count)
    w=np.array(source['w']);J=old.jacobian(w,target);forcing=np.zeros(3)
    try:
        for parameter,delta,h in [('alpha',alpha-source['alpha'],1e-6),('m0',m0-source['m0'],1e-7)]:
            if not delta:continue
            vals=[]
            for s in [-1,1]:
                kwargs=dict(alpha=source['alpha'],m0=source['m0']);kwargs[parameter]+=s*h
                vals.append(EventMap(Chua(**kwargs),count).residual(w,target))
            forcing+=delta*(vals[1]-vals[0])/(2*h)
        prediction=w-np.linalg.solve(J,forcing)
        f=EventMap(Chua(alpha=alpha,m0=m0),count);new,aux=f.solve(prediction,target)
        scale=np.array([.1,.1,1.]);ratio=float(np.linalg.norm((new-prediction)/scale)/max(np.linalg.norm((prediction-w)/scale),1e-8))
        if ratio>.3 or abs(aux['period']/source['period']-1)>.01:raise RuntimeError('Nonlocal event-map correction')
        if 'half_sign' in source and (aux['half_inversion_error']>2e-7 or source['half_sign']*aux['half_trace']<0):
            raise RuntimeError('Symmetry or half-map branch changed')
        result=dict(alpha=float(alpha),m0=float(m0),w=new.tolist(),period=aux['period'],correction_ratio=ratio,
                    half_inversion_error=aux['half_inversion_error'],half_trace=aux['half_trace'])
        if 'half_sign' in source:result['half_sign']=source['half_sign']
        return result
    except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:
        if depth>=8:raise RuntimeError(str(e))
        mid=step(source,(alpha+source['alpha'])/2,(m0+source['m0'])/2,count,target,depth+1)
        return step(mid,alpha,m0,count,target,depth+1)


def continue_branch(source,alpha,m0,count,target):
    n=max(1,int(np.ceil(abs(alpha-source['alpha'])/.002)),int(np.ceil(abs(m0-source['m0'])/.0002)))
    aa=np.linspace(source['alpha'],alpha,n+1)[1:];mm=np.linspace(source['m0'],m0,n+1)[1:];path=[]
    for a,m in zip(aa,mm):source=step(source,float(a),float(m),count,target);path.append(source)
    return source,path


def record(source,symmetric):
    model=Chua(alpha=source['alpha'],m0=source['m0'],rtol=2e-12,step=.02);y,z,b=source['w']
    family,q=from_orbit(model,b,[y,z,source['period']],symmetric);family=FastFlights(model,family.regions,family.targets,symmetric)
    r=family.record(b,q,True)
    if r['region_violation']>1e-7 or r['spectral_radius']>=.999:raise RuntimeError('Invalid bridged orbit')
    return dict(alpha=source['alpha'],m0=source['m0'],representative=r,regions=family.regions,targets=family.targets,
                status='continued_with_automatic_itinerary')


def local_step(template,source,alpha,m0,depth=0):
    """Small fixed-itinerary correction with an explicit locality check."""
    if time.monotonic()>deadline:raise ContinuationBudget('Per-node local continuation time limit')
    target=template['target_multiplier'];sym=template['symmetric'];r=source['representative']
    old=FastFlights(Chua(alpha=source['alpha'],m0=source['m0']),template['regions'],template['targets'],sym)
    w=np.r_[r['q'],r['beta']];forcing=np.zeros(len(w))
    try:
        for parameter,delta,h in [('alpha',alpha-source['alpha'],1e-6),('m0',m0-source['m0'],1e-7)]:
            if not delta:continue
            vals=[]
            for sign in [-1,1]:
                kw=dict(alpha=source['alpha'],m0=source['m0']);kw[parameter]+=sign*h
                f=FastFlights(Chua(**kw),template['regions'],template['targets'],sym)
                vals.append(f.boundary_system(w,target))
            forcing+=delta*(vals[1]-vals[0])/(2*h)
        prediction=w-np.linalg.solve(old.boundary_system(w,target,True),forcing)
        f=FastFlights(Chua(alpha=alpha,m0=m0),template['regions'],template['targets'],sym)
        beta,q=f.boundary(prediction[-1],prediction[:-1],target);new=np.r_[q,beta]
        scale=np.r_[.1,.1,np.ones(len(w)-2)]
        ratio=float(np.linalg.norm((new-prediction)/scale)/max(np.linalg.norm((prediction-w)/scale),1e-8))
        if ratio>.3:raise RuntimeError('Nonlocal affine correction')
        result=f.record(beta,q,True)
        if abs(result['period']/r['period']-1)>.01 or result['region_violation']>1e-7 or result['spectral_radius']>=.999:
            raise RuntimeError('Invalid local affine orbit')
        return dict(alpha=float(alpha),m0=float(m0),representative=result,regions=f.regions,targets=f.targets,
                    correction_ratio=ratio,status='continued_with_bounded_local_affine_steps')
    except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:
        if depth>=8:raise RuntimeError(str(e))
        mid=local_step(template,source,(source['alpha']+alpha)/2,(source['m0']+m0)/2,depth+1)
        return local_step(template,mid,alpha,m0,depth+1)


def continue_reference_branch(template,source,alpha,m0):
    n=max(1,int(np.ceil(abs(alpha-source['alpha'])/.002)),int(np.ceil(abs(m0-source['m0'])/.0002)))
    aa=np.linspace(source['alpha'],alpha,n+1)[1:];mm=np.linspace(source['m0'],m0,n+1)[1:];history=[]
    for a,m in zip(aa,mm):
        source=local_step(template,source,float(a),float(m))
        history.append(dict(alpha=source['alpha'],m0=source['m0'],beta=source['representative']['beta'],
                            period=source['representative']['period'],correction_ratio=source['correction_ratio']))
    return source,history


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--recover-off21',action='store_true');p.add_argument('--families',nargs='*',default=['E09','E13'])
    p.add_argument('--rebuild',action='store_true')
    p.add_argument('--seconds-per-node',type=float,default=40.)
    p.add_argument('--fast-local',action='store_true',help='Fill remaining holes with bounded small affine continuation steps.')
    args=p.parse_args();out=args.output
    if args.recover_off21:
        local=json.loads((out/'off21-local-branch-relations.json').read_text());r=local['root']
        start=dict(alpha=r['alpha'],m0=r['m0'],w=[*r['q'][:2],r['beta']],period=r['period'],half_sign=1)
        endpoint,path=continue_branch(start,9.,-8/7,3,.25);node=record(endpoint,True)
        records=json.loads((out/'offcenter-cycles.json').read_text());r=records[21]
        r['rejected_continuation']=dict(previous=r.get('canonical'),reason='A single large parameter jump did not establish branch continuity.')
        r['canonical']=dict(representative=node['representative'],regions=node['regions'],targets=node['targets'],symmetric=True,target=.5)
        r['canonical_recovery']='Automatic event itinerary, small parameter steps, and bounded predictor corrections.'
        (out/'off21-continuation.json').write_text(json.dumps(path,indent=2)+'\n',encoding='utf-8',newline='\n')
        (out/'offcenter-cycles.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8',newline='\n')
        print('off21 verified continuation',node['representative']['beta'],node['representative']['period'],node['regions'],flush=True);return
    registry={r['id']:r for r in json.loads((out/'family-registry.json').read_text())['families']}
    for id in args.families:
        path=out/(id+'-sheet.json');sheet=json.loads(path.read_text());nodes={tuple(n['index']):n for n in sheet['nodes']};failures=[];bridges=sheet.get('bridges',[])
        if args.fast_local:sheet.setdefault('prior_bridge_failures',[]).extend(sheet.get('bridge_failures',[]))
        if args.rebuild:nodes={(4,4):nodes[4,4]};bridges=[]
        count=registry[id]['section_returns'];target=sheet['target_multiplier']**2 if sheet['symmetric'] else sheet['target_multiplier']
        missing=sorted([(i,j) for i in range(9) for j in range(9) if (i,j) not in nodes],key=lambda ij:abs(ij[0]-4)+abs(ij[1]-4))
        for ij in missing:
            global deadline
            deadline=time.monotonic()+args.seconds_per_node
            near=min(nodes,key=lambda k:abs(k[0]-ij[0])+abs(k[1]-ij[1]));n=nodes[near];r=n['representative']
            source=dict(alpha=n['alpha'],m0=n['m0'],w=[*r['q'][:2],r['beta']],period=r['period'])
            if sheet['symmetric']:source['half_sign']=1 if sheet['target_multiplier']>0 else -1
            try:
                if args.fast_local:
                    template=dict(symmetric=sheet['symmetric'],target_multiplier=sheet['target_multiplier'],regions=n.get('regions',sheet['regions']),targets=n.get('targets',sheet['targets']))
                    node,history=continue_reference_branch(template,n,sheet['alpha'][ij[0]],sheet['m0'][ij[1]])
                else:
                    endpoint,history=continue_branch(source,sheet['alpha'][ij[0]],sheet['m0'][ij[1]],count,target)
                    node=record(endpoint,sheet['symmetric'])
                node['index']=list(ij);nodes[ij]=node
                bridges.append(dict(index=list(ij),source=list(near),method='bounded_local_affine' if args.fast_local else 'automatic_event_map',steps=history));print(id,ij,'recovered',node['representative']['beta'],len(node['regions']),flush=True)
            except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError,ContinuationBudget) as e:failures.append(dict(index=list(ij),error=str(e)));print(id,ij,'unresolved',str(e)[:120],flush=True)
            sheet['nodes']=list(nodes.values());sheet['complete']=len(nodes)==81;sheet['bridge_failures']=failures;sheet['bridges']=bridges
            path.write_text(json.dumps(sheet,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
