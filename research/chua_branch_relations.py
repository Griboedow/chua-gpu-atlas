"""Trace attracting branches back through doubling and symmetry transitions.

The +1 edge is approached along the multiplier, rather than found by a single
unconstrained jump. Half-period closure distinguishes a doubled daughter from
an inversion-symmetry branch. Failures remain unresolved identities.
"""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_layer_bifurcations import from_orbit
from chua_fast_flights import FastFlights

ROOT=Path(__file__).resolve().parents[1]


def canonical(f,b,q,target):
    b,q=f.boundary(b,q,target);r=f.record(b,q,True)
    if r['region_violation']>1e-7 or r['spectral_radius']>=.9999:raise RuntimeError('Canonical orbit invalid')
    return b,q,r


def branch_relations(template):
    origin=template['representative'];model=Chua(alpha=origin['alpha'],m0=origin['m0'],m1=origin['m1'])
    symmetric=template['symmetric'];f=FastFlights(model,template['regions'],template['targets'],symmetric)
    r=template['representative'];q=np.array(r['q']);b=r['beta'];chain=[]
    for depth in range(6):
        if symmetric:
            b,q,r=canonical(f,b,q,.5)
            return dict(status='identified',symmetric=True,root=r,regions=f.regions,targets=f.targets,chain=chain)
        current=f.record(b,q);eigen=current['map_multipliers'][0][0]
        prevperiod=current['period']
        edge_approach=False
        for target in np.linspace(eigen,1.,41)[1:]:
            try:bb,qq=f.boundary(b,q,float(target))
            except (RuntimeError,ValueError,np.linalg.LinAlgError):
                if target!=1.:raise
                bb,qq=f.boundary(b,q,.99999);edge_approach=True
            rr=f.record(bb,qq,True)
            if rr['region_violation']>1e-7 or abs(rr['period']/prevperiod-1)>.03:raise RuntimeError('Edge continuation jumped')
            b,q=bb,qq;prevperiod=rr['period']
        edge=f.record(b,q,True);u=np.r_[0.,q[:2]];T=edge['period'];half=model.integrate(b,u,T/2)
        returning=float(max(abs(half.y[:,-1]-u)));inverting=float(max(abs(half.y[:,-1]+u)))
        # Daughter amplitudes vanish at the edge, but the degenerate closing
        # root can retain a small amplitude after its multiplier has converged.
        # A near half-return is only a candidate: the independently closed
        # parent must also have its multiplier near -1 below.
        tolerance=5e-4
        if returning<tolerance:
            parent,pq=from_orbit(model,b,[*q[:2],T/2],False);parent=FastFlights(model,parent.regions,parent.targets,False)
            pr=parent.record(b,pq,True)
            if abs(pr['map_multipliers'][0][0]+1)>1e-3:raise RuntimeError('Half-period parent is not at a doubling')
            exact_b,exact_q=parent.boundary(b,pq,-1)
            chain.append(dict(type='period_doubling',beta=exact_b,daughter_period=T,parent_period=T/2,half_return_error=returning,near_edge_fallback=edge_approach))
            f=parent;q=pq
            for target in [-.99,-.9,-.8,-.7,-.6,-.5]:b,q,r=canonical(f,b,q,target)
            # A stable parent may itself be inversion symmetric.
            sol=model.integrate(b,[0.,*q[:2]],r['period']/2)
            symmetric=max(abs(sol.y[:,-1]+np.r_[0.,q[:2]]))<2e-5
            if symmetric:
                half,hq=from_orbit(model,b,[*q[:2],r['period']],True)
                f=FastFlights(model,half.regions,half.targets,True);q=hq
        elif inverting<tolerance:
            chain.append(dict(type='symmetry_breaking',beta=b,period=T,half_inversion_error=inverting))
            half,hq=from_orbit(model,b,[*q[:2],T],True)
            f=FastFlights(model,half.regions,half.targets,True);q=hq;symmetric=True
        else:
            # A fold-born asymmetric pair: return to an interior root on the
            # same full-period branch and use it as the identity representative.
            chain.append(dict(type='fold_candidate',beta=b,period=T,half_return_error=returning,half_inversion_error=inverting))
            source=current;b,q,r=canonical(f,source['beta'],np.array(source['q']),-.5)
            return dict(status='identified',symmetric=False,root=r,regions=f.regions,targets=f.targets,chain=chain)
    raise RuntimeError('Branch-relation search depth limit')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--offcenter',action='store_true')
    p.add_argument('--local-offcenter-index',type=int)
    p.add_argument('--families',nargs='*')
    args=p.parse_args();out=args.output;items=[]
    if args.local_offcenter_index is not None:
        r=json.loads((out/'offcenter-cycles.json').read_text())[args.local_offcenter_index]
        item=dict(representative=r['cycle'],regions=r['regions'],targets=r['targets'],symmetric=False)
        result=branch_relations(item)
        (out/f'off{args.local_offcenter_index}-local-branch-relations.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(result['status'],result['root']['beta'],flush=True);return
    if args.offcenter:
        for r in json.loads((out/'offcenter-cycles.json').read_text()):
            if 'canonical' in r:items.append(dict(id='off'+str(r['candidate_index']),**r['canonical']))
        filename='offcenter-branch-relations.json'
    else:
        for path in sorted(out.glob('*-sheet.json')):
            r=json.loads(path.read_text());center=next(x for x in r['nodes'] if x['index']==[4,4])
            if args.families and r['id'] not in args.families:continue
            items.append(dict(id=r['id'],regions=r['regions'],targets=r['targets'],symmetric=r.get('symmetric',True),representative=center['representative']))
        filename='sheet-branch-relations.json'
    results=[]
    for item in items:
        try:result=dict(id=item['id'],**branch_relations(item))
        except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as e:result=dict(id=item['id'],status='unresolved',error=str(e))
        results.append(result);(out/filename).write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(result['id'],result['status'],result.get('root',{}).get('beta'),[x['type'] for x in result.get('chain',[])],result.get('error'),flush=True)


if __name__=='__main__':main()
