"""Finite-iterate homoclinic candidates for the boundary saddle of a hidden lobe."""
import argparse
import json
import numpy as np
from scipy.optimize import brentq,root
from chua_lobe_map import LobeMap
from chua_research_followup import OUT,dump


def parent(f,beta,guess):
    p,J,aux=f.fixed(guess,beta);ev,V=np.linalg.eig(J);order=np.argsort(-abs(ev));ev,V=ev[order].real,V[:,order].real
    if not abs(ev[0])>1 or not abs(ev[1])<1:raise ValueError('Expected real saddle')
    if V[0,0]<0:V[:,0]*=-1
    return p,ev,V,np.linalg.inv(V),aux


def connection(f,beta,guess,a,n):
    p,ev,V,W,aux=parent(f,beta,guess);state=p+a*V[:,0];tangent=V[:,0].copy();duration=0.
    for _ in range(n):
        state,J,aux=f.apply(state,beta);tangent=J@tangent;duration+=aux['period']
    return float(W[0]@(state-p)),float(W[0]@tangent),float(W[1]@(state-p)),state


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--index',type=int,default=818);ap.add_argument('--solve',action='store_true');ap.add_argument('--embedded',action='store_true');a=ap.parse_args()
    data=json.loads((OUT/'hidden-search/lobe-saddles.json').read_text())['records']
    row=next(r for r in data if r['parameter_index']==a.index)
    saddle=next(r for r in row['roots'] if (r['multipliers'][0][0]<-1 if a.embedded else r['multipliers'][0][0]>1))
    tag=str(a.index)+('-embedded' if a.embedded else '')
    p=row['parameters'];f=LobeMap(p);q,ev,V,W,aux=parent(f,p['beta'],saddle['point'])
    if a.solve:
        screen=json.loads((OUT/'hidden-search'/f'connection-screen-{tag}.json').read_text())
        candidates=[(r,z) for r in screen['records'] for z in r['folds']]
        first,fold=min(candidates,key=lambda rz:abs(rz[1]['gap'])+abs(rz[1]['stable'])) if a.embedded else candidates[0]
        w=np.array([fold['a']/1e-4,p['beta']]);n=first['iterations'];history=[];status='unresolved'
        for iteration in range(160):
            def fun(v):return np.array(connection(f,v[1],q,v[0]*1e-4,n)[:2])*[1,1e-4]
            try:
                residual=fun(w)
                history.append(dict(iteration=iteration,coordinate=(w*[1e-4,1]).tolist(),residual=residual.tolist()))
                if max(abs(residual))<1e-8:status='finite_iterate_root';break
                J=np.column_stack([(fun(w+np.eye(2)[j]*h)-fun(w-np.eye(2)[j]*h))/(2*h) for j,h in enumerate([1e-4,1e-5])])
                change=np.linalg.solve(J,-residual)
                change*=min(1.,.15/max(abs(change[1]),1e-99),.25/max(abs(change[0]),1e-99))
                for factor in [1.,.5,.25,.125,.0625,.03125]:
                    trial=w+factor*change
                    try:improved=max(abs(fun(trial)))<max(abs(residual))
                    except (ValueError,RuntimeError,np.linalg.LinAlgError):improved=False
                    if improved:w=trial;break
                else:raise ValueError('Connection corrector stalled')
            except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
                history.append(dict(error=str(exc)));break
            if iteration%10==0:print('connection solve',iteration,w,residual,flush=True)
            dump(OUT/'hidden-search'/f'connection-newton-{tag}.json',dict(status=status,iterations=n,history=history))
        dump(OUT/'hidden-search'/f'connection-newton-{tag}.json',dict(status=status,iterations=n,history=history,
            scope='Finite-iterate line tangency only. Exact invariant manifolds, admissibility and hiddenness need separate certificates.'))
        print(status,history[-1],flush=True);return
    amin=1e-4;amax=amin*(ev[0]**2 if ev[0]<0 else ev[0]);nbase=int(np.ceil(np.log(.3/amin)/np.log(abs(ev[0]))))
    records=[];failures=[]
    for n in range(max(2,nbase-2),nbase+9):
        sample=[]
        for amplitude in np.linspace(amin,amax,81):
            try:
                g,d,s,state=connection(f,p['beta'],q,amplitude,n)
                sample.append(dict(a=float(amplitude),gap=g,derivative=d,stable=s))
            except (ValueError,RuntimeError,np.linalg.LinAlgError):sample.append(dict(a=float(amplitude),failed=True))
        folds=[]
        for l,r in zip(sample[:-1],sample[1:]):
            if l.get('failed') or r.get('failed') or l['derivative']*r['derivative']>=0:continue
            try:
                aa=brentq(lambda x:connection(f,p['beta'],q,x,n)[1],l['a'],r['a'],xtol=1e-14)
                g,d,s,state=connection(f,p['beta'],q,aa,n)
                folds.append(dict(a=aa,gap=g,stable=s,point=state.tolist()))
            except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:failures.append(str(exc))
        records.append(dict(iterations=n,samples=sample,folds=folds))
        dump(OUT/'hidden-search'/f'connection-screen-{tag}.json',dict(parameters=p,saddle=saddle,records=records,failures=failures,
            scope='Finite-iterate line images and stationary folds, not exact invariant-manifold intersections.'))
        print('folds',a.index,n,[(z['a'],z['gap'],z['stable']) for z in folds],flush=True)


if __name__=='__main__':main()
