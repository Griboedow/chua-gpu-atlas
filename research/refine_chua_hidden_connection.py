"""High-precision numerical centers for the separate exact-manifold certificate."""
import argparse
import json
import numpy as np
from flint import arb,arb_mat,ctx
from chua_research_followup import OUT,dump
from chua_interval_lobe import IntervalLobe,saddle_center
from chua_interval_map import col,midpoint,bounds
from certify_chua_layer_organizer import enclose_saddle,submatrix
from chua_interval_connection import connect


def main(pre,post,digits):
    data=json.loads((OUT/'hidden-search/lobe-saddles.json').read_text())['records'];row=next(r for r in data if r['parameter_index']==818)
    root=next(r for r in row['roots'] if r['multipliers'][0][0]<-1)
    seed=json.loads((OUT/'hidden-search/connection-newton-818-embedded.json').read_text())
    a0,b0=seed['history'][-1]['coordinate'];f=IntervalLobe(row['parameters'],digits);b=arb(str(b0))
    p=saddle_center(f,b,root['point']);_,J,_,_=f.evaluate(p.entries(),b)
    eig,Vn=np.linalg.eig(np.array([[float(J[i,j].mid()) for j in range(2)] for i in range(2)]))
    order=np.argsort(-abs(eig));eig,Vn=eig[order].real,Vn[:,order].real
    if Vn[0,0]<0:Vn[:,0]*=-1
    V=arb_mat([[arb(str(v)) for v in r] for r in Vn]);W=V.inv()
    a=arb(str(a0))/arb(str(eig[0]))**pre;n=seed['iterations']+pre+post;history=[]
    for iteration in range(15):
        center=saddle_center(f,b,root['point']);p,pb,pbb,cert=enclose_saddle(f,center,b,'1e-65',inflation_floor='1e-75')
        # This stage computes a numerical seed. Midpoints prevent a deliberately
        # inflated root enclosure from dominating many subsequent iterations.
        # The separate certificate will restore the full fixed-point enclosure.
        p,pb,pbb=midpoint(p),midpoint(pb),midpoint(pbb)
        g,J,H,records=connect(f,a,arb(0),b,b,p,pb,pb,pbb,V,W,n)
        residual=col([g[0,0],J[0,0]])
        history.append(dict(iteration=iteration,a=bounds(a),beta=bounds(b),gap=bounds(g[0,0]),slope=bounds(J[0,0]),stable=bounds(g[1,0])))
        print('center',pre,post,iteration,'gap',g[0,0].str(10),'slope',J[0,0].str(10),flush=True)
        if abs(g[0,0])<arb('1e-28') and abs(J[0,0])<arb('1e-18'):break
        jac=arb_mat([[J[0,0],J[0,2]],[H[0][0,0],H[0][0,2]]])
        step=midpoint(jac.inv()*residual);a=(a-step[0,0]).mid();b=(b-step[1,0]).mid()
    else:
        dump(OUT/'hidden-search'/f'connection-center-{pre}-{post}-failed.json',dict(history=history))
        raise ValueError('Connection center refinement exhausted')
    dump(OUT/'hidden-search'/f'connection-center-{pre}-{post}.json',dict(parameters=row['parameters'],
        a=a.mid().str(digits,radius=False),beta=b.mid().str(digits,radius=False),iterations=n,
        pre=pre,post=post,coordinates=[[V[i,j].mid().str(digits,radius=False) for j in range(2)] for i in range(2)],
        parent_guess=root['point'],history=history,stable_return=bounds(g[1,0]),
        scope='Numerical tangent-line root. Does not assert an exact invariant-manifold tangency.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pre',type=int,default=14);ap.add_argument('--post',type=int,default=2);ap.add_argument('--digits',type=int,default=140);a=ap.parse_args();main(a.pre,a.post,a.digits)
