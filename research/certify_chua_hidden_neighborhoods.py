"""Exclude non-equilibrium attraction near all equilibria of the new candidates.

This certificate does not assert existence of a chaotic attractor. It proves
    that every attracting periodic orbit at these parameters must be hidden.
"""
import argparse
import json
import numpy as np
from flint import arb,arb_mat,acb,ctx
from chua_research_followup import ROOT,OUT,dump
from chua_interval_map import col,bounds
from certify_chua_coexisting_cycles import ReturnMap
from certify_chua_hidden_cycle import lyapunov,outer_path


class Affine(ReturnMap):
    def __init__(self,p,digits):
        ctx.dps=digits
        self.alpha,self.beta,self.m0,self.m1=[arb(str(p[k])) for k in ['alpha','beta','m0','m1']]


def modal_basis(A,a,b,m,root_radius='1e-50'):
    a1=1+a*(1+m);a2=b+a*m;a3=a*b*(1+m)
    poly=lambda x:x**3+a1*x*x+a2*x+a3
    deriv=lambda x:3*x*x+2*a1*x+a2
    An=np.array([[float(A[i,j].mid()) for j in range(3)] for i in range(3)])
    guess=max(np.linalg.eigvals(An).real);t=arb(str(guess))
    for _ in range(20):t=(t-poly(t)/deriv(t)).mid()
    box=arb(t,arb(root_radius));root=t-poly(t)/deriv(box)
    assert box.contains_interior(root) and root>0
    lam=root;sigma=(-a1-lam)/2;omega=(-a3/lam-sigma*sigma).sqrt()
    assert sigma<0 and omega>0
    mu=acb(sigma,omega);yc=1+m+mu/a;zc=-b*yc/mu
    v=col([1,1+m+lam/a,-b*(1+m+lam/a)/lam])
    S=arb_mat([[1,1,0],[v[1,0],yc.real,yc.imag],[v[2,0],zc.real,zc.imag]])
    assert not S.det().contains(0)
    return S,lam,dict(unstable=bounds(lam),stable_real=bounds(sigma),stable_imaginary=bounds(omega))


def main(digits,uniform=False):
    candidates=json.loads((OUT/'hidden-search/outer-screen.json').read_text())['candidates'];rows=[]
    if uniform:
        center=json.loads((OUT/'hidden-search/connection-center-14-2.json').read_text())
        candidates=[dict(parameter_index=818,parameters={**center['parameters'],'beta':center['beta']})]
    for row in candidates:
        p=row['parameters'];f=Affine(p,digits);result=dict(parameter_index=row['parameter_index'],parameters=p,passed=False)
        if uniform:
            f.beta=arb(f.beta,arb('1e-8'))
            result['uniform_beta_interval']=bounds(f.beta)
        try:
            A0,c0,_=f.matrices(None,0);A1,c1,_=f.matrices(None,1)
            assert f.m0<-1 and f.m1>-1 and c1[0,0]>1
            S,lam,eig=modal_basis(A0,f.alpha,f.beta,f.m0,'1e-6' if uniform else '1e-50');P,lyap=lyapunov(A1)
            level=(c1[0,0]-1)**2/(2*P.inv()[0,0]);delta=arb('1e-6')
            start=S*col([arb('.5'),arb(0,delta),arb(0,delta)])
            guess=float((arb(2).log()/lam).mid())
            T=f.crossing(A0,c0,start,1,guess);hit,_=f.flow(A0,c0,start,T);hit[0,0]=1
            central=f.verify_arc(A0,c0,dict(duration=T,start=start,region=0,target=1))
            for duration in [10,20,40,80,160,320,640]:
                end,_=f.flow(A1,c1,hit,arb(duration));offset=end-c1
                energy=(offset.transpose()*P*offset)[0,0]
                if energy<level:break
            else:raise ValueError('Outer affine flight not captured')
            cells=outer_path(f,A1,c1,hit,arb(duration))
            W=S.inv();norm=sum(v*v for v in W.entries()).sqrt()
            origin_radius=(delta/(2*norm)).lower()
            Pnorm=max(sum(abs(P[i,j]) for j in range(3)).abs_upper() for i in range(3))
            outer_radius=(level/(2*Pnorm)).sqrt().lower()
            result.update(passed=True,spectrum=eig,outer_lyapunov=lyap,central_path=central,
                modal_basis=[[bounds(S[i,j]) for j in range(3)] for i in range(3)],
                exit_state=[bounds(v) for v in hit.entries()],outer_equilibrium=[bounds(v) for v in c1.entries()],
                central_exit_time=bounds(T),outer_path_cells=cells,outer_duration=duration,
                captured_energy=bounds(energy),ellipsoid_level=bounds(level),
                origin_ball_radius=bounds(origin_radius),outer_ball_radius=bounds(outer_radius),
                unstable_exit_coordinate='.5',stable_coordinate_radius='1e-6')
        except (ValueError,AssertionError,ZeroDivisionError) as exc:result['error']=str(exc)
        rows.append(result);print('neighborhood exclusion',row['parameter_index'],result['passed'],result.get('error'),flush=True)
        prefix='uniform-equilibrium-exclusion' if uniform else 'equilibrium-exclusion'
        dump(OUT/'hidden-search'/f'{prefix}-{digits}.json',dict(working_digits=digits,records=rows,
            passed=all(r['passed'] for r in rows),complete=len(rows)==len(candidates),
            scope='Certified equilibrium neighborhoods exclude the basin of every attracting periodic orbit, or compact attractor disjoint from all equilibria, if one exists. No chaotic-attractor or homoclinic-tangency existence is inferred.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--digits',type=int,default=80);ap.add_argument('--uniform',action='store_true');a=ap.parse_args();main(a.digits,a.uniform)
