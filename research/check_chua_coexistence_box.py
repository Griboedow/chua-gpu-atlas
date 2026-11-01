"""Independent mpmath finite-difference check of all shooting derivatives."""
import json
import mpmath as mp
from flint import ctx
from certify_chua_coexistence_box import MultipleShooting,nominal
from chua_research_followup import ROOT,OUT,dump


def residual(q,p,schedule,start_x):
    a,b,m0,m1,g=p;out=[];n=len(schedule)
    for k,(reg,target,_,sign) in enumerate(schedule):
        slope=m0 if reg==0 else m1
        H=mp.matrix([[-a*(1+slope),a,0,-a*(m0-m1)*reg],
                     [1,-1,1,0],[0,-b,-g,0],[0,0,0,0]])
        X=mp.expm(H*q[3*k+2])*mp.matrix([start_x[k],q[3*k],q[3*k+1],1])
        nxt=(k+1)%n
        out.extend([X[0]-target,X[1]-q[3*nxt],X[2]-q[3*nxt+1]])
    return mp.matrix(out)


def main():
    ctx.dps=100;mp.mp.dps=80
    records=json.loads((ROOT/'results/chua-layer-atlas/four-cycles.json').read_text())['records']
    rows=[];h=mp.mpf('1e-20')
    for rec in records:
        f=MultipleShooting(rec,100);q=[mp.mpf(v.mid().str(90,radius=False)) for v in f.center.entries()]
        p=[mp.mpf(v.mid().str(90,radius=False)) for v in nominal()]
        B=f.parameter_jacobian(f.center,nominal());_,J,_=f.evaluate(f.center,nominal())
        errors=[]
        for kind,values,M in [('parameters',p,B),('shooting',q,J)]:
            maximum=mp.mpf(0)
            for j in range(len(values)):
                plus=values.copy();minus=values.copy();plus[j]+=h;minus[j]-=h
                if kind=='parameters':d=(residual(q,plus,f.schedule,f.start_x)-residual(q,minus,f.schedule,f.start_x))/(2*h)
                else:d=(residual(plus,p,f.schedule,f.start_x)-residual(minus,p,f.schedule,f.start_x))/(2*h)
                for i,value in enumerate(d):
                    reference=mp.mpf(M[i,j].mid().str(90,radius=False))
                    maximum=max(maximum,abs(value-reference))
            assert maximum<mp.mpf('1e-32'),(rec['id'],kind,maximum)
            errors.append(dict(variables=kind,max_absolute_difference=mp.nstr(maximum,16)))
        rows.append(dict(id=rec['id'],checks=errors,passed=True))
        print(rec['id'],errors,flush=True)
    dump(OUT/'coexistence-box/derivative-checks.json',dict(records=rows,working_digits=80,
        centered_difference_step='1e-20',method='Independent mpmath 4x4 affine exponentials',
        scope='Numerical implementation cross-check, not an independent interval proof.',passed=True))


if __name__=='__main__':main()
