"""Independent derivative and continuous-flow controls for the basin certificate.

mpmath finite differences and SciPy ODE solves are implementation controls;
they do not replace the outward-rounded inclusions in the certificate.
"""
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from flint import arb,ctx
from certify_chua_basin_separator import (ROOT,OUT,refine_root,local_block,
    propagate,continuous_origin_capture,origin_ellipsoid,col,npmat)


def mp_return(point,schedule):
    state=mp.matrix([0,*point])
    for reg,target,guess,sign in schedule:
        slope=-mp.mpf(19)/20 if reg==0 else -mp.mpf(1)/20
        A=mp.matrix([[-mp.mpf(42)/5*(1+slope),mp.mpf(42)/5,0],[1,-1,1],[0,-12,0]])
        x=reg*(mp.mpf(18)/20)/(mp.mpf(19)/20)
        c=mp.matrix([x,0,-x]);orbit=lambda t:c+mp.expm(A*t)*(state-c)
        g=mp.mpf(str(guess));t=mp.findroot(lambda t:orbit(t)[0]-target,(g*mp.mpf('.9999'),g*mp.mpf('1.0001')),tol=mp.mpf('1e-75'))
        state=orbit(t);state[0]=target
    return -state[1],-state[2]


def main():
    ctx.dps=100;mp.mp.dps=90
    seeds=json.loads((ROOT/'results/chua-hidden-boundary/cycle-seeds.json').read_text())['full_records']
    rows=[]
    for record in seeds:
        f,p,pbox,J0=refine_root(record,100)
        value,J,arcs=f.evaluate(p.entries());H=arcs[-1]['DD']
        point=[mp.mpf(v.mid().str(90,radius=False)) for v in p.entries()]
        base=mp_return(point,f.schedule_rows);errors=[]
        for hs in ['1e-12','3e-13']:
            h=mp.mpf(hs);first=[];second=[]
            for i in range(2):
                left=point.copy();right=point.copy();left[i]-=h;right[i]+=h
                yl=mp_return(left,f.schedule_rows);yr=mp_return(right,f.schedule_rows)
                for k in range(2):
                    d1=(yr[k]-yl[k])/(2*h);d2=(yr[k]-2*base[k]+yl[k])/(h*h)
                    first.append(abs(d1-mp.mpf(J[k,i].mid().str(90,radius=False)))/(1+abs(d1)))
                    second.append(abs(d2+mp.mpf(H[i][i][k+1,0].mid().str(90,radius=False)))/(1+abs(d2)))
            corners=[mp_return([point[0]+sx*h,point[1]+sy*h],f.schedule_rows) for sx,sy in [(1,1),(1,-1),(-1,1),(-1,-1)]]
            for k in range(2):
                d2=(corners[0][k]-corners[1][k]-corners[2][k]+corners[3][k])/(4*h*h)
                second.append(abs(d2+mp.mpf(H[0][1][k+1,0].mid().str(90,radius=False)))/(1+abs(d2)))
            errors.append(dict(step=hs,jacobian=mp.nstr(max(first),12),hessian=mp.nstr(max(second),12)))
        value_error=max(abs(base[k]-mp.mpf(value[k,0].mid().str(90,radius=False))) for k in range(2))
        assert value_error<mp.mpf('1e-65') and all(mp.mpf(r['jacobian'])<mp.mpf('1e-18') and mp.mpf(r['hessian'])<mp.mpf('1e-18') for r in errors)
        rows.append(dict(section=[str(x) for x in point],value_error=mp.nstr(value_error,12),differences=errors))
    # Reject a wrong fixed-point box, not merely a malformed input.
    rejected=False
    try:refine_root(seeds[0],100,shift='1e-12')
    except AssertionError:rejected=True
    assert rejected
    # A whole box whose members can have different diode itineraries.
    f,p,pbox,J=refine_root(seeds[0],100);V,W,d,local,umax,smax=local_block(f,p,pbox,J)
    ua=(-(d+(umax-d)*arb(23)/64)).union(-(d+(umax-d)*arb(24)/64))
    centre=p+V*col([ua.mid(),arb(0)]);r=[ua.rad()+arb('1e-38'),smax+arb('1e-38')];basis=V
    for _ in range(21):centre,basis,r,_,_=propagate(centre,basis,r,100)
    A,_,_=f.matrices(None,0);P,pc=origin_ellipsoid(A);level=arb('.99')/P.inv()[0,0]
    flow=continuous_origin_capture(centre,basis,r,P,level)
    end=np.array([float(arb(v['enclosure']).mid()) for v in flow['centre']]);rho=float(arb(flow['norm_radius']['enclosure']).lower())
    duration=float(arb(flow['time']['enclosure']).mid());Pn=npmat(P)
    def rhs(t,z):
        x,y,w=z;h=-.05*x+(-.95+.05)*(abs(x+1)-abs(x-1))/2
        return [8.4*(y-x-h),x-y+w,-12*y]
    checks=[]
    for a,b in [(0,0),(-1,-1),(-1,1),(1,-1),(1,1)]:
        initial=np.r_[0.,npmat(centre).ravel()+npmat(basis)@np.array([a*float(r[0]),b*float(r[1])])]
        sol=solve_ivp(rhs,[0,duration],initial,method='DOP853',rtol=2e-12,atol=2e-14,max_step=.002)
        error=sol.y[:,-1]-end;distance=float(np.sqrt(error@Pn@error));energy=float(sol.y[:,-1]@Pn@sol.y[:,-1])
        assert sol.success and distance<rho and energy<float(level.lower())
        checks.append(dict(corner=[a,b],distance_in_P_norm=distance,enclosure_radius=rho,energy=energy))
    report=dict(passed=True,derivatives=rows,shifted_root_rejected=rejected,
                grazing_flow_certificate=flow,independent_flow_checks=checks,
                scope='Independent numerical controls, not replacements for interval or analytic proofs.')
    (OUT/'independent-checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Derivative, shifted-root and grazing-flow controls passed',flush=True)


if __name__=='__main__':main()
