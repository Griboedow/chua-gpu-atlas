"""Verify four attracting cycles in an overlap of asymmetric layer branches."""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_layer_bifurcations import from_orbit
from chua_fast_flights import FastFlights
from audit_chua_census import independent

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output;out.mkdir(parents=True,exist_ok=True)
    alpha=9.06;m0=-8/7+.006;beta=14.3458;model=Chua(alpha=alpha,m0=m0,rtol=2e-12,step=.02)
    first=json.loads((ROOT/'results/chua-layers/extended.json').read_text())['surfaces']
    second=json.loads((ROOT/'results/chua-layer-census/second-surfaces.json').read_text())
    records=[];arrays={}
    for j,rows in enumerate([first,second]):
        edge=min(rows,key=lambda r:abs(r['alpha']-alpha)+abs(r['m0']-m0))['period_doubling']
        family=FastFlights(model,[a['region'] for a in edge['arcs']],[round(a['end'][0]) for a in edge['arcs'][:-1]],False)
        q=np.array(edge['q'])
        for b in np.linspace(edge['beta'],beta,41)[1:]:q=family.solve(b,q)
        for partner in [0,1]:
            if partner:
                inv=model.integrate(beta,-np.r_[0.,q[:2]],sum(q[2:])*1.01)
                state=inv.y_events[2][0];f,qq=from_orbit(model,beta,[*state[1:],sum(q[2:])],False)
                f=FastFlights(model,f.regions,f.targets,False)
            else:f,qq=family,q
            r=f.record(beta,qq,True);u=np.r_[0.,qq[:2]];T=r['period'];orbit=model.integrate(beta,u,T)
            radau=Chua(alpha=alpha,m0=m0,rtol=2e-12,step=.02,method='Radau').integrate(beta,u,T)
            var=independent(r);closure=float(max(abs(radau.y[:,-1]-u)))
            name=f'F{j+1}{"b" if partner else "a"}'
            records.append(dict(id=name,cycle=r,regions=f.regions,targets=f.targets,peak_count=len(orbit.t_events[3]),
                                radau_closure=closure,variational=var))
            t=np.linspace(0,T,4001);arrays[name+'_t']=t;arrays[name+'_xyz']=orbit.sol(t).T
            print(name,T,r['spectral_radius'],closure,var,flush=True)
    points=np.array([r['cycle']['q'][:2] for r in records]);dist=np.linalg.norm(points[:,None]-points[None,:],axis=2)
    passed=all(r['cycle']['spectral_radius']<.999 and r['cycle']['region_violation']<1e-7 and r['radau_closure']<2e-7
        and r['variational']['success'] and abs(r['variational']['dominant'])<.999 and r['variational']['multiplier_error']<2e-5 for r in records)
    report=dict(alpha=alpha,beta=beta,m0=m0,m1=-5/7,records=records,section_distances=dist.tolist(),passed=passed)
    (out/'four-cycles.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    np.savez_compressed(out/'four-cycles.npz',**arrays)
    print('passed',passed,'minimum section separation',dist[np.triu_indices(4,1)].min(),flush=True)


if __name__=='__main__':main()
