"""Reproduce exploratory seeds and continuations for the accessibility study."""
import json
import numpy as np
from trace_chua_window import Chua
from chua_layer_bifurcations import from_orbit
from chua_fast_flights import FastFlights
from certify_chua_coexisting_cycles import OUT


def main():
    OUT.mkdir(exist_ok=True)
    model=Chua(alpha=8.4,m0=-1.2,m1=-.05,rtol=2e-12,step=.02)
    sol=model.integrate(12,[0,1,0],150)
    t=sol.t_events[2]; states=sol.y_events[2]
    guess=[*states[-1,1:],t[-1]-t[-2]]
    f,q=from_orbit(model,12,guess,False)
    f=FastFlights(model,f.regions,f.targets,False); q=f.solve(12,q)
    seed=dict(id='H1',regions=f.regions,targets=f.targets,cycle=f.record(12,q,True))
    (OUT/'hidden-cycle-seed.json').write_text(json.dumps(seed,indent=2)+'\n',encoding='utf-8')
    for parameter,values,name in [('m0',np.linspace(-1.2,-.5,141),'hidden-inner-continuation'),
                                  ('m1',np.arange(-.05,.801,.005),'hidden-slope-continuation')]:
        rows=[]; current=q.copy()
        for value in values:
            settings=dict(alpha=8.4,m0=-1.2,m1=-.05); settings[parameter]=float(value)
            branch=FastFlights(Chua(**settings),f.regions,f.targets,False)
            current=branch.solve(12,current); record=branch.record(12,current,True)
            assert record['region_violation']<1e-7 and record['spectral_radius']<1
            rows.append(record)
        (OUT/(name+'.json')).write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
        print(name,len(rows),'numerically continued cycles',flush=True)


if __name__=='__main__':
    main()
