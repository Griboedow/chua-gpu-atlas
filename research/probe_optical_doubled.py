"""Probe a genuine doubled orbit near the first optical cycle flip candidate."""
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from continue_optical_cycle import close,shoot,ROOT,OUT
from smooth_operators import rhs_and_jac


def run():
    source=json.loads((OUT/'from-narrow.json').read_text())['records']
    records=[];wall=time.perf_counter()
    for width in [.0064,.0067]:
        anchor=min(source,key=lambda r:abs(r['parameters'][4]-width))
        p=np.array(anchor['parameters']);p[4]=width
        base=close(p,np.array(anchor['section_state']),anchor['period'],rtol=3e-12,max_step=.005)
        y=np.array(base['section_state']);T=base['period'];end,M=shoot(p,y,T,rtol=3e-12,max_step=.005)
        f=rhs_and_jac(end,p)[0];S=((np.eye(6)-np.outer(f,np.eye(6)[0])/f[0])@M)[1:,1:]
        vals,vec=np.linalg.eig(S);j=np.argmin(abs(vals+1))
        assert abs(vals[j].imag)<1e-9
        scale=np.array([1,np.sqrt(p[0]),1,1,np.sqrt(p[0]),1]);v=vec[:,j].real/scale[1:];v/=np.linalg.norm(v)
        for amplitude in [.005,.02,.05,.1,.002]:
            if time.perf_counter()-wall>150:break
            seed=y.copy();seed[1:]+=amplitude*v*scale[1:]
            row=dict(width=width,seed_amplitude=amplitude,base_period=T)
            try:
                orbit=close(p,seed,2*T,rtol=3e-12,max_step=.005)
                z=np.array(orbit['section_state']);period=orbit['period']
                half,_=shoot(p,z,period/2,rtol=3e-12,max_step=.005)
                error=float(max(abs((half-z)/scale)))
                row.update(orbit=orbit,half_period_nonclosure=error,
                           classification='nontrivial doubled-orbit candidate' if error>1e-5 else 'repeated parent cycle')
                if error>1e-5:
                    def ev(t,q):return q[0]-p[2]
                    ev.direction=1
                    sol=solve_ivp(lambda t,q:rhs_and_jac(q,p)[0],(0,period),z,
                        method='Radau',rtol=1e-12,atol=1e-14,max_step=.0025,events=ev)
                    events=sol.t_events[0];inside=events[(events>1e-6)&(events<period-1e-6)]
                    row.update(independent_radau_success=bool(sol.success),
                               independent_radau_closure=float(max(abs((sol.y[:,-1]-z)/scale))),
                               upward_crossings_per_period=1+len(inside))
            except Exception as exc:row.update(classification='unresolved shooting attempt',error=f'{type(exc).__name__}: {exc}')
            records.append(row)
            (OUT/'doubled-probes.json').write_text(json.dumps(dict(attempts=records,wall_seconds=time.perf_counter()-wall),indent=2)+'\n',encoding='utf8')
            print(json.dumps({k:v for k,v in row.items() if k!='orbit'}),flush=True)
            if row['classification']=='nontrivial doubled-orbit candidate':break


if __name__=='__main__':run()
