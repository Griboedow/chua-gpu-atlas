"""Ideal two-op-amp realization and a reproducible component-tolerance screen.

Equations (3), (6), (8) of Kuznetsov et al. (2023) supply the component map.
The saturated diode below is independently reconstructed by continuity. Finite
op-amp bandwidth, slew rate and noise are not modeled.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp
from trace_chua_window import Chua
from chua_fast_flights import FastFlights
from certify_chua_coexisting_cycles import OUT


class LossyChua(Chua):
    def __init__(self,gamma=0,**kwargs):
        super().__init__(**kwargs); self.gamma=gamma

    def field(self,beta,u):
        v=super().field(beta,u); v[2]-=self.gamma*u[2]; return v

    def jacobian(self,beta,inner):
        A=super().jacobian(beta,inner); A[2,2]=-self.gamma; return A


class LossyFlights(FastFlights):
    def matrices(self,beta):
        if getattr(self,'_beta',None)!=beta:
            self._beta=beta; self._matrices={}
            for region in set(self.regions):
                A=self.model.jacobian(beta,region==0)
                b=np.array([-self.model.alpha*(self.model.m0-self.model.m1)*region,0,0])
                c=np.linalg.solve(A,-b)
                eig,V=np.linalg.eig(A); W=np.linalg.inv(V)
                self._matrices[region]=(A,c,eig,V,W)
        return self._matrices


def parameters(c):
    R,R1,R2,R3,R4,R5,R6,C1,C2,L,RL=[c[k] for k in ['R','R1','R2','R3','R4','R5','R6','C1','C2','L','RL']]
    K1=1+R3/R4; K2=1+R5/R6
    return dict(alpha=C2/C1,beta=C2*R*R/L,gamma=C2*R*RL/L,
                m0=R*(-R3/(R4*R1)-R5/(R6*R2)),m1=R*(1/R1-R5/(R6*R2))),K1/K2


def main():
    R=1000.; K1=40.; K2=2.; m0=-1.2; m1=-.05
    delta=m1-m0
    components=dict(R=R,R1=R*K1/delta,R2=R*(K2-1)/(delta/K1-m1),
                    R3=(K1-1)*1000,R4=1000.,R5=(K2-1)*1000,R6=1000.,
                    C1=1e-7/8.4,C2=1e-7,L=1e-7*R*R/12,RL=.5)
    nominal,limit=parameters(components)
    seed=json.loads((OUT/'hidden-cycle-seed.json').read_text())
    rng=np.random.default_rng(20261101)
    records=[]
    for i in range(129):
        c={k:v*(1 if i==0 else rng.uniform(.99,1.01)) for k,v in components.items()}
        pars,saturation=parameters(c); beta=pars.pop('beta')
        model=LossyChua(**pars)
        f=LossyFlights(model,seed['regions'],seed['targets'],False)
        try:
            q=f.solve(beta,np.array(seed['cycle']['q']))
            rec=f.record(beta,q,True)
            # Dense samples locate maxima for a numerical saturation screen.
            xmax=0.
            for arc in rec['arcs']:
                A,eq,eig,V,W=f.matrices(beta)[arc['region']]
                ts=np.linspace(0,arc['duration'],401)
                x=eq[0]+np.real(np.exp(ts[:,None]*eig[None,:])*(V[0]*(W@(np.array(arc['start'])-eq)))[None,:]).sum(axis=1)
                xmax=max(xmax,float(max(abs(x))))
            passed=rec['region_violation']<1e-7 and rec['spectral_radius']<1 and xmax<saturation
            records.append(dict(index=i,parameters=dict(pars,beta=beta),components=c,period=rec['period'],
                spectral_radius=rec['spectral_radius'],xmax_sampled=xmax,outer_breakpoint=saturation,
                passed=bool(passed),q=q.tolist()))
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
            records.append(dict(index=i,parameters=dict(pars,beta=beta),components=c,passed=False,error=str(error)))
    selected={0,min((r for r in records if r['passed']),key=lambda r:r['outer_breakpoint']-r['xmax_sampled'])['index'],
              max((r for r in records if r['passed']),key=lambda r:r['spectral_radius'])['index']}
    independent=[]
    for index in sorted(selected):
        rec=records[index]; p=rec['parameters']; comp=rec['components']
        m2=comp['R']*(1/comp['R1']+1/comp['R2'])
        def physical_field(t,u):
            x,y,z=u
            h=p['m1']*x+(p['m0']-p['m1'])*np.clip(x,-1,1)
            h+=(m2-p['m1'])*np.sign(x)*max(abs(x)-rec['outer_breakpoint'],0)
            return [p['alpha']*(y-x-h),x-y+z,-p['beta']*y-p['gamma']*z]
        u=np.r_[0.,rec['q'][:2]]
        sol=solve_ivp(physical_field,[0,rec['period']],u,method='Radau',rtol=2e-12,atol=2e-13,max_step=.008)
        error=float(max(abs(sol.y[:,-1]-u)))
        assert sol.success and error<2e-8
        independent.append(dict(index=index,method='Radau',five_segment_diode=True,closure=error,passed=True))
    result=dict(nominal_components_SI=components,nominal_parameters=dict(nominal,beta=12),independent_checks=independent,
        supply_voltage=9,nominal_inner_breakpoint_voltage=9/K1,nominal_outer_breakpoint=limit,
        tolerance_fraction=.01,random_seed=20261101,samples=records,
        passed_count=sum(r['passed'] for r in records),
        source='https://doi.org/10.1007/s11071-022-08078-y',
        necessary_constraints=['m1-m0 = R*K1/R1 > 0','m0 < 0','gamma >= 0 for nonnegative inductor resistance'],
        scope='An ideal component construction and finite numerical perturbation screen. This does not certify the entire tolerance box, hiddenness after perturbation, or operation of a real device.')
    (OUT/'circuit-constraints.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Component realizations passing the numerical cycle/saturation checks:',result['passed_count'],'/',len(records),flush=True)


if __name__=='__main__':
    main()
