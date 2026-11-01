"""Shared numerical helpers for the three follow-up studies.

These routines select and check numerical candidates. They do not replace the
separate interval certificates or imply completeness of a basin/parameter scan.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import root
import mpmath as mp
from trace_chua_window import Chua
from trace_chua_layer_organizer import HalfMap
from chua_fast_flights import FastFlights
from chua_precision_flights import PrecisionFlights

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/chua-followup'
BASE = dict(alpha=9., beta=14.30635909996146, m0=-8/7, m1=-5/7)


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def half_map(parameters, step=.06):
    f = HalfMap()
    f.model = Chua(**{k: parameters[k] for k in ['alpha', 'm0', 'm1']})
    f.bracket_step = step
    return f


def equilibrium_test(p):
    a, b, m = p['alpha'], p['beta'], p['m0']
    coeff = [1+a*(1+m), b+a*m, a*b*(1+m)]
    guards = coeff+[coeff[0]*coeff[1]-coeff[2]]
    return dict(unique=bool(m > -1 and p['m1'] > -1),
                stable=bool(all(v > 0 for v in guards)), hurwitz_guards=guards,
                eigenvalues=[[float(v.real), float(v.imag)] for v in
                             np.linalg.eigvals(Chua(a,m,p['m1']).jacobian(b,True))])


class ParameterPrecisionFlights(PrecisionFlights):
    """Allow one physical parameter, rather than only beta, to vary."""
    def __init__(self, regions, targets, parameters, free='beta', digits=50, spectral=True):
        super().__init__(regions, targets, True, digits)
        self.parameters = {k: mp.mpf(str(v)) for k,v in parameters.items()}
        self.free = free
        self.spectral=spectral
        self.eigen_cache={}

    def exponential(self,A,duration):
        if not self.spectral:return mp.expm(A*duration)
        key=tuple(A)
        if key not in self.eigen_cache:
            ev,V=mp.eig(A);self.eigen_cache[key]=(ev,V,V**-1)
            if len(self.eigen_cache)>30:
                del self.eigen_cache[next(iter(self.eigen_cache))]
        ev,V,W=self.eigen_cache[key]
        E=V*mp.diag([mp.exp(v*duration) for v in ev])*W
        return E.apply(mp.re)

    def evaluate(self, value, q, gradients=True):
        self.parameters[self.free] = value
        self.alpha, self.m0, self.m1 = [self.parameters[k] for k in ['alpha','m0','m1']]
        return super().evaluate(self.parameters['beta'], q, gradients)

    def boundary(self, value, q, target):
        rec = super().boundary(value, q, target)
        rec['coordinate'] = rec.pop('beta')
        rec['free_parameter'] = self.free
        rec['parameters'] = {k: mp.nstr(v,self.dps-5) for k,v in self.parameters.items()}
        rec['parameters'][self.free] = rec['coordinate']
        return rec


def classical_templates():
    registry = json.loads((ROOT/'results/chua-layer-prediction/extended-family-registry.json').read_text())
    peaks = {'S04':15,'S05':25,'S06':35,'S07':45,'P55':55,'P65':65,'P75':75,'P85':85}
    return {peaks[r['id']]:r for r in registry['families'] if r['id'] in peaks}


def parent(parameters,guess=None):
    seed = json.loads((ROOT/'results/chua-layers/affine-central.json').read_text())
    f = FastFlights(Chua(parameters['alpha'],parameters['m0'],parameters['m1']),
                    seed['regions'],seed['targets'],True)
    q = f.solve(parameters['beta'],np.array(seed['boundaries']['-1']['q'] if guess is None else guess))
    r = f.record(parameters['beta'],q,True)
    if r['region_violation'] > 1e-7:
        raise ValueError('Parent itinerary invalid')
    J = np.array(r['section_matrix'])
    ev,V = np.linalg.eig(J)
    order = np.argsort(-abs(ev)); ev,V = ev[order],V[:,order]
    if np.max(abs(ev.imag)) > 1e-10 or not abs(ev[0])>1>abs(ev[1]):
        raise ValueError('Parent is not a real saddle')
    V=V.real
    if V[0,0]<0: V[:,0]*=-1
    return q[:2],ev.real,V,np.linalg.inv(V),r


def approximate_tangency(parameters, free='beta', guess=None, iterations=8, parent_seed=None):
    """Finite-iterate tangent-line approximation; NOT an exact tangency proof."""
    p0=parameters.copy()
    if guess is None: guess=[-.0006075330794225644,p0[free]]
    def quantities(a,t,detail=False):
        p=p0.copy();p[free]=float(t)
        point,ev,V,W,r=parent(p,parent_seed)
        h=half_map(p);state=point+a*V[:,0];tangent=V[:,0].copy()
        words=[];minimum=1e30
        for _ in range(iterations):
            state,J,aux=h.apply(state,p['beta']);tangent=J@tangent
            words.append(aux['itinerary']);minimum=min(minimum,aux['minimum_speed'])
        g=float(W[0]@(state-point));slope=float(W[0]@tangent)
        return g,slope,dict(parameters=p,expansion=float(ev[0]**2),
            period=r['period'],minimum_speed=minimum,words=words,
            stable_displacement=float(W[1]@(state-point)),parent_q=r['q'])
    scale=abs(guess[0])
    fun=lambda w:np.array(quantities(w[0]*scale,w[1])[:2])*[1,scale]
    w=np.array([guess[0]/scale,guess[1]],float)
    for iteration in range(15):
        residual=fun(w)
        if max(abs(residual))<2e-11: break
        steps=[1e-4,2e-7]
        J=np.column_stack([(fun(w+np.eye(2)[j]*h)-fun(w-np.eye(2)[j]*h))/(2*h)
                           for j,h in enumerate(steps)])
        change=np.linalg.solve(J,-residual)
        for factor in [1.,.5,.25,.125,.0625,.03125]:
            trial=w+factor*change
            try: improved=max(abs(fun(trial)))<max(abs(residual))
            except (ValueError,RuntimeError,np.linalg.LinAlgError): improved=False
            if improved: w=trial;break
        else: raise ValueError('Tangency Newton could not improve')
    else: raise ValueError('Tangency Newton exhausted')
    a,t=w[0]*scale,w[1]
    g,slope,details=quantities(a,t)
    if abs(g)>1e-8 or abs(slope*scale)>1e-8:
        raise ValueError(f'Tangency candidate residual {g}, {slope*scale}')
    details.update(starting_coordinate=float(a),coordinate=float(t),free_parameter=free,
        normal_residual=g,scaled_tangent_residual=slope*scale,iterations=iterations,
        solver_success=True,equilibrium=equilibrium_test(details['parameters']))
    return details


def continue_tangency(parameters, free='beta', source=None):
    if source is None: source=approximate_tangency(BASE)
    source=dict(source,coordinate=source['parameters'][free],free_parameter=free)
    start=source['parameters'].copy();target=parameters.copy()
    done=0.;step=.02;record=source;previous=None;path=[]
    while done<1:
        nxt=min(1.,done+step)
        p={k:start[k]+nxt*(target[k]-start[k]) for k in start}
        guess=[record['starting_coordinate'],record['parameters'][free]]
        if previous is not None:
            oldt,oldrec=previous
            for j,key in enumerate(['starting_coordinate','coordinate']):
                guess[j]+=(nxt-done)/(done-oldt)*(record[key]-oldrec[key])
        try: new=approximate_tangency(p,free,guess,parent_seed=record.get('parent_q'))
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            step/=2
            if step<1e-6: raise ValueError(f'Continuation stopped at {done}: {exc}')
            continue
        path.append(dict(fraction=nxt,coordinate=new['coordinate']))
        previous=(done,record);record=new;done=nxt;step=min(.3,step*1.7)
    record['continuation']=path
    return record
