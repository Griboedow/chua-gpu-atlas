"""Short rectangular-pulse law from the already known basin separator.

The coefficient uses a normal to the continued saddle stable manifold, not a
fit to pulse data. Finite-iterate normals are numerical approximations.
"""
import json
import numpy as np
from study_chua_trigger import setup,PARS
from chua_research_followup import half_map,OUT,dump


def main():
    p,pc,P,level,e,V,W,*_=setup();f=half_map(PARS)
    Astar=.34061011644046296
    A=np.array([[-8.4*.05,8.4,0],[1,-1,1],[0,-12,0]])
    ey=np.array([0.,1.,0.]);rows=[]
    for n in [4,6,8]:
        state=np.array([Astar,0.]);D=np.eye(2)
        for _ in range(n):state,J,_=f.apply(state,12.);D=J@D
        g=W[0]@D;g=g/g[0]
        normal=np.r_[(g[0]+12*g[1])/8.4,g]
        coefficient=-Astar*(normal@(A@A@ey))/24
        rows.append(dict(half_returns=n,normal=normal.tolist(),quadratic_coefficient=float(coefficient)))
    c=rows[-1]['quadratic_coefficient']
    pulses=json.loads((OUT/'excitation-plane/pulses.json').read_text());comparisons=[]
    for row in pulses['rows']:
        if row['duration']>.15:continue
        tr=next(t for t in row['transitions'] if t['left']==1)
        tau=row['duration'];prediction=Astar+c*tau*tau
        comparisons.append(dict(duration=tau,measured_area=tr['impulse'],predicted_area=prediction,
            measured_quadratic_coefficient=(tr['impulse']-Astar)/(tau*tau),
            quadratic_increment_relative_error=abs((tr['impulse']-Astar)/(c*tau*tau)-1)))
    result=dict(voltage_ray_threshold=Astar,normal_approximations=rows,comparisons=comparisons,
        area_law='a*tau = Astar + c*tau^2 + O(tau^4)',
        cost_law='J = Astar^2/tau + 2*Astar*c*tau + O(tau^3)',
        coefficient_formula='c = -Astar * n dot (A0^2 e_y) / (24 * n dot e_y)',
        scope='Analytic short-pulse expansion on the smooth local separator; numerical coefficient from finite-iterate normals. Comparison with existing pulse data, not a newly frozen prospective experiment.')
    dump(OUT/'excitation-plane/short-pulse-law.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
