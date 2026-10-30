"""Independent Radau departure controls; no trajectories or source results overwritten."""
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from smooth_operators import rhs_and_jac

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/light-operators/saddle-departure'


def run():
    if (OUT/'independent-controls.json').exists():
        raise RuntimeError('Independent control already exists')
    source = json.loads((OUT/'probe.json').read_text())
    p = np.array(source['parameters']); scale = np.array(source['coordinate_scale'])
    stable = np.array(source['stable_section_balanced'])
    begin = time.perf_counter(); rows = []
    def field(t, y):
        if time.perf_counter()-begin > 420:
            raise TimeoutError('Independent CPU control wall budget')
        return rhs_and_jac(y*scale, p)[0]/scale
    def jac(t, y):
        return rhs_and_jac(y*scale, p)[1]/scale[:,None]*scale[None,:]
    def section(t, y):
        return y[0]-p[2]
    section.direction = 1
    for i, original in enumerate(source['rows']):
        row = dict(amplitude=original['amplitude'], sign=original['sign'], completed=False)
        try:
            start = np.array(original['start_balanced'])
            ts = np.linspace(0, 20, 1001)
            short = {}
            for method in ['DOP853', 'Radau']:
                opts = dict(jac=jac) if method=='Radau' else {}
                sol=solve_ivp(field,(0,20),start,method=method,rtol=2e-12,atol=2e-14,
                              max_step=.01,t_eval=ts,**opts)
                if not sol.success: raise RuntimeError(sol.message)
                short[method] = sol.y
            row['short_dop_radau_max_difference'] = float(np.max(abs(short['DOP853']-short['Radau'])))
            sol=solve_ivp(field,(0,500),start,method='Radau',jac=jac,rtol=2e-10,atol=2e-12,
                          max_step=.03,events=section,t_eval=np.linspace(0,500,2501))
            if not sol.success: raise RuntimeError(sol.message)
            distances = np.linalg.norm(sol.y_events[0]-stable,axis=1)
            periods = np.diff(sol.t_events[0])
            close = bool(len(distances)>=20 and np.all(distances[-20:]<1e-6))
            row.update(completed=True, completed_time=float(sol.t[-1]),
                final_20_section_distances=distances[-20:].tolist(),
                final_20_return_intervals=periods[-20:].tolist(),
                last_20_sections_near_known_cycle=close,
                original_cycle_classification_agrees=close==original['last_20_sections_near_known_cycle'],
                thermal_min=float((sol.y*scale[:,None])[[2,5]].min()),
                thermal_max=float((sol.y*scale[:,None])[[2,5]].max()))
            np.savez_compressed(OUT/f'independent-radau-{i}.npz',times=sol.t,states=sol.y.T*scale,
                section_times=sol.t_events[0],section_states_balanced=sol.y_events[0],
                section_distance_to_stable=distances)
        except (TimeoutError, RuntimeError) as error:
            row['error']=str(error)
        rows.append(row)
        record=dict(rows=rows,wall_seconds=time.perf_counter()-begin,
            settings=dict(method='Radau',duration=500,rtol=2e-10,atol=2e-12,max_step=.03,
                          short_both_solvers_rtol=2e-12,short_duration=20,short_max_step=.01),
            scope='Independent finite departure classification, not proof of an invariant basin boundary or chaotic attractor.')
        (OUT/'independent-controls.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(row),flush=True)
        if not row['completed']: break


if __name__=='__main__': run()
