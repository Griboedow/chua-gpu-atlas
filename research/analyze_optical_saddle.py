"""Summarize archived controls and check periodic finite-time phase growth."""
import json
from pathlib import Path
import numpy as np
from smooth_operators import rhs_and_jac

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/light-operators/saddle-departure'


def main():
    source=json.loads((OUT/'probe.json').read_text())
    independent=json.loads((OUT/'independent-controls.json').read_text())
    scale=np.array(source['coordinate_scale']);p=np.array(source['parameters'])
    phase=[]
    for row in source['rows']:
        if not row['last_20_sections_near_known_cycle']:continue
        data=np.load(OUT/row['array_file'])
        at100=np.flatnonzero(np.isclose(data['times'],100,rtol=0,atol=1e-10))
        assert len(at100)==1
        f0=rhs_and_jac(data['states'][at100[0]],p)[0]/scale
        f1=rhs_and_jac(data['states'][-1],p)[0]/scale
        tangent=data['final_balanced_state_tangent'][6:]
        phase.append(dict(amplitude=row['amplitude'],sign=row['sign'],
            tangent_mean=row['leading_finite_time_mean'],
            flow_endpoint_mean=float(np.log(np.linalg.norm(f1)/np.linalg.norm(f0))/400),
            final_tangent_absolute_cosine_to_flow=float(abs(np.dot(tangent,f1))/(np.linalg.norm(tangent)*np.linalg.norm(f1)))))
    count=len(independent['rows']);completed=sum(row['completed'] for row in independent['rows'])
    record=dict(independent_plan=4,independent_attempted=count,
        independent_completed=completed,independent_incomplete=count-completed,
        independent_not_started=4-count,independent_wall_seconds=independent['wall_seconds'],
        completed_proximity_classifications_agree=all(row['original_cycle_classification_agrees'] for row in independent['rows'] if row['completed']),
        note='Radau checks known-cycle proximity; it does not independently calculate a Lyapunov exponent. Failed callbacks are not accepted endpoints.',
        periodic_phase_growth=phase)
    (OUT/'analysis.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
