"""Adaptive full-field checks for the unstable-origin prediction samples.

This resolves outcomes left ambiguous by fixed-step breakpoint integration.
The original fixed-step results remain available as numerical controls.
"""
import json
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

OUT=Path(__file__).resolve().parents[1]/'results/chua-accessibility'


def integrate(task):
    m0,initial,root=task
    def field(t,u):
        x,y,z=u
        h=-.05*x+(m0+.05)*(abs(x+1)-abs(x-1))/2
        return [8.4*(y-x-h),x-y+z,-12*y]
    def section(t,u):return u[0]
    section.direction=1
    sol=solve_ivp(field,[0,1000],initial,method='DOP853',events=section,
                  rtol=2e-10,atol=2e-13,max_step=.12)
    last=sol.y_events[0][-1,1:] if len(sol.y_events[0]) else np.array([np.nan,np.nan])
    distance=float(np.linalg.norm(last-root))
    return dict(m0=m0,initial=list(initial),section=last.tolist(),distance=distance,
                reached_cycle=bool(sol.success and distance<1e-6),method='DOP853',duration=1000)


def refine(task):
    record,root=task
    if record['reached_cycle']:return record
    m0=record['m0']
    def field(t,u):
        x,y,z=u; h=-.05*x+(m0+.05)*(abs(x+1)-abs(x-1))/2
        return [8.4*(y-x-h),x-y+z,-12*y]
    def section(t,u):return u[0]
    section.direction=1
    sol=solve_ivp(field,[0,400],[0,*record['section']],method='DOP853',events=section,
                  rtol=2e-12,atol=2e-14,max_step=.04)
    last=sol.y_events[0][-1,1:]
    distance=float(np.linalg.norm(last-root))
    return dict(record,first_horizon_distance=record['distance'],first_horizon_reached_cycle=False,
                section=last.tolist(),distance=distance,reached_cycle=bool(sol.success and distance<1e-6),
                additional_time_from_last_section=400,refined_rtol=2e-12)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refine-existing',action='store_true')
    args=parser.parse_args()
    prior=json.loads((OUT/'stability-basin-check.json').read_text())
    rows=json.loads((OUT/'hidden-inner-continuation.json').read_text())
    tasks=[]
    for pred in prior['records']:
        if pred['prediction']['origin_stable']:continue
        m0=float(pred['prediction']['m0']); row=min(rows,key=lambda r:abs(r['m0']-m0))
        for radius in prior['radii']:
            for direction in prior['directions']:
                tasks.append((m0,radius*np.array(direction),np.array(row['q'][:2])))
    with ProcessPoolExecutor(max_workers=4) as pool:
        if args.refine_existing:
            results=json.loads((OUT/'adaptive-prediction-check.json').read_text())['records']
        else:
            results=list(pool.map(integrate,tasks))
        first_horizon_count=sum(r['reached_cycle'] for r in results)
        refined_tasks=[(r,np.array(min(rows,key=lambda row:abs(row['m0']-r['m0']))['q'][:2])) for r in results]
        results=list(pool.map(refine,refined_tasks))
    report=dict(records=results,reached_cycle_count=sum(r['reached_cycle'] for r in results),sample_count=len(results),
                first_horizon_reached_cycle_count=first_horizon_count,
                scope='Adaptive checks of these finite-radius starts, not a proof of self-excitation at every parameter on the unstable side.')
    (OUT/'adaptive-prediction-check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(report['reached_cycle_count'],'/',report['sample_count'],'adaptive trajectories reached the continued cycle',flush=True)


if __name__=='__main__':
    main()
