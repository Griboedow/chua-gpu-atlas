"""Summarize retained NIFTE evidence without launching new trajectories."""
import json
from pathlib import Path
import numpy as np
from nifte import ROOT


def main():
    records=[]
    for stem in ['scout-5d','scout-9d','focus-5d','random-5d','random-9d']:
        for suffix in ['', '-refined']:
            name=stem+suffix;path=ROOT/(name+'.npz')
            if not path.exists():continue
            a=np.load(path);meta=json.loads(path.with_suffix('.json').read_text());ok=a['status']==0
            counts={str(k):int(np.sum(a['status']==k)) for k in range(5)}
            weak=ok&np.all(a['lle_windows']>0,axis=1)
            record=dict(dataset=name,lanes=len(ok),status_counts=counts,all_four_windows_positive=int(weak.sum()),
                largest_completed_mean=float(np.max(a['lle'][ok])) if np.any(ok) else None,
                largest_completed_minimum_window=float(np.max(a['lle_windows'][ok].min(axis=1))) if np.any(ok) else None,
                passed_declared_threshold=int(np.sum(a['retained_candidate' if suffix else 'preliminary_candidate'])),
                completed_with_nonpositive_reconstructed_absolute_pressure=int(np.sum(ok&(a['stats'][:,0]<=-1))),
                pressure_diagnostic_scope='Observation minima only; necessary condition, not a full physical-admissibility test',
                elapsed_min=float(np.min(a['elapsed'])),elapsed_max=float(np.max(a['elapsed'])),settings=meta['settings'],wall_seconds=meta['wall_seconds'])
            records.append(record)
    controls=[]
    for name in ['source-5d-cycle','source-9d-cycle','source-5d-unstable-cycle','cpu-row215-cycle']:
        a=np.load(ROOT/(name+'.npz'));p=a['parameters'];s=a['states'];meta=json.loads((ROOT/(name+'.json')).read_text())
        volume=1.09*p[-2]*p[9]-p[-2]*(p[10]*s[:,1]+p[11]*s[:,2])
        controls.append(dict(dataset=name,period_seconds=meta.get('period_seconds',meta.get('period_s')),minimum_absolute_pressure_Pa=float(np.min(p[-2]*(1+s[:,0]))),
            conditional_volume_min_m3=float(volume.min()),conditional_volume_max_m3=float(volume.max()),
            displacer_displacement_range_m=(np.array([s[:,1].min(),s[:,1].max()])*p[-2]/(621*9.81)).tolist(),
            physical_scope='Assumes source fluid density, gravity and cylindrical-volume reconstruction; not a calibrated vessel-operating envelope'))
    out=dict(surveys=records,cycle_diagnostics=controls,
        scout_trajectories=sum(r['lanes'] for r in records if not r['dataset'].endswith('-refined')),
        followup_trajectories=sum(r['lanes'] for r in records if r['dataset'].endswith('-refined')),
        interpretation='No validated chaos candidate in this finite survey; unresolved budgets and small finite-window positives are not a global no-chaos result')
    (ROOT/'summary.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf-8')
    for r in records:print({k:r[k] for k in ['dataset','lanes','status_counts','largest_completed_mean','all_four_windows_positive','passed_declared_threshold']})


if __name__=='__main__':main()
