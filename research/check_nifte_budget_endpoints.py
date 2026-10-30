"""Bounded CPU inspection of six saved step-budget NIFTE endpoints.

This samples unresolved long runs; it does not classify the other budget cases.
The existing cycle_check performs DOP853 section checks/Newton shooting and
independent dimensional mass-balance/Radau closure. No GPU work is repeated.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import queue
import time
import numpy as np
from nifte import ROOT, COLUMNS, matrix
from check_nifte_spectrum import cycle_check


def worker(p,endpoint,seconds,result_queue):
    started=time.perf_counter()
    def remaining():
        if time.perf_counter()-started>seconds:
            raise TimeoutError('Per-case CPU deadline reached')
    try:
        result=cycle_check(p,endpoint,matrix(p,len(endpoint)==9),len(endpoint)==9,remaining)
        result['worker_wall_seconds']=time.perf_counter()-started
    except Exception as exc:
        result={'status':'CPU check unresolved','error':f'{type(exc).__name__}: {exc}',
                'worker_wall_seconds':time.perf_counter()-started}
    result_queue.put(result)


def choose():
    chosen=[]
    for name,count in [('scout-5d',None),('random-5d',2),('scout-9d',1)]:
        refined_path=ROOT/(name+'-refined.npz');scout_path=ROOT/(name+'.npz')
        refined=np.load(refined_path);scout=np.load(scout_path)
        unresolved=np.flatnonzero(refined['status']==3)
        if name=='scout-5d' and len(unresolved)!=3:
            raise ValueError('Expected exactly three scout-5d step-budget rows')
        scores=scout['lle_windows'][refined['source_indices'][unresolved]].min(axis=1)
        order=np.argsort(-scores,kind='stable')
        selected=unresolved[order] if count is None else unresolved[order[:count]]
        for row in selected:
            source=int(refined['source_indices'][row])
            if not np.array_equal(refined['parameters'][row],scout['parameters'][source]):
                raise ValueError('Source-index parameter map mismatch')
            chosen.append(dict(survey=name,refined_row=int(row),source_scout_row=source,
                               refined_file=refined_path.name,scout_file=scout_path.name,
                               refined_sha256=hashlib.sha256(refined_path.read_bytes()).hexdigest(),
                               scout_sha256=hashlib.sha256(scout_path.read_bytes()).hexdigest(),
                               original_scout_windows=scout['lle_windows'][source].copy(),
                               original_scout_minimum_window=float(scout['lle_windows'][source].min()),
                               gpu_elapsed_seconds=float(refined['elapsed'][row]),
                               gpu_step_counts=refined['step_counts'][row].copy(),
                               parameters=refined['parameters'][row].copy(),
                               initial_saved_endpoint=refined['final_states'][row].copy()))
    return chosen


def write(records,seconds,started,finished):
    arrays={};rows=[]
    for i,record in enumerate(records):
        row={}
        for key,value in record.items():
            if isinstance(value,np.ndarray):
                target=f'case{i}_{key}';arrays[target]=value
                row[key]={'npz_key':target,'shape':list(value.shape)}
            else:row[key]=value
        rows.append(row)
    out=ROOT/'budget-endpoint-controls'
    np.savez_compressed(str(out)+'.npz',**arrays)
    summary={
        'created_at_utc':datetime.now(timezone.utc).isoformat(),
        'status':'completed requested six-case inspection' if finished else 'partial inspection',
        'selection':'All three scout-5d step-budget endpoints; two random-5d and one scout-9d step-budget endpoints ranked by largest original scout minimum window, with source_indices mapping.',
        'cases_requested':6,'cases_recorded':len(records),'per_case_hard_wall_limit_seconds':seconds,
        'wall_seconds':time.perf_counter()-started,'parameter_columns':COLUMNS,
        'time_unit':'physical seconds','input_status_code':3,
        'method':'Existing check_nifte_spectrum.cycle_check: DOP853 60-second section inspection, Newton Poincare shooting, independent mass_rhs/Radau closure.',
        'scope':'Only these six saved endpoints are sampled. No blanket classification of remaining step-budget cases and no physical vessel-domain certification.',
        'nine_dimensional_limit':'Only resolved dominant transverse multipliers support stability classification. Tiny contracting multipliers and direct monodromy determinants are not certified.',
        'stable_closed_samples':sum(r.get('classification')=='stable periodic endpoint control' for r in records),
        'unresolved_samples':sum(r.get('classification')!='stable periodic endpoint control' for r in records),
        'cases':rows}
    Path(str(out)+'.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(seconds=38.):
    started=time.perf_counter();ctx=mp.get_context('spawn');records=[]
    selected=choose()
    print('selected',[(r['survey'],r['refined_row'],r['source_scout_row']) for r in selected],flush=True)
    for selected_case in selected:
        result_queue=ctx.Queue();begin=time.perf_counter()
        process=ctx.Process(target=worker,args=(selected_case['parameters'],selected_case['initial_saved_endpoint'],seconds-2,result_queue))
        process.start()
        try:
            result=result_queue.get(timeout=seconds)
        except queue.Empty:
            result={'status':'CPU check unresolved','error':'Per-case hard wall-clock deadline; no cycle classification'}
        finally:
            if process.is_alive():process.join(timeout=.2)
            if process.is_alive():process.terminate()
            process.join(timeout=1.)
            result_queue.close()
        row={**selected_case,**result,'case_wall_seconds':time.perf_counter()-begin}
        amplitude=float(np.max(np.ptp(result['states'],axis=0))) if 'states' in result else None
        row['nontrivial_peak_to_peak_max']=amplitude
        stable=(result.get('status')=='closed periodic orbit' and result.get('maximum_transverse_multiplier',2)<1
                and amplitude is not None and amplitude>1e-7)
        row['classification']='stable periodic endpoint control' if stable else 'unresolved by this bounded CPU check'
        if 'section_multipliers' in result:
            row['dominant_transverse_multiplier_modulus']=float(max(abs(result['section_multipliers'])))
            row['multiplier_accuracy_scope']='Dominant multipliers only; tiny modes and determinant may be numerically unresolved.' if len(selected_case['initial_saved_endpoint'])==9 else 'FP64 periodic-orbit multipliers; no rigorous interval enclosure.'
        records.append(row)
        write(records,seconds,started,len(records)==6)
        print(json.dumps({k:row.get(k) for k in ['survey','refined_row','source_scout_row','status','classification','period_seconds','section_closure','independent_mass_radau_closure','dominant_transverse_multiplier_modulus','case_wall_seconds']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--per-case-wall',type=float,default=38.)
    args=parser.parse_args()
    if not 1<args.per_case_wall<=40:raise ValueError('Per-case bound must be between 1 and 40 seconds')
    main(args.per_case_wall)
