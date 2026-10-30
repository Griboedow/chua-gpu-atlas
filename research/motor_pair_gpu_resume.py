"""Exact continuation of a wall-censored motor-pair GPU checkpoint.

No state, tangent, accumulated logarithm, observation boundary or step count is
reset. The input NPZ stays immutable. Completed lanes are copied unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pyopencl as cl

from motor_pair_gpu import MotorPairGPU, from_balanced, save


def resume(source, output, wall=600.):
    source=Path(source);metadata=json.loads(source.with_suffix('.json').read_text())
    original=np.load(source);result={k:original[k].copy() for k in original.files}
    result.update({k:v for k,v in metadata.items() if k not in result})
    settings=metadata['settings'];transient=settings['transient_scaled'];duration=settings['observation_scaled']
    gpu=MotorPairGPU();params=result['parameters'];n=len(params)
    p=np.zeros((n,4));p[:,:2]=params;p[:,2]=np.minimum(1.,params[:,0])
    state=result['raw_state'].copy();clock=result['raw_clock'].copy();logs=result['raw_log_sums'].copy()
    stats=result['stats'].copy();status=result['status'].copy();counts=result['step_counts'].copy()
    selected=status==4
    if not selected.any():raise RuntimeError('No wall-censored lane to continue')
    completed_before=(status==0).copy();status[selected]=0
    arrays=[state,p,clock,logs,stats,status,counts]
    buffers=[cl.Buffer(gpu.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
    began=time.perf_counter();last=began;calls=0
    while True:
        gpu.kernel(gpu.queue,(n,),None,*buffers,np.int32(n),np.int32(128),np.float64(transient),np.float64(duration),
            np.float64(settings['rtol']),np.float64(settings['atol']),np.float64(settings['step_factor']),
            np.int32(settings['max_steps_per_lane'])).wait();calls+=1
        if calls%8==0:
            cl.enqueue_copy(gpu.queue,clock,buffers[2]).wait();cl.enqueue_copy(gpu.queue,status,buffers[5]).wait()
            active=(status==0)&(clock[:,0]<transient+duration)
            if time.perf_counter()-last>25:
                print(json.dumps(dict(complete=int(np.sum((status==0)&~active)),remaining=int(active.sum()),
                    elapsed_scaled=clock[:,0].tolist(),wall=round(time.perf_counter()-began,1))),flush=True);last=time.perf_counter()
            if not active.any():break
            if time.perf_counter()-began>wall:
                status[active]=4;cl.enqueue_copy(gpu.queue,buffers[5],status).wait();break
    for array,buffer in zip(arrays,buffers):cl.enqueue_copy(gpu.queue,array,buffer).wait()
    for key,array in [('raw_state',state),('raw_clock',clock),('raw_log_sums',logs),('stats',stats),('step_counts',counts)]:
        if not np.array_equal(array[completed_before],original[key][completed_before]):
            raise RuntimeError('A previously completed lane changed')
    rates=p[:,2];windows=logs/(duration/4)*rates[:,None];windows[status!=0]=np.nan
    result.update(raw_state=state,raw_clock=clock,raw_log_sums=logs,stats=stats,status=status,step_counts=counts,
        balanced_states=state[:,:3].copy(),balanced_tangents=state[:,3:6].copy(),
        final_states=np.array([from_balanced(x,pa) for x,pa in zip(state[:,:3],params)]),
        elapsed_scaled=clock[:,0].copy(),elapsed_original=clock[:,0]/rates,lle_windows=windows,lle=windows.mean(axis=1),
        mean_divergence_original=stats[:,8]/duration*rates,
        continuation=dict(input_npz=source.as_posix(),input_npz_exact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            input_metadata_canonical_lf_sha256=hashlib.sha256(source.with_suffix('.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
            source_ids_resumed=result['source_ids'][selected].tolist(),previously_completed_lanes_bitwise_unchanged=True,
            no_reset='States, tangents, clocks, steps, logarithms and integration settings carried exactly',
            script_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
            wall_seconds=time.perf_counter()-began,wall_budget_seconds=wall),
        wall_seconds=metadata['wall_seconds']+time.perf_counter()-began)
    save(output,result)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--wall',type=float,default=600.)
    args=parser.parse_args();resume(args.source,args.output,args.wall)
