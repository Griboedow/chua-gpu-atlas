"""GPU control matched to the aligned independent CPU interval at alpha2048."""
import json
import numpy as np
from chua_scaling_gpu import OUT,ScalingGPU,save


if __name__=='__main__':
    source=np.load(OUT/'scout.npz');row=2198
    result=ScalingGPU().scan(source['parameters'][row:row+1],source['final_states'][row:row+1],
                             transient=5,duration=30,rtol=5e-12,max_step=.075,
                             max_attempts=12000000,max_wall_seconds=180.)
    result['source_index']=row
    result['scope']='Matched requested start, original-time window, tolerance and maximum fast step to cpu-2198-aligned. Distinct solvers/control norms can separate chaotic long trajectories.'
    save(OUT/'gpu-2198-matched',result)
    print(dict(status=result['status'].tolist(),lle=result['lle'].tolist(),windows=result['lle_windows'].tolist(),wall_seconds=result['wall_seconds']))
