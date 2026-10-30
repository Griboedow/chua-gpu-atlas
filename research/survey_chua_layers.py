"""Resolution-aware census of periodic-window candidates on a fixed transect.

Counts describe sampled intervals for prescribed starts and peak tolerances.
They are not a count of all connected components of the exact 3D parameter set.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from chua_gpu import GPU, DEFAULT_SEEDS
from detail_analysis import collect_peaks


def recurrence(values,tolerance=5e-4):
    v=values[np.isfinite(values)][-384:]
    if len(v)<24:return 0
    scale=max(1.,np.max(abs(v)))
    for lag in range(1,min(96,len(v)//4)+1):
        if np.max(abs(v[lag:]-v[:-lag]))<tolerance*scale:return lag
    return 0


def runs(axis,mask):
    changes=np.diff(np.r_[False,mask,False].astype(int))
    return [(int(l),int(r)) for l,r in zip(np.flatnonzero(changes==1),np.flatnonzero(changes==-1))]


def summarize(data):
    beta=data['beta'];period=data['period'];classes=data['classes']
    recurrent=(period>0)&(data['peak_status']==0)
    strict=recurrent&np.isin(classes,[1,2]).all(axis=1)
    windows=[]
    for l,r in runs(beta,recurrent):
        rows=period[l:r]
        windows.append(dict(first_index=l,last_index=r-1,start=float(beta[l]),end=float(beta[r-1]),
                            nodes=r-l,periods=np.unique(rows).tolist(),all_starts_near_zero_nodes=int(strict[l:r].sum()),
                            representative_index=int(l+np.argmax(np.bincount(rows)[rows]))))
    levels=[]
    for stride in [256,64,16,4,1]:
        x=beta[::stride];mask=recurrent[::stride]
        levels.append(dict(stride=stride,nodes=len(x),step=float(x[1]-x[0]),candidate_intervals=len(runs(x,mask)),
                           strict_intervals=len(runs(x,strict[::stride])),periodic_nodes=int(mask.sum())))
    return dict(windows=windows,resolution_levels=levels,
                definitions=dict(candidate='Seed zero has a repeated peak pattern (four repeats, relative tolerance 5e-4) in the later peak collection.',
                                 strict='Additionally all prescribed starts are near-zero/nonpositive or stationary in the preceding finite-time screen.',
                                 interval='Consecutive successful nodes; a one-node interval is retained, not silently removed.'),
                limitations=['Finite samples miss narrower windows.','Long transients can split or create candidate intervals.',
                             'A transect interval is not automatically a distinct globally connected layer.',
                             'Peak period is a diagnostic, not a proof of an attracting cycle.'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/chua-layer-census')
    p.add_argument('--analyze-only',action='store_true')
    p.add_argument('--tight',action='store_true')
    p.add_argument('--label')
    p.add_argument('--beta-low',type=float,default=14.24)
    p.add_argument('--beta-high',type=float,default=14.40)
    p.add_argument('--nodes',type=int,default=10241)
    p.add_argument('--one-seed',action='store_true')
    p.add_argument('--dt',type=float)
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    name=args.label or ('tight' if args.tight else 'line')
    if args.analyze_only:
        raw=np.load(out/(name+'.npz'));d={k:raw[k] for k in raw.files}
        d['period']=np.array([recurrence(v) if s==0 else 0 for v,s in zip(d['peaks'],d['peak_status'])],dtype=np.int16)
        summary=summarize(d)
        np.savez_compressed(out/(name+'-recurrence.npz'),beta=d['beta'],period=d['period'])
        (out/(name+'-recurrence.json')).write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
        print('levels',json.dumps(summary['resolution_levels']),flush=True)
        print('windows',json.dumps(summary['windows']),flush=True);return
    gpu=GPU();beta=np.linspace(args.beta_low,args.beta_high,args.nodes)
    parameters=np.tile([9.,0.,-8/7,-5/7],(len(beta),1));parameters[:,1]=beta
    dt=args.dt or (.00125 if args.tight else .0025)
    transient=2000 if args.tight else 1000
    duration=2400 if args.tight else 1600
    last=[-1]
    def progress(fraction,seconds):
        tenth=int(fraction*10)
        if tenth>last[0]:print('scan',tenth*10,'percent',round(seconds,1),'seconds',flush=True);last[0]=tenth
    started=time.monotonic()
    result=gpu.scan(parameters,seeds=DEFAULT_SEEDS[:1] if args.one_seed else DEFAULT_SEEDS,
                    dt=dt,transient=transient,duration=duration,progress=progress)
    peaks,counts,status=collect_peaks(gpu,parameters,result['final_states'][:,0],duration=1200,dt=dt,capacity=512)
    periods=np.array([recurrence(v) if s==0 else 0 for v,s in zip(peaks,status)],dtype=np.int16)
    arrays=dict(beta=beta,parameters=parameters,period=periods,peaks=peaks,peak_counts=counts,peak_status=status,
                **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    np.savez_compressed(out/(name+'.npz'),**arrays)
    report=dict(device=gpu.info,settings=result['settings'],timing=result['timing'],elapsed_seconds=time.monotonic()-started,
                peak_settings=dict(duration=1200,capacity=512,max_period=96,tolerance=5e-4,min_repeats=4),
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),**summarize(arrays))
    (out/(name+'-summary.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('levels',json.dumps(report['resolution_levels']),flush=True)
    print('windows',json.dumps(report['windows']),flush=True)


if __name__=='__main__':main()
