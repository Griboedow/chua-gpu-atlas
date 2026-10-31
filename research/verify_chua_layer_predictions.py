"""Search the frozen target itineraries, retaining double and high precision trials."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from trace_chua_window import Chua
from chua_fast_flights import FastFlights
from chua_precision_flights import PrecisionFlights
from predict_chua_layers import ROOT, OUT, parent_at


def insert_return(edge, regions):
    beta=float(edge['beta']);q=np.array(edge['q'],dtype=float)
    parent=parent_at(beta)
    central=json.loads((ROOT/'results/chua-layers/affine-central.json').read_text())
    hreg=central['regions'];hq=np.array(parent['q']);hd=hq[2:]
    preg=hreg[:-1]+[0]+[-r for r in hreg[1:]]
    pd=np.r_[hd[:-1],hd[-1]+hd[0],hd[1:]]
    ptarg=[preg[j] or preg[j+1] for j in range(len(preg)-1)]
    pf=FastFlights(Chua(),preg,ptarg,False)
    pq=pf.solve(beta,np.r_[hq[:2],pd])
    model=Chua(rtol=2e-12,step=.02)
    sol=model.integrate(beta,[0.,*q[:2]],float(edge['period'])/2)
    ts=sol.t_events[2];us=sol.y_events[2]
    inds=np.flatnonzero((ts>1e-6)&(ts<float(edge['period'])/2-1e-6))
    i=inds[np.argmin(np.linalg.norm(us[inds,1:]-pq[:2],axis=1))]
    t=ts[i];ends=np.cumsum(q[2:]);starts=np.r_[0.,ends[:-1]];k=int(np.searchsorted(ends,t))
    assert regions[k]==0 and preg[0]==preg[-1]==0
    durations=np.r_[q[2:2+k],t-starts[k]+pq[2],pq[3:-1],pq[-1]+ends[k]-t,q[3+k:]]
    newregions=regions[:k]+preg+regions[k+1:]
    targets=[newregions[j] or newregions[j+1] for j in range(len(newregions)-1)]
    return newregions,targets,np.r_[q[:2],durations],float(t)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=OUT)
    p.add_argument('--digits',type=int,default=55)
    p.add_argument('--through',type=int,default=85)
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    frozen=out/'predictions-frozen.json';forecast=json.loads(frozen.read_text())
    old=json.loads((ROOT/'results/chua-layer-atlas/extended-adding.json').read_text())[-1]
    edge=old['fold'];regions=old['regions'];records=[]
    result_path=out/'search-results.json'
    if result_path.exists():
        raise SystemExit('Search results already exist; use a new output directory for a separate trial.')
    start=time.monotonic()
    def save():
        result_path.write_text(json.dumps(dict(forecast_sha256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
            records=records,elapsed_seconds=time.monotonic()-start),indent=2)+'\n',encoding='utf-8',newline='\n')
    for prediction in forecast['predictions']:
        if prediction['peaks']>args.through:break
        rec=dict(peaks_requested=prediction['peaks'],status='unresolved',trials=[])
        records.append(rec);save()
        try:
            newregions,targets,guess,t=insert_return(edge,regions)
            rec.update(regions=newregions,targets=targets,insertion_time=t)
            word=''.join('R' if r==1 else 'L' for r in newregions if r)
            assert word==prediction['half_outer_word']
            f=FastFlights(Chua(),newregions,targets,True)
            beta=prediction['fold_beta_corrected'];q=guess
            try:
                beta,q=f.boundary(beta,guess,1)
                rr=f.record(beta,q,True)
                rec['trials'].append(dict(method='float64 affine eigendecomposition',fold=rr))
                print(prediction['peaks'],'float fold',beta,'residual',rr['residual'],flush=True)
            except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError) as exc:
                rec['trials'].append(dict(method='float64 affine eigendecomposition',error=str(exc)))
                print(prediction['peaks'],'float trial failed',exc,flush=True)
            save()
            high=PrecisionFlights(newregions,targets,True,args.digits)
            fold=high.boundary(beta,q,1);rec['fold']=fold;save()
            print(prediction['peaks'],'precision fold',fold['beta'],'residual',fold['residual'],flush=True)
            symmetry=high.boundary(fold['beta'],fold['q'],-1);rec['symmetry']=symmetry;save()
            core=high.boundary(fold['beta'],fold['q'],.5);rec['core']=core
            import mpmath as mp
            rec['symmetric_core_width']=mp.nstr(mp.mpf(fold['beta'])-mp.mpf(symmetry['beta']),args.digits-5)
            rec['status']='high_precision_candidate'
            rec['forecast_errors']=dict(
                fold_absolute=float(abs(mp.mpf(fold['beta'])-mp.mpf(str(prediction['fold_beta_corrected'])))),
                fold_in_core_widths=float(abs(mp.mpf(fold['beta'])-mp.mpf(str(prediction['fold_beta_corrected'])))/mp.mpf(rec['symmetric_core_width'])),
                width_relative=float(abs(mp.mpf(rec['symmetric_core_width'])/mp.mpf(str(prediction['symmetric_core_width']))-1)),
                period_absolute=abs(float(fold['period'])-prediction['fold_period']))
            edge=fold;regions=newregions;save()
            print(prediction['peaks'],'width',rec['symmetric_core_width'],rec['forecast_errors'],flush=True)
        except (RuntimeError,ValueError,OverflowError,np.linalg.LinAlgError,ZeroDivisionError) as exc:
            rec['error']=str(exc);save();print('STOP',exc,flush=True);break
    print('Saved search-results.json',flush=True)


if __name__=='__main__':main()
