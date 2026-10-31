"""Check arbitrary-precision shooting against known roots and finite differences."""
import json
import mpmath as mp
import numpy as np
from chua_precision_flights import PrecisionFlights
from chua_fast_flights import FastFlights
from trace_chua_window import Chua
from predict_chua_layers import ROOT,OUT


def main():
    rows=json.loads((ROOT/'results/chua-layer-atlas/family-registry.json').read_text())['families']
    result=[]
    for name in ['S08','S04']:
        row=next(r for r in rows if r['id']==name);core=row['core']
        f=PrecisionFlights(row['regions'],row['targets'],True,55)
        w=mp.matrix([mp.mpf(str(v)) for v in core['q']]+[mp.mpf(str(core['beta']))])
        residual,J=f.system(w,mp.mpf('.5'),True)
        errors=[];h=mp.mpf('1e-15')
        for j in [0,1,2,len(w)//2,len(w)-1]:
            a=w.copy();b=w.copy();a[j]+=h;b[j]-=h
            col=(f.system(a,mp.mpf('.5'))-f.system(b,mp.mpf('.5')))/(2*h)
            errors.append(float(max(abs(col[i]-J[i,j]) for i in range(len(w)))/max(1,max(abs(v) for v in col))))
        refined=f.boundary(core['beta'],core['q'],.5)
        fast=FastFlights(Chua(),row['regions'],row['targets'],True)
        r=fast.record(float(refined['beta']),np.array(refined['q'],float),True)
        record=dict(id=name,maximum_scaled_jacobian_error=max(errors),
            beta_discrepancy=abs(float(refined['beta'])-core['beta']),
            high_precision_residual=refined['residual'],float_closure=r['residual'],
            float_multiplier=r['map_multipliers'][0][0])
        assert max(errors)<1e-20
        assert record['beta_discrepancy']<1e-10
        assert r['residual']<1e-9 and abs(r['map_multipliers'][0][0]-.5)<1e-7
        result.append(record);print(record,flush=True)
    (OUT/'precision-implementation-checks.json').write_text(json.dumps(dict(passed=True,checks=result),indent=2)+'\n')


if __name__=='__main__':main()
