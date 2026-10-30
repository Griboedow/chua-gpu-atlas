"""Compare the faster affine implementation and its boundary derivatives."""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_layer_bifurcations import Flights
from chua_fast_flights import FastFlights

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output;out.mkdir(parents=True,exist_ok=True);checks=[]
    edges=json.loads((ROOT/'results/chua-layer-census/family-edges.json').read_text())
    for item in edges:
        r=item['symmetry'];model=Chua();args=(model,item['regions'],item['targets'],True)
        slow=Flights(*args);fast=FastFlights(*args);q=np.array(r['q']);beta=r['beta'];w=np.r_[q,beta]
        a=slow.evaluate(beta,q,True);b=fast.evaluate(beta,q,True)
        J=fast.boundary_system(w,-1,True);h=1e-6
        numeric=np.column_stack([(fast.boundary_system(w+h*e,-1)-fast.boundary_system(w-h*e,-1))/(2*h) for e in np.eye(len(w))])
        err=abs(J-numeric);relative=float(np.max(err/np.maximum(1.,abs(numeric))))
        checks.append(dict(beta=beta,peak_count=item['peak_count'],closure_difference=float(max(abs(a[0]-b[0]))),
            state_jacobian_difference=float(np.max(abs(a[1]-b[1]))),section_matrix_difference=float(np.max(abs(a[2]-b[2]))),
            state_jacobian_scaled_error=float(np.max(abs(a[1]-b[1])/np.maximum(1.,abs(a[1])))),
            boundary_jacobian_scaled_error=relative,region_violation=b[3]['violation']))
    report=dict(checks=checks,metric_note='State Jacobians are compared relative to each entry with an absolute floor of one because long-flight sensitivities grow. Absolute discrepancies are retained.',passed=all(r['closure_difference']<1e-10 and r['state_jacobian_scaled_error']<1e-9 and
                 r['section_matrix_difference']<2e-6 and r['boundary_jacobian_scaled_error']<5e-4 and r['region_violation']<1e-7 for r in checks))
    (out/'fast-flight-checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(2)


if __name__=='__main__':main()
