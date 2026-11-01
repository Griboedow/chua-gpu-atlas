"""Resolve the 55/65-maximum symmetric representatives on their actual beta scale."""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua
from chua_fast_flights import FastFlights
from chua_layer_bifurcations import Flights
from audit_chua_census import independent
from refine_readme_layers import slim

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output;out.mkdir(parents=True,exist_ok=True)
    templates=json.loads((ROOT/'results/chua-layer-atlas/extended-adding.json').read_text());rows=[]
    for template in templates:
        f=FastFlights(Chua(),template['regions'],template['targets'],True)
        slow=Flights(Chua(),template['regions'],template['targets'],True)
        low=template['symmetry']['beta'];high=template['fold']['beta'];q=np.array(template['symmetry']['q']);nodes=[]
        for fraction in np.linspace(.02,.98,25):
            beta=low+float(fraction)*(high-low);q=f.solve(beta,q);r=f.record(beta,q,True);check=slow.record(beta,q)
            error=abs(check['spectral_radius']-r['spectral_radius'])
            assert r['region_violation']<1e-7 and check['residual']<2e-8 and error<2e-5 and check['spectral_radius']<.999
            nodes.append(dict(fraction=float(fraction),representative=slim(r),expm_closure=check['residual'],expm_rho_error=error))
        control=independent(nodes[len(nodes)//2]['representative'])
        assert control['success'] and control['closure']<2e-7 and control['multiplier_error']<2e-5 and abs(control['dominant'])<.999
        rows.append(dict(id='P'+str(template['peak_count']),maxima=template['peak_count'],regions=template['regions'],targets=template['targets'],
            symmetry_beta=low,fold_beta=high,width=high-low,beta_step=float((high-low)*.04),nodes=nodes,independent_control=control))
        print(rows[-1]['id'],'25 attracting cycles','beta spacing',rows[-1]['beta_step'],flush=True)
    result=dict(status='passed',families=rows,scope='Interior samples of the symmetric stable representatives only, at alpha=9, m0=-8/7, m1=-5/7. Boundary estimates come from the published continuation.')
    (out/'thin-intervals.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
