"""Tighten an unresolved long-cycle Radau replay without discarding its trial."""
import argparse,json
from pathlib import Path
import numpy as np
from trace_chua_window import Chua


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output;path=out/'extended-adding.json';records=json.loads(path.read_text())
    for r in records:
        if r['status']!='unresolved' or 'core' not in r:continue
        core=r['core'];u=np.r_[0.,core['q'][:2]];model=Chua(rtol=2.5e-14,step=.005,method='Radau')
        sol=model.integrate(core['beta'],u,core['period']);closure=float(max(abs(sol.y[:,-1]-u)))
        r.setdefault('radau_refinements',[]).append(dict(rtol=2.5e-14,atol=2.5e-15,max_step=.005,closure=closure))
        v=r['independent_variational']
        if closure<2e-7 and v['success'] and abs(v['dominant'])<.999 and v['closure']<2e-7 and v['multiplier_error']<2e-5 and r['peak_count']==r['maxima_requested']:
            r['prior_failure']=r.pop('error',None);r['status']='verified_attracting_cycle'
        print(r['maxima_requested'],closure,r['status'],flush=True)
    path.write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
