"""Compare interpolated sheets with fresh solves between every mesh node."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.interpolate import RectBivariateSpline
from refine_chua_sheets import advance as advance_core
from refine_chua_satellites import advance as advance_stable


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--families',nargs='*')
    args=p.parse_args();out=args.output;summary=[]
    for path in sorted(out.glob('*-sheet.json')):
        r=json.loads(path.read_text())
        if args.families and r['id'] not in args.families:continue
        if not r['complete']:continue
        aa=np.array(r['alpha']);mm=np.array(r['m0']);nodes={tuple(n['index']):n for n in r['nodes']};bounded='fold' in next(iter(nodes.values()))
        keys=['core','fold','symmetry'] if bounded else ['core']
        splines={key:RectBivariateSpline(aa,mm,np.array([[nodes[i,j][key]['beta'] for j in range(len(mm))] for i in range(len(aa))]),kx=3,ky=3) for key in keys}
        template=dict(regions=r['regions'],targets=r['targets'],symmetric=r.get('symmetric',True),target=r.get('target_multiplier'))
        checks=[]
        for i in range(len(aa)-1):
            for j in range(len(mm)-1):
                alpha=(aa[i]+aa[i+1])/2;m0=(mm[j]+mm[j+1])/2
                predicted=float(splines['core'](alpha,m0)[0,0])
                if not 14.23<=predicted<=14.41:continue
                exact=(advance_core if bounded else advance_stable)(template,nodes[i,j],alpha,m0)
                errors={key:float(splines[key](alpha,m0)[0,0]-exact[key]['beta']) for key in keys}
                checks.append(dict(alpha=float(alpha),m0=float(m0),errors=errors,
                                   width=exact.get('width'),core_beta=exact['core']['beta']))
        maximum=max(abs(e) for c in checks for e in c['errors'].values())
        item=dict(id=r['id'],nodes_per_axis=len(aa),checks=checks,max_absolute_error=maximum,
                  max_error_over_core_width=max(abs(c['errors']['core'])/c['width'] for c in checks) if bounded else None)
        summary.append(item);(out/'interpolation-checks.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(r['id'],'checks',len(checks),'max error',maximum,'width ratio',item['max_error_over_core_width'],flush=True)


if __name__=='__main__':main()
