"""Map the two diode slopes at fixed alpha and beta using the existing GPU screen."""
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from chua_gpu import GPU,CLASS_NAMES

CONFIGS={
    'wide':dict(m0=(-1.5,-.9),m1=(-.98,-.10),grid=241,dt=.005,transient=600,duration=1600),
    'zoom':dict(m0=(-8/7-.006,-8/7+.006),m1=(-5/7-.006,-5/7+.006),grid=241,dt=.005,transient=1000,duration=2000),
}


def categories(classes):
    """Keep disagreements between tested starts separate from uniform outcomes."""
    c=np.asarray(classes)
    out=np.full(len(c),3,dtype=np.uint8)
    for value in [0,1,2,4]:out[(c==value).all(axis=1)]=value
    out[(c==5).any(axis=1)]=5
    # Mixed means differing screen labels, not proven distinct attractors.
    out[(c!=c[:,:1]).any(axis=1)&~(c==5).any(axis=1)]=6
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--kind',choices=['wide','zoom','all'],default='all')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True);gpu=GPU()
    for name in CONFIGS if args.kind=='all' else [args.kind]:
        cfg=CONFIGS[name];m0=np.linspace(*cfg['m0'],cfg['grid']);m1=np.linspace(*cfg['m1'],cfg['grid'])
        x,y=np.meshgrid(m0,m1,indexing='xy');params=np.column_stack([np.full(x.size,9.),np.full(x.size,100/7),x.ravel(),y.ravel()])
        last=[-1]
        def progress(f,s):
            k=int(f*10)
            if k>last[0]:print(name,k*10,'percent',round(s,1),'seconds',flush=True);last[0]=k
        result=gpu.scan(params,dt=cfg['dt'],transient=cfg['transient'],duration=cfg['duration'],progress=progress)
        labels=categories(result['classes']);counts={str(k):int((labels==k).sum()) for k in range(7)}
        np.savez_compressed(args.output/(name+'.npz'),parameters=params,m0=m0,m1=m1,labels=labels,
            **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        report=dict(alpha=9.,beta=100/7,grid_shape=[len(m1),len(m0)],settings=result['settings'],device=gpu.info,timing=result['timing'],
            grid_counts=counts,per_start_classes=CLASS_NAMES,labels={**CLASS_NAMES,6:'different classifications across starts'},
            scope='A new two-dimensional parameter section with m0 and m1 varied independently. Finite-time screen at three starts; no proven periodicity, complete basin coverage, or exhaustive family count.')
        (args.output/(name+'.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(name,'complete',counts,flush=True)


if __name__=='__main__':main()
