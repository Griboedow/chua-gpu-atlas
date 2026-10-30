"""Extend the m0-m1 section along alpha or beta with genuine 3D GPU scans."""
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from chua_gpu import GPU,CLASS_NAMES
from scan_chua_slope_plane import categories


def parameters(m0,m1,control,kind):
    z,y,x=np.meshgrid(control,m1,m0,indexing='ij')
    p=np.column_stack([np.full(x.size,9.),np.full(x.size,100/7),x.ravel(),y.ravel()])
    p[:,0 if kind=='alpha' else 1]=z.ravel()
    return p


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--kind',choices=['alpha','beta','all'],default='all');p.add_argument('--resume',action='store_true')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);gpu=GPU()
    source=ROOT/'results/chua-slope-plane/wide.npz';plane=np.load(source)
    # The central layer uses exactly the same coordinates as every fourth
    # row and column of the earlier 241-square section.
    m0=plane['m0'][::4];m1=plane['m1'][::4]
    for kind in ['alpha','beta'] if a.kind=='all' else [a.kind]:
        if a.resume and (a.output/(kind+'.npz')).exists() and (a.output/(kind+'.json')).exists():continue
        control=np.linspace(7,11,41) if kind=='alpha' else 100/7+np.linspace(-4,4,41)
        params=parameters(m0,m1,control,kind);last=[-1]
        def progress(f,s):
            k=int(f*20)
            if k>last[0]:print(kind,k*5,'percent',round(s,1),'seconds',flush=True);last[0]=k
        result=gpu.scan(params,dt=.005,transient=600,duration=1600,progress=progress)
        labels=categories(result['classes']);shape=(len(control),len(m1),len(m0))
        arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray) and k!='stats'}
        np.savez_compressed(a.output/(kind+'.npz'),m0=m0,m1=m1,control=control,labels=labels,**arrays)
        center=result['classes'].reshape(*shape,3)[20]
        previous=plane['classes'].reshape(241,241,3)[::4,::4]
        c_lle=result['lle'].reshape(*shape,3)[20];p_lle=plane['lle'].reshape(241,241,3)[::4,::4]
        finite=np.isfinite(c_lle)&np.isfinite(p_lle)
        report=dict(varied_control=kind,fixed_parameter={'beta':100/7} if kind=='alpha' else {'alpha':9.},
            axes_order=['control','m1','m0'],shape=list(shape),nodes=len(params),initial_conditions_per_node=3,
            ranges=dict(m0=[float(m0[0]),float(m0[-1])],m1=[float(m1[0]),float(m1[-1])],control=[float(control[0]),float(control[-1])]),
            settings=result['settings'],device=gpu.info,timing=result['timing'],
            counts={str(k):int((labels==k).sum()) for k in range(7)},labels={**CLASS_NAMES,6:'different classifications across starts'},
            central_plane_reproduction=dict(nodes=int(len(m0)*len(m1)),per_start_class_agreement=float(np.mean(center==previous)),
                max_finite_lle_difference=float(np.max(abs(c_lle[finite]-p_lle[finite]))),source=source.relative_to(ROOT).as_posix()),
            compacting='Raw accumulation statistics are omitted; classes, four-window exponents, statuses, outer-region flags and final states are retained. Parameter tuples are reconstructed from the stored axes.',
            scope='Three independent coordinates varied at each node. Finite-time classification of three starts; no exhaustive layer count or proven continuous boundary.')
        (a.output/(kind+'.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(kind,'complete',report['counts'],'central agreement',report['central_plane_reproduction']['per_start_class_agreement'],flush=True)


if __name__=='__main__':main()
