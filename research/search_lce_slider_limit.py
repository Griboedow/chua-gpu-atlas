"""Exact normalized slider and its small-intensity relay limit.

At I>0 this is a coordinate change with q=I/(1+I). The I=0 normalized
equations are a singular rescaled limit, NOT finite oscillation in darkness.
Holding delta=d/q fixed requires retuning the physical illumination threshold.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from lce_slider import SliderGPU,return_map,affine_flow
from check_lce_slider import dop853_crossing


def reference():
    old=np.load('results/lce-slider/wide/scan.npz');p=old['parameters'][2222].copy();s=old['starts'][2222].copy()
    q=p[2]/(1+p[2]);pn=p.copy();pn[3]/=q
    return p,s,pn,s/q,q


def validate(path):
    p,s,pn,sn,q=reference();normal=SliderGPU(normalized=True);original=SliderGPU()
    old=original.scan(p[None],s[None],transient=0,cycles=4,subdivisions=256)
    new=normal.scan(pn[None],sn[None],transient=0,cycles=4,subdivisions=256)
    result=dict(normalized_is_coordinate_change_for_positive_I=True,
        gpu_scaled_four_return_error=float(np.max(abs(old['final_sections'][0]/q-new['final_sections'][0]))),
        gpu_scaled_log_stretch_error=float(np.max(abs(old['log_stretch']-new['log_stretch']))),
        gpu_status=[int(old['status'][0]),int(new['status'][0])],cases=[])
    for intensity in [0.,pn[2],1.]:
        test=pn.copy();test[2]=intensity
        exact,dt,J=return_map(sn,test,256,True)
        middle,t1,J1=dop853_crossing(sn,test,True,1e-12,True)
        end,t2,J2=dop853_crossing(middle,test,False,1e-12,True)
        h=1e-6
        FD=np.column_stack([(return_map(sn+h*e,test,256,True)[0]-return_map(sn-h*e,test,256,True)[0])/(2*h) for e in np.eye(2)])
        cp=sn.copy()
        for _ in range(4):cp=return_map(cp,test,256,True)[0]
        gp=normal.scan(test[None],sn[None],transient=0,cycles=4,subdivisions=256)
        result['cases'].append(dict(intensity=float(intensity),analytic_dop_state_error=float(np.max(abs(end-exact))),
            analytic_dop_time_error=float(abs(dt-t1-t2)),analytic_dop_jacobian_error=float(np.max(abs(J-J2@J1))),
            jacobian_fd_error=float(np.max(abs(J-FD))),gpu_cpu_four_return_error=float(np.max(abs(cp-gp['final_sections'][0]))),
            gpu_status=int(gp['status'][0])))
    result['passed']=(max(result['gpu_status'])==0 and result['gpu_scaled_four_return_error']<1e-8
        and result['gpu_scaled_log_stretch_error']<1e-8
        and all(r['gpu_status']==0 and r['analytic_dop_state_error']<1e-8 and r['analytic_dop_jacobian_error']<1e-7
                and r['jacobian_fd_error']<1e-5 and r['gpu_cpu_four_return_error']<1e-8 for r in result['cases']))
    path.mkdir(parents=True,exist_ok=True);(path/'validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2),flush=True)
    if not result['passed']:raise RuntimeError('Normalized event-map validation failed')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--validate',action='store_true')
    parser.add_argument('--seeds',type=int,default=16);parser.add_argument('--transient',type=int,default=2000)
    parser.add_argument('--cycles',type=int,default=8000);parser.add_argument('--subdivisions',type=int,default=256)
    parser.add_argument('--output',type=Path,default=Path('results/lce-slider/intensity-limit'))
    args=parser.parse_args()
    if args.validate:return validate(args.output)
    p,s,pn,sn,q=reference();levels=np.array([0.,1e-6,1e-5,1e-4,.001,.003,p[2],.01,.03,.1,.3,1.,3.,10.])
    rng=np.random.default_rng(2026102952)
    seeds=np.column_stack((np.sqrt(p[0])*10**rng.uniform(-1,1,args.seeds),rng.uniform(.01,.99,args.seeds)));seeds[0]=sn
    parameters=np.tile(pn,(len(levels)*args.seeds,1));parameters[:,2]=np.repeat(levels,args.seeds)
    starts=np.tile(seeds,(len(levels),1))
    result=SliderGPU(normalized=True).scan(parameters,starts,args.transient,args.cycles,args.subdivisions,
        progress=lambda f,t:print(f'Normalized intensity scan {f:.1%}, {t:.1f}s',flush=True))
    args.output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output/'scan.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    groups=[]
    for j,intensity in enumerate(levels):
        ids=slice(j*args.seeds,(j+1)*args.seeds);valid=result['status'][ids]==0
        retained=valid&(result['lle'][ids]>1e-4)&(result['lle_windows'][ids].min(axis=1)>1e-5)
        groups.append(dict(intensity=float(intensity),count=args.seeds,completed=int(valid.sum()),retained=int(retained.sum()),
            lle=result['lle'][ids].tolist(),windows=result['lle_windows'][ids].tolist(),status=result['status'][ids].tolist()))
    meta={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    meta.update(groups=groups,seed=2026102952,original_reference_parameters=p.tolist(),normalized_reference_parameters=pn.tolist(),
        source_index=2222,physical_caveat='I>0 exact coordinate change; fixed delta requires threshold retuning. I=0 is a rescaled limiting relay model, not nonzero motion at zero light. No device calibration.')
    (args.output/'summary.json').write_text(json.dumps(meta,indent=2))
    print([(g['intensity'],g['retained'],g['completed']) for g in groups],flush=True)


if __name__=='__main__':main()
