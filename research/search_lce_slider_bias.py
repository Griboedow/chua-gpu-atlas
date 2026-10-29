"""Offset dependence of the normalized I=0 relay limit.

At I=0 reflection (U,V,Y,delta)->(-U,-V,1-Y,1-delta) is exact.
Only delta<=1/2 is surveyed. No claim of chaos at zero physical illumination.
"""
import json
from pathlib import Path
import numpy as np
from lce_slider import SliderGPU,flow
from search_lce_slider_limit import reference


def main():
    p,s,pn,sn,q=reference();rng=np.random.default_rng(2026102953)
    errors=[]
    for _ in range(128):
        state=rng.uniform(-1,1,3);state[2]=(state[2]+1)/2;delta=rng.uniform(.001,.999)
        parameter=np.array([pn[0],pn[1],0.,delta]);mirror=parameter.copy();mirror[3]=1-delta
        reflected=np.array([-state[0],-state[1],1-state[2]])
        f=flow(state,parameter,state[0]>0,True);g=flow(reflected,mirror,reflected[0]>0,True)
        errors.append(float(max(abs(f+g))))
    levels=np.unique(np.r_[np.linspace(.005,.5,101),pn[3]])
    seeds=np.column_stack((np.sqrt(pn[0])*10**rng.uniform(-1,1,8),rng.uniform(.01,.99,8)));seeds[0]=sn
    parameters=np.tile(pn,(len(levels)*8,1));parameters[:,2]=0;parameters[:,3]=np.repeat(levels,8)
    starts=np.tile(seeds,(len(levels),1))
    result=SliderGPU(normalized=True).scan(parameters,starts,transient=2000,cycles=8000,subdivisions=256,
        progress=lambda f,t:print(f'Relay bias scan {f:.1%}, {t:.1f}s',flush=True))
    out=Path('results/lce-slider/relay-bias');out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'scan.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    groups=[]
    for j,delta in enumerate(levels):
        ids=slice(8*j,8*(j+1));valid=result['status'][ids]==0
        retained=valid&(result['lle'][ids]>1e-4)&(result['lle_windows'][ids].min(axis=1)>1e-5)
        groups.append(dict(delta=float(delta),completed=int(valid.sum()),retained=int(retained.sum()),
            lle=result['lle'][ids].tolist(),windows=result['lle_windows'][ids].tolist(),status=result['status'][ids].tolist()))
    meta={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    meta.update(seed=2026102953,groups=groups,maximum_reflection_identity_error=max(errors),
        symmetry='At I=0, reflection of U,V,Y about (0,0,1/2) exchanges delta and 1-delta; unsurveyed half not labeled as data.',
        caveat='Normalized limiting relay model; finite parameter/start ensemble, not proof of asymptotic chaos or basin-volume estimates.')
    (out/'summary.json').write_text(json.dumps(meta,indent=2))
    print('Positive groups',[(g['delta'],g['retained']) for g in groups if g['retained']],flush=True)


if __name__=='__main__':main()
