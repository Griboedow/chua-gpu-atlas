"""Refine sampled slope-plane labels and independently inspect representative outcomes."""
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from chua_gpu import GPU
from detail_analysis import collect_peaks
from independent_reference import reference_lle
from scan_chua_slope_plane import categories
from survey_chua_layers import recurrence


def dump(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')


def refine(out):
    gpu=GPU();rng=np.random.default_rng(20261030);reports=[]
    for name in ['wide','zoom']:
        d=np.load(out/(name+'.npz'));meta=json.loads((out/(name+'.json')).read_text());labels=d['labels']
        selected=[]
        for label in np.unique(labels):
            ids=np.flatnonzero(labels==label);selected.extend(rng.choice(ids,min(64,len(ids)),replace=False).tolist())
        rest=np.setdiff1d(np.arange(len(labels)),selected)
        selected=np.sort(np.r_[selected,rng.choice(rest,512-len(selected),replace=False)]).astype(int)
        cfg=meta['settings'];last=[-1]
        def progress(f,s):
            k=int(f*4)
            if k>last[0]:print(name,'refinement',k*25,'percent',round(s,1),'seconds',flush=True);last[0]=k
        result=gpu.scan(d['parameters'][selected],dt=cfg['dt']/2,transient=cfg['transient']*2,duration=cfg['duration']*2,progress=progress)
        new=categories(result['classes']);old=labels[selected]
        table=np.zeros((7,7),dtype=int);np.add.at(table,(old,new),1)
        np.savez_compressed(out/(name+'-refinement.npz'),indices=selected,parameters=d['parameters'][selected],old_labels=old,new_labels=new,
            **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        r=dict(name=name,samples=len(selected),label_agreement=float(np.mean(old==new)),
            per_start_agreement=float(np.mean(d['classes'][selected]==result['classes'])),
            chaos_any_agreement=float(np.mean((d['classes'][selected]==4).any(axis=1)==(result['classes']==4).any(axis=1))),
            confusion_old_rows_new_columns=table.tolist(),numerical_failures=int((result['classes']==5).sum()),settings=result['settings'],timing=result['timing'])
        reports.append(r);print(name,'agreement',r['label_agreement'],flush=True)
    dump(out/'refinement.json',dict(device=gpu.info,selection='At most 64 samples per original displayed category, then uniform fill to 512, RNG seed 20261030.',
        scope='Stratified sensitivity checks, not a confidence interval or a complete finer map. Main images retain the original grids.',planes=reports))


def peaks(out):
    gpu=GPU();d=np.load(out/'zoom.npz');ids=np.flatnonzero(d['labels']==2)
    values,counts,status=collect_peaks(gpu,d['parameters'][ids],d['final_states'][ids,0],duration=1200,dt=.0025,capacity=512)
    periods=np.array([recurrence(v) if s==0 else 0 for v,s in zip(values,status)],dtype=np.int16)
    np.savez_compressed(out/'zoom-peaks.npz',indices=ids,periods=periods,peaks=values,counts=counts,status=status)
    unique,n=np.unique(periods,return_counts=True)
    dump(out/'zoom-peaks.json',dict(nodes=len(ids),peak_pattern_counts={str(int(k)):int(v) for k,v in zip(unique,n)},
        dt=.0025,duration=1200,capacity=512,last_peaks_checked=384,maximum_lag=96,minimum_repeats=4,relative_tolerance=.0005,
        scope='Seed zero at grid nodes with all three preceding screen labels equal to regular candidate. Repeated peak counts are diagnostics, not closed-cycle proofs or family identities. Zero means unresolved.'))
    print('peak patterns',dict(zip(unique,n)),flush=True)


def cpu(out):
    wide=np.load(out/'wide.npz');zoom=np.load(out/'zoom.npz');jobs=[]
    # Interior examples chosen by proximity, without altering the source grids.
    for name,d,label,where in [('chaotic',zoom,4,[-8/7,-5/7]),('regular',wide,2,[-1.29,-.70]),
                              ('equilibrium',wide,1,[-1.30,-.95]),('mixed',wide,6,[-1.18,-.14])]:
        eligible=np.flatnonzero(d['labels']==label)
        distance=((d['parameters'][eligible,2:]-where)/np.ptp(d['parameters'][:,2:],axis=0))**2
        i=int(eligible[np.argmin(distance.sum(axis=1))]);seeds=range(3) if name=='mixed' else [0]
        for seed in seeds:
            p=d['parameters'][i];initial=d['final_states'][i,seed]
            ref=reference_lle(p,initial,transient=40,duration=160)
            jobs.append(dict(name=name,source='zoom' if name=='chaotic' else 'wide',index=i,seed=seed,parameters=p.tolist(),initial=initial.tolist(),
                original_class=int(d['classes'][i,seed]),original_lle=float(d['lle'][i,seed]),reference=ref))
            dump(out/'cpu-controls.json',dict(controls=jobs,scope='Full-field DOP853 state and variational checks started from saved late GPU states. These do not recheck basin capture from the original initial conditions or the whole grid.'))


def cycles(out):
    from trace_chua_window import Chua
    from audit_chua_census import independent
    d=np.load(out/'zoom.npz');peaks=np.load(out/'zoom-peaks.npz');records=[]
    for count in [5,7]:
        ids=peaks['indices'][(peaks['periods']==count)&np.isclose(d['parameters'][peaks['indices'],3],-5/7,atol=1e-12,rtol=0)]
        i=int(ids[len(ids)//2]);p=d['parameters'][i];model=Chua(alpha=p[0],m0=p[2],m1=p[3])
        late=model.integrate(p[1],d['final_states'][i,0],120)
        times,states=late.t_events[2],late.y_events[2]
        for lag in range(1,min(20,len(times)//2)):
            if np.max(abs(states[-1]-states[-1-lag]))<2e-5:
                guess=[*states[-1,1:],times[-1]-times[-1-lag]];break
        else:raise RuntimeError('No recurrent section for the selected peak pattern')
        v=model.shoot(p[1],guess);r=model.record(p[1],v)
        rho=max(abs(complex(*m)) for m in r['multipliers'])
        assert r['residual']<2e-8 and rho<.999 and len(r['maxima'])==count
        control=independent({**r,'q':v[:2].tolist(),'full_multipliers':r['multipliers']})
        assert control['success'] and control['closure']<2e-7 and control['multiplier_error']<2e-5 and abs(control['dominant'])<.999
        records.append(dict(index=i,peak_count=count,spectral_radius=rho,cycle=r,independent_control=control))
        print('closed cycle',count,'period',r['period'],'rho',rho,'residual',r['residual'],flush=True)
    dump(out/'closed-cycles.json',dict(status='passed',cycles=records,scope='Two selected periodic orbits closed by full-field shooting; affine transverse multipliers cross-checked with full-field variational integration. This does not assign identities from the earlier 32-family catalogue.'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--stage',choices=['refine','peaks','cpu','cycles'],required=True);a=p.parse_args()
    {'refine':refine,'peaks':peaks,'cpu':cpu,'cycles':cycles}[a.stage](a.output)


if __name__=='__main__':main()
