"""Sampled step/time refinement and full-field controls for the two 3D sections."""
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from chua_gpu import GPU
from independent_reference import reference_lle
from scan_chua_slope_plane import categories
from scan_chua_slope_volumes import parameters


def dump(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')


def refine(out):
    gpu=GPU();rng=np.random.default_rng(20261030);reports=[]
    for kind in ['alpha','beta']:
        d=np.load(out/(kind+'.npz'));labels=d['labels'];p=parameters(d['m0'],d['m1'],d['control'],kind)
        selected=[]
        for label in np.unique(labels):
            ids=np.flatnonzero(labels==label);selected.extend(rng.choice(ids,min(48,len(ids)),replace=False).tolist())
        remaining=np.setdiff1d(np.arange(len(labels)),selected)
        selected=np.sort(np.r_[selected,rng.choice(remaining,512-len(selected),replace=False)]).astype(int)
        last=[-1]
        def progress(f,s):
            k=int(f*4)
            if k>last[0]:print(kind,'refinement',k*25,'percent',round(s,1),'seconds',flush=True);last[0]=k
        result=gpu.scan(p[selected],dt=.0025,transient=1200,duration=3200,progress=progress)
        new=categories(result['classes']);old=labels[selected];table=np.zeros((7,7),dtype=int);np.add.at(table,(old,new),1)
        np.savez_compressed(out/(kind+'-refinement.npz'),indices=selected,parameters=p[selected],old_labels=old,new_labels=new,
            **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        r=dict(kind=kind,samples=len(selected),label_agreement=float(np.mean(old==new)),per_start_agreement=float(np.mean(d['classes'][selected]==result['classes'])),
            chaos_any_agreement=float(np.mean((d['classes'][selected]==4).any(axis=1)==(result['classes']==4).any(axis=1))),
            confusion_old_rows_new_columns=table.tolist(),numerical_failures=int((result['classes']==5).sum()),settings=result['settings'],timing=result['timing'])
        reports.append(r);print(kind,'label agreement',r['label_agreement'],flush=True)
    dump(out/'refinement.json',dict(device=gpu.info,selection='Up to 48 nodes per original displayed category, uniform fill to 512; RNG seed 20261030.',
        scope='Stratified sensitivity checks, not whole-volume accuracy or a confidence interval. Original complete grids remain the plotted inputs.',volumes=reports))


def cpu(out):
    controls=[]
    for kind in ['alpha','beta']:
        d=np.load(out/(kind+'.npz'));p=parameters(d['m0'],d['m1'],d['control'],kind);column=0 if kind=='alpha' else 1
        for label,side in [(4,'low'),(4,'high'),(2,'middle')]:
            ids=np.flatnonzero(d['labels']==label)
            level=p[ids,column].min() if side=='low' else p[ids,column].max() if side=='high' else d['control'][len(d['control'])//2]
            ids=ids[np.isclose(p[ids,column],level,rtol=0,atol=1e-12)]
            i=int(ids[np.argmin(np.sum(((p[ids,2:]-[-8/7,-5/7])/[.6,.88])**2,axis=1))])
            check=reference_lle(p[i],d['final_states'][i,0],transient=40,duration=160)
            controls.append(dict(kind=kind,selection=side,index=i,seed=0,parameters=p[i].tolist(),initial=d['final_states'][i,0].tolist(),
                original_label=label,original_lle=float(d['lle'][i,0]),reference=check))
            dump(out/'cpu-controls.json',dict(controls=controls,scope='Independent full-field state and variational integration at selected late GPU states. Does not recheck original basin capture or every volume node.'))


def low_alpha(out):
    r=next(r for r in json.loads((out/'cpu-controls.json').read_text())['controls'] if r['kind']=='alpha' and r['selection']=='low')
    p=np.array(r['parameters']);initial=np.array(r['initial']);gpu=GPU()
    result=gpu.scan(p[None],initial_states=initial.reshape(1,1,3),dt=.00125,transient=200,duration=800)
    ref=reference_lle(p,initial,transient=200,duration=800)
    dump(out/'low-alpha-recheck.json',dict(parameters=p.tolist(),initial=initial.tolist(),original_control=r,
        gpu=dict(settings=result['settings'],classification=int(result['classes'][0,0]),lle=float(result['lle'][0,0]),windows=result['lle_windows'][0,0].tolist()),
        reference=ref,scope='Longer matched-start check of the low-alpha chaos candidate whose first independent observation had a nonpositive window. The original volume is retained.'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--stage',choices=['refine','cpu','low-alpha'],required=True)
    a=p.parse_args();{'refine':refine,'cpu':cpu,'low-alpha':low_alpha}[a.stage](a.output)


if __name__=='__main__':main()
