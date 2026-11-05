"""Search symmetry partners and catalogue overlaps for non-split joint links."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
for key,folder in [('MPLCONFIGDIR','matplotlib'),('NUMBA_CACHE_DIR','numba')]:
    path=ROOT/'.publication'/folder;path.mkdir(parents=True,exist_ok=True);os.environ[key]=str(path)
TASK_TEMP=ROOT.parent/'.agent-tmp';TASK_TEMP.mkdir(parents=True,exist_ok=True)
os.environ['TEMP']=os.environ['TMP']=str(TASK_TEMP)
os.environ['PIP_CACHE_DIR']='D:/_TMP/python-cache/pip'
import argparse,json,itertools,random,time,hashlib
import numpy as np
from spherogram import Link
from chua_link_diagrams import diagram,clean,jones_A
import mpmath as mp
from chua_precision_flights import PrecisionFlights
from reconstruct_chua_knot_orbits import sample
OUT=ROOT/'results/chua-link-search'


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def save(path,value):path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')
def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n') if path.suffix in ['.json','.md','.py'] else path.read_bytes()).hexdigest()


def split_screen(curves,seed=0):
    link,r=diagram(curves,seed);random.seed(43061+seed)
    stages=[]
    for mode in ['basic','global']:
        link.simplify(mode,type_III_limit=500)
        stages.append(dict(mode=mode,crossings=len(link.crossings),free=link.unlinked_unknot_components,components=len(link.link_components)))
    pieces=link.split_link_diagram() if link.crossings else []
    sizes=[len(p.link_components) for p in pieces]+[1]*link.unlinked_unknot_components
    assert sum(sizes)==len(curves)
    split=all(n==1 for n in sizes)
    r.update(reductions=stages,simplified_pd_code=link.PD_code(),simplified_crossings=len(link.crossings),
        final_free_components=link.unlinked_unknot_components,split_piece_component_counts=sizes,
        status='split_union' if split else 'unresolved',scope='Numerical polygons; split components may themselves be knots.')
    if not split:
        from chua_link_invariants import fox_nonsplit_test
        values=fox_nonsplit_test(link)['evaluations']
        r['alexander_minors']=values
        if any(v!=0 for v in values.values()):r['status']='non_split_by_alexander'
        elif np.any(r['oriented_linking_matrix']):r['status']='non_split_by_linking_number'
    return r


def symmetry(digits,subdivisions,seeds,ids):
    registry=read(ROOT/'results/chua-layer-atlas/family-registry.json')['families']
    chosen=[r for r in registry if not r['symmetric'] and (not ids or r['id'] in ids)]
    for row in chosen:
        id=row['id'];path=ROOT/'results/chua-knot-atlas'/f'{id}-d{digits}-n{subdivisions}.npz'
        if not path.exists():print('Missing input',path.name,flush=True);continue
        p=np.load(path)['xyz']
        for seed in seeds:
            target=OUT/'symmetry'/f'{id}-d{digits}-n{subdivisions}-s{seed}.json'
            target.parent.mkdir(exist_ok=True)
            if target.exists():continue
            start=time.monotonic();r=split_screen([p,-p],seed)
            r.update(id=id,source=path.relative_to(ROOT).as_posix(),source_sha256=digest(path),
                digits=digits,subdivisions=subdivisions,elapsed_seconds=time.monotonic()-start)
            save(target,r);print(id,r['crossings'],'->',r['simplified_crossings'],r['status'],r['split_piece_component_counts'],r.get('alexander_minors'),flush=True)


def node_key(r):return tuple(round(r[k],10) for k in ['alpha','m0','m1'])


def catalogue():
    sheets={p.stem[:-6]:read(p) for p in (ROOT/'results/chua-layer-atlas').glob('*-sheet.json')}
    groups={};intervals={}
    for id,s in sheets.items():
        for n in s['nodes']:
            if n.get('status')!='continued':continue
            key=node_key(n['representative']);groups.setdefault(key,{})[id]=n
            if 'fold' in n and 'symmetry' in n:
                lo,hi=sorted([n['fold']['beta'],n['symmetry']['beta']]);intervals.setdefault(key,{})[id]=(lo,hi)
    candidates=[]
    for key,nodes in groups.items():
        spans=intervals.get(key,{})
        trials=[(id,n['representative']['beta']) for id,n in nodes.items()]
        for (a,span),(b,other) in itertools.combinations(spans.items(),2):
            lo=max(span[0],other[0]);hi=min(span[1],other[1])
            if hi>lo:trials.append((a,(lo+hi)/2))
        for source,beta in trials:
            included=sorted(set([source]+[id for id,(lo,hi) in spans.items() if lo<beta<hi]))
            if len(included)<2:continue
            depth=min([min(beta-spans[id][0],spans[id][1]-beta)/(spans[id][1]-spans[id][0]) for id in included if id in spans] or [0])
            witness=dict(ids=included,parameters=dict({k:nodes[source]['representative'][k] for k in ['alpha','m0','m1']},beta=beta),score=depth,
                inputs={id:dict(regions=sheets[id]['regions'],targets=sheets[id]['targets'],symmetric=sheets[id].get('symmetric',True),q=nodes[id]['representative']['q'],source_beta=nodes[id]['representative']['beta']) for id in included})
            candidates.append(witness)
    best={}
    for r in candidates:
        key=tuple(r['ids'])
        if key not in best or r['score']>best[key]['score']:best[key]=r
    result=dict(candidate_count=len(candidates),combinations=len(best),selected=sorted(best.values(),key=lambda r:(-len(r['ids']),-r['score'])),all_candidates=candidates)
    save(OUT/'catalogue-candidates.json',result)
    for r in result['selected']:print(r['ids'],r['parameters'],'depth',r['score'],flush=True)


def reconstruct_case(case,digits,subdivisions,seeds):
    candidates=read(OUT/'catalogue-candidates.json')['selected'];row=candidates[case]
    folder=OUT/'overlaps'/('-'.join(row['ids']));folder.mkdir(parents=True,exist_ok=True)
    curves=[];names=[]
    for id in row['ids']:
        entry=row['inputs'][id];stem=f'{id}-d{digits}-n{subdivisions}';target=folder/(stem+'.json')
        if not target.exists():
            f=PrecisionFlights(entry['regions'],entry['targets'],entry['symmetric'],digits)
            params={k:mp.mpf(str(v)) for k,v in row['parameters'].items()}
            f.alpha=params['alpha'];f.m0=params['m0'];f.m1=params['m1'];beta=params['beta']
            q=mp.matrix([mp.mpf(str(v)) for v in entry['q']]);errors=[]
            for _ in range(25):
                ev=f.evaluate(beta,list(q),False);err=max(abs(v) for v in ev['residual']);errors.append(mp.nstr(err,15))
                if err<mp.mpf(10)**(-(digits-20)):break
                q+=mp.lu_solve(ev['jacobian'],-ev['residual'])
                assert min(q[2:])>0
            else:raise RuntimeError('Shooting failed at common parameters')
            eigen=mp.polyroots([1,-ev['trace'],ev['determinant']],maxsteps=1000)
            radius=max(abs(v) for v in eigen)**(2 if entry['symmetric'] else 1)
            assert radius<1
            factor=2 if entry['symmetric'] else 1
            root=dict(id=id,**{k:mp.nstr(v,digits) for k,v in params.items()},gamma='0',q=[mp.nstr(v,digits-5) for v in q],
                period=mp.nstr(factor*sum(q[2:]),digits-5),regions=f.regions,targets=f.targets,symmetric=entry['symmetric'],
                trace=mp.nstr(ev['trace'],30),determinant=mp.nstr(ev['determinant'],30),spectral_radius=float(radius),
                residual=mp.nstr(err,15),working_digits=digits,newton_residuals=errors,source_case=case)
            xyz,t,diagnostic=sample(f,root,ev,subdivisions);root['sampling']=diagnostic
            assert diagnostic['maximum_sampled_region_violation']<mp.mpf('1e-25')
            np.savez_compressed(folder/(stem+'.npz'),xyz=xyz,t=t);save(target,root)
            print(folder.name,id,'closed','rho',float(radius),flush=True)
        xyz=np.load(folder/(stem+'.npz'))['xyz'];curves.append(xyz);names.append(id)
        if not entry['symmetric']:curves.append(-xyz);names.append(id+'-inverted')
    for seed in seeds:
        target=folder/f'joint-d{digits}-n{subdivisions}-s{seed}.json'
        if target.exists():continue
        r=split_screen(curves,seed);r.update(ids=names,parameters=row['parameters'],digits=digits,subdivisions=subdivisions)
        save(target,r);print(folder.name,'joint',r['crossings'],'->',r['simplified_crossings'],r['status'],r.get('alexander_minors'),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['symmetry','catalogue','overlaps']);p.add_argument('--digits',type=int,default=60);p.add_argument('--subdivisions',type=int,default=96);p.add_argument('--seeds',type=int,nargs='+',default=[0]);p.add_argument('--ids',nargs='*');p.add_argument('--cases',type=int,nargs='*');a=p.parse_args();OUT.mkdir(exist_ok=True)
    if a.stage=='symmetry':symmetry(a.digits,a.subdivisions,a.seeds,a.ids)
    elif a.stage=='catalogue':catalogue()
    else:
        for case in a.cases if a.cases is not None else range(len(read(OUT/'catalogue-candidates.json')['selected'])):
            try:reconstruct_case(case,a.digits,a.subdivisions,a.seeds)
            except Exception as e:
                print('FAILED case',case,type(e).__name__,str(e),flush=True)
                save(OUT/f'failed-case{case}-d{a.digits}-n{a.subdivisions}.json',dict(case=case,error_type=type(e).__name__,error=str(e)))


if __name__=='__main__':main()
