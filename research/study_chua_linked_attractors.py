"""Fixed-parameter reconstruction and joint topology of four Chua attractors."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.publication/matplotlib'))
os.environ['TEMP']=os.environ['TMP']=str(ROOT/'.publication/tmp')
Path(os.environ['TEMP']).mkdir(parents=True,exist_ok=True)
import argparse,json,hashlib,itertools,time
import mpmath as mp
import numpy as np
from chua_precision_flights import PrecisionFlights
from reconstruct_chua_knot_orbits import sample
from chua_link_diagrams import classify,clean,polygon_gauss,canonical_unoriented
from spherogram import Link
OUT=ROOT/'results/chua-linked-attractors'


def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,r):p.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n') if p.suffix in ['.json','.md','.py'] else p.read_bytes()).hexdigest()


def reconstruct(digits,subdivisions):
    source=read(ROOT/'results/chua-layer-atlas/four-cycles.json')
    for row in source['records']:
        id=row['id'];stem=f'{id}-d{digits}-n{subdivisions}';path=OUT/(stem+'.json')
        if path.exists():continue
        f=PrecisionFlights(row['regions'],row['targets'],False,digits)
        f.alpha=mp.mpf(453)/50;f.m0=-mp.mpf(8)/7+mp.mpf(3)/500;f.m1=-mp.mpf(5)/7;beta=mp.mpf(71729)/5000
        q=mp.matrix([mp.mpf(str(v)) for v in row['cycle']['q']]);history=[]
        for _ in range(20):
            ev=f.evaluate(beta,list(q),False);error=max(abs(v) for v in ev['residual']);history.append(mp.nstr(error,15))
            if error<mp.mpf(10)**(-(digits-20)):break
            change=mp.lu_solve(ev['jacobian'],-ev['residual']);q+=change
            assert min(q[2:])>0
        else:raise RuntimeError('Fixed-parameter root failed')
        eigen=mp.polyroots([1,-ev['trace'],ev['determinant']],maxsteps=1000)
        r=dict(id=id,beta=mp.nstr(beta,digits),q=[mp.nstr(v,digits-5) for v in q],period=mp.nstr(sum(q[2:]),digits-5),
            alpha='453/50',m0='-8/7+3/500',m1='-5/7',gamma='0',regions=row['regions'],targets=row['targets'],symmetric=False,
            trace=mp.nstr(ev['trace'],30),determinant=mp.nstr(ev['determinant'],30),spectral_radius=float(max(abs(v) for v in eigen)),
            residual=mp.nstr(error,15),newton_residuals=history,working_digits=digits)
        xyz,t,diag=sample(f,r,ev,subdivisions);r['sampling']=diag;r['minimum_duration']=float(min(q[2:]));r['source_sha256']=digest(ROOT/'results/chua-layer-atlas/four-cycles.json')
        assert r['spectral_radius']<1 and diag['maximum_sampled_region_violation']<1e-30
        np.savez_compressed(OUT/(stem+'.npz'),xyz=xyz,t=t);save(path,r)
        print(stem,'samples',len(xyz),'residual',r['residual'],'rho',r['spectral_radius'],flush=True)


def control_curves(n=800):
    t=np.linspace(0,2*np.pi,n,endpoint=False);c=np.cos(t);s=np.sin(t);zero=np.zeros(n)
    borromean=[np.c_[2*c,s,zero],np.c_[zero,2*c,s],np.c_[s,zero,2*c]]
    hopf=[np.c_[c,s,zero],np.c_[1+c,zero,s]]
    unlink=[np.c_[c+4*i,s,zero] for i in range(3)]
    return dict(borromean=borromean,hopf=hopf,unlink=unlink)


def controls():
    rows=[];standard=canonical_unoriented(Link('L6a4'))
    for name,curves in control_curves().items():
        for size in range(1,len(curves)+1):
            for subset in itertools.combinations(range(len(curves)),size):
                for seed in [0,2]:
                    r,link=classify([curves[i] for i in subset],seed)
                    expected='borromean' if name=='borromean' and size==3 else 'hopf' if name=='hopf' and size==2 else 'unlink'
                    passed=(r['status']=='unlink') if expected=='unlink' else (abs(r['oriented_linking_matrix'][0][1])==1) if expected=='hopf' else canonical_unoriented(link)==standard
                    r.update(control=name,subset=list(subset),expected=expected,passed=passed)
                    if expected=='borromean' and passed:r['status']='borromean_by_standard_diagram'
                    assert passed,(name,subset,seed,r['simplified_crossings']);rows.append(r)
    curves=control_curves()['hopf'];tests=[]
    for name,cc in [('original',curves),('reflected',[c*np.array([-1,1,1]) for c in curves]),('one_reversed',[curves[0][::-1],curves[1]])]:
        r,_=classify(cc);value=polygon_gauss(clean(cc[0]),clean(cc[1]));tests.append(dict(name=name,diagram=r['oriented_linking_matrix'][0][1],solid_angle=float(value)))
    sign=tests[0]['diagram'];assert sign!=0 and tests[1]['diagram']==tests[2]['diagram']==-sign
    # With the projection convention above, the solid-angle integral has the opposite sign.
    convention=round(tests[0]['solid_angle'])/sign
    assert all(abs(t['solid_angle']-convention*t['diagram'])<1e-8 for t in tests)
    save(OUT/'controls.json',dict(rows=rows,orientation_controls=tests,gauss_to_diagram_sign=convention,passed=all(r['passed'] for r in rows)))
    np.savez_compressed(OUT/'standard-link-curves.npz',**{f'{k}_{i}':c for k,v in control_curves().items() for i,c in enumerate(v)})
    print('Controls passed',len(rows),'orientation',tests,flush=True)


def analyze(digits,subdivisions):
    ids=['F1a','F1b','F2a','F2b'];curves={id:np.load(OUT/f'{id}-d{digits}-n{subdivisions}.npz')['xyz'] for id in ids};records=[]
    dest=OUT/'diagrams';dest.mkdir(exist_ok=True)
    for size in [1,2,3,4]:
        for subset in itertools.combinations(ids,size):
            for seed in [0,2]:
                label='-'.join(subset)+f'-d{digits}-n{subdivisions}-s{seed}';path=dest/(label+'.json')
                if path.exists():records.append(read(path));continue
                started=time.monotonic();r,link=classify([curves[id] for id in subset],seed)
                r.update(ids=list(subset),digits=digits,subdivisions=subdivisions,source_sha256={id:digest(OUT/f'{id}-d{digits}-n{subdivisions}.npz') for id in subset},elapsed_seconds=time.monotonic()-started)
                save(path,r);records.append(r);print(label,r['crossings'],'->',r['simplified_crossings'],r['status'],flush=True)
    # Separate algorithm, directly in physical coordinates and independent of projection.
    gauss=[]
    for a,b in itertools.combinations(ids,2):
        v=float(polygon_gauss(clean(curves[a]),clean(curves[b])));gauss.append(dict(ids=[a,b],value=v,nearest_integer=round(v),integer_error=abs(v-round(v))))
    save(OUT/f'analysis-d{digits}-n{subdivisions}.json',dict(records=records,gauss_integrals=gauss))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['reconstruct','controls','analyze']);p.add_argument('--digits',type=int,default=60);p.add_argument('--subdivisions',type=int,default=96);args=p.parse_args();OUT.mkdir(exist_ok=True)
    if args.stage=='controls':controls()
    elif args.stage=='reconstruct':reconstruct(args.digits,args.subdivisions)
    else:analyze(args.digits,args.subdivisions)


if __name__=='__main__':main()
