"""Long GPU basin probe at eta=1.082, where a stable orbit also exists."""
import json
from pathlib import Path
import numpy as np
from continue_synchronematic_orbit import Orbit
from synchronematic import SynchronematicGPU,COEFFICIENTS,gradient_coefficients
from search_synchronematic import save


def main():
    root=Path('results/synchronematic');branch=json.loads((root/'n4-periodic-branch/branch.json').read_text())
    row=min(branch['records'],key=lambda r:abs(r['eta']-1.082));model=Orbit(branch['positions'])
    base=model.E@np.array(row['state']);rng=np.random.default_rng(2026102948)
    starts=[];labels=[]
    for scale in [.001,.01,.1]:
        for _ in range(8):
            starts.append(base+rng.normal(0,scale,8));labels.append(f'near_periodic_sigma_{scale}')
    for _ in range(32):starts.append(rng.uniform(-np.pi,np.pi,8));labels.append('global_uniform_angles')
    old=np.load(root/'n4-transition/scan.npz');ids=np.flatnonzero(abs(old['experiment_labels'][:,0]-1.082)<1e-8)
    for i in ids:starts.append(old['starts'][i]);labels.append(f'previous_start_{int(i)}')
    starts=np.array(starts);positions=np.broadcast_to(branch['positions'],(len(starts),4,2)).copy()
    g=gradient_coefficients();c=g+1.082*(COEFFICIENTS-g)
    result=SynchronematicGPU(4).scan(positions,starts,dt=.0025,transient=4000,duration=8000,
        coefficients=c,progress=lambda f,t:print(f'Basin probe {f:.1%}, {t:.1f}s',flush=True))
    result['experiment_labels']=np.array(labels)
    result['experiment']=dict(eta=1.082,seed=2026102948,changes_equations=True,
        description='Local neighborhoods of the attracting periodic orbit, random global angles, and eight previously sampled starts; fractions are not basin-volume estimates.')
    out=root/'n4-coexistence-gpu';save(out,result)
    groups=[]
    for label in sorted(set(labels)):
        mask=np.array(labels)==label;valid=result['status'][mask]==0
        positive=valid&(result['lle'][mask]>.002)&(result['lle_windows'][mask].min(axis=1)>.001)
        groups.append(dict(label=label,count=int(mask.sum()),completed=int(valid.sum()),retained=int(positive.sum()),
                           lle=result['lle'][mask].tolist(),windows=result['lle_windows'][mask].tolist()))
    (out/'groups.json').write_text(json.dumps(groups,indent=2));print('Saved coexistence probe',flush=True)


if __name__=='__main__':main()
