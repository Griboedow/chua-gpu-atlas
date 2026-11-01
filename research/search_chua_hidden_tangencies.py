"""Candidate searches for hidden coexistence, retaining negative outcomes.

Finite-time Lyapunov classifications are screens, never homoclinic certificates.
The direct continuation follows only the already known organizer branch.
"""
import argparse
import json
import os
import sys
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from chua_research_followup import ROOT,OUT,BASE,dump,equilibrium_test,continue_tangency,approximate_tangency

DIR=OUT/'hidden-search'


def continuation():
    source=approximate_tangency(BASE);rows=[source];failures=[]
    # A fixed-alpha, fixed-outer-slope path tests whether the existing organizer
    # reaches the unique-equilibrium region. Failure is not a nonexistence proof.
    for m0 in np.linspace(BASE['m0'], -.94, 82)[1:]:
        p=source['parameters'].copy();p['m0']=float(m0)
        try:
            source=continue_tangency(p,source=source);rows.append(source)
            print('organizer m0',m0,'beta',source['coordinate'],'stable',source['equilibrium']['stable'],flush=True)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            failures.append(dict(target=p,error=str(exc)));break
        dump(DIR/'classical-continuation.json',dict(records=rows,failures=failures))
    dump(DIR/'classical-continuation.json',dict(records=rows,failures=failures,
        scope='Numerical continuation of one fixed itinerary organizer. A stopped branch does not exclude other connections.'))


def screen(kind='stable'):
    os.environ['PYOPENCL_CACHE_DIR']=str(ROOT/'.publication/opencl-cache')
    os.environ['TEMP']=str(ROOT.parent/'.agent-tmp');os.environ['TMP']=os.environ['TEMP']
    sys.path.insert(0,str(ROOT));from chua_gpu import GPU
    exponent=11 if kind=='stable' else 14 if kind=='broad' else 12
    u=qmc.Sobol(5,scramble=True,seed=20261101).random_base2(exponent)
    params=[];regions=[]
    for a,b,c,d,e in u:
        alpha=(3+57*a) if kind=='broad' else 6+24*a
        beta=1.05*alpha+b*(.245*(alpha+1)**2-1.05*alpha)
        disc=np.sqrt((alpha+1)**2-4*beta)
        low=(-(alpha+1)-disc)/(2*alpha);high=(-(alpha+1)+disc)/(2*alpha)
        inner=-1+(.05+.9*c)*(low+1)
        outer=high+.005+d*(.4-high-.005)
        if e<.5:m0,m1=inner,outer;region='low_inner'
        else:m0,m1=outer,inner;region='high_inner'
        if kind=='broad':m1=-.999+1.499*d
        if kind=='outer':m0=-1.005-1.495*c;m1=outer;region='stable_outer'
        p=dict(alpha=alpha,beta=beta,m0=m0,m1=m1)
        if kind!='outer':assert equilibrium_test(p)['unique'] and equilibrium_test(p)['stable']
        params.append([alpha,beta,m0,m1]);regions.append(region)
    params=np.array(params)
    # Multiple finite-amplitude starts are needed because small starts decay.
    seeds=np.array([[0,1,0],[0,3,20],[5,0,-5],[0,8,40]],float)
    protocol=dict(parameters=len(params),seeds=seeds.tolist(),sampling='Scrambled Sobol, seed 20261101, '+kind+' sector screen',
        scope='Stable unique equilibrium is checked at every sample. Four-window positive exponents select candidates only.')
    if kind=='outer':protocol['scope']='Two stable outer equilibria and an unstable origin. Hiddenness still requires basin exclusion near the origin.'
    dump(DIR/(kind+'-screen-protocol.json'),protocol)
    gpu=GPU();last=[-1]
    def progress(frac,seconds):
        pc=int(frac*10)
        if pc>last[0]:print('hidden screen',pc*10,'percent',round(seconds,1),flush=True);last[0]=pc
    result=gpu.scan(params,seeds=seeds,dt=.005,transient=400,duration=1200,progress=progress)
    stem='stable-origin' if kind=='stable' else kind
    np.savez_compressed(DIR/(stem+'-screen.npz'),parameters=params,regions=np.array(regions),
        **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    ids=np.argwhere(result['classes']==4)
    candidates=[dict(parameter_index=int(i),seed_index=int(j),parameters=dict(zip(['alpha','beta','m0','m1'],params[i].tolist())),
        lle=float(result['lle'][i,j]),initial=seeds[j].tolist(),final=result['final_states'][i,j].tolist()) for i,j in ids]
    dump(DIR/(stem+'-screen.json'),dict(protocol=protocol,settings=result['settings'],device=gpu.info,timing=result['timing'],
        counts={str(k):int((result['classes']==k).sum()) for k in range(6)},candidates=candidates,
        scope='No finite-time classification proves a chaotic attractor, a tangency, infinite coexistence, or absence of hidden regimes.'))
    print('candidate trajectories',len(candidates),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['continuation','screen','broad','outer']);a=ap.parse_args()
    DIR.mkdir(parents=True,exist_ok=True)
    if a.phase=='continuation':continuation()
    else:screen('stable' if a.phase=='screen' else a.phase)
