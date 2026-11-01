"""Local stability predictions, a global sector certificate and a hidden cycle.

The branch seeds are exploratory continuation data. The basin experiment uses
fixed new directions and radii and is reported separately from the interval
statements. No finite-time classification is promoted to a proof.
"""
import argparse
import json
import numpy as np
from flint import arb, arb_mat, ctx
from certify_chua_coexisting_cycles import ROOT, OUT, ReturnMap, certify
from certify_chua_hidden_cycle import lyapunov, positive_definite
from chua_interval_map import bounds
from probe_chua_accessibility import simulate


def matrix(m):
    return arb_mat([[-arb(42)/5*(1+m),arb(42)/5,0],[1,-1,1],[0,-12,0]])


def run(digits=80):
    ctx.dps=digits
    # A common rational quadratic Lyapunov function for an entire slope sector.
    P, pc=lyapunov(matrix(-arb(1)/20))
    endpoints=[-arb(11)/200,-arb(9)/200]
    sector=[]
    for m in endpoints:
        Q=-(matrix(m).transpose()*P+P*matrix(m))
        sector.append(dict(slope=str(m),Q_minors=positive_definite(Q)))
    rejected=False
    wrong=matrix(-arb(1)/2)
    try:
        positive_definite(-(wrong.transpose()*P+P*wrong))
    except AssertionError:
        rejected=True
    assert rejected
    # A second benchmark has a UNIQUE stable equilibrium, but also a large
    # periodic sink. Its origin ellipsoid suffices to establish hiddenness.
    rows=json.loads((OUT/'hidden-inner-continuation.json').read_text())
    row=min(rows,key=lambda r:abs(r['m0']+.95))
    seed=json.loads((OUT/'hidden-cycle-seed.json').read_text())
    seed=dict(seed,id='H2',cycle=row)
    pars=(arb(42)/5,arb(12),-arb(19)/20,-arb(1)/20)
    orbit=certify(seed,digits,parameters=pars)
    P0, local=lyapunov(matrix(pars[2]))
    level=1/(2*P0.inv()[0,0])
    pnorm=max(sum(abs(P0[i,j]) for j in range(3)).abs_upper() for i in range(3))
    radius=(level/(2*pnorm)).sqrt().lower()
    # Routh--Hurwitz equalities are computed in Arb; admissibility of the
    # coefficients must also be checked before interpreting the two roots.
    alpha=arb(42)/5
    b=alpha+alpha*alpha
    disc=b*b-4*alpha*alpha*12
    thresholds=[(-b-disc.sqrt())/(2*alpha*alpha),(-b+disc.sqrt())/(2*alpha*alpha)]
    predictions=[]
    for value in ['-0.97','-0.95','-0.93','-0.90']:
        m=arb(value); a1=1+alpha*(1+m); a2=12+alpha*m; a3=12*alpha*(1+m)
        margin=a1*a2-a3
        stable=all(x>0 for x in [a1,a2,a3,margin])
        predictions.append(dict(m0=value,origin_stable=stable,routh_margin=bounds(margin),
            basin_prediction='small starts approach origin' if stable else 'small starts approach the continued cycle (hypothesis)'))
    result=dict(passed=True,working_digits=digits,alpha='42/5',beta='12',gamma='0',
        global_sector=dict(m0_and_m1_interval=['-11/200','-9/200'],lyapunov=pc,endpoint_checks=sector,unstable_slope_control_rejected=rejected,
            conclusion='Global exponential convergence to the unique origin for every pair of slopes in this closed sector.'),
        unique_equilibrium_hidden_cycle=dict(parameters=dict(alpha='42/5',beta='12',m0='-19/20',m1='-1/20',gamma='0'),
            periodic_orbit=orbit,origin_lyapunov=local,ellipsoid_level=bounds(level),excluded_origin_radius=bounds(radius),
            conclusion='A stable origin and a hidden attracting periodic orbit coexist. Global convergence to the equilibrium is false.'),
        routh_thresholds=[bounds(v) for v in thresholds],basin_predictions=predictions,
        scope='The global sector and coexistence point are separated in parameter space. They do not bracket a unique or first global-stability boundary.')
    return result


def basin_experiment():
    predictions=json.loads((OUT/'stability-certificate-80.json').read_text())['basin_predictions']
    directions=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,-1,1],[-1,1,1],[1,1,-1],[-1,-1,-1]],float)
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    radii=[1e-3,1e-6,1e-9]
    initial=np.array([r*d for r in radii for d in directions])
    rows=json.loads((OUT/'hidden-inner-continuation.json').read_text())
    records=[]
    for pred in predictions:
        m0=float(pred['m0']); row=min(rows,key=lambda r:abs(r['m0']-m0))
        roots=np.array([row['q'][:2]])
        controls=[]
        for step in [.01,.005]:
            labels,errors,sections,final=simulate(initial,dict(alpha=8.4,beta=12,m0=m0,m1=-.05),step,1000,roots)
            eq=np.linalg.norm(final,axis=1)<1e-8
            labels[eq]=1
            # 0 cycle, 1 equilibrium, -1 unresolved.
            controls.append(dict(step=step,labels=labels.tolist(),counts=[int(np.sum(labels==k)) for k in [-1,0,1]],
                                 final_norms=np.linalg.norm(final,axis=1).tolist()))
        disagreement=int(np.sum(np.array(controls[0]['labels'])!=np.array(controls[1]['labels'])))
        records.append(dict(prediction=pred,controls=controls,refinement_disagreements=disagreement))
        print(m0,'labels',controls[-1]['counts'],'disagreements',disagreement,flush=True)
    result=dict(duration=1000,radii=radii,directions=directions.tolist(),records=records,section_interpolation='cubic Hermite',
                label_order=['unresolved','cycle','origin'],
                scope='Finite-time validation of a local-stability-based prediction. Self-excitation and the exact hiddenness boundary are not certified by these samples.')
    (OUT/'stability-basin-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--digits',type=int,default=80)
    p.add_argument('--basins',action='store_true')
    a=p.parse_args()
    if a.basins:
        basin_experiment()
    else:
        report=run(a.digits)
        (OUT/f'stability-certificate-{a.digits}.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Certified global sector and hidden cycle with a unique equilibrium.',flush=True)
        print('Origin stability thresholds:',[v['ball'] for v in report['routh_thresholds']],flush=True)


if __name__=='__main__':
    main()
