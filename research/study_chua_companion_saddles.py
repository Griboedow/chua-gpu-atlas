"""Solve saddle companions at the attracting references' identical parameters."""
import argparse
import json
import time
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import expm
from chua_precision_flights import PrecisionFlights
from study_chua_return_geometry import ROOT, OUT, IDS, read, digest


def fixed_beta(f, beta, guess):
    q = mp.matrix(guess)
    history = []
    for iteration in range(24):
        ev = f.evaluate(beta, list(q), False)
        err = max(abs(v) for v in ev['residual'])
        history.append(mp.nstr(err, 12))
        if err < mp.mpf(10)**(-(f.dps-20)):
            break
        change = mp.lu_solve(ev['jacobian'], -ev['residual'])
        for scale in [1, mp.mpf('.5'), mp.mpf('.25'), mp.mpf('.125')]:
            trial = q+scale*change
            if min(trial[2:]) <= 0:
                continue
            if max(abs(v) for v in f.evaluate(beta, list(trial), False)['residual']) < err:
                q = trial
                break
        else:
            raise RuntimeError('Fixed-parameter Newton stalled: '+str(history))
    else:
        raise RuntimeError('Fixed-parameter Newton exhausted: '+str(history))
    tr, det = ev['trace'], ev['determinant']
    large = (tr+mp.sqrt(tr*tr-4*det))/2
    return dict(beta=mp.nstr(beta, f.dps-5), q=[mp.nstr(v, f.dps-5) for v in q],
                period=mp.nstr(2*sum(q[2:]), f.dps-5),
                half_multiplier=mp.nstr(large, 24), full_multiplier=mp.nstr(large*large, 24),
                residual=mp.nstr(err, 12), newton_residuals=history), ev


def path_checks(f, beta, ev):
    violation = 0.
    speed = 1e100
    for arc in ev['arcs']:
        reg = arc['region']
        slope = f.m0 if reg == 0 else f.m1
        A = np.array([[-float(f.alpha*(1+slope)), float(f.alpha), 0], [1,-1,1], [0,-float(beta),0]])
        c = np.array([reg*float((f.m1-f.m0)/(1+f.m1)), 0, -reg*float((f.m1-f.m0)/(1+f.m1))])
        u = np.array(list(arc['start']), float)
        for t in np.linspace(0, float(arc['duration']), 129):
            state = c+expm(A*t)@(u-c)
            violation = max(violation, abs(state[0])-1 if reg == 0 else 1-reg*state[0])
        speed = min(speed, abs((A@(np.array(list(arc['end']), float)-c))[0]))
    return dict(sampled_region_violation=violation, minimum_end_crossing_speed=speed)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ids', nargs='*', default=IDS)
    p.add_argument('--digits', type=int, default=70)
    p.add_argument('--root-files',nargs='*',type=Path)
    args=p.parse_args()
    OUT.mkdir(exist_ok=True)
    sources=args.root_files or [ROOT/f'results/chua-knot-atlas/{id}-d70-n192.json' for id in args.ids]
    for source in sources:
        source=source.resolve()
        row=read(source);id=row['id']
        target=OUT/f'{id}-companion-d{args.digits}.json'
        if target.exists():
            print('Retained', target.name, flush=True)
            continue
        started=time.monotonic()
        f=PrecisionFlights(row['regions'], row['targets'], True, args.digits)
        beta=mp.mpf(row['beta'])
        stable=mp.matrix([mp.mpf(x) for x in row['q']])
        fold=f.boundary(row['beta'], row['q'], 1)
        fold_q=mp.matrix([mp.mpf(x) for x in fold['q']])
        saddle,ev=fixed_beta(f, beta, 2*fold_q-stable)
        saddle_q=mp.matrix([mp.mpf(x) for x in saddle['q']])
        delta=saddle_q-stable
        distance=mp.sqrt(delta[0]**2+delta[1]**2)
        path=path_checks(f, beta, ev)
        result=dict(id=id, digits=args.digits, source=source.relative_to(ROOT).as_posix(),
                    source_sha256=digest(source), fold=fold, companion=saddle,
                    section_displacement=[mp.nstr(delta[i],24) for i in range(2)],
                    section_distance=mp.nstr(distance,24),
                    fold_to_reference=mp.nstr(mp.mpf(fold['beta'])-beta,24),
                    path_checks=path, elapsed_seconds=time.monotonic()-started)
        result['passed']=bool(float(saddle['residual'])<1e-30 and distance>mp.mpf('1e-12') and mp.mpf(saddle['half_multiplier'])>1 and path['sampled_region_violation']<1e-9 and path['minimum_end_crossing_speed']>0.01)
        target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(id,'companion',saddle['half_multiplier'],'distance',result['section_distance'],
              'seconds',round(result['elapsed_seconds'],2),'passed',result['passed'],flush=True)
        if not result['passed']:
            raise SystemExit('Companion failed its numerical checks')


if __name__=='__main__':main()
