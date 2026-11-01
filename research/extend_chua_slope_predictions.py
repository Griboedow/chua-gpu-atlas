"""Out-of-sample replication of the fixed slope-window prediction rule.

Calibration uses only 35-, 45- and 55-peak cycles at the new nodes. The previous
central bias correction is reused without refitting. Freeze and publish before
running search. All outputs are separate from the original experiment.
"""
import argparse
import concurrent.futures
import json
import subprocess
import sys
from pathlib import Path
import predict_chua_slope_sheets as study
from chua_research_followup import ROOT, OUT, BASE, dump

OLD = study.DIR
DIR = OUT / 'extended-slopes'


def nodes():
    for free in ['beta', 'alpha']:
        for i, j in [(-1,-1), (-1,1), (1,-1), (1,1)]:
            p = BASE.copy()
            p['m0'] += i*.0008
            p['m1'] += j*.0008
            yield f'{free}_{i+1}{j+1}', free, i, j, p


def initialize():
    if (DIR/'protocol.json').exists():
        return
    original = json.loads((OLD/'protocol.json').read_text())
    original.update(slope_offsets=[-.0008,.0008],
        nodes=[name for name,*_ in nodes()],
        purpose='Test the unchanged model outside the previous slope square; twice its half-width.',
        central_bias_correction='Reuse the original central residuals; no fitting to new target endpoints.',
        target_count=24,
        publication_rule='Commit and push forecasts-frozen.json before any new target solve.',
        failure_policy='Keep all misses and unresolved targets; do not refit this experiment.')
    dump(DIR/'protocol.json', original)
    central=json.loads((OLD/'central-residual-calibration.json').read_text())
    for free in ['alpha','beta']:
        for row in json.loads((OLD/(free+'_11-calibration.json')).read_text())['training']:
            central[free][str(row['peaks'])]=dict(fold=row['fold'],
                role='Previously published central training fold; shooting initialization only.')
    dump(DIR/'central-residual-calibration.json',central)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('phase', choices=['init','calibrate','freeze','search'])
    ap.add_argument('--only')
    ap.add_argument('--workers',type=int,default=2)
    args=ap.parse_args()
    initialize()
    study.DIR=DIR
    study.nodes=nodes
    study.initialize=initialize
    # Follow each training fold from its known central value. A direct coarse
    # shooting guess can converge to another root with the same formal word.
    original_seed=study.boundary_seed
    templates=study.classical_templates()
    def training_seed(seed,p,free,guess):
        peaks=next(k for k,v in templates.items() if v['id']==seed['id'])
        if peaks in [35,45,55]:
            return study.continued_target_seed(peaks,p,free)
        return original_seed(seed,p,free,guess)
    study.boundary_seed=training_seed
    if args.phase=='init': return
    if args.phase=='freeze': study.freeze(); return
    if args.only:
        getattr(study,args.phase)(args.only)
        return
    def run(name):
        log=DIR/(name+'-'+args.phase+'.log')
        with log.open('w',encoding='utf-8') as stream:
            result=subprocess.run([sys.executable,'-X','utf8',__file__,args.phase,'--only',name],
                                  stdout=stream,stderr=subprocess.STDOUT,cwd=ROOT)
        print(name,args.phase,'exit',result.returncode,flush=True)
        return result.returncode
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        codes=list(pool.map(run,[name for name,*_ in nodes()]))
    if any(codes): raise SystemExit('Some nodes failed; see retained logs.')


if __name__=='__main__': main()
