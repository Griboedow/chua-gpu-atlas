"""Freeze a theory-guided forecast using only previously published layers.

This script does not solve the target long cycles. Its only new orbit solve is
the short organizing saddle, used to calibrate the expansion factor and period.
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import root_scalar

from trace_chua_window import Chua
from chua_fast_flights import FastFlights

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/chua-layer-prediction'
INPUTS = ['results/chua-layer-census/census-analysis.json',
          'results/chua-layer-atlas/extended-adding.json',
          'results/chua-layers/affine-central.json',
          'results/chua-layer-atlas/family-registry.json']


def old_layers():
    old = json.loads((ROOT / INPUTS[0]).read_text())['adding_sequence']
    extra = json.loads((ROOT / INPUTS[1]).read_text())
    return [dict(n=(r['maxima'] - 5)//10, peaks=r['maxima'],
                 fold=r['fold'], symmetry=r['symmetry'],
                 width=r['symmetric_core_width']) for r in old] + [
        dict(n=(r['peak_count']-5)//10, peaks=r['peak_count'],
             fold=r['fold'], symmetry=r['symmetry'],
             width=r['symmetric_core_width']) for r in extra
        if r['status'] == 'verified_attracting_cycle']


def parent_at(beta):
    central = json.loads((ROOT / INPUTS[2]).read_text())
    f = FastFlights(Chua(), central['regions'], central['targets'], True)
    q = f.solve(beta, np.array(central['boundaries']['-1']['q']))
    r = f.record(beta, q, True)
    assert r['region_violation'] < 1e-8
    assert r['full_multipliers'][0][0] > 1
    return r


def corrected_fit(rows):
    """Four known folds fix limit,A,B,C; the parent fixes r self-consistently."""
    ns = np.array([r['n'] for r in rows], dtype=float)
    bs = np.array([r['fold']['beta'] for r in rows])
    base = bs[-1]

    def fit(limit):
        parent = parent_at(limit)
        ratio = 1 / parent['full_multipliers'][0][0]
        # Scale by the last training n to avoid unnecessarily tiny columns.
        t = ratio**(ns-ns[-1])
        design = np.column_stack([np.ones(len(t)), t, t*t, (ns-ns[-1])*t*t])
        coef = np.linalg.solve(design, bs-base)
        return base+coef[0], ratio, coef, parent

    sol = root_scalar(lambda b: fit(b)[0]-b,
                      x0=14.30636, x1=14.30637, xtol=2e-14)
    limit, ratio, coef, parent = fit(sol.root)
    return dict(limit=limit, ratio=ratio, coefficients=coef.tolist(),
                beta_origin=base, n_origin=int(ns[-1]), parent=parent,
                training_n=ns.astype(int).tolist(), root_converged=sol.converged)


def predict_beta(fit, n):
    k = n-fit['n_origin']; r = fit['ratio']; t = r**k
    c = fit['coefficients']
    return fit['beta_origin']+c[0]+c[1]*t+c[2]*t*t+c[3]*k*t*t


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    path = args.output/'predictions-frozen.json'
    if path.exists():
        raise SystemExit('Frozen predictions already exist; refusing to overwrite.')
    rows = old_layers(); fit = corrected_fit(rows[-4:])
    previous_fit = corrected_fit(rows[-5:-1])
    r = fit['ratio']; parent = fit['parent']; last=rows[-1]; prev=rows[-2]
    targets=[]
    for n in [7, 8]:
        k=n-last['n']
        corrected=predict_beta(fit,n)
        simple=last['fold']['beta']+(last['fold']['beta']-prev['fold']['beta'])*sum(r**j for j in range(1,k+1))
        width=last['width']*r**(2*k)
        period=last['fold']['period']+2*k*parent['period']
        targets.append(dict(n=n,peaks=10*n+5,section_returns=2*n+1,
            symmetric=True,half_outer_word='RRLL'*n+'RR',
            fold_beta_corrected=corrected,fold_beta_leading=simple,
            symmetric_core_width=width,fold_period=period,
            acceptance=dict(corrected_location_absolute_error=2e-8,
                width_relative_error=.20,period_absolute_error=.002,
                rationale='Pre-search engineering gates, not statistical confidence intervals; location is also scored in measured core widths.')))
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),
        stage='Frozen before any new 75/85-peak cycle solve in this experiment',
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in INPUTS},
        known_layers=[dict(n=x['n'],peaks=x['peaks'],fold_beta=x['fold']['beta'],
            width=x['width'],fold_period=x['fold']['period']) for x in rows],
        theory=dict(position_ratio='1/Lambda',width_ratio='1/Lambda^2',
            period_increment='2*Tp',position_correction='limit + A*r^n + (B+C*n)*r^(2*n)',
            mechanism_status='Quadratic global-return/tangency hypothesis; not a proved global connection'),
        calibration=fit,
        retrospective_control=dict(note='P65 was already known; this is NOT a blind test.',
            fit=previous_fit,predicted_fold=predict_beta(previous_fit,6),
            known_fold=last['fold']['beta'],absolute_error=abs(predict_beta(previous_fit,6)-last['fold']['beta'])),
        predictions=targets)
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(limit=fit['limit'],Lambda=1/r,r=r,width_ratio=r*r,
        parent_period=parent['period'],retrospective_error=report['retrospective_control']['absolute_error'],
        predictions=targets),indent=2))


if __name__=='__main__':
    main()
