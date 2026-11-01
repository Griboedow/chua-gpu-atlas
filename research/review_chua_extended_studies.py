"""Check stored inequalities and provenance without trusting passed flags."""
import copy
import hashlib
import importlib.metadata
import json
import subprocess
from flint import arb,ctx
from chua_research_followup import ROOT,OUT,dump
from study_chua_circuit_constraints import parameters
from certify_chua_coexistence_box import enclosure


def check_certificate(data):
    assert len(data['records'])==4
    for rec in data['records']:
        Q=[enclosure(v,True) for v in rec['shooting_box']]
        K=[enclosure(v) for v in rec['krawczyk_image']]
        assert len(Q)==len(K)==30
        assert all(q.contains_interior(k) for q,k in zip(Q,K))
        assert arb(rec['contraction']['ball'])<1
        assert arb(rec['determinant']['ball'])<1
        assert {r['multiplier'] for r in rec['multiplier_exclusions']}=={-1,1}
        assert all(arb(r['inverse_defect']['ball'])<1 for r in rec['multiplier_exclusions'])
        assert all(arb(v['ball'])>0 for v in rec['nominal_stability_certificate']['schur_margins'])
        assert len(rec['path_cells'])==10 and min(rec['path_cells'])>0
        assert all(not arb(v['ball']).contains(0) for v in rec['crossing_speeds'])
        assert all(arb(rec['section_box'][i]['ball']).contains(arb(rec['nominal_stability_certificate']['root'][i]['ball'])) for i in range(2))
    gaps=[]
    for i,a in enumerate(data['records']):
        for b in data['records'][i+1:]:
            gap=abs(arb(a['section_box'][0]['ball'])-arb(b['section_box'][0]['ball']))
            assert gap>0;gaps.append(True)
    assert len(gaps)==6


def main():
    ctx.dps=100
    cert=json.loads((OUT/'coexistence-box/certificate.json').read_text());check_certificate(cert)
    rejected=[]
    for kind in ['contraction','inclusion']:
        bad=copy.deepcopy(cert)
        if kind=='contraction':bad['records'][0]['contraction']['ball']='2'
        else:bad['records'][0]['krawczyk_image'][0].update(ball='100',lower='100',upper='100')
        try:check_certificate(bad)
        except AssertionError:rejected.append(kind)
    assert rejected==['contraction','inclusion']
    comp=json.loads((OUT/'coexistence-box/component-image.json').read_text())
    image,limit=parameters({k:enclosure(v) for k,v in comp['component_box'].items()})
    assert all(enclosure(cert['parameter_box'][k],True).contains(v) for k,v in image.items() if k!='gamma')
    assert comp['gamma_physical_lower']=='0' and comp['inductor_resistance_ohm']=='[0, 3e-7]'
    assert all(enclosure(comp['component_box'][k])>0 for k in ['C2','R','L'])
    assert enclosure(cert['parameter_box']['gamma'],True).contains(0)
    assert image['gamma'].upper()<arb(cert['parameter_box']['gamma']['upper']).lower()
    assert all(arb(r['x_abs_upper']['ball'])<limit for r in cert['records'])
    frozen=OUT/'extended-slopes/forecasts-frozen.json'
    published=subprocess.check_output(['git','show','0c39117:results/chua-followup/extended-slopes/forecasts-frozen.json'],cwd=ROOT)
    assert json.loads(published)==json.loads(frozen.read_text())
    forecast=json.loads(published)
    for name,digest in forecast['calibration_sha256'].items():
        assert hashlib.sha256((OUT/'extended-slopes'/name).read_bytes()).hexdigest()==digest
    score=json.loads((OUT/'extended-slopes/scores.json').read_text())
    assert score['total']==24 and score['corrected_hits']==24 and score['raw_hits']==16
    assert len(json.loads((OUT/'extended-slopes/independent-checks.json').read_text())['records'])==12
    derivatives=json.loads((OUT/'coexistence-box/derivative-checks.json').read_text())
    assert len(derivatives['records'])==4
    assert all(float(c['max_absolute_difference'])<1e-32 for r in derivatives['records'] for c in r['checks'])
    check=json.loads((OUT/'excitation-plane/independent-checks.json').read_text())
    assert len(check['records'])==14,len(check['records'])
    assert all(len(r['checks'])==2 and all(c['label']==r['expected'] for c in r['checks']) for r in check['records'])
    scripts=['certify_chua_coexistence_box','study_chua_excitation_plane','extend_chua_slope_predictions',
             'check_chua_extended_predictions','study_chua_coexistence_tolerances','plot_chua_extended_studies','review_chua_extended_studies','check_chua_coexistence_box','derive_chua_pulse_limit']
    result=dict(passed=True,cycle_certificates=4,pair_separations=6,negative_controls_rejected=rejected,
        independent_excitation_cases=len(check['records']),independent_prediction_cycles=12,
        frozen_forecast_commit='0c39117',versions={p:importlib.metadata.version(p) for p in ['python-flint','numpy','scipy','mpmath','numba']},
        source_sha256={f'research/{name}.py':hashlib.sha256((ROOT/'research'/(name+'.py')).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in scripts},
        scope='Stored inequality and provenance review; does not replace re-executing the interval program or the written mathematical argument.')
    dump(OUT/'extended-review.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
