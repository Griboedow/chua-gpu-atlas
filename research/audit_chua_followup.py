"""Check finite evidence, prediction provenance and theorem-domain containment."""
import json
import subprocess
import hashlib
from decimal import Decimal
from flint import arb,ctx
from chua_research_followup import ROOT,OUT,dump
from publication_schema import current_names


def read(p):return json.loads((OUT/p).read_text())


def main():
    ctx.dps=100;checks={}
    frozen=OUT/'sheets/forecasts-frozen.json';original=subprocess.check_output(['git','show','112b01f:results/chua-followup/sheets/forecasts-frozen.json'],cwd=ROOT)
    checks['published_forecasts_unchanged']=current_names(original.decode('utf-8')).replace('\r\n','\n')==frozen.read_text(encoding='utf-8')
    # Existing search records identify the original published file. New runs
    # identify its renamed version; exact contents are compared above.
    digests={hashlib.sha256(b.replace(b'\r\n',b'\n')).hexdigest() for b in [original,frozen.read_bytes()]}
    searches=[json.loads(p.read_text()) for p in (OUT/'sheets').glob('*-search.json')]
    checks['all_searches_use_frozen_forecasts']=len(searches)==18 and all(d['forecast_sha256'] in digests for d in searches)
    checks['all_54_targets_solved']=all(len(d['records'])==3 and all(r['status']=='numerical_root' for r in d['records']) for d in searches)
    a=read('analysis.json');checks['48_held_out_corrected_hits']=a['sheets']['held_out']==48 and a['sheets']['corrected_midpoints_inside']==48
    checks['16_uncorrected_misses_retained']=a['sheets']['raw_midpoints_inside']==32
    for part in ['sheets','trigger']:
        c=read(part+'/independent-checks.json');checks[part+'_12_independent_checks']=c['passed'] and len(c['records'])==12
    for digits in [80,100]:
        local=read(f'trigger/local-certificate-{digits}.json');checks[f'pulse_local_{digits}']=local['passed'] and len(local['records'])==7 and all(all(r['gates'].values()) for r in local['records'])
        e=read(f'hidden-search/uniform-equilibrium-exclusion-{digits}.json');checks[f'hidden_neighborhoods_{digits}']=e['passed'] and e['complete']
    below=read('trigger/below-certificate-80.json');rows=below['records']
    checks['whole_lower_amplitude_coverage']=below['passed'] and not below['partial'] and below['central']['passed'] and all(Decimal(x['upper'])==Decimal(y['lower']) for x,y in zip(rows[:-1],rows[1:]))
    center=read('hidden-search/connection-center-14-2.json');b=arb(center['beta']);J=arb('29.516618475').union(arb('29.516618493'))
    checks['J_inside_uniform_exclusion']=arb(b,arb('1e-8')).contains(J)
    checks['tangency_box_inside_interior_J']=J.contains_interior(arb(b,arb('5e-12')))
    for digits in [140,160]:
        c=read(f'hidden-search/hidden-tangency-{digits}.json')
        checks[f'tangency_{digits}']=c['status']=='passed' and all(c['gates'].values()) and all(c['local_gates'].values()) and all(c['graph_gates'].values())
        checks[f'complete_path_{digits}']=len(c['itinerary_cells'])==36 and all(len(row)==2 and all(n>0 for n in row) for row in c['itinerary_cells'])
    for sign in ['minus','plus']:
        c=read(f'hidden-search/shifted-tangency-{sign}.json');checks['shifted_'+sign+'_rejected']=c['status']=='failed' and not all(c['gates'].values())
    checks['independent_hidden_controls']=read('hidden-search/independent-tangency-checks.json')['passed']
    checks['classical_regression']=read('hidden-search/classical-regression-80.json')['status']=='passed'
    dump(OUT/'audit.json',dict(passed=all(checks.values()),checks=checks,
        scope='Audit of finite records and parameter containment. The analytic localization, graph-transform implication and imported Newhouse–Robinson theorem are stated in the proof report; this audit is not their formal verification.'))
    assert all(checks.values()),checks
    print('Passed',len(checks),'evidence and provenance checks')


if __name__=='__main__':main()
