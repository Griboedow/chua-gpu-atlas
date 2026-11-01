"""Score the frozen predictions and collect explicitly scoped evidence."""
import json
from decimal import Decimal
from chua_research_followup import OUT,dump


def read(p):return json.loads((OUT/p).read_text())


def main():
    rows=[]
    for source in sorted((OUT/'sheets').glob('*-search.json')):
        d=json.loads(source.read_text());cal=read('sheets/'+d['node']+'-calibration.json')
        for r in d['records']:
            assert r['status']=='numerical_root'
            row={k:r[k] for k in ['peaks','midpoint_inside','transport_midpoint_inside','signed_fold_error_in_widths','transport_fold_error_in_widths','width_relative_error']}
            low=min(Decimal(r['fold']['coordinate']),Decimal(r['symmetry']['coordinate']));width=abs(Decimal(r['signed_width']))
            row.update(node=d['node'],section=d['free_parameter'],held_out=not d['node'].endswith('_11'),
                raw_midpoint_in_core_units=float((Decimal(r['forecast']['midpoint'])-low)/width),
                transport_midpoint_in_core_units=float((Decimal(r['forecast']['transport_midpoint'])-low)/width),
                measured_width=float(width),organizer_control_difference=cal['organizer_difference'])
            rows.append(row)
    held=[r for r in rows if r['held_out']];assert len(held)==48
    sheet=dict(held_out=48,calibration_targets=6,raw_midpoints_inside=sum(r['midpoint_inside'] for r in held),
        transported_midpoints_inside=sum(r['transport_midpoint_inside'] for r in held),
        max_transport_fold_error_in_widths=max(abs(r['transport_fold_error_in_widths']) for r in held),
        max_absolute_width_relative_error=max(abs(r['width_relative_error']) for r in held),
        independent_checks=read('sheets/independent-checks.json'),records=rows)
    hidden=read('hidden-search/hidden-tangency-160.json');uniform=read('hidden-search/uniform-equilibrium-exclusion-100.json')
    assert hidden['status']=='passed' and uniform['passed']
    assert read('hidden-search/independent-tangency-checks.json')['passed']
    report=dict(sheets=sheet,trigger_local=read('trigger/local-certificate-100.json')['passed'],
        trigger_lower_amplitude=dict(passed=read('trigger/below-certificate-80.json')['passed'],
            interval_count=len(read('trigger/below-certificate-80.json')['records'])),
        trigger_independent=read('trigger/independent-checks.json')['passed'],
        hidden=dict(tangency_passed=True,uniform_hiddenness_passed=True,
            flight_count=2*hidden['return_iterations'],time_cells=sum(sum(r) for r in hidden['itinerary_cells']),
            beta_interval=hidden['parameter_interval'],squared_multipliers=hidden['full_return_multipliers'],
            squared_determinant=hidden['full_return_determinant'],
            scope='Computer-assisted hypotheses plus a localized application of the parameter Newhouse–Robinson theorem. No decimal witness of infinite coexistence, chaotic-attractor certificate, complete Lean proof, or priority claim.'))
    dump(OUT/'analysis.json',report)
    print({k:v for k,v in sheet.items() if k not in ['records','independent_checks']})
    print('hidden cells',report['hidden']['time_cells'])


if __name__=='__main__':main()
