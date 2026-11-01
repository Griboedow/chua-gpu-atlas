"""Enclose a saddle cycle coexisting with the certified hidden periodic sink."""
import json
from flint import arb,ctx
from certify_chua_coexisting_cycles import ReturnMap,certify
from chua_interval_map import col,eye,midpoint,bounds
from check_chua_slope_refinement import independent_segmented
from study_chua_hidden_boundary import OUT,dump


def saddle(record,digits,shifted=False):
    ctx.dps=digits
    parameters=(arb(42)/5,arb(12),-arb(19)/20,-arb(1)/20)
    f=ReturnMap(record,digits,parameters)
    centre=col([arb(str(v)) for v in record['cycle']['q'][:2]])
    for _ in range(10):
        value,J,_=f.evaluate(centre.entries())
        centre=midpoint(centre-(J-eye(2)).inv()*(value-centre))
    if shifted:centre[0,0]+=arb('1e-12')
    box=col([arb(v,arb('1e-35')) for v in centre.entries()])
    value,J0,_=f.evaluate(centre.entries())
    _,J,arcs=f.evaluate(box.entries(),True)
    C=midpoint((J0-eye(2)).inv());R=eye(2)-C*(J-eye(2))
    K=centre-C*(value-centre)+R*(box-centre)
    contraction=max(sum(abs(R[i,j]) for j in range(2)).abs_upper() for i in range(2))
    if not all(box[i,0].contains_interior(K[i,0]) for i in range(2)) or not contraction<1:
        raise ValueError('Krawczyk inclusion rejected')
    trace,det=J.trace(),J.det();disc=(trace*trace-4*det).sqrt()
    large=(trace+disc)/2;small=det/large
    assert large>1 and small>0 and small<1
    assert 1-trace+det<0  # The sink criterion must reject this saddle.
    return dict(passed=True,id=record['id'],root=[bounds(v) for v in box.entries()],
        krawczyk_image=[bounds(v) for v in K.entries()],contraction=bounds(contraction),
        period=bounds(sum(a['duration'] for a in arcs)),multipliers=[bounds(large),bounds(small)],
        path_cells=[a['path']['accepted_cells'] for a in arcs],first_upward_return=True,
        sink_classification_rejected=True)


def main():
    records=json.loads((OUT/'cycle-seeds.json').read_text())['full_records']
    for digits in [80,100]:
        s=saddle(records[0],digits)
        a=certify(records[1],digits,parameters=(arb(42)/5,arb(12),-arb(19)/20,-arb(1)/20))
        rejected=False
        try:saddle(records[0],digits,True)
        except ValueError as exc:rejected='Krawczyk' in str(exc)
        assert rejected
        result=dict(passed=True,working_digits=digits,
            parameters=dict(alpha='42/5',beta='12',m0='-19/20',m1='-1/20',gamma='0'),
            saddle=s,attractor=a,shifted_root_rejected=rejected,
            scope='Two distinct first-return periodic orbits, one saddle and one sink. '
                  'The stable origin and hiddenness of the sink were certified in the accessibility study. '
                  'The global stable manifold and its basin-separation role are not interval-certified here.')
        dump(OUT/f'saddle-certificate-{digits}.json',result)
        print(digits,'digits saddle and sink certified',flush=True)
    controls=[]
    for record in records:
        for method in ['DOP853','Radau']:
            check=independent_segmented(record['cycle'],method)
            assert check['success'] and check['closure']<2e-8 and check['multiplier_error']<2e-6
            controls.append(dict(id=record['id'],**check))
    dump(OUT/'cycle-independent-controls.json',dict(passed=True,controls=controls))


if __name__=='__main__':main()
