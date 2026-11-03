"""Compare frozen forecasts, short saddles and transported perturbations."""
import mpmath as mp
import numpy as np
from test_chua_return_families import OUT,ROOT,read,save,parent_at,digest
from analyze_chua_return_geometry import sections
from study_chua_return_geometry import mesh


def parent_mesh(parent):
    data,diag=mesh(parent,128,85)
    # The shared mesh routine defaults to the attracting reference multiplier
    # 1/4. A continued short saddle must instead be compared with its own root.
    diag['distance_from_attracting_reference_multiplier']=diag.pop('multiplier_error')
    diag['expected_multiplier']=parent['full_multiplier']
    diag['multiplier_error']=float(abs(mp.mpf(diag['dominant_multiplier'])-mp.mpf(parent['full_multiplier'])))
    diag['reference_type']='continued short saddle'
    return data,diag


def direction_curves(data,record):
    # Start at the published root, where the product is well conditioned.
    eig,V=np.linalg.eig(np.array(record['diagnostics']['full_product']))
    v=np.real(V[:,np.argmax(abs(eig))]);v/=np.linalg.norm(v);floquet=[1.]
    for L in data['local']:v=L@v;floquet.append(float(np.linalg.norm(v)))
    # Also transport the actual sampled maximizing input without reoptimizing
    # its direction at each elapsed time.
    start=record['maximizing_start_index'];N=len(data['local'])
    v=data['basis'][start].T@np.array(record['input_direction']);growth=[float(np.linalg.norm(v))]
    for k in range(N):v=data['local'][(start+k)%N]@v;growth.append(float(np.linalg.norm(v)))
    times=np.r_[data['t'],data['t']+float(data['period']),2*float(data['period'])]
    return dict(floquet_time=np.r_[data['t'],float(data['period'])].tolist(),floquet_amplitude=floquet,
        maximizing_time=(times[start:start+N+1]-times[start]).tolist(),maximizing_amplitude=growth,
        floquet_final=floquet[-1],maximizing_peak=max(growth),maximizing_final=growth[-1],
        scope='Each curve transports one fixed initial direction in the Euclidean transverse quotient norm')


def main():
    forecast=read(OUT/'forecasts-frozen.json');scores=[]
    for t in forecast['targets']:
        id=t['id'];r=read(OUT/f'{id}-reference.json');e=read(OUT/f'{id}-edges.json');g=read(OUT/f'{id}-gain-s8.json')
        measured=dict(width=float(e['width']),gain=g['refined']['maximum']['gain'],period=float(r['period']))
        err={k:measured[k]/t['corrected'][k]-1 for k in ['width','gain']}
        leading={k:measured[k]/t['leading'][k]-1 for k in ['width','gain']}
        err['period']=measured['period']-t['period']
        scores.append(dict(id=id,peaks=t['peaks'],measured=measured,corrected_error=err,leading_error=leading,
            passed=dict(width=abs(err['width'])<=t['gates']['width_relative'],gain=abs(err['gain'])<=t['gates']['gain_relative'],period=abs(err['period'])<=t['gates']['period_absolute'])))
    save(OUT/'scores.json',dict(forecast_sha256=digest(OUT/'forecasts-frozen.json'),rows=scores,
        passed_by_quantity={k:sum(s['passed'][k] for s in scores) for k in ['width','gain','period']},
        error_convention='measured / forecast - 1 for width and gain; measured - forecast for period'))
    calibration=read(OUT/'calibration.json');geometry=[]
    for seq in calibration['sequences']:
        refs=[dict(r['representative'],id=r['id'],n=i+1,organizer=seq['id']) for i,r in enumerate(seq['calibration'])]
        refs += [read(OUT/f"{seq['id']}-n{n}-reference.json") for n in [3,4,5]]
        for r in refs:
            path=OUT/f"{r['id']}-geometry.json"
            if path.exists():geometry.append(read(path));continue
            parent=parent_at(seq['id'],r['beta']);p,t=sections(r,r)
            d=np.linalg.norm(p-np.array(parent['q'][:2],float),axis=1);anchor=int(np.argmax(d));order=np.roll(np.arange(len(p)),-anchor)
            result=dict(id=r['id'],n=r['n'],organizer=seq['id'],parent=parent,
                section_points=p.tolist(),section_times=t.tolist(),distance_to_parent=d.tolist(),
                excursion_anchor=anchor,aligned_distance_to_parent=d[order].tolist(),
                minimum_distance=float(d.min()),section_count=len(p))
            if r['n']==5:
                data,diag=parent_mesh(parent);np.savez_compressed(OUT/f"{r['id']}-parent.npz",**data)
                result['parent_diagnostics']=diag
            save(path,result);geometry.append(result)
            print(r['id'],'parent expansion',parent['full_multiplier'],'nearest',d.min(),flush=True)
    save(OUT/'geometry-summary.json',dict(rows=geometry))
    if not (OUT/'P55-parent.npz').exists():
        old=read(ROOT/'results/chua-knot-atlas/P55-d70-n192.json')
        parent=parent_at('S08',old['beta']);data,diag=parent_mesh(parent)
        np.savez_compressed(OUT/'P55-parent.npz',**data)
        save(OUT/'P55-parent.json',dict(parent=parent,diagnostics=diag))
    for id,stem in [('P55',ROOT/'results/chua-return-geometry/P55-d85-s128'),
                    ('S11-n5',OUT/'S11-n5-gain-s128'),('E05-n5',OUT/'E05-n5-gain-s128')]:
        r=read(stem.with_suffix('.json'))
        with np.load(stem.with_suffix('.npz')) as z:data={k:z[k] for k in z.files}
        curves=direction_curves(data,r);save(OUT/f'{id}-directions.json',curves)
        print(id,'fixed Floquet direction',curves['floquet_final'],'actual maximum',curves['maximizing_peak'],flush=True)


if __name__=='__main__':main()
