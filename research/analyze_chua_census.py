"""Connect verified windows, the period-adding sequence, and the original mesh."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import label,generate_binary_structure

from trace_chua_window import Chua
from chua_layer_bifurcations import Flights,from_orbit
from continue_chua_families import boundary


ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/chua-layer-census')
    args=p.parse_args();out=args.output
    records=[]
    for source in ['cycles','cycles-extra','adding-cycles']:
        for r in json.loads((out/(source+'.json')).read_text())['records']:
            r['source']=source;records.append(r)
    verified=[r for r in records if r['status']=='verified_attracting_cycle']
    # Two independently checked representatives belong to the same continued
    # period-nine window. Preserve both raw records, count that window once.
    anchors=sorted([r for r in verified if abs(r['beta']-14.259953125)>1e-10],key=lambda r:r['beta'])
    tight=np.load(out/'tight.npz');beta=tight['beta'];period=tight['period'];gaps=[]
    for left,right in zip(anchors[:-1],anchors[1:]):
        available=np.flatnonzero((beta>left['beta'])&(beta<right['beta'])&(period==0)&(tight['classes']==4).all(axis=1))
        if len(available):
            i=available[np.argmax(tight['lle_windows'][available].min(axis=(1,2)))]
            gaps.append(dict(left_beta=left['beta'],right_beta=right['beta'],witness_beta=float(beta[i]),
                             four_window_lle=tight['lle_windows'][i].tolist(),
                             scope='A finite-time three-start chaotic screen between verified windows; not exclusion of hidden periodic basins.'))
        else:gaps.append(dict(left_beta=left['beta'],right_beta=right['beta'],status='no_saved_gap_witness'))
    model=Chua(rtol=2e-12,step=.02)
    central=json.loads((ROOT/'results/chua-layers/affine-central.json').read_text())
    parent=Flights(model,central['regions'],central['targets'],True)
    parent_q=np.array(central['boundaries']['-1']['q'])
    adding=[]
    for b in [14.3036875,14.30525,14.305921875,14.30619]:
        r=next(r for r in verified if abs(r['beta']-b)<1e-10)
        T=r['cycle']['period'];v=[*r['cycle']['q'][:2],T]
        family,q=from_orbit(model,b,v,True)
        row=dict(beta=b,maxima=r['peak_count'],period=T,section_period=r['section_period'])
        for name,target in [('fold',1),('symmetry',-1)]:
            bb,qq=boundary(family,b,q,target);rr=family.record(bb,qq,True)
            if rr['region_violation']>1e-7:raise RuntimeError('Invalid adding-family itinerary')
            row[name]=rr
        row['symmetric_core_width']=abs(row['fold']['beta']-row['symmetry']['beta'])
        pp=parent.solve(b,parent_q);pr=parent.record(b,pp)
        sol=model.integrate(b,[0.,*r['cycle']['q'][:2]],T)
        returns=sol.y_events[2][sol.t_events[2]<T-1e-6]
        distances=np.linalg.norm(returns[:,1:]-pp[:2],axis=1)
        row.update(parent_five_cycle=pr,section_returns=returns.tolist(),
                   distances_to_parent_section=distances.tolist(),minimum_parent_distance=float(min(distances)))
        adding.append(row)
        print('adding',b,r['peak_count'],T,row['symmetric_core_width'],row['minimum_parent_distance'],flush=True)
    # Assign only the attracting five-maximum cores; doubled daughters and
    # satellites are deliberately left outside this narrow comparison.
    atlas_path=ROOT/'outputs/detailed/window_zoom_scan.npz'
    atlas=np.load(atlas_path);params=atlas['parameters'];classes=atlas['classes']
    regular=np.isin(classes,[1,2]).all(axis=1)&(classes==2).any(axis=1)
    old_shape=(len(atlas['alpha']),len(atlas['beta']),len(atlas['m0']))
    topology=[]
    for connectivity in [1,2,3]:
        labels,n=label(regular.reshape(old_shape),generate_binary_structure(3,connectivity))
        sizes=np.bincount(labels.ravel())[1:]
        topology.append(dict(neighbors=[6,18,26][connectivity-1],sampled_components=int(n),largest_component_nodes=int(sizes.max())))
    first=json.loads((ROOT/'results/chua-layers/extended.json').read_text())['surfaces']
    second=json.loads((out/'second-surfaces.json').read_text())
    aa=np.linspace(8.94,9.06,7);mm=np.linspace(-8/7-.006,-8/7+.006,7);masks=[]
    coverage=[]
    for name,rows in [('first_five',first),('second_five',second)]:
        folds=np.array([r['fold']['beta'] for r in rows]).reshape(7,7)
        pds=np.array([r['period_doubling']['beta'] for r in rows]).reshape(7,7)
        coords=params[:,[0,2]]
        upper=RegularGridInterpolator((aa,mm),folds,bounds_error=False,fill_value=None)(coords)
        lower=RegularGridInterpolator((aa,mm),pds,bounds_error=False,fill_value=None)(coords)
        inside=(params[:,1]>lower)&(params[:,1]<upper);masks.append(inside)
        coverage.append(dict(family=name,regular_nodes_inside=int((inside&regular).sum()),total_nodes_inside=int(inside.sum())))
    union=np.logical_or.reduce(masks)
    summaries=[dict(beta=r['beta'],maxima=r['peak_count'],section_period=r['section_period'],period=r['cycle']['period'],
                    rho=r['independent_variational_rho'],closure=r.get('radau_refinement',{}).get('closure',r['radau_closure']),source=r['source']) for r in anchors]
    report=dict(verified_anchor_windows=len(anchors),anchors=summaries,gap_witnesses=gaps,adding_sequence=adding,
                old_grid=dict(shape=old_shape,regular_nodes=int(regular.sum()),sampled_connectivity=topology,
                              source_sha256=hashlib.sha256(atlas_path.read_bytes()).hexdigest(),five_maximum_core_coverage=coverage,
                              union_regular_nodes=int((regular&union).sum()),union_fraction=float((regular&union).sum()/regular.sum()),
                              interpolation_note='Surface boundaries are bilinearly interpolated from 7x7 continuation points; the comparison is approximate.'),
                scope='Verified local periodic windows along one transect, not the number of globally connected 3D sheets.',
                infinity_status='No finite total or infinitude is proved. Four period-adding members and shrinking cores support an accumulating hierarchy hypothesis.')
    (out/'census-analysis.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    np.savez_compressed(out/'atlas-comparison.npz',alpha=atlas['alpha'],beta=atlas['beta'],m0=atlas['m0'],
                        regular=regular.reshape(old_shape),first_core=masks[0].reshape(old_shape),second_core=masks[1].reshape(old_shape))
    print('count',len(anchors),'gap witnesses',sum('witness_beta' in g for g in gaps),'old grid',report['old_grid'],flush=True)


if __name__=='__main__':main()
