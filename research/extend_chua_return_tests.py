"""Freeze forecasts, then solve two new cycles in a separate evidence directory."""
import argparse
from datetime import datetime, timezone
import json
import time
import mpmath as mp
from study_chua_return_geometry import ROOT, OUT, read, digest
from chua_precision_flights import PrecisionFlights
from verify_chua_layer_predictions import insert_return


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--freeze',action='store_true')
    args=p.parse_args()
    path=OUT/'forecasts-frozen.json'
    if args.freeze:
        if path.exists():raise SystemExit('Frozen forecasts already exist')
        mp.mp.dps=60
        base=read(ROOT/'results/chua-layer-prediction/search-results.json')['records'][-1]
        gain=read(OUT/'P85-optimized.json')['maximum']['gain']
        distance=mp.mpf(read(OUT/'P85-companion-d70.json')['section_distance'])
        width=mp.mpf(base['symmetric_stability_width'])
        expansion=mp.mpf('2.7134317829')
        targets=[]
        for k in [1,2]:
            targets.append(dict(id=f'P{85+10*k}',n=8+k,peaks=85+10*k,
                width=mp.nstr(width/expansion**(2*k),30),
                maximum_gain=float(gain*expansion**k),
                companion_section_distance=mp.nstr(distance/expansion**(2*k),30),
                companion_half_multiplier=1.5,
                relative_gates=dict(width=.10,maximum_gain=.10,companion_section_distance=.15),
                multiplier_absolute_gate=.001))
        sources=[ROOT/'results/chua-layer-prediction/search-results.json',
                 OUT/'P85-optimized.json',OUT/'P85-companion-d70.json',OUT/'protocol.md']
        result=dict(created_utc=datetime.now(timezone.utc).isoformat(),targets=targets,
                    expansion=str(expansion),source_sha256={s.relative_to(ROOT).as_posix():digest(s) for s in sources},
                    scope='New prospective tests of gain, width and companion separation at 95/105 peaks; existing 85-peak data calibrate the constants.',
                    separation_rule='Same initial normalized x=0 crossing inherited during return insertion from P85; predicts Lambda^(-2k), not the originally conjectured uniform Lambda^(-k) law at every phase.',
                    limitation='Root locations are corrected by shooting; this is not a prospective point-location or basin-escape forecast.')
        path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(json.dumps(result,indent=2));return
    forecast=read(path)
    previous=read(ROOT/'results/chua-layer-prediction/search-results.json')['records'][-1]
    edge=previous['fold'];regions=previous['regions']
    for target in forecast['targets']:
        id=target['id'];root_path=OUT/f'{id}-reference.json';edge_path=OUT/f'{id}-edges.json'
        if root_path.exists() and edge_path.exists():
            edge=read(edge_path)['fold'];regions=read(root_path)['regions'];continue
        started=time.monotonic()
        newregions,targets,guess,insertion=insert_return(edge,regions)
        f=PrecisionFlights(newregions,targets,True,75)
        # The exact tangency interval supplies an initial location only; Newton
        # correction is allowed and the location is explicitly not a forecast.
        limit=mp.mpf('14.30635909996146')
        beta=limit+(mp.mpf(edge['beta'])-limit)/mp.mpf(forecast['expansion'])
        fold=f.boundary(beta,guess,1,max_iterations=36)
        symmetry=f.boundary(fold['beta'],fold['q'],-1)
        representative=f.boundary(fold['beta'],fold['q'],'.5')
        width=mp.mpf(fold['beta'])-mp.mpf(symmetry['beta'])
        representative.update(id=id,symmetric=True,regions=newregions,targets=targets,
                              alpha='9',m0='-8/7',m1='-5/7',gamma='0',
                              scope='Targeted extension for return-geometry tests; not added to the 34-reference display catalogue')
        edges=dict(id=id,fold=fold,symmetry=symmetry,width=mp.nstr(width,65),
                   insertion_time=insertion,forecast_sha256=digest(path),
                   elapsed_seconds=time.monotonic()-started)
        root_path.write_text(json.dumps(representative,indent=2)+'\n',encoding='utf-8',newline='\n')
        edge_path.write_text(json.dumps(edges,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(id,'width',edges['width'],'period',representative['period'],'seconds',edges['elapsed_seconds'],flush=True)
        edge=fold;regions=newregions


if __name__=='__main__':main()
