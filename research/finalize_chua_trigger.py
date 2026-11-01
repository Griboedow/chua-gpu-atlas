"""Join completed, independently run amplitude bands without hiding a gap."""
import json
import hashlib
import argparse
from pathlib import Path
from decimal import Decimal
from chua_research_followup import OUT,dump


def main(directory=None,local_path=None):
    directory=directory or OUT/'trigger';prefix=directory/'below-prefix-80.json';sources=[prefix]
    data=json.loads(prefix.read_text());rows=list(data['records']);failures=list(data['rejected_enclosures'])
    assert data['central']['passed']
    if 'first_band' in data:
        assert data['first_band']==0 and data['total_bands']==98 and data['working_digits']==80
        next_band=data['last_completed_band']+1
    else:
        # The original prefix predates band metadata. Accept that one retained
        # file by its published digest; do not silently infer missing metadata
        # for an arbitrary incomplete run. Endpoint coverage is checked below.
        published=json.loads((OUT/'trigger/below-certificate-80.json').read_text())
        assert hashlib.sha256(prefix.read_bytes()).hexdigest()==published['source_sha256'][prefix.name]
        next_band=31
    parts=[(path,json.loads(path.read_text())) for path in directory.glob('below-part-80-*.json')]
    for path,chunk in sorted(parts,key=lambda item:item[1]['first_band']):
        first,last=map(int,path.stem.split('-')[-2:])
        assert chunk.get('chunk_complete') and chunk['last_completed_band']==last-1
        assert chunk['first_band']==first==next_band and chunk['total_bands']==98
        assert chunk['central']==data['central'] and chunk['working_digits']==80
        next_band=last
        sources.append(path);rows.extend(chunk['records']);failures.extend(chunk['rejected_enclosures'])
    assert next_band==98
    rows.sort(key=lambda r:Decimal(r['lower']))
    local=json.loads((local_path or directory/'local-certificate-80.json').read_text());crossing=next(r for r in local['records'] if r['angle']==0)
    assert local['passed'] and crossing['passed'] and all(crossing['gates'].values())
    assert Decimal(rows[0]['lower'])==Decimal('.299')
    assert Decimal(rows[-1]['upper'])==Decimal(crossing['bracket'][0])
    assert all(Decimal(a['upper'])==Decimal(b['lower']) for a,b in zip(rows[:-1],rows[1:]))
    assert all(Decimal(r['lower'])<Decimal(r['upper']) and r['label'].startswith('origin') for r in rows)
    report=dict(passed=True,partial=False,working_digits=80,central=data['central'],records=rows,
        rejected_enclosures=failures,threshold_bracket=crossing['bracket'],
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        coverage=dict(adjoining_intervals=len(rows),first=rows[0]['lower'],last=rows[-1]['upper'],
            convention='Decimal strings label shared subdivision endpoints; each original interval is their Arb enclosure.'),
        scope='Every amplitude from zero to the lower bracket endpoint converges to the origin. Together with the local crossing theorem, Astar is the infimum excitation amplitude along (0,A,0).')
    dump(directory/'below-certificate-80.json',report);print('complete lower-amplitude coverage:',len(rows),'intervals')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path)
    parser.add_argument('--local-certificate',type=Path)
    args=parser.parse_args();main(args.directory,args.local_certificate)
