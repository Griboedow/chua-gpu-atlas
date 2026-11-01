"""Join completed, independently run amplitude bands without hiding a gap."""
import json
import hashlib
from decimal import Decimal
from chua_research_followup import OUT,dump


def main():
    directory=OUT/'trigger';prefix=directory/'below-prefix-80.json';sources=[prefix]
    data=json.loads(prefix.read_text());rows=list(data['records']);failures=list(data['rejected_enclosures'])
    # The prefix was saved only after completing original bands 0 through 30.
    for first,last in [(31,43),(43,55),(55,67),(67,79),(79,89),(89,98)]:
        path=directory/f'below-part-80-{first}-{last}.json';chunk=json.loads(path.read_text())
        assert chunk.get('chunk_complete') and chunk['last_completed_band']==last-1
        assert chunk['first_band']==first and chunk['total_bands']==98
        sources.append(path);rows.extend(chunk['records']);failures.extend(chunk['rejected_enclosures'])
    rows.sort(key=lambda r:Decimal(r['lower']))
    local=json.loads((directory/'local-certificate-80.json').read_text());crossing=next(r for r in local['records'] if r['angle']==0)
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


if __name__=='__main__':main()
