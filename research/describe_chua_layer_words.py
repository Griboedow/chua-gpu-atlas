"""Describe the existing 32 roots and expose collisions of coarse symbolic codes."""
import csv
import json
from collections import defaultdict
from predict_chua_layers import ROOT,OUT


def canonical(word):
    inverted=word.translate(str.maketrans('LR','RL'))
    return min(w[i:]+w[:i] for w in [word,inverted] for i in range(len(w)))


def primitive(word):
    for k in range(1,len(word)+1):
        if len(word)%k==0 and word==word[:k]*(len(word)//k):return word[:k],len(word)//k


def main():
    registry=json.loads((ROOT/'results/chua-layer-atlas/family-registry.json').read_text())['families']
    rows=[];groups=defaultdict(list);primitives=defaultdict(list)
    for r in sorted(registry,key=lambda x:x['representative']['beta']):
        short=''.join('R' if reg==1 else 'L' for reg in r['regions'] if reg)
        full=short+short.translate(str.maketrans('LR','RL')) if r['symmetric'] else short
        word=canonical(full);prim,reps=primitive(word)
        representative=r['representative'];row=dict(id=r['id'],beta=representative['beta'],period=representative['period'],peaks=r['peak_count'],
            section_returns=r['section_returns'],symmetric=r['symmetric'],
            half_or_full_outer_word=short,canonical_full_outer_word=word,
            primitive_outer_word=canonical(prim),symbolic_repetitions=reps,
            regions=r['regions'],residence_times=representative['q'][2:],section_state=representative['q'][:2],
            full_multipliers=representative['full_multipliers'],branch_relations=r['sources'])
        rows.append(row);groups[word].append(r['id']);primitives[canonical(prim)].append(r['id'])
    summary=dict(reference_families=len(rows),distinct_full_outer_words=len(groups),
        distinct_primitive_outer_words=len(primitives),
        full_word_collisions=[ids for ids in groups.values() if len(ids)>1],
        primitive_word_groups=[dict(word=k,families=v) for k,v in primitives.items()],
        convention='Cyclic phase rotations and global inversion are identified; time reversal is not. Counts concern this finite catalogue only.',
        interpretation='The L/C/R region partition is not a generating code here. Residence times, section coordinates and branch identity must be retained.',families=rows)
    (OUT/'motion-descriptions.json').write_text(json.dumps(summary,indent=2)+'\n')
    fields=['id','beta','period','peaks','section_returns','symmetric','canonical_full_outer_word','primitive_outer_word','symbolic_repetitions']
    with (OUT/'motion-descriptions.csv').open('w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print({k:v for k,v in summary.items() if k not in ['families','primitive_word_groups']})


if __name__=='__main__':main()
