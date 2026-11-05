"""Compare a joint link with the split union of its actual component knots."""
from search_chua_links import ROOT,OUT,read,save,split_screen
from chua_link_diagrams import diagram
from chua_link_bracket import frontier_jones
from chua_link_invariants import split_union_jones
from spherogram import Link
import argparse,random,numpy as np


def compare(record,curves):
    link=Link(record['simplified_pd_code']);link.unlinked_unknot_components=record['final_free_components']
    actual,diagnostic=frontier_jones(link);components=[]
    for xyz in curves:
        k,r=diagram([xyz],record['seed']);random.seed(7682);k.simplify('global',type_III_limit=500)
        polynomial,diag=frontier_jones(k)
        components.append(dict(polynomial=polynomial,pd_code=k.PD_code(),free_components=k.unlinked_unknot_components,diagnostic=diag))
    expected=split_union_jones(c['polynomial'] for c in components)
    return dict(joint_polynomial=actual,split_union_polynomial=expected,distinct=actual!=expected,components=components,
        joint_diagnostic=diagnostic,conclusion=('non_split_pair' if len(curves)==2 else 'not_completely_split') if actual!=expected else 'inconclusive',
        scope='Exact polynomial comparison conditional on numerical diagram extraction; equality is inconclusive.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--folder',default='S02-S08');p.add_argument('--digits',type=int,default=60);p.add_argument('--subdivisions',type=int,default=96);p.add_argument('--seed',type=int,default=0);a=p.parse_args()
    folder=OUT/'overlaps'/a.folder;stem=f'd{a.digits}-n{a.subdivisions}-s{a.seed}';record=read(folder/f'joint-{stem}.json');curves=[]
    for id in record['ids']:
        inverted=id.endswith('-inverted');base=id.removesuffix('-inverted')
        xyz=np.load(folder/f'{base}-d{a.digits}-n{a.subdivisions}.npz')['xyz'];curves.append(-xyz if inverted else xyz)
    result=compare(record,curves);result.update(ids=record['ids'],joint_record=f'joint-{stem}.json')
    save(folder/f'jones-{stem}.json',result);print(a.folder,result['conclusion'],'joint',result['joint_polynomial'],'split',result['split_union_polynomial'],flush=True)


if __name__=='__main__':main()
