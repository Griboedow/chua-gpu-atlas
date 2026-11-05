"""Exact Kauffman bracket by frontier connectivity, with explicit size limits."""
from collections import Counter,defaultdict
import random,time


def crossing_order(pd):
    candidates=[]
    for seed in range(12):
        rng=random.Random(seed);pending=list(range(len(pd)));rng.shuffle(pending);active=set();order=[];width=0;cost=0
        while pending:
            options=[]
            for i in pending:
                # Each edge occurs twice; a loop at this crossing toggles twice.
                new=active.symmetric_difference(k for k,v in Counter(pd[i]).items() if v%2)
                options.append((len(new),-len(active.intersection(pd[i])),pending.index(i),i,new))
            _,_,_,chosen,active=min(options,key=lambda r:r[:3]);order.append(chosen);pending.remove(chosen)
            width=max(width,len(active));cost+=2**(len(active)//2)
        candidates.append((width,cost,order))
    return min(candidates)[:3]


def multiply_delta(poly,k):
    for _ in range(k):
        new=Counter()
        for e,v in poly.items():new[e+2]-=v;new[e-2]-=v
        poly={e:v for e,v in new.items() if v}
    return poly


def frontier_jones(link,max_states=250000,seconds=180):
    pd=link.PD_code();free=int(link.unlinked_unknot_components)
    if not pd:
        from chua_link_diagrams import unlink_jones
        return unlink_jones(free),dict(crossings=0,maximum_frontier=0,maximum_states=1)
    width,cost,order=crossing_order(pd);active=set();start=time.monotonic()
    states={((),free>0):multiply_delta({0:1},max(free-1,0))};counts=[]
    for index in order:
        ports=pd[index];after=active.symmetric_difference(k for k,v in Counter(ports).items() if v%2)
        newstates=defaultdict(Counter)
        for (pairs,normalized),poly in states.items():
            for smooth in [0,1]:
                labels=active.union(ports);parent={k:k for k in labels}
                def find(x):
                    while parent[x]!=x:x=parent[x]
                    return x
                def join(a,b):parent[find(a)]=find(b)
                for a,b in pairs:join(a,b)
                a,b,c,d=ports
                for u,v in ([(a,b),(c,d)] if smooth==0 else [(a,d),(b,c)]):join(u,v)
                groups=defaultdict(list)
                for label in labels:
                    key=find(label)
                    if key not in groups:groups[key]=[]
                    if label in after:groups[key].append(label)
                assert all(len(g) in [0,2] for g in groups.values())
                loops=sum(not g for g in groups.values());newpairs=tuple(sorted(tuple(sorted(g)) for g in groups.values() if g))
                exponent=loops-(1 if loops and not normalized else 0)
                updated=multiply_delta(poly,exponent);target=newstates[newpairs,normalized or loops>0];shift=1 if smooth==0 else -1
                for e,v in updated.items():target[e+shift]+=v
        states={k:{e:v for e,v in p.items() if v} for k,p in newstates.items() if any(p.values())};active=after;counts.append(len(states))
        if len(states)>max_states or time.monotonic()-start>seconds:raise RuntimeError('Frontier calculation resource limit; no invariant conclusion')
    assert not active and set(states)=={((),True)}
    w=link.writhe();sign=-1 if w%2 else 1
    result={str(e-3*w):v*sign for e,v in sorted(states[(),True].items()) if v}
    return result,dict(crossings=len(pd),order=order,maximum_frontier=width,maximum_states=max(counts),states_per_step=counts,elapsed_seconds=time.monotonic()-start)


def controls():
    from spherogram import Link
    from chua_link_diagrams import jones_A
    rows=[]
    for name in ['3_1','4_1','5_1','L2a1','L5a1','L6a4','T(3,5)']:
        for mirror in [False,True]:
            link=Link(name);link=link.mirror() if mirror else link
            value,diagnostic=frontier_jones(link);expected=jones_A(link);assert value==expected,(name,mirror,value,expected)
            rows.append(dict(name=name,mirror=mirror,polynomial=value,diagnostic=diagnostic))
    return rows


if __name__=='__main__':
    from search_chua_links import OUT,save
    save(OUT/'bracket-controls.json',dict(passed=True,rows=controls(),scope='Agreement of frontier connectivity and exhaustive states on standard diagrams'))
    print('14 frontier/exhaustive Jones comparisons passed',flush=True)
