"""Locate numerical saddle cycles in the controlled hidden-chaos candidates."""
import json
import numpy as np
from chua_research_followup import ROOT,OUT,dump
from chua_lobe_map import LobeMap
from trace_chua_window import Chua


def main():
    source=json.loads((OUT/'hidden-search/outer-screen.json').read_text())['candidates'];results=[]
    for row in source:
        p=row['parameters'];f=LobeMap(p);initial=np.array(row['final']);initial*=np.sign(initial[0])
        sol=f.model.integrate(p['beta'],initial,100)
        crossings=sol.y_events[1];points=np.array([q[1:] for q in crossings if f.model.field(p['beta'],q)[0]>0])
        guesses=list(points[::max(1,len(points)//20)]);roots=[];failures=0
        if len(points):
            for y in np.linspace(points[:,0].min(),points[:,0].max(),7):
                for z in np.linspace(points[:,1].min(),points[:,1].max(),7):guesses.append([y,z])
        for guess in guesses:
            try:
                q,J,aux=f.fixed(guess,p['beta'])
                if all(np.linalg.norm(q-np.array(r['point']))>1e-6 for r in roots):
                    e,V=np.linalg.eig(J);order=np.argsort(-abs(e));e,V=e[order],V[:,order]
                    roots.append(dict(point=q.tolist(),jacobian=J.tolist(),
                        multipliers=[[float(v.real),float(v.imag)] for v in e],**aux))
            except (ValueError,RuntimeError,np.linalg.LinAlgError):failures+=1
        item=dict(parameter_index=row['parameter_index'],parameters=p,roots=roots,
            section_points=points.tolist(),failed_starts=failures,total_starts=len(guesses))
        results.append(item);dump(OUT/'hidden-search/lobe-saddles.json',dict(records=results,
            scope='Numerical positive-lobe fixed points; no completeness, interval certification, or homoclinic conclusion.'))
        print(row['parameter_index'],'roots',[(r['point'],r['multipliers']) for r in roots],flush=True)


if __name__=='__main__':main()
