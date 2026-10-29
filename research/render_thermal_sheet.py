"""Illustrate transient false positives and the domain of a mathematical candidate."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from scipy.integrate import solve_ivp

from thermomechanical_sheet import equilibrium,field_jacobian


def main():
    base=Path('results/thermomechanical-sheet')
    cpu=json.loads((base/'cpu-reference.json').read_text());p=np.array(cpu['parameters'])
    times=np.linspace(0,50,25001)
    sol=solve_ivp(lambda t,y:field_jacobian(y,p)[0],(0,50),cpu['final_state'],method='DOP853',
                  rtol=1e-11,atol=1e-13,max_step=.01,t_eval=times)
    if not sol.success:raise RuntimeError(sol.message)
    np.savez_compressed(base/'reference-trace.npz',time=times,state=sol.y.T,parameters=p)
    trace_info=dict(source='cpu-reference.json final_state',duration=50,samples=len(times),
                    method='DOP853',rtol=1e-11,atol=1e-13,max_step=.01,
                    interpretation='A finite mathematical orbit outside the physical sheet domain; no physical attractor claim.')
    (base/'reference-trace.json').write_text(json.dumps(trace_info,indent=2))
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    a=np.load(base/'scout/scan.npz');b=np.load(base/'refined/scan.npz')
    sm=json.loads((base/'scout/summary.json').read_text());bm=json.loads((base/'refined/summary.json').read_text())
    reference=5906;j=bm['parent_indices'].index(reference)
    for result,meta,index,label,color in [(a,sm,reference,'Scout','#d98118'),(b,bm,j,'Longer, tighter follow-up','#136a9b')]:
        setting=meta['settings'];duration=setting['duration']
        centers=setting['transient']+duration*(np.arange(4)+.5)/4
        axes[0].plot(centers,result['lle_windows'][index],'o-',label=label,color=color)
    q=a['parameters'][reference];rate=np.max(np.linalg.eigvals(field_jacobian(equilibrium(q),q)[1]).real)
    axes[0].axhline(0,color='.4',linewidth=.7)
    axes[0].axhline(rate,color='.5',ls='--',label='Linear growth at equilibrium')
    axes[0].set(xlabel='Dimensionless time',ylabel='Largest finite-time exponent',title='Small-start candidate: longer observation matters')
    axes[0].legend(fontsize=8,loc='best');axes[0].grid(alpha=.15)
    x,y=sol.y[0],sol.y[2];inside=(x>=0)&(y>=0)&(y<=1)&(1-p[3]*y>0)
    axes[1].add_patch(Rectangle((0,0),x.max()*1.1,min(1,1/p[3]),facecolor='#d5efd9',alpha=.8,label='Admissible projection'))
    axes[1].plot(x,y,color='.65',linewidth=.3,alpha=.7)
    for flag,color,label in [(True,'#136a9b','Inside domain'),(False,'#b43336','Outside domain')]:
        mask=(inside==flag)&(np.arange(len(x))%4==0)
        axes[1].scatter(x[mask],y[mask],s=1.5,color=color,label=label,alpha=.6,rasterized=True)
    spectrum=', '.join(f'{v:.4f}' for v in cpu['spectrum'])
    axes[1].text(.02,.98,'CPU spectrum: '+spectrum,transform=axes[1].transAxes,va='top',fontsize=8,
                 bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    axes[1].set(xlabel='Displacement X',ylabel='Temperature rise Y',title='Strong candidate: outside the physical sheet domain')
    axes[1].legend(fontsize=8,loc='lower right',markerscale=3);axes[1].grid(alpha=.15)
    fig.suptitle('Thermomechanical sheet: numerical evidence and its limits',fontsize=13)
    fig.savefig('docs/assets/thermal_sheet_screening.png',dpi=160)
    plt.close(fig)


if __name__=='__main__':main()
