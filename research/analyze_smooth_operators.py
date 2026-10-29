"""Summarize finite-width scans and independently inspect a regular window.

The logistic profile is an explicit modification, not an apparatus calibration.
No interpolation of QR knots is used to locate the smooth Poincare crossings.
"""
import json
import argparse
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from smooth_operators import rhs_and_jac
from audit_coupled_operators import smooth_equilibrium

OUT=Path('results/light-operators')


def sample(width,name):
    source=json.loads((OUT/f'smooth-cpu-1026-{name}.json').read_text())
    p=np.r_[source['parameters'],width];start=np.array(source['final_state'])
    def event(t,y):return y[0]-p[2]
    event.direction=1
    sol=solve_ivp(lambda t,y:rhs_and_jac(y,p)[0],(0,200),start,
        method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01,
        events=event,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    times=np.arange(0,200,.005);states=sol.sol(times).T
    sections=sol.y_events[0];crossings=sol.t_events[0]
    np.savez_compressed(OUT/f'smooth-section-{name}.npz',times=times,states=states,
                        sections=sections,crossing_times=crossings)
    result=dict(width=width,source=f'smooth-cpu-1026-{name}.json',duration=200,
        rtol=1e-12,atol=1e-14,max_step=.01,sample_interval=.005,
        section='x_minus=delta_minus, v_minus>0; five-dimensional section',
        crossings=len(crossings),short_return_distances={str(j):float(np.max(abs(sections[-1]-sections[-1-j])))
            for j in range(1,min(17,len(sections)))})
    candidates=[int(j) for j,error in result['short_return_distances'].items() if error<1e-7]
    if candidates:
        order=min(candidates);period=crossings[-1]-crossings[-1-order];y=sections[-1-order]
        def tangent(t,z):
            f,j=rhs_and_jac(z[:6],p);return np.r_[f,(j@z[6:].reshape(6,6)).ravel()]
        check=solve_ivp(tangent,(0,period),np.r_[y,np.eye(6).ravel()],
            method='DOP853',rtol=1e-12,atol=1e-14,max_step=.005)
        if not check.success:raise RuntimeError(check.message)
        last=check.y[:6,-1];m=check.y[6:,-1].reshape(6,6);f=rhs_and_jac(last,p)[0]
        if abs(f[0])<1e-7:raise RuntimeError('Near-grazing regular control')
        section_matrix=((np.eye(6)-np.outer(f,np.eye(6)[0])/f[0])@m)[1:,1:]
        eig=np.linalg.eigvals(section_matrix)
        result['periodic_control']=dict(section_return_order=order,period=float(period),
            closure_error=float(np.max(abs(last-y))),
            poincare_eigenvalues=[[float(z.real),float(z.imag)] for z in eig],
            spectral_radius=float(max(abs(eig))),
            caveat='Numerical periodic control, not a computer-assisted existence proof; tiny multipliers have limited relative precision')
    return result


def main(plot_only=False):
    summaries=[]
    for stage in ['scout','refined']:
        data=np.load(OUT/f'smooth-width-{stage}.npz');p=data['parameters'];r=data['reference_index']
        good=(data['status']==0)&(data['lle']>.005)&np.all(data['lle_windows']>.002,axis=1)
        for ref in np.unique(r):
            for width in np.unique(p[:,4]):
                mask=(r==ref)&(p[:,4]==width);ok=mask&(data['status']==0)
                summaries.append(dict(stage=stage,reference_index=int(ref),width=float(width),
                    runs=int(mask.sum()),completed=int(ok.sum()),retained=int(good[mask].sum()),
                    leading_exponents=data['lle'][ok].tolist(),window_exponents=data['lle_windows'][ok].tolist()))
    controls=json.loads((OUT/'smooth-width-summary.json').read_text())['section_controls'] if plot_only else [sample(.01,'w001'),sample(.04,'w004')]
    result=dict(groups=summaries,section_controls=controls,
        scout_seed=2026102967,refinement_selection='Two largest completed scout exponents per reference and width, resumed from saved final states',
        caveat='Finite-time numerical screening; unselected scout trajectories were not all refined. No probability, exact boundary, novelty or optical-calibration claim.')
    (OUT/'smooth-width-summary.json').write_text(json.dumps(result,indent=2))
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for ax,ref,label in zip(axes[0],[1026,3054],['A','B']):
        for stage,color,marker in [('scout','#a1b4c7','.'),('refined','#176b91','o')]:
            rows=[row for row in summaries if row['stage']==stage and row['reference_index']==ref]
            for row in rows:ax.scatter([row['width']]*row['completed'],row['leading_exponents'],s=16 if stage=='refined' else 10,color=color,marker=marker)
            ax.scatter([],[],s=16,color=color,marker=marker,label=stage)
        unresolved=[row['width'] for row in summaries if row['stage']=='refined' and row['reference_index']==ref and row['completed']<row['runs']]
        if unresolved:ax.scatter(unresolved,[.02]*len(unresolved),transform=ax.get_xaxis_transform(),marker='x',color='#ae7700',label='unresolved (bottom markers)')
        parameters=np.load(OUT/'coupled-scout.npz')['parameters'][ref]
        widths=np.geomspace(.0005,.1,160)
        eig=[max(z[0] for z in smooth_equilibrium(parameters,w)['eigenvalues']) for w in widths]
        ax.plot(widths,eig,color='#be6948',lw=1,label='equilibrium max Re(eigenvalue)')
        ax.axhline(0,color='.5',lw=.7);ax.set_xscale('log');ax.set_ylim((-.9,.7) if ref==1026 else (-1.15,6))
        ax.set(title=f'Reference {label}: finite optical width',xlabel='Logistic width w',ylabel='Rate per thermal time')
        ax.legend(fontsize=8);ax.grid(alpha=.15)
    for ax,name,title,color in zip(axes[1],['w001','w004'],['w = 0.01: regular-window control','w = 0.04: positive-exponent candidate'],['#386b91','#b65069']):
        data=np.load(OUT/f'smooth-section-{name}.npz');s=data['states'];s=s[-12000:]
        ax.plot(s[:,0],s[:,3],lw=.45,alpha=.75,color=color)
        ax.set(title=title,xlabel='Minus-operator displacement',ylabel='Plus-operator displacement');ax.grid(alpha=.15)
    fig.suptitle('Smooth optical feedback can lose and regain chaotic dynamics\nExplicit logistic modification of a six-state model; lower panels are physical projections',fontsize=13)
    fig.savefig('docs/assets/coupled_operators_smooth_width.png',dpi=170);plt.close(fig)
    print(json.dumps(controls,indent=2))


def regular_control():
    result=sample(.04,'w004-regular')
    (OUT/'smooth-coexisting-periodic-control.json').write_text(json.dumps(result,indent=2))
    chaotic=np.load(OUT/'smooth-section-w004.npz');regular=np.load(OUT/'smooth-section-w004-regular.npz')
    fig=plt.figure(figsize=(12,8),layout='constrained')
    ax=fig.add_subplot(221,projection='3d')
    for data,color,label in [(chaotic,'#b65069','Positive-exponent candidate'),(regular,'#386b91','Regular control')]:
        s=data['states'][-6000:]
        ax.plot(s[:,0],s[:,3],s[:,2]-s[:,5],lw=.6,color=color,label=label)
    ax.set(xlabel='Minus displacement',ylabel='Plus displacement',zlabel='Temperature difference',title='Physical projection of six states')
    ax.tick_params(labelsize=8);ax.legend(fontsize=8)
    ax=fig.add_subplot(222)
    for data,color,label in [(chaotic,'#b65069','Candidate'),(regular,'#386b91','Regular')]:
        s=data['sections'];ax.scatter(s[:,3],s[:,2],s=10,color=color,alpha=.65,label=label)
    ax.set(xlabel='Plus displacement',ylabel='Minus temperature',title='Projection of five-dimensional Poincare section');ax.legend()
    ax=fig.add_subplot(223)
    for data,color in [(chaotic,'#b65069'),(regular,'#386b91')]:
        intervals=np.diff(data['crossing_times']);ax.scatter(intervals[:-1],intervals[1:],s=10,color=color,alpha=.65)
    ax.set(xlabel='Return interval n',ylabel='Return interval n+1',title='Return-time map')
    ax=fig.add_subplot(224)
    for data,color in [(chaotic,'#b65069'),(regular,'#386b91')]:
        mask=data['times']<25;ax.plot(data['times'][mask],data['states'][mask,0],lw=.7,color=color)
    ax.set(xlabel='Thermal time from saved final state',ylabel='Minus displacement',title='Two starts at identical parameters')
    fig.suptitle('Coexisting numerical regimes at reference A, logistic width w = 0.04\nFinite-time evidence; not a proof of an asymptotic chaotic attractor',fontsize=13)
    fig.savefig('docs/assets/coupled_operators_coexistence.png',dpi=170);plt.close(fig)
    print(json.dumps(result,indent=2))


def rotation():
    from matplotlib.animation import FuncAnimation, PillowWriter
    chaotic=np.load(OUT/'smooth-section-w004.npz')['states'][-12000::2]
    regular=np.load(OUT/'smooth-section-w004-regular.npz')['states'][-760:]
    fig=plt.figure(figsize=(6.4,6.4));ax=fig.add_subplot(projection='3d')
    for s,color,label in [(chaotic,'#b65069','Chaotic candidate'),(regular,'#176b91','Stable periodic orbit')]:
        ax.plot(s[:,0],s[:,3],s[:,2]-s[:,5],lw=.65,color=color,label=label)
    ax.set(xlabel='Minus displacement',ylabel='Plus displacement',zlabel='Temperature difference')
    ax.tick_params(labelsize=8);ax.legend(loc='upper left',fontsize=8)
    fig.suptitle('Same parameters, different initial states\nReference A, smooth optical width 0.04',fontsize=12)
    fig.text(.5,.035,'3D projection of six states; numerical chaos candidate\nOnly the camera rotates; the data stay fixed.',ha='center',fontsize=9)
    def frame(i):
        angle=2*np.pi*i/72;ax.view_init(elev=22+14*np.sin(angle),azim=360*i/72,roll=8*np.sin(2*angle))
        return []
    animation=FuncAnimation(fig,frame,frames=72)
    animation.save('docs/assets/coupled_operators_smooth_coexistence.gif',writer=PillowWriter(fps=12),dpi=90)
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--plot-only',action='store_true');parser.add_argument('--regular-control',action='store_true');parser.add_argument('--gif',action='store_true');args=parser.parse_args()
    if args.gif:rotation()
    elif args.regular_control:regular_control()
    else:main(args.plot_only)
