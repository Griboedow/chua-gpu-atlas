"""Periodic control and physical-coordinate figures for the two-strain study."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from two_strain import ROOT,decode,encode,log_rhs,log_jac,save


def cycle_control():
    data=np.load(ROOT/'cpu-small-cycle-log.npz');p=data['parameters']
    start=encode(data['sections'][-2],p);period=float(np.diff(data['section_times'])[-1])
    def integrate(q,t):
        def field(t,z):return np.r_[log_rhs(z[:5],p),(log_jac(z[:5],p)@z[5:].reshape(5,5)).ravel()]
        sol=solve_ivp(field,(0,t),np.r_[q,np.eye(5).ravel()],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.25)
        if not sol.success:raise RuntimeError(sol.message)
        return sol.y[:5,-1],sol.y[5:,-1].reshape(5,5)
    def residual(v):
        q=np.r_[start[0],v[:4]];end,m=integrate(q,v[4]);return end-q
    shot=root(residual,np.r_[start[1:],period],tol=1e-9)
    q=np.r_[start[0],shot.x[:4]];period=float(shot.x[4]);end,m=integrate(q,period)
    f=log_rhs(end,p);normal=np.eye(5)[0]
    section_map=((np.eye(5)-np.outer(f,normal)/f[0])@m)[1:,1:]
    eig=np.linalg.eigvals(section_map);closure=float(np.max(abs(end-q)))
    result=dict(parameters=p.tolist(),period=period,section_log_state=q.tolist(),
        closure_log_state=closure,poincare_eigenvalues=[[float(z.real),float(z.imag)] for z in eig],
        spectral_radius=float(max(abs(eig))),stable_numerically=bool(closure<1e-8 and max(abs(eig))<1),
        scope='Numerical periodic orbit and Poincare stability; smallest strongly contracting multiplier is limited by double precision',
        root_solver_success=bool(shot.success))
    (ROOT/'small-cycle-control.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
    if not result['stable_numerically']:raise RuntimeError('Periodic control failed')


def densify():
    data=np.load(ROOT/'cpu-refined32-log.npz');p=data['parameters']
    sol=solve_ivp(lambda t,q:log_rhs(q,p),(0,15000),data['final_log_state'],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.5,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    times=np.arange(0,15000.01,.5);q=sol.sol(times).T;states=np.array([decode(v,p) for v in q])
    save(ROOT/'candidate-render-trace',dict(parameters=p,times=times,log_states=q,states=states,
        scope='Additional CPU trace for physical projections; not an independent full-spectrum run',rtol=1e-12,atol=1e-14,max_step=.5))


def plot(gif=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation,PillowWriter
    from matplotlib.ticker import MaxNLocator
    def ticks(ax):
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_major_locator(MaxNLocator(4))
    assets=Path('docs/assets');assets.mkdir(exist_ok=True)
    scout=np.load(ROOT/'slice-scout.npz');ref=np.load(ROOT/'slice-refined.npz')
    cpu=np.load(ROOT/'cpu-refined32-log.npz');trace=np.load(ROOT/'candidate-render-trace.npz')
    fig=plt.figure(figsize=(12,10),layout='constrained');grid=fig.add_gridspec(3,2)
    ax=fig.add_subplot(grid[0,:]);r2=scout['parameters'][:,2]/(scout['parameters'][:,4]+scout['parameters'][:,0])
    for kind,color in [('coexistence','#347fa5'),('invasion','#d47337')]:
        indices=scout['start_kinds']==kind;ax.scatter(r2[indices],scout['lle'][indices]*1e3,s=12,color=color,alpha=.5,label=f'Scout, near {kind}')
    x=ref['parameters'][:,2]/(ref['parameters'][:,4]+ref['parameters'][:,0]);y=ref['lle']*1e3
    ax.scatter(x,y,marker='x',c='black',s=25,label='400,000-unit refinement')
    ax.axhline(0,color='gray',lw=.7);ax.set(xlabel='Strain-2 reproduction number',ylabel='Largest FT exponent × 1,000',title='Original compartment tangent; fixed diagonal scaling')
    ax.legend(ncol=3,fontsize=8);ax.grid(alpha=.2)
    ax=fig.add_subplot(grid[1,0],projection='3d');s=trace['states']
    ax.plot(s[:,0],1e3*s[:,1],1e3*s[:,2],lw=.45,color='#227b91',alpha=.7)
    ax.set(xlabel='Susceptible S',ylabel='1,000 × I₁',zlabel='1,000 × I₂',title='Physical fractions: 3D projection of 5D flow')
    ticks(ax)
    ax=fig.add_subplot(grid[1,1]);sec=cpu['sections'];ax.scatter(sec[:,1],sec[:,2],c=np.arange(len(sec)),s=13,cmap='viridis')
    ax.set(xlabel='I₁ on section',ylabel='I₂ on section',title='Projection of the four-dimensional section S=0.4');ax.ticklabel_format(style='sci',scilimits=(-3,3));ax.grid(alpha=.2)
    ax=fig.add_subplot(grid[2,0]);ax.plot(trace['times'],s[:,1],label='I₁');ax.plot(trace['times'],s[:,2],label='I₂',alpha=.8)
    ax.set(xlabel='Time / primary recovery time',ylabel='Population fraction',title='Unequal outbreak intervals and amplitudes');ax.legend();ax.grid(alpha=.2)
    ax=fig.add_subplot(grid[2,1]);ax.plot(trace['times'],np.log10(trace['parameters'][0])+trace['log_states'][:,1]/np.log(10),label='log₁₀ I₁');ax.plot(trace['times'],np.log10(trace['parameters'][0])+trace['log_states'][:,2]/np.log(10),label='log₁₀ (I₂+I₁₂)')
    ax.set(xlabel='Time / primary recovery time',ylabel='Log₁₀ population fraction',title='Very small prevalence limits physical interpretation');ax.legend();ax.grid(alpha=.2)
    fig.suptitle('Two-strain outbreak model: a checked numerical chaos candidate\nRelated epidemic chaos and possibly aperiodic source regimes are already published',fontsize=14)
    fig.savefig(assets/'two_strain_candidate.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(10,6),layout='constrained',sharex=True)
    for ax,kind in zip(axes,['coexistence','invasion']):
        d=np.load(ROOT/f'source-r2-1.9-{kind}.npz');ax.plot(d['times'],d['states'][:,1],label='I₁');ax.plot(d['times'],d['states'][:,2],label='I₂')
        ax.set(ylabel='Population fraction',title=f'Identical parameters, initial state near {kind}');ax.legend();ax.grid(alpha=.2)
    axes[-1].set_xlabel('Time / primary recovery time');fig.suptitle('Source-parameter reproduction at R₂=1.9: small and large oscillations')
    fig.savefig(assets/'two_strain_source_controls.png',dpi=160);plt.close(fig)
    fig=plt.figure(figsize=(12,5),layout='constrained')
    ax=fig.add_subplot(121,projection='3d');ax.plot(s[:,0],1e3*s[:,1],1e3*s[:,2],lw=.5,color='#227b91',alpha=.7)
    ax.set(xlabel='S',ylabel='1,000 × I₁',zlabel='1,000 × I₂',title='Two primary infection fractions')
    ticks(ax)
    ax=fig.add_subplot(122,projection='3d');ax.plot(s[:,0],s[:,4],1e3*s[:,1],lw=.5,color='#855c9b',alpha=.7)
    ax.set(xlabel='S',ylabel='Recovered R₁',zlabel='1,000 × I₁',title='Susceptibility, immunity and infection')
    ticks(ax)
    fig.suptitle('Two physical projections of the same five-state trajectory\nAxis multipliers only change displayed units; no new topology is claimed')
    fig.savefig(assets/'two_strain_projections.png',dpi=150);plt.close(fig)
    if gif:
        fig=plt.figure(figsize=(6,6),facecolor='white');ax=fig.add_subplot(projection='3d')
        ax.plot(s[::2,0],1e3*s[::2,1],1e3*s[::2,2],lw=.55,color='#227b91',alpha=.7)
        ax.set(xlabel='S',ylabel='1,000 × I₁',zlabel='1,000 × I₂',title='Two-strain candidate: physical fractions\n3D projection of a five-state model')
        ticks(ax)
        def frame(i):ax.view_init(elev=25+18*np.sin(2*np.pi*i/72),azim=5*i);return []
        animation=FuncAnimation(fig,frame,frames=72,interval=80)
        animation.save(assets/'two_strain_rotation.gif',writer=PillowWriter(fps=12),dpi=100);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--cycle',action='store_true');parser.add_argument('--densify',action='store_true');parser.add_argument('--gif',action='store_true');args=parser.parse_args()
    if args.cycle:cycle_control()
    elif args.densify:densify()
    else:plot(args.gif)
