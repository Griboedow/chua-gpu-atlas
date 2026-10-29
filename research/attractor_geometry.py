"""Geometry baseline for existing candidates, reusing saved verified trajectories.

No new numerical orbit is integrated here. Observables are labeled explicitly;
projection appearance and finite-spectrum dimension estimates are not novelty tests.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from synchronematic_observables import collective_order,pair_signal


def kaplan_yorke(spectrum):
    values=np.sort(np.asarray(spectrum,float))[::-1];total=np.cumsum(values)
    good=np.flatnonzero(total>=0)
    if not len(good):return 0.
    j=int(good[-1])+1
    return float(j if j==len(values) else j+total[j-1]/abs(values[j]))


def main():
    root=Path('results');out=root/'attractor-geometry';out.mkdir(parents=True,exist_ok=True)
    sync=np.load(root/'synchronematic/n4-cpu-extended/trace.npz')
    sm=json.loads((root/'synchronematic/n4-cpu-extended/spectrum.json').read_text())
    y=sync['state'];C,Q=collective_order(y)
    xyz=np.column_stack([Q.real,Q.imag,C]);alternative=np.column_stack([C,abs(Q),pair_signal(y,0,1)])
    slider=np.load(root/'lce-slider/geometry-check/trace.npz');p=slider['parameters'];q=p[2]/(1+p[2]);dense=slider['trajectory']
    geometric=dense[:,1:4]/np.array([q,np.sqrt(p[0])*q,q])
    sec=np.load(root/'lce-slider/analytic-2222/trace.npz')['section'];section=sec[:,1:3]/np.array([np.sqrt(p[0])*q,q])
    lm=json.loads((root/'lce-slider/analytic-2222/spectrum.json').read_text())
    full_lce=sorted([*lm['transverse_spectrum'],0.],reverse=True)
    summary=dict(purpose='Baseline morphology audit; no newly discovered attractor claimed',
                 synchronematic=dict(source='results/synchronematic/n4-cpu-extended',state_dimension_after_symmetry_reduction=7,
                     spectrum=sm['spectrum'],kaplan_yorke_estimate=kaplan_yorke(sm['spectrum']),
                     estimate_by_window=[kaplan_yorke(v) for v in sm['spectrum_windows']],
                     projection=['Re(mean(exp(2i*chi)))','Im(mean(exp(2i*chi)))','sqrt((abs(mean(exp(i*(theta+chi))))**2+abs(mean(exp(i*(theta-chi))))**2)/2)'],
                     alternative_color='cos(theta_2-theta_1)*cos(chi_2-chi_1)',
                     invariance='Independent (theta_i,chi_i)->(theta_i+pi,chi_i+pi), and uniform phase shift',
                     warning='Three physical collective observables of a seven-dimensional angular flow; not an embedding or a three-state ODE. Corrects earlier gauge-dependent R=abs(mean(exp(i*theta))) axis.'),
                 slider=dict(source='results/lce-slider/analytic-2222',state_dimension=3,flow_spectrum=full_lce,
                     kaplan_yorke_estimate=kaplan_yorke(full_lce),neutral_exponent='Inserted zero flow exponent; two measured values are transverse exponents per elapsed time',
                     state_view=['u/q','v/(sqrt(k)*q)','phi/q'],section='u=0, upward crossings; coordinates v/(sqrt(k)*q),phi/q'),
                 limitations='Kaplan-Yorke estimates use finite-time spectra and are not measured or rigorously proved Hausdorff dimensions. Different spectra or pictures alone do not establish a new topological type.')
    (out/'baseline.json').write_text(json.dumps(summary,indent=2))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(12,9),layout='constrained')
    ax=fig.add_subplot(221,projection='3d');points=xyz[::8]
    ax.scatter(*points.T,s=.6,c='#227d9c',alpha=.35,rasterized=True)
    ax.set(xlabel='Re Q',ylabel='Im Q',zlabel='Coherence C',title='7D candidate: three gauge-invariant observables')
    ax.view_init(24,38)
    ax=fig.add_subplot(222,projection='3d');ax.plot(*geometric.T,lw=.35,c='#b56835',alpha=.6)
    ax.set(xlabel='u / q',ylabel='v / (sqrt(k) q)',zlabel='Cis fraction / q',title='3D slider candidate: full state, rescaled units');ax.view_init(24,38)
    ax=fig.add_subplot(223);ax.scatter(alternative[::4,0],alternative[::4,1],s=.5,c=alternative[::4,2],cmap='coolwarm',vmin=-1,vmax=1,alpha=.4,rasterized=True)
    ax.set(xlabel='Synchronematic coherence C',ylabel='Nematic order |Q|',title='Same 7D trajectory, different physical observables')
    ax.text(.02,.03,'Color: cos(delta theta) cos(delta chi), pair 1-2\nThis is not a Poincare section.',transform=ax.transAxes,fontsize=9)
    ax=fig.add_subplot(224);ax.scatter(section[:,0],section[:,1],s=.5,c='#b56835',alpha=.4)
    ax.set(xlabel='Upward crossing v / (sqrt(k) q)',ylabel='Cis fraction / q',title='Slider Poincare section: 8,000 event returns')
    fig.suptitle('Existing chaos candidates: geometry baseline, scientific novelty unconfirmed',fontsize=14)
    fig.savefig('docs/assets/attractor_geometry_baseline.png',dpi=160);plt.close(fig)
    fig=plt.figure(figsize=(6.5,5.8));ax=fig.add_subplot(projection='3d')
    ax.scatter(*points.T,s=.8,c='#227d9c',alpha=.35)
    ax.set(xlabel='Re Q',ylabel='Im Q',zlabel='Coherence C',xlim=(-1,1),ylim=(-1,1),zlim=(0,1))
    ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1]);ax.set_zticks([0,.5,1])
    fig.suptitle('Existing synchronematic candidate\n3 gauge-invariant observables of a 7D flow',fontsize=12)
    fig.text(.5,.02,'Camera rotation only | Novelty unconfirmed',ha='center',fontsize=10)
    def frame(i):
        angle=2*np.pi*i/72
        ax.view_init(elev=25+20*np.sin(2*angle),azim=360*i/72)
        return []
    animation=FuncAnimation(fig,frame,frames=72,interval=80)
    animation.save('docs/assets/synchronematic_projection_rotation.gif',writer=PillowWriter(fps=12),dpi=95)
    plt.close(fig);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
