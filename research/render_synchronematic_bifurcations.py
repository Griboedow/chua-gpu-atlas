"""Scientific figure for periodic branches and the finite basin probe."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from synchronematic_observables import pair_signal


def main():
    root=Path('results/synchronematic');p=root/'n4-periodic-branch'
    branch=json.loads((p/'branch.json').read_text());bif=json.loads((p/'bifurcations.json').read_text())
    doubled=json.loads((p/'doubled-branch.json').read_text());folds=json.loads((p/'folds.json').read_text())
    # The fixed-parameter Newton pilot jumped branches below 1.0735.
    # Do not connect its unrelated solutions to the regular branch.
    rows=[r for r in branch['records'] if r['eta']>=1.0735-1e-10]
    x=np.array([r['eta'] for r in rows]);T=np.array([r['period'] for r in rows]);rho=np.array([r['spectral_radius'] for r in rows])
    two=doubled['records'];x2=np.array([r['eta'] for r in two]);T2=np.array([r['period'] for r in two]);r2=np.array([r['spectral_radius'] for r in two])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12,8.3),layout='constrained')
    a,b,c,d=axs.ravel();blue='#2864ac';orange='#cf691b';green='#26836c'
    a.plot(x,T,color='#9ca3af',lw=1)
    for stable,color,label in [(True,green,'Parent: attracting'),(False,orange,'Parent: unstable')]:
        mask=(rho<1)==stable;a.scatter(x[mask],T[mask],c=color,s=18,label=label,zorder=3)
    a.plot(x2,T2/2,color=blue,lw=1.8,label='Doubled orbit: period / 2')
    arc=bif['lower_fold_arc'];a.plot([r['eta'] for r in arc],[r['period'] for r in arc],color=orange,lw=1)
    for f in folds['folds']:a.scatter(f['eta'],f['period'],marker='*',c='black',s=65,zorder=4)
    a.set(xlabel=r'Interaction parameter $\eta$',ylabel='Time (doubled period divided by 2)',title='A  Periodic branches and fold candidates')
    a.legend(fontsize=8,loc='lower right');a.ticklabel_format(useOffset=False,axis='x')
    b.plot(x,rho,color=orange,label='Parent orbit');b.plot(x2,r2,color=blue,label='Doubled orbit')
    b.axhline(1,color='black',ls='--',lw=.9,label='Stability boundary')
    for q in bif['flip_crossings']:b.axvline(q['eta'],color='#999999',ls=':',lw=.9)
    b.set(xlabel=r'Interaction parameter $\eta$',ylabel='Largest transverse multiplier magnitude',title='B  Loss and recovery of parent stability',ylim=(.48,1.25))
    b.legend(fontsize=8,loc='lower right');b.ticklabel_format(useOffset=False,axis='x')
    orbit=np.load(p/'orbit-1.078.npz');q=orbit['state'];mid=(len(q)-1)//2
    full=np.column_stack([np.zeros(len(q)),q[:,:3],q[:,3:]])
    pair=pair_signal(full,0,1);axis=np.sin(2*full[:,4])
    c.plot(pair[:mid+1],axis[:mid+1],color=blue,label='First half of the doubled cycle')
    c.plot(pair[mid:],axis[mid:],color=orange,label='Second half')
    c.scatter(pair[[0,mid]],axis[[0,mid]],s=32,c=[blue,orange],zorder=4)
    c.set(xlabel=r'$\cos(\theta_2-\theta_1)\cos(\chi_2-\chi_1)$',ylabel=r'$\sin(2\chi_1)$',title=r'C  Physical half-cycles at $\eta=1.078$')
    c.legend(fontsize=8,loc='best')
    gpu=np.load(root/'n4-coexistence-gpu/scan.npz');rng=np.random.default_rng(2026102951)
    for center,ids,color,label in [(0,np.arange(24),green,'Near periodic orbit'),(1,np.arange(24,56),blue,'Global angle starts'),(2,np.arange(56,64),orange,'Eight earlier starts')]:
        jitter=rng.uniform(-.16,.16,len(ids));d.scatter(center+jitter,gpu['lle'][ids],c=color,s=22,alpha=.8)
        for xx,i in zip(center+jitter,ids):d.plot([xx,xx],[gpu['lle_windows'][i].min(),gpu['lle_windows'][i].max()],c=color,lw=.6,alpha=.5)
    cpu=root/'n4-coexistence-chaos-cpu/spectrum.json'
    if cpu.exists():
        check=json.loads(cpu.read_text());d.scatter([2.35],[check['spectrum'][0]],marker='*',c='black',s=100,label='Independent CPU example');d.legend(fontsize=8)
    d.axhline(0,color='black',lw=.8);d.set_xticks([0,1,2],['Near cycle\n24 starts','Global\n32 starts','Earlier\n8 starts'])
    d.set(ylabel='Largest finite-time Lyapunov estimate',title=r'D  Different outcomes at $\eta=1.082$',ylim=(-.005,.1),xlim=(-.4,2.65))
    fig.suptitle('A periodic window in the four-oscillator Fourier model\n'
                 r'$F_\eta=G+\eta(F-G)$ changes the interaction law; numerical evidence, not device calibration',fontsize=13)
    out=Path('docs/assets/synchronematic_bifurcations.png');fig.savefig(out,dpi=180);plt.close(fig);print(out)


if __name__=='__main__':main()
