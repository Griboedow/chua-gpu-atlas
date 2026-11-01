"""Figures from the global storage argument and recorded periodic-orbit checks."""
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.publication/matplotlib'))
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from plot_chua_slope_refinement import sample_orbit

OUT = ROOT/'results/chua-hidden-boundary'
ASSETS = ROOT/'docs/assets'
GREEN, BLUE, RED, GRAY = '#258b68', '#197da1', '#bd536b', '#77838c'
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'axes.titleweight': 'semibold',
    'text.color': '#243747', 'axes.labelcolor': '#243747', 'axes.spines.top': False,
    'axes.spines.right': False, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})


def read(name):
    return json.loads((OUT/name).read_text())


def save(fig, name):
    fig.savefig(ASSETS/(name+'.png'), dpi=185)
    fig.savefig(OUT/(name+'.pdf'))
    plt.close(fig)


def threshold(a, b):
    return (-(a+1)+np.sqrt((a+1)**2-4*b))/(2*a)


def sector():
    fig = plt.figure(figsize=(16, 7.5))
    fig.subplots_adjust(left=.065, right=.96, bottom=.23, top=.78, wspace=.32)
    ax = fig.add_subplot(1, 3, 1)
    m = threshold(8.4, 12)
    ax.set_facecolor('#f4f6f8')
    ax.add_patch(Rectangle((m,m), .15-m, .15-m, facecolor=GREEN, alpha=.24))
    ax.plot([m,m],[m,.15], color=GREEN, lw=2.5)
    ax.plot([m,.15],[m,m], color=GREEN, lw=1.5, ls='--')
    ax.axhline(-.05, color=GRAY, lw=1, ls=':')
    ax.scatter([-.95],[-.05],color=RED,s=50,zorder=10)
    ax.annotate('Hidden cycle + stable origin\n+ saddle cycle',(-.95,-.05),
                xytext=(-.98,.035), fontsize=9.5, arrowprops=dict(arrowstyle='-', color=RED))
    ax.text(-.63,-.25,'Not classified by this theorem',ha='center',fontsize=9.5,color=GRAY)
    ax.text((m+.15)/2,.085,'All states\nconverge to 0',ha='center',fontsize=10,color=GREEN,weight='semibold')
    ax.add_patch(Rectangle((-.055,-.055),.01,.01,facecolor=GREEN,edgecolor='#163b2c'))
    ax.annotate('Earlier small sector',(-.05,-.05),xytext=(-.61,-.14),fontsize=9,
                arrowprops=dict(arrowstyle='-',color=GRAY))
    ax.set(xlim=(-1.03,.15),ylim=(-.3,.15),xlabel=r'$m_0$',ylabel=r'$m_1$',
           title='A  Slope plane: α = 8.4, β = 12')
    for i,kind in enumerate(['alpha','beta'],2):
        ax = fig.add_subplot(1,3,i,projection='3d')
        t = np.linspace(6,11,75) if kind=='alpha' else np.linspace(8,17,75)
        ms = threshold(t,12) if kind=='alpha' else threshold(8.4,t)
        tt, ss = np.meshgrid(t,np.linspace(0,1,30),indexing='ij')
        mm = np.broadcast_to(ms[:,None],tt.shape)
        free = mm+ss*(.08-mm)
        ax.plot_surface(mm,free,tt,color=GREEN,alpha=.56,linewidth=0,shade=False)
        ax.plot_surface(free,mm,tt,color=GREEN,alpha=.22,linewidth=0,shade=False)
        ax.plot(ms,ms,t,color=GREEN,lw=2)
        ax.plot(ms,np.full_like(ms,.08),t,color=GREEN,lw=1.5)
        ax.set(xlabel=r'$m_0$',ylabel=r'$m_1$',zlabel='α' if kind=='alpha' else 'β',
               title='B  β = 12' if kind=='alpha' else 'C  α = 8.4')
        ax.view_init(elev=23,azim=-57)
        ax.set_proj_type('ortho')
        ax.set_box_aspect((1,1,1.1))
        ax.tick_params(labelsize=9,pad=0)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.set_major_locator(MaxNLocator(4))
            axis.set_pane_color((.97,.98,.985,1))
        ax.text2D(.5,-.15,'Convergence region lies toward\nlarger m₀ and m₁ from these faces.',
                  transform=ax.transAxes,ha='center',fontsize=9.5,color=GREEN)
    fig.suptitle('A region where no oscillatory layer can exist',y=.96,fontsize=23,weight='semibold')
    fig.text(.5,.865,r'$m_0>m_*(\alpha,\beta),\quad m_1\geq m_*(\alpha,\beta),\qquad '
             r'm_*=[-(\alpha+1)+\sqrt{(\alpha+1)^2-4\beta}]/(2\alpha)$',ha='center',fontsize=14)
    fig.text(.5,.065,'Analytical storage proof, γ = 0. Assumptions: α > 0, β > 0, positive discriminant and r > 1.\n'
             'The m₀ face contains neutral periodic orbits. The m₁ face is a sufficient-sector edge, not an exact global boundary.\n'
             'The figures show a bounded part of an unbounded slope sector; no absence claim is made outside it.',
             ha='center',fontsize=10.5,color='#536573')
    save(fig,'chua-global-sector')


def separator():
    seeds = read('cycle-seeds.json'); basin = read('basin-separator.json')
    fig, axes = plt.subplots(1,2,figsize=(14,7.8))
    fig.subplots_adjust(left=.075,right=.96,bottom=.22,top=.79,wspace=.36)
    ax = axes[0]
    for name,color,label in [('attractor',BLUE,'Hidden attracting cycle'),('saddle',RED,'Saddle cycle')]:
        r=seeds[name]; states=sample_orbit(r,np.linspace(0,r['period'],2400))
        ax.plot(states[:,0],states[:,2],color=color,lw=1.8,label=label)
    ax.scatter([0],[0],color=GREEN,s=34,zorder=10,label='Stable equilibrium')
    ax.set(xlabel='x',ylabel='z',title='A  Coexisting invariant objects')
    ax.legend(frameon=False,fontsize=10,loc='upper right')
    inset=ax.inset_axes([.06,.09,.37,.30])
    r=seeds['saddle'];states=sample_orbit(r,np.linspace(0,r['period'],1000))
    inset.plot(states[:,0],states[:,2],color=RED,lw=1.6)
    inset.scatter([0],[0],color=GREEN,s=16)
    inset.set_title('Saddle: central detail',fontsize=9)
    inset.tick_params(labelsize=8)
    ax=axes[1]
    ss=np.array([c['stable_coordinate'] for c in basin['curve']])
    aa=np.array([c['approximations'][-1]['unstable_coordinate'] for c in basin['curve']])
    ax.plot(ss,aa,color=RED,lw=2,label='Computed local stable manifold')
    W=np.linalg.inv(basin['basis']);centre=np.array(basin['cycle']['q'][:2])
    for label,color,mark,title in [('origin',GREEN,'o','8 trajectories → origin'),
                                  ('outer_cycle',BLUE,'^','8 trajectories → hidden cycle')]:
        points=np.array([W@(np.array(c['initial'][1:])-centre) for c in basin['controls'] if c['label']==label])
        ax.scatter(points[:,1],points[:,0],color=color,marker=mark,s=38,label=title,zorder=10)
    ax.scatter([0],[0],color=RED,marker='x',s=60,zorder=11)
    ax.set(xlabel='Stable eigenvector coordinate s',ylabel='Transverse eigenvector coordinate a',
           title='B  Local section x = 0, dx/dt > 0')
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax.legend(frameon=False,fontsize=9.5,loc='lower center')
    fig.suptitle('A saddle between two attracting motions',y=.96,fontsize=24,weight='semibold')
    fig.text(.5,.875,'α = 8.4    β = 12    m₀ = −0.95    m₁ = −0.05    γ = 0',ha='center',fontsize=13)
    fig.text(.5,.095,'Cycle existence and stability: interval certificates at 80 and 100 digits.\n'
        'Saddle multipliers: 3.35201949 and 0.0000231261. Two attractors are exhibited; the saddle is not an attractor.\n'
        'The separator curve and 16 two-sided outcomes are numerical. The complete basin boundary is not certified.',
        ha='center',fontsize=10.5,color='#536573')
    save(fig,'chua-basin-saddle')


def contact():
    rows=read('linear-contact-results.json')['records']; branches=read('branches.json')['records']
    predictions=read('linear-contact-prediction.json')['forecasts']
    fig,axes=plt.subplots(2,2,figsize=(14,10))
    fig.subplots_adjust(left=.075,right=.96,bottom=.15,top=.85,wspace=.27,hspace=.43)
    for branch,color,label in [('saddle',RED,'Saddle'),('attractor',BLUE,'Attracting cycle')]:
        curve=[r['cycle'] for r in branches if r['branch']==branch]
        data=[r for r in rows if r['branch']==branch]
        ax=axes[0,0]
        ax.plot([r['m0'] for r in curve],[r['amplitude'] for r in curve],color=color,lw=2,label=label)
        ax.plot([r['cycle']['m0'] for r in data],[r['cycle']['amplitude'] for r in data],color=color,lw=2)
        delta=np.array([r['delta'] for r in data]); eta=np.array([r['cycle']['amplitude']-1 for r in data])
        axes[0,1].loglog(delta,eta,'o',ms=3.3,color=color,label=label+' (computed)')
        axes[0,1].loglog(delta,[r['predicted_excess'] for r in data],color=color,lw=1.5)
        axes[1,0].loglog(delta,[r['cycle']['q'][3] for r in data],'o',ms=3.3,color=color)
        axes[1,0].loglog(delta,[r['predicted_outer_duration'] for r in data],color=color,lw=1.5)
        axes[1,1].semilogx(delta,[r['amplitude_ratio'] for r in data],color=color,lw=2,label=label+': amplitude')
        axes[1,1].semilogx(delta,[r['duration_ratio'] for r in data],color=color,lw=1.5,ls='--',label=label+': duration')
    ax=axes[0,0]
    ax.axhline(1,color=GRAY,ls=':',lw=1)
    for p in predictions:ax.scatter([p['mstar']],[1],facecolor='white',edgecolor=GRAY,s=45,zorder=10)
    lo,hi=[p['mstar'] for p in predictions]
    ax.plot([-.999,lo],[0,0],color=GREEN,lw=3)
    ax.plot([lo,hi],[0,0],color=GREEN,lw=1.3,ls='--')
    ax.plot([hi,-.12],[0,0],color=GREEN,lw=3)
    ax.set(xlabel=r'$m_0$',ylabel=r'$A=\max |x|$',title='A  Two distinct centre contacts',xlim=(-1.015,-.12))
    ax.legend(frameon=False,fontsize=9.5)
    ax.text(-.55,.55,'A → 1, not 0',fontsize=10,color=GRAY)
    axes[0,1].set(xlabel=r'$\delta=m_* - m_0$',ylabel=r'$A-1$',title=r'B  Amplitude excess: $c_A\delta^{2/3}$')
    axes[0,1].legend(frameon=False,fontsize=9.5)
    axes[1,0].set(xlabel=r'$\delta=m_* - m_0$',ylabel='Time in x > 1',title=r'C  Outer flight: $c_T\delta^{1/3}$')
    axes[1,1].axhline(1,color=GRAY,lw=1)
    axes[1,1].set(xlabel=r'$\delta=m_* - m_0$',ylabel='Computed / predicted',title='D  Test of the coefficients')
    axes[1,1].legend(frameon=False,fontsize=9)
    fig.suptitle('How the cycles meet the linear central region',y=.965,fontsize=23,weight='semibold')
    fig.text(.5,.905,'α = 8.4, β = 12, m₁ = −0.05     ·     31 distances per branch, from 10⁻² to 10⁻⁷',ha='center',fontsize=12)
    fig.text(.5,.045,'Markers: piecewise-linear closed orbits. Lines in B/C: coefficients fixed before the fine-distance sweep.\n'
        'A: the green origin branch is solid where locally stable and dashed where unstable.\n'
        'At δ = 10⁻⁷, all four coefficient errors are below 0.0061%. This is an asymptotic numerical check, not a remainder proof.',
        ha='center',fontsize=10.5,color='#536573')
    save(fig,'chua-linear-contact')


if __name__=='__main__':
    sector(); separator(); contact()
