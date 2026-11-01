"""Figures from the accessibility certificates and their numerical controls."""
import json
import os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'.publication/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from certify_chua_coexisting_cycles import ROOT, OUT

ASSETS=ROOT/'docs/assets'
COLORS=['#256e9e','#cc6246','#448f75','#9963a6']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.titleweight':'bold','axes.labelcolor':'#293642',
                     'text.color':'#293642','axes.edgecolor':'#8e9aa5','grid.alpha':.18,
                     'figure.facecolor':'#ffffff','savefig.facecolor':'#ffffff'})


def load(name):return json.loads((OUT/name).read_text())
def val(x):return float(x['midpoint'])


def trace(cycle,pars,count=2200):
    alpha,beta,m0,m1=pars
    def fun(t,u):
        x,y,z=u; h=m1*x+(m0-m1)*np.clip(x,-1,1)
        return [alpha*(y-x-h),x-y+z,-beta*y]
    t=np.linspace(0,val(cycle['period']),count)
    u=[0,*[val(v) for v in cycle['root']]]
    sol=solve_ivp(fun,[0,t[-1]],u,t_eval=t,rtol=2e-11,atol=2e-12,max_step=.008)
    assert sol.success
    return t,sol.y,fun


def save(fig,name):
    fig.savefig(ASSETS/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)


def main():
    four=load('four-cycle-certificate-80.json')
    hidden=load('hidden-cycle-certificate-80.json')
    stable=load('stability-certificate-80.json')
    h2=stable['unique_equilibrium_hidden_cycle']
    basin=load('stability-basin-check.json')
    adaptive=load('adaptive-prediction-check.json')
    assert basin.get('section_interpolation')=='cubic Hermite','Wait for the corrected basin run'
    fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    t,u,fun=trace(h2['periodic_orbit'],[8.4,12,-.95,-.05])
    ax=axs[0,0]; ax.plot(u[0],u[2],color=COLORS[0],lw=1.8,label='Attracting periodic orbit')
    ax.scatter([0],[0],s=35,color=COLORS[1],zorder=5,label='Stable equilibrium')
    ax.set(xlabel='x',ylabel='z',title='A  A stable equilibrium coexists with a hidden cycle')
    ax.text(.03,.94,'m₀ = −0.95,  m₁ = −0.05\nα = 8.4,  β = 12,  γ = 0',transform=ax.transAxes,va='top',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    ax.legend(loc='lower left',fontsize=9,frameon=False); ax.grid()
    P=np.array(h2['origin_lyapunov']['P'],float); level=val(h2['ellipsoid_level'])
    ax=axs[0,1]
    x,z=np.meshgrid(np.linspace(-.8,.8,401),np.linspace(-.5,.5,401))
    V=P[0,0]*x*x+2*P[0,2]*x*z+P[2,2]*z*z
    ax.contourf(x,z,V,levels=[0,level],colors=['#d9e9f1'])
    ax.contour(x,z,V,levels=[level],colors=[COLORS[0]],linewidths=1)
    radius=val(h2['excluded_origin_radius'])
    ax.add_patch(Circle((0,0),radius,fill=False,color=COLORS[1],lw=2))
    for initial in [[.02,0,.01],[-.02,0,-.01]]:
        sol=solve_ivp(fun,[0,80],initial,rtol=1e-10,atol=1e-12,max_step=.08)
        ax.plot(sol.y[0],sol.y[2],color=COLORS[1],lw=.8)
    ax.scatter([0],[0],color=COLORS[1],s=25)
    ax.set(xlabel='x',ylabel='z',title='B  Certified region returning to the equilibrium',xlim=(-.22,.22),ylim=(-.18,.18))
    ax.text(.03,.97,'Shading: ellipsoid section y = 0; traces: x–z projections\nEuclidean ball radius ≥ 0.0305',transform=ax.transAxes,va='top',fontsize=8,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    ax.grid()
    ax=axs[1,0]
    grid=np.linspace(-.99,-.89,401)
    spectral=[]
    for m in grid:
        A=np.array([[-8.4*(1+m),8.4,0],[1,-1,1],[0,-12,0]])
        spectral.append(max(np.linalg.eigvals(A).real))
    threshold=val(stable['routh_thresholds'][0])
    ax.axvspan(-.99,threshold,color='#e4eff4'); ax.axhline(0,color='#89929a',lw=1)
    ax.plot(grid,spectral,color='#293642',lw=1.6,label='Largest equilibrium eigenvalue: real part')
    ax.axvline(threshold,ls='--',color=COLORS[1],lw=1.3)
    ax.text(threshold+.002,.02,'−0.9376757184',rotation=90,color=COLORS[1],fontsize=9)
    for row in basin['records']:
        m=float(row['prediction']['m0']); count=row['controls'][-1]['counts']
        if not row['prediction']['origin_stable']:
            passed=sum(r['reached_cycle'] for r in adaptive['records'] if r['m0']==m)
            count=[24-passed,passed,0]
        color=COLORS[0] if count[2]==24 else COLORS[2] if count[1]==24 else '#777777'
        text='24 → equilibrium' if count[2]==24 else '24 → cycle' if count[1]==24 else f'{count[1]} → cycle; {count[0]} unresolved'
        ax.text(m,.75,text,transform=ax.get_xaxis_transform(),rotation=35,ha='left',fontsize=9,color=color)
    ax.set(xlabel='Inner slope m₀',ylabel='Real part of equilibrium eigenvalue',title='C  Local prediction and new initial-condition checks',ylim=(-.18,.38))
    ax.text(.03,.03,'24 starts per parameter; radii 10⁻³, 10⁻⁶, 10⁻⁹\nAdaptive follow-up on the unstable side.\nThreshold is local stability, not a certified basin boundary.',transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    ax.grid()
    ax=axs[1,1]
    _,u,fun1=trace(hidden['periodic_orbit'],[8.4,12,-1.2,-.05])
    ax.plot(u[0],u[2],color=COLORS[0],lw=1.7)
    lam=val(hidden['origin_spectrum']['unstable']); v=np.array([1,-.2+lam/8.4,-12*(-.2+lam/8.4)/lam])
    for sign in [-1,1]:
        sol=solve_ivp(fun1,[0,120],sign*.5*v,rtol=2e-11,atol=2e-12,max_step=.05)
        ax.plot(sol.y[0],sol.y[2],color=COLORS[1],lw=1.5)
    e=23/19; ax.scatter([e,-e],[-e,e],color=COLORS[1],s=32,zorder=5)
    ax.scatter([0],[0],facecolor='white',edgecolor='#293642',s=32,zorder=5)
    ax.set(xlabel='x',ylabel='z',title='D  Hiddenness also certified with three equilibria')
    ax.text(.03,.97,'m₀ = −1.20\nCentral exits enter the two outer basins.\nEquilibrium neighborhoods are excluded by bounds.',transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    ax.grid()
    fig.suptitle('Attractor existence and accessibility are different questions',fontsize=17,fontweight='bold')
    save(fig,'chua-accessibility')

    fig,axs=plt.subplots(2,4,figsize=(15,7),layout='constrained',sharex='row',sharey='row')
    for j,c in enumerate(four['records']):
        t,u,_=trace(c,[9.06,14.3458,-8/7+.006,-5/7])
        rho=max(abs(val(v)) for v in c['multipliers'])
        axs[0,j].plot(u[0],u[2],lw=1.2,color=COLORS[j])
        axs[0,j].set(title=f"{c['id']}   T = {val(c['period']):.6f}",xlabel='x')
        axs[0,j].text(.03,.04,f'|μ|max ≈ {rho:.6f}',transform=axs[0,j].transAxes,fontsize=9)
        axs[1,j].plot(t,u[0],lw=1.2,color=COLORS[j]); axs[1,j].set(xlabel='Time',xlim=(0,15))
        for i in [0,1]:axs[i,j].grid()
    axs[0,0].set_ylabel('z'); axs[1,0].set_ylabel('x(t)')
    fig.suptitle('Four distinct attracting cycles at one exact parameter point',fontsize=17,fontweight='bold')
    fig.supxlabel('α = 453/50,  β = 71729/5000,  m₀ = −8/7 + 3/500,  m₁ = −5/7,  γ = 0\nFirst-return roots, complete itineraries and attraction certified with interval arithmetic.',fontsize=10)
    save(fig,'chua-four-cycle-certificate')

    circuit=load('circuit-constraints.json'); rows=circuit['samples']; c=circuit['nominal_components_SI']
    fig,axs=plt.subplots(1,2,figsize=(13,5.6),layout='constrained')
    ax=axs[0]; x=np.linspace(-24,24,1201); m2=c['R']*(1/c['R1']+1/c['R2'])
    h=-.05*x+(-1.2+.05)*np.clip(x,-1,1)+(m2+.05)*np.sign(x)*np.maximum(abs(x)-20,0)
    ax.plot(x,h,color=COLORS[0],lw=2)
    for v in [-20,-1,1,20]:ax.axvline(v,color='#a4adb5',ls=':',lw=1)
    nominal=rows[0]['xmax_sampled']; ax.axvspan(-nominal,nominal,color='#e4eff4',zorder=-1)
    ax.set(xlabel='Normalized voltage x',ylabel='Diode characteristic h(x)',title='A  Full five-segment diode and the orbit range')
    ax.text(.03,.04,'Shading: nominal orbit range\nOuter breakpoints: ±20\nPositive coil resistance: 0.5 Ω',transform=ax.transAxes,fontsize=9); ax.grid()
    ax=axs[1]; indices=np.arange(len(rows))
    ax.plot(indices,[r['outer_breakpoint'] for r in rows],color='#293642',lw=1,label='Second breakpoint')
    ax.scatter(indices,[r['xmax_sampled'] for r in rows],s=12,color=COLORS[2],label='Sampled orbit amplitude')
    ax.set(xlabel='Component sample (0 = nominal)',ylabel='|x|',title='B  Nominal circuit and 128 component perturbations',ylim=(0,23))
    ax.text(.03,.12,'All 129 numerical cycles remain attracting.\nEach component varies independently within ±1%.\nThis is a finite sample, not a tolerance-box proof.',transform=ax.transAxes,fontsize=9)
    ax.legend(frameon=False,loc='lower right',bbox_to_anchor=(1,.3),fontsize=9); ax.grid()
    fig.suptitle('Component constraints and numerical robustness',fontsize=17,fontweight='bold')
    save(fig,'chua-circuit-accessibility')


if __name__=='__main__':
    main()
