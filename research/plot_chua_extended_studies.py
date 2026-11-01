"""Figures from the excitation, coexistence-box and held-out prediction data."""
import os
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path(__file__).resolve().parents[1]/'.publication/matplotlib'))
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from scipy.integrate import solve_ivp
from chua_research_followup import ROOT,OUT,dump
from study_chua_excitation_plane import field

ASSETS=ROOT/'docs/assets'
COLORS=['#34699a','#de7c35','#368b72','#9c5c96']


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':13,
        'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,
        'figure.facecolor':'#fcfbf8','axes.facecolor':'#fcfbf8','savefig.facecolor':'#fcfbf8',
        'grid.alpha':.2,'axes.labelcolor':'#243647','text.color':'#243647'})


def save(fig,name):
    fig.savefig(ASSETS/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(ASSETS/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)


def excitation():
    directory=OUT/'excitation-plane';d=np.load(directory/'plane.npz')
    meta=json.loads((directory/'plane.json').read_text());pulses=json.loads((directory/'pulses.json').read_text())
    opt=json.loads((directory/'pulse-optimum.json').read_text())
    fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    ax=axs[0,0];aa,rr=d['angles'],d['radii']
    ae=np.r_[aa[0]-(aa[1]-aa[0])/2,(aa[:-1]+aa[1:])/2,aa[-1]+(aa[1]-aa[0])/2]
    re=np.r_[0,(rr[:-1]+rr[1:])/2,rr[-1]+(rr[1]-rr[0])/2]
    a,r=np.meshgrid(ae,re,indexing='ij');y=r*np.cos(a);z=r*np.sin(a)
    cmap=ListedColormap(['#bcb9b1','#eadfcf','#86b9b0']);norm=BoundaryNorm([-.5,.5,1.5,2.5],3)
    for sign in [-1,1]:
        ax.pcolormesh(sign*y,sign*z,d['labels'],cmap=cmap,norm=norm,rasterized=True)
        ax.plot(sign*d['thresholds']*np.cos(aa),sign*d['thresholds']*np.sin(aa),color='#243647',lw=1.4)
    theta=np.linspace(0,2*np.pi,500);rmin=meta['minima'][0]['radius']
    ax.plot(rmin*np.cos(theta),rmin*np.sin(theta),'--',color='#a94040',lw=1.5,label='Equal-energy circle at lowest sampled threshold')
    ax.scatter([.34061011644,0],[0,-rmin],color=['#34699a','#a94040'],s=45,zorder=5)
    ax.text(.03,.07,'Origin',fontsize=12);ax.text(.24,.32,'Hidden cycle',fontsize=12)
    ax.set(xlim=(-.43,.43),ylim=(-.43,.43),xlabel='$y$ (capacitor voltage)',ylabel=r'$z/\sqrt{12}$ (weighted current)',title='A  Initial states: same energy is a circle',aspect='equal')
    ax=axs[0,1]
    ax.plot(np.degrees(aa),d['thresholds']**2/2,color='#243647',lw=2)
    ax.axhline(.34061011644**2/2,color='#34699a',ls='--',label='Voltage-only excitation infimum')
    ax.scatter([-90],[rmin*rmin/2],s=65,color='#a94040',zorder=5)
    ax.annotate('Near current-only start\n~27% of voltage-only energy',(-90,rmin*rmin/2),(-72,.034),arrowprops={'arrowstyle':'->','color':'#a94040'})
    ax.set(xlabel='Direction in energy coordinates (degrees)',ylabel='$E/(C_2 B_p^2)$',title='B  Direction changes the excitation cost')
    ax.legend(loc='upper right',fontsize=9);ax.grid()
    pulse=np.load(directory/'pulses.npz');ax=axs[1,0]
    ax.pcolormesh(pulse['durations'],pulse['impulses'],pulse['labels'].T,cmap=cmap,norm=norm,shading='nearest')
    ax.axvline(pulses['rows'][23]['duration'],color='#a94040',lw=1.2,ls='--')
    ax.set(xscale='log',xlabel=r'Pulse duration $\tau$',ylabel=r'Pulse area $a\tau$',title='C  Rectangular forcing in the voltage equation')
    ax.text(.025,2.6,'Cycle',fontsize=12);ax.text(3,.5,'Origin',fontsize=12)
    ax=axs[1,1]
    entries=[(row['duration'],t['squared_input_cost']) for row in pulses['rows'] for t in row['transitions'] if t['left']==1]
    ax.plot(*np.array(entries).T,color='#243647',marker='o',ms=4)
    ax.scatter(opt['duration'],opt['squared_input_cost'],s=90,marker='*',color='#a94040',zorder=5)
    ax.annotate('Local minimum\n'+fr"$\tau\approx${opt['duration']:.3f}, $a\approx${opt['amplitude']:.3f}",
        (opt['duration'],opt['squared_input_cost']),(.07,.32),arrowprops={'arrowstyle':'->','color':'#a94040'})
    ax.set(xscale='log',yscale='log',xlabel=r'Pulse duration $\tau$',ylabel=r'Input cost $\int u^2 dt=a^2\tau$',title='D  A finite-duration optimum within this pulse class')
    ax.grid()
    fig.suptitle('Excitation of a hidden periodic orbit',fontsize=19,fontweight='bold')
    save(fig,'chua-extended-excitation')


def pulse_demo():
    directory=OUT/'excitation-plane';tau=2.345526436848252
    ts=np.linspace(0,600,18001);traces=[];records=[]
    from study_chua_trigger import setup
    pc=setup()[1]
    for amplitude in [.7,1.2]:
        def rhs(t,X,u=0):return field.py_func(*X,u)
        on=solve_ivp(lambda t,X:rhs(t,X,amplitude),[0,tau],[0,0,0],method='DOP853',rtol=2e-12,atol=2e-13,dense_output=True,max_step=.02)
        def event(t,X):return X[0]
        event.direction=1
        off=solve_ivp(rhs,[tau,600],on.y[:,-1],method='DOP853',rtol=2e-12,atol=2e-13,dense_output=True,max_step=.03,events=event)
        assert on.success and off.success
        a=ts<=tau;X=np.empty((3,len(ts)));X[:,a]=on.sol(ts[a]);X[:,~a]=off.sol(ts[~a]);traces.append(X)
        distance=float(np.linalg.norm(off.y_events[0][-1,1:]-pc))
        finalnorm=float(np.linalg.norm(off.y[:,-1]))
        expected='cycle' if amplitude==.7 else 'origin'
        assert (distance<1e-6 if expected=='cycle' else finalnorm<1e-8)
        records.append(dict(amplitude=amplitude,duration=tau,destination=expected,
                            last_section_distance=distance,final_norm=finalnorm))
    np.savez_compressed(directory/'pulse-demonstration.npz',time=ts,states=np.array(traces),amplitudes=[.7,1.2])
    dump(directory/'pulse-demonstration.json',dict(records=records,method='DOP853',passed=True))
    fig,axs=plt.subplots(2,2,figsize=(13,6.8),gridspec_kw={'width_ratios':[2,1]},layout='constrained')
    for i,(X,rec) in enumerate(zip(traces,records)):
        color=['#368b72','#a94040'][i]
        axs[i,0].plot(ts,X[0],color=color,lw=.85)
        axs[i,0].axvspan(0,tau,color='#de7c35',alpha=.3)
        axs[i,0].set(xlim=(0,110),ylabel='$x(t)$',xlabel='Time after pulse onset',title=f"Amplitude {rec['amplitude']}: {'sustained oscillation' if i==0 else 'return to rest'}")
        axs[i,1].plot(ts,X[0],color=color,lw=1)
        axs[i,1].set(xlim=(570,580),ylim=(-9,9),xlabel='Late time',ylabel='$x(t)$')
    fig.suptitle(f'A stronger pulse can fail to excite the cycle — same duration {tau:.3f}',fontsize=17,fontweight='bold')
    save(fig,'chua-extended-pulse-return')


def coexistence():
    directory=OUT/'coexistence-box';d=np.load(directory/'slope-section.npz')
    cert=json.loads((directory/'certificate.json').read_text())
    fig,axs=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    ax=axs[0];off=d['offsets']*1000
    counts=d['distinct_attracting_counts']
    cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#bcb9b1')
    masked=np.ma.masked_where(np.any(d['status']==0,axis=2),counts)
    im=ax.pcolormesh(off,off,masked.T,vmin=0,vmax=4,cmap=cmap,shading='nearest')
    ax.scatter([0],[0],s=70,facecolors='none',edgecolors='#fff',lw=1.6)
    ax.set(xlabel='$10^3(m_0-m_{0,c})$',ylabel='$10^3(m_1-m_{1,c})$',title='A  Numerically continued attracting cycles')
    fig.colorbar(im,ax=ax,ticks=range(5),label='Attracting cycles; grey = unresolved')
    ax=axs[1];ids=[r['id'] for r in cert['records']]
    con=[float(__import__('flint').arb(r['contraction']['ball']).upper()) for r in cert['records']]
    exc=[max(float(__import__('flint').arb(v['inverse_defect']['ball']).upper()) for v in r['multiplier_exclusions']) for r in cert['records']]
    x=np.arange(4);ax.bar(x-.18,con,width=.35,color='#34699a',label='Root contraction bound')
    ax.bar(x+.18,exc,width=.35,color='#de7c35',label='Variational inverse defect')
    ax.axhline(1,color='#a94040',ls='--');ax.set(xticks=x,xticklabels=ids,ylim=(0,1.12),ylabel='Required to be below 1',title='B  Uniform interval tests for all four cycles');ax.legend(fontsize=8,loc='upper left')
    ax=axs[2];ax.axis('off')
    text='Certified simultaneously\n\n'+r'$|\alpha-9.06|\leq10^{-8}$'+'\n'+r'$|\beta-14.3458|\leq10^{-8}$'+'\n'+r'$|m_0-(-8/7+0.006)|\leq10^{-8}$'+'\n'+r'$|m_1+5/7|\leq10^{-8}$'+'\n'+r'$0\leq\gamma\leq10^{-8}$'+'\n\nThe certified box is much smaller\nthan one cell of the numerical map.\nIt is a conservative proof region,\nnot the measured stability boundary.'
    ax.text(.03,.96,text,va='top',fontsize=12,linespacing=1.65)
    fig.suptitle('Four coexisting periodic orbits persist on a parameter box',fontsize=18,fontweight='bold')
    save(fig,'chua-extended-coexistence')


def predictions():
    directory=OUT/'extended-slopes';scores=json.loads((directory/'scores.json').read_text())
    fig,axs=plt.subplots(1,2,figsize=(14,7),gridspec_kw={'width_ratios':[1.8,1]},layout='constrained')
    ax=axs[0]
    for i,rec in enumerate(scores['records']):
        node=rec['id'].rsplit('_',1)[0];peaks=int(rec['id'].rsplit('_',1)[1])
        data=json.loads((directory/(node+'-search.json')).read_text());row=next(r for r in data['records'] if r['peaks']==peaks)
        import mpmath as mp
        mp.mp.dps=50;a=mp.mpf(row['fold']['coordinate']);b=mp.mpf(row['symmetry']['coordinate']);w=a-b
        raw=float((mp.mpf(row['forecast']['midpoint'])-a)/w)
        corrected=float((mp.mpf(row['forecast']['corrected_midpoint'])-a)/w)
        ax.plot([-1,0],[i,i],color='#b6d4ca',lw=7,solid_capstyle='butt',zorder=1)
        ax.scatter(raw,i,color='#de7c35',marker='x',s=30,zorder=3,label='Uncorrected forecast' if i==0 else None)
        ax.scatter(corrected,i,color='#243647',s=18,zorder=3,label='Previously calibrated forecast' if i==0 else None)
    ax.set(yticks=np.arange(24),yticklabels=[r['id'].replace('_',' · ') for r in scores['records']],
        xlabel='Position relative to measured fold, divided by measured interval width',title='A  Predictions published before target searches')
    ax.tick_params(axis='y',labelsize=8);ax.invert_yaxis();ax.axvline(-1,color='#368b72',lw=.8);ax.axvline(0,color='#368b72',lw=.8);ax.legend(fontsize=9,loc='lower left')
    ax=axs[1]
    for free,color in [('beta','#34699a'),('alpha','#9c5c96')]:
        recs=[r for r in scores['records'] if r['id'].startswith(free)]
        ax.scatter([abs(r['corrected_fold_error'])*100 for r in recs],[r['width_relative_error']*100 for r in recs],label=f'Predict {free}',color=color,s=45,alpha=.8)
    ax.axhline(10,color='#a94040',ls='--',label='Predeclared width-error limit')
    ax.set(xlabel='Fold error (% of measured width)',ylabel='Relative width error (%)',title='B  Accuracy on the enlarged slope square',ylim=(0,11));ax.grid();ax.legend(fontsize=9)
    fig.suptitle('24 / 24 new corrected midpoints fall inside the measured windows',fontsize=18,fontweight='bold')
    save(fig,'chua-extended-predictions')


def short_pulse():
    data=json.loads((OUT/'excitation-plane/short-pulse-law.json').read_text())
    rows=data['comparisons'];tau=np.array([r['duration'] for r in rows])
    actual=np.array([r['measured_area']-data['voltage_ray_threshold'] for r in rows])
    c=data['normal_approximations'][-1]['quadratic_coefficient']
    fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    ax=axs[0];ax.plot(tau*tau,c*tau*tau,color='#34699a',lw=2,label=r'Prediction $c\tau^2$; no pulse-data fit')
    ax.scatter(tau*tau,actual,color='#de7c35',s=40,zorder=3,label='Integrated pulse thresholds')
    ax.set(xlabel=r'Squared pulse duration $\tau^2$',ylabel=r'Excess area $I_*(\tau)-A_*$',title='A  The linear duration term cancels');ax.legend(fontsize=9);ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2));ax.grid()
    ax=axs[1];ax.plot(tau,[100*r['quadratic_increment_relative_error'] for r in rows],color='#243647',marker='o')
    ax.set(xlabel=r'Pulse duration $\tau$',ylabel='Error relative to the area correction (%)',title='B  Error is small but visible');ax.grid()
    fig.suptitle(r'Short-pulse law: $a\tau=A_*+0.0535732\,\tau^2+O(\tau^4)$',fontsize=18,fontweight='bold')
    save(fig,'chua-extended-short-pulse')


if __name__=='__main__':
    style();excitation();coexistence();predictions();short_pulse();pulse_demo()
