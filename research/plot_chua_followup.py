"""Scientific figures drawn from the frozen tests and interval evidence."""
import json
from decimal import Decimal
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.integrate import solve_ivp
from flint import arb,arb_mat
from chua_research_followup import ROOT,OUT
from chua_lobe_map import LobeMap
from chua_interval_lobe import IntervalLobe,saddle_center
from certify_chua_layer_organizer import enclose_saddle
from chua_interval_connection import connect

INK='#213647';BLUE='#246ca0';RED='#b85247';TEAL='#237e70';GOLD='#b47b2a'
COLORS={65:BLUE,75:TEAL,85:RED}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelcolor':INK,'text.color':INK,'axes.titlecolor':INK,'axes.spines.top':False,'axes.spines.right':False})
def read(p):return json.loads((OUT/p).read_text())
def save(fig,name):
    fig.savefig(ROOT/'docs/assets'/f'{name}.png',dpi=180,bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
    plt.close(fig)


def sheets():
    d=read('analysis.json')['sheets'];rows=[r for r in d['records'] if r['held_out']]
    fig=plt.figure(figsize=(15,11));gs=fig.add_gridspec(2,2,hspace=.38,wspace=.24,left=.07,right=.96,bottom=.10,top=.85)
    fig.suptitle('Predicting thin periodic layers across the two slope sections',x=.06,y=.98,ha='left',fontsize=21,weight='bold')
    fig.text(.06,.925,'48 held-out windows  |  65, 75 and 85 peaks  |  Forecasts published before endpoint searches',fontsize=12)
    ax=fig.add_subplot(gs[0,0]);ax.axvspan(0,1,color=TEAL,alpha=.13,label='Measured symmetric core')
    for k,r in enumerate(rows):
        y=k+(2 if k>=24 else 0)
        ax.plot([r['raw_midpoint_in_core_units'],r['transport_midpoint_in_core_units']],[y,y],color='#b7c2c8',lw=.7)
        ax.scatter(r['raw_midpoint_in_core_units'],y,c=RED,s=16,marker='x',label='Uncorrected forecast' if k==0 else None)
        ax.scatter(r['transport_midpoint_in_core_units'],y,c=BLUE,s=12,label='Transported central residual' if k==0 else None)
    ax.set(xlabel='Predicted midpoint / measured core width',yticks=[11.5,37.5],yticklabels=['Fixed β','Fixed α'])
    ax.set_title('A  Every forecast at the scale of its own window',loc='left',fontsize=12,weight='bold')
    ax.legend(frameon=False,fontsize=9,loc='upper right');ax.grid(axis='x',alpha=.15)
    ax=fig.add_subplot(gs[0,1])
    for section,marker,offset in [('beta','o',-.08),('alpha','s',.08)]:
        for n in [65,75,85]:
            group=[r for r in rows if r['section']==section and r['peaks']==n]
            ax.scatter(np.arange(len(group))*.018+n+offset,[100*r['width_relative_error'] for r in group],s=24,c=COLORS[n],marker=marker,
                label=('Fixed α; predict β' if section=='beta' else 'Fixed β; predict α') if n==65 else None)
    ax.axhline(10,color=INK,ls='--',lw=1,label='Predeclared 10% limit');ax.axhline(0,color=INK,lw=.7)
    ax.set(xlabel='Peaks in one full period',ylabel='Relative width error (%)',xticks=[65,75,85],ylim=(-1,11))
    ax.set_title('B  Width prediction retained without correction',loc='left',fontsize=12,weight='bold')
    ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.15)
    for j,free in enumerate(['beta','alpha']):
        ax=fig.add_subplot(gs[1,j],projection='3d')
        for n in [65,75,85]:
            x=[];y=[];z=[]
            for path in sorted((OUT/'sheets').glob(f'{free}_??-search.json')):
                data=json.loads(path.read_text());cal=read('sheets/'+data['node']+'-calibration.json');r=next(r for r in data['records'] if r['peaks']==n)
                x.append((data['parameters']['m0']+8/7)*1e4);y.append((data['parameters']['m1']+5/7)*1e4)
                z.append(float(Decimal(r['fold']['coordinate'])-Decimal(str(cal['organizer']['parameters'][free])))*1e6)
            ax.plot_trisurf(x,y,z,color=COLORS[n],alpha=.58,linewidth=.5,edgecolor='white')
            ax.scatter(x,y,z,color=COLORS[n],s=16)
            ax.text(x[-1]+.5,y[-1],z[-1],f'{n} peaks',fontsize=9,color=COLORS[n])
        ax.view_init(24,-48);ax.set_xlabel(r'$10^4\,\Delta m_0$');ax.set_ylabel(r'$10^4\,\Delta m_1$')
        symbol=r'\beta' if free=='beta' else r'\alpha'
        ax.set_zlabel(r'$10^6('+symbol+'-'+symbol+'_*)$',labelpad=9)
        ax.set_title(('C  Fixed α = 9; predict β' if free=='beta' else 'D  Fixed β ≈ 14.30636; predict α'),loc='left',fontsize=12,weight='bold',pad=4)
    fig.text(.06,.025,'Surfaces join nine solved folds per family. The vertical coordinate subtracts the numerical organizing surface to expose the layer spacing.',fontsize=10)
    save(fig,'chua-followup-sheets')


def field(p):
    a,b,m0,m1=p
    def rhs(t,v):
        x,y,z=v;h=m1*x+(m0-m1)*(abs(x+1)-abs(x-1))/2
        return [a*(y-x-h),x-y+z,-b*y]
    return rhs


def pulse():
    d=read('trigger/delay.json');thresholds=read('trigger/thresholds.json')['records'];fig,axes=plt.subplots(1,3,figsize=(16,6))
    fig.subplots_adjust(left=.06,right=.97,bottom=.22,top=.75,wspace=.35)
    fig.suptitle('A finite impulse reaches a hidden periodic attractor',x=.06,y=.97,ha='left',fontsize=22,weight='bold')
    fig.text(.06,.88,r'$\alpha=8.4,\ \beta=12,\ m_0=-0.95,\ m_1=-0.05$'+'  |  Seven certified local basin-boundary crossings',fontsize=12)
    ax=axes[0]
    for r in thresholds:
        A=r['approximants'][-1]['amplitude'];v=np.array(r['direction'],float)*A
        ax.plot([0,v[0]],[0,v[1]],color='#b7c7d0',lw=1);ax.scatter(v[0],v[1],c=BLUE,s=32)
        ax.annotate(f'{A:.4f}',(v[0],v[1]),xytext=(5,4),textcoords='offset points',fontsize=9)
    ax.scatter(0,0,color=INK,s=20);ax.set(xlabel='Initial y',ylabel='Initial z',xlim=(-.02,.56));ax.set_aspect('equal',adjustable='datalim');ax.grid(alpha=.15)
    ax.set_title('A  Threshold depends on impulse direction',loc='left',fontsize=12,weight='bold')
    ax=axes[1];A=thresholds[3]['approximants'][-1]['amplitude'];rhs=field([8.4,12,-.95,-.05]);t=np.linspace(0,90,9001)
    for side,color,label in [(-1,BLUE,'Below: origin'),(1,RED,'Above: hidden cycle')]:
        sol=solve_ivp(rhs,[0,90],[0,A+side*1e-6,0],method='DOP853',rtol=2e-11,atol=2e-13,t_eval=t,max_step=.03)
        ax.plot(t,sol.y[0],color=color,lw=.8,label=label)
    ax.set(xlabel='Time',ylabel='x');ax.legend(frameon=False,fontsize=9,loc='upper left')
    ax.set_title('B  Two impulses separated by 2 × 10⁻⁶',loc='left',fontsize=12,weight='bold');ax.grid(alpha=.15)
    ax=axes[2];slope=d['prediction']['predicted_time_per_decimal_digit']
    for angle,color in [(0.,BLUE),(.8,RED)]:
        rows=[r for r in d['records'] if r['angle']==angle];x=np.array([-np.log10(r['distance']) for r in rows]);y=np.array([r['residence'] for r in rows])
        ax.scatter(x,y,s=12,color=color,alpha=.65,label=f'θ = {angle:g}')
        xx=np.array([4,10]);ax.plot(xx,slope*xx+np.mean(y-slope*x),color=color,ls='--',lw=1.5)
    ax.set(xlabel=r'$-\log_{10}|A-A_*|$',ylabel='Time near the saddle cycle',xticks=[4,6,8,10]);ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=9)
    ax.set_title('C  Predicted delay: 5.85255 per decade',loc='left',fontsize=12,weight='bold')
    fig.text(.06,.07,'A: brackets are narrower than the markers (2 × 10⁻¹¹).  B: independent numerical illustration.  C: dashed slopes are predicted; only their offsets are fitted.',fontsize=10)
    save(fig,'chua-followup-trigger')


def hidden():
    center=read('hidden-search/connection-center-14-2.json');cert=read('hidden-search/hidden-tangency-160.json');f=IntervalLobe(center['parameters'],140)
    beta=arb(center['beta']);p=saddle_center(f,beta,center['parent_guess']);p,pb,pbb,_=enclose_saddle(f,p,beta,'1e-100',inflation_floor='1e-115')
    V=arb_mat([[arb(v) for v in row] for row in center['coordinates']]);_,_,_,records=connect(f,arb(center['a']),arb(0),beta,beta,p,pb,pb,pbb,V,V.inv(),center['iterations'])
    starts=[p+V*arb_mat([[arb(center['a'])],[0]])]+[r['state'] for r in records[:-1]];model=LobeMap(center['parameters']);curve=[]
    for start in starts:
        q=np.array([float(v.mid()) for v in start.entries()]);_,_,aux=model.apply(q,float(beta.mid()));mats=model.matrices(float(beta.mid()))
        for arc in aux['arcs']:
            A,c,e,U,Wi=mats[arc['region']];state=np.array(arc['start']);ts=np.linspace(0,arc['duration'],80)
            curve.extend(c+np.real(np.exp(ts[:,None]*e)@(np.diag(Wi@(state-c))@U.T)))
    curve=np.array(curve);fig=plt.figure(figsize=(16,7));gs=fig.add_gridspec(1,3,width_ratios=[1.1,1,1],left=.05,right=.97,bottom=.20,top=.78,wspace=.32)
    fig.suptitle('A homoclinic tangency in a region of hidden periodic attractors',x=.05,y=.97,ha='left',fontsize=21,weight='bold')
    fig.text(.05,.89,'72 transverse affine flights  |  Interval checks at 140 and 160 digits  |  Equilibrium exclusion uniform in β',fontsize=12)
    ax=fig.add_subplot(gs[0,0]);ax.plot(curve[:,0],curve[:,2],color=BLUE,lw=.7,alpha=.72)
    eq=(center['parameters']['m1']-center['parameters']['m0'])/(1+center['parameters']['m1']);ax.scatter(eq,-eq,c=RED,s=35,zorder=5,label='Stable outer equilibrium')
    ax.set(xlabel='x',ylabel='z');ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.15)
    ax.set_title('A  Numerical connection near β*',loc='left',fontsize=12,weight='bold')
    ax=fig.add_subplot(gs[0,1]);labels=['β lower face','β upper face','a lower face','a upper face']
    for y,key,scale in zip(range(4),['beta_lower','beta_upper','a_lower','a_upper'],[1e-11,1e-11,30,30]):
        ball=arb(cert['face_signs'][key]['ball']);lo=float(ball.lower())/scale;hi=float(ball.upper())/scale
        ax.plot([lo,hi],[y,y],color=BLUE,lw=5,solid_capstyle='round');ax.scatter((lo+hi)/2,y,color=BLUE,s=14)
    ax.axvline(0,color=RED,lw=1,ls='--');ax.set(yticks=range(4),yticklabels=labels,xlabel='Signed interval / displayed scale',xlim=(-1.5,1.5));ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    ax.set_title('B  Four strict signs for exact manifolds',loc='left',fontsize=12,weight='bold')
    ax.text(0,-.19,'Parameter faces: scale 10⁻¹¹\nDerivative faces: scale 30',transform=ax.transAxes,fontsize=9,linespacing=1.5)
    ax=fig.add_subplot(gs[0,2]);ax.axis('off')
    steps=[('Verified local dynamics','Quadratic tangency; nonzero unfolding;\ndissipative saddle; complete flight paths.'),('Classical parameter theorem','Some nearby β values support infinitely\nmany simultaneous periodic sinks.'),('Verified equilibrium exclusion','Every one of these attracting cycles\nhas a hidden basin.')]
    for k,(title,body) in enumerate(steps):
        y=.98-k*.34;ax.text(0,y,title,fontsize=12,weight='bold',va='top',color=TEAL if k==2 else INK)
        ax.text(0,y-.09,body,fontsize=10,va='top',linespacing=1.6)
        if k<2:ax.annotate('',xy=(.48,y-.30),xytext=(.48,y-.22),arrowprops=dict(arrowstyle='->',color=INK))
    ax.set_title('C  The logical implication',loc='left',fontsize=12,weight='bold',pad=15)
    fig.text(.05,.055,'The numerical curve illustrates the connection. Infinity follows from the certified hypotheses and the cited theorem; no decimal β∞ or infinite orbit list is computed.',fontsize=10)
    save(fig,'chua-followup-hidden')


if __name__=='__main__':sheets();pulse();hidden()
