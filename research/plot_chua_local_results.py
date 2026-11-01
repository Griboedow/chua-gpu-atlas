"""Scientific figures from the strip certificates and a labelled theorem schematic."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyBboxPatch,FancyArrowPatch
from flint import arb
from audit_chua_local_results import ROOT,OUT,leaves

BLUE,ORANGE,INK,TEAL='#276a9c','#bd5735','#233442','#237765'

def read(path):return json.loads((ROOT/path).read_text())
def number(v):return float(v['midpoint'])
def save(fig,name):
    fig.savefig(ROOT/'docs/assets'/f'{name}.png',dpi=180,bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight',facecolor='white',metadata={'CreationDate':None,'ModDate':None})
    plt.close(fig)


def basin():
    r=read('results/chua-basin-proof/separator-100.json');assert r['passed']
    controls=read('results/chua-basin-proof/independent-checks.json');assert controls['passed']
    fig,axes=plt.subplots(1,3,figsize=(16,6.0),gridspec_kw={'width_ratios':[1,1.2,1.15]})
    fig.subplots_adjust(left=.06,right=.98,bottom=.27,top=.77,wspace=.40)
    fig.suptitle('A local boundary between two Chua attraction basins',x=.055,ha='left',y=.97,fontsize=22,weight='bold',color=INK)
    fig.text(.055,.89,r'$\alpha=8.4,\ \beta=12,\ m_0=-0.95,\ m_1=-0.05$'+'   |   Interval coverage at 80 and 100 digits',fontsize=12,color=INK)
    ax=axes[0];s=np.linspace(-1,1,401);w=.02*abs(s)
    ax.fill_betweenx(s,-1,-w,color=BLUE,alpha=.20);ax.fill_betweenx(s,w,1,color=ORANGE,alpha=.20)
    ax.fill_betweenx(s,-w,w,color='#8b8795',alpha=.65)
    ax.plot(-w,s,color='#635767',lw=1);ax.plot(w,s,color='#635767',lw=1)
    ax.scatter([0],[0],c=INK,s=30,zorder=4)
    ax.text(-.55,.43,'Origin',ha='center',color=BLUE,weight='bold',fontsize=12)
    ax.text(.55,.43,'Hidden\ncycle',ha='center',color=ORANGE,weight='bold',fontsize=12)
    ax.annotate('Exact graph lies\ninside this wedge',xy=(.01,-.68),xytext=(.32,-.6),fontsize=10,
                arrowprops={'arrowstyle':'-','color':INK},ha='left')
    ax.set(xlim=(-1,1),ylim=(-1,1),xlabel=r'Transverse coordinate $u/\delta$',ylabel=r'Stable coordinate $s/\delta$')
    ax.set_title('A  Certified local square',loc='left',fontsize=12,weight='bold',pad=12)
    ax.text(.02,-.34,r'$\delta=10^{-7}$; stable-graph slope $\leq 0.02$',transform=ax.transAxes,fontsize=10)
    ax=axes[1];maximum=number(r['exit_u_max'])/1e-7
    for side,y,color,label in [(-1,1,BLUE,'Origin: 64 / 64'),(1,0,ORANGE,'Cycle: 64 / 64')]:
        for e in (e for e in r['exits'] if e['side']==side):
            left=1+(maximum-1)*e['index']/64;width=(maximum-1)/64
            ax.add_patch(Rectangle((left,y-.15),width,.30,facecolor=color,edgecolor='white',lw=.4,
                                    hatch='///' if e['continuous'] else None))
        ax.text(1,y+.30,label,color=color,weight='bold',fontsize=12)
    ax.set(xlim=(.99,maximum+.01),ylim=(-.45,1.60),xlabel=r'Exit distance $|u|/\delta$',yticks=[])
    ax.set_title('B  Whole strips, not sampled starts',loc='left',fontsize=12,weight='bold',pad=12)
    ax.text(0,-.34,'Hatching: continuous flow across diode contact.\nEach rectangle includes the entire stable-coordinate range.',transform=ax.transAxes,fontsize=10,linespacing=1.5)
    for spine in ['top','right','left']:ax.spines[spine].set_visible(False)
    ax=axes[2]
    # Use a quadratic-norm trace: its radius is recorded at the displayed time.
    choices=[f for e in r['exits'] for f in leaves(e['continuous']) if f['method']=='quadratic_norm']
    f=choices[len(choices)//2];trace=f['trace'];t=np.array([q['time'] for q in trace]);x=np.array([q['point'][0] for q in trace])
    P=np.array(r['origin_capture']['lyapunov']['P'],dtype=float);scale=np.sqrt(np.linalg.inv(P)[0,0])
    rad=np.array([q['radius'] for q in trace])*scale;select=(t>=.45)&(t<=1.1)
    ax.errorbar(t[select],x[select],yerr=rad[select],fmt='o-',color=BLUE,ecolor=TEAL,capsize=3,ms=4,lw=1.3,label='Enclosed states at recorded times')
    ax.axhline(1,color=ORANGE,lw=1.3,ls='--',label='Diode switching surface')
    ax.set(xlabel='Time after the last section return',ylabel=r'$x$',ylim=(.87,1.06))
    ax.set_title('C  Passing through the contact',loc='left',fontsize=12,weight='bold',pad=12)
    ax.legend(frameon=False,fontsize=9,loc='lower right');ax.grid(alpha=.18)
    ax.text(0,-.34,'Error bars enclose whole sets at the marked times.\nPicard tubes certify the intervals between them.',transform=ax.transAxes,fontsize=10,linespacing=1.5)
    for ax in axes:
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    save(fig,'chua-basin-proof')


def coexistence():
    fig,axes=plt.subplots(1,2,figsize=(14,6.8),gridspec_kw={'width_ratios':[1.25,1]})
    fig.subplots_adjust(left=.04,right=.98,bottom=.12,top=.78,wspace=.12)
    fig.suptitle('Infinite coexistence: why the return-map sinks belong to Chua',x=.04,ha='left',y=.96,fontsize=20,weight='bold',color=INK)
    fig.text(.04,.88,r'$\alpha=9,\ m_0=-8/7,\ m_1=-5/7$'+'   |   Theorem schematic; parameter intervals are not measured',fontsize=12,color=INK)
    ax=axes[0];ax.set(xlim=(0,10),ylim=(0,6));ax.axis('off')
    ax.add_patch(FancyBboxPatch((.2,.9),9.1,4.8,boxstyle='round,pad=.1',facecolor='#f1f6f6',edgecolor=TEAL,ls='--',lw=1.7))
    ax.text(.45,5.3,'Agreement region Ω: original Chua return map',color=TEAL,fontsize=12,weight='bold')
    for x,y,w,h,text in [(1,2.4,3,1.3,'Local saddle\nrectangle U'),(5.7,2.4,2.8,1.3,'Outgoing\nrectangle V')]:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.08',facecolor='white',edgecolor=BLUE,lw=1.5));ax.text(x+w/2,y+h/2,text,ha='center',va='center',color=INK,fontsize=13)
    ax.add_patch(FancyArrowPatch((4.1,3.1),(5.6,3.1),arrowstyle='-|>',mutation_scale=16,color=BLUE,lw=2))
    ax.text(4.8,3.55,r'$T_0^n$',ha='center',fontsize=15,color=BLUE)
    ax.add_patch(FancyArrowPatch((6.5,2.3),(2.5,2.3),connectionstyle='arc3,rad=-.45',arrowstyle='-|>',mutation_scale=16,color=ORANGE,lw=2))
    ax.text(4.5,1.23,r'$T_1$: finite excursion',ha='center',color=ORANGE,fontsize=12)
    ax.text(.3,.15,'Every intermediate iterate stays in Ω.\nThe arbitrary compact extension is never used by these orbits.',fontsize=11,color=INK,linespacing=1.5)
    ax=axes[1];ax.set(xlim=(0,1),ylim=(0,6));ax.axis('off')
    for n,(l,r,y) in enumerate([(.07,.93,4.8),(.27,.78,3.6),(.39,.68,2.4),(.47,.61,1.2)],1):
        ax.plot([l,r],[y,y],color=BLUE,lw=10,solid_capstyle='butt',alpha=.32+.12*n)
        ax.text(.03,y+.34,'One attracting cycle persists' if n==1 else f'At least {n} attracting cycles persist',fontsize=12,color=INK)
    ax.axvline(.54,ymin=.12,ymax=.84,color=ORANGE,ls='--',lw=1.5)
    ax.text(.58,.39,r'$\beta_\infty$: infinitely many sinks',fontsize=13,color=ORANGE,weight='bold')
    ax.text(.03,5.4,'Nested parameter construction',fontsize=13,weight='bold',color=INK)
    fig.text(.04,.035,'Existence near the certified tangency; no decimal witness or positive sampling probability is asserted.',fontsize=11,color=INK)
    save(fig,'chua-coexistence-localization')


if __name__=='__main__':
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
    basin();coexistence()
