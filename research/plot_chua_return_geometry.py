"""Publication figures from recorded Chua return-geometry measurements."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.publication/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
import numpy as np
from study_chua_return_geometry import OUT, IDS, read
FIGURES = ROOT/'docs/assets' if OUT == ROOT/'results/chua-return-geometry' else OUT

BG='#fafaf7';INK='#182d39';BLUE='#167c9d';ORANGE='#c66b2b';GREY='#70818a'


def setup():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'text.color':INK,
        'axes.labelcolor':INK,'axes.edgecolor':GREY,'axes.facecolor':BG,
        'figure.facecolor':BG,'axes.spines.top':False,'axes.spines.right':False,
        'xtick.color':INK,'ytick.color':INK,'axes.titleweight':'bold'})


def save(fig,name):
    fig.savefig(FIGURES/f'{name}.png',dpi=150,facecolor=BG)
    fig.savefig(OUT/f'{name}.pdf',facecolor=BG)
    plt.close(fig)


def sequence():
    rows=read(OUT/'sequence-analysis.json')['rows']
    for r in rows:r['gain']=read(OUT/f"{r['id']}-optimized.json")['maximum']['gain']
    for id in ['P95','P105']:
        r=read(OUT/f'{id}-reference.json');e=read(OUT/f'{id}-edges.json')
        g=read(OUT/f'{id}-optimized.json');c=read(OUT/f'{id}-companion-d70.json')
        rows.append(dict(id=id,n=int(id[1:])//10,peaks=int(id[1:]),period=float(r['period']),
            width=float(e['width']),gain=g['maximum']['gain'],
            anchor_pair_separation=float(c['section_distance']),
            companion_half_multiplier=float(c['companion']['half_multiplier'])))
    return rows


def overview(rows):
    fig,axs=plt.subplots(2,2,figsize=(14.5,9.3))
    fig.subplots_adjust(left=.08,right=.97,bottom=.14,top=.86,hspace=.48,wspace=.28)
    fig.suptitle('Why this sequence has thinner and more sensitive layers',fontsize=21,y=.97,x=.5)
    fig.text(.5,.916,'Extra returns near one saddle predict both parameter precision and temporary amplification.',ha='center',fontsize=13)
    ax=axs[0,0]
    for id,col in [('S04',GREY),('S07',BLUE),('P85',ORANGE)]:
        r=next(x for x in rows if x['id']==id)
        d=np.array(r['section_distance_to_parent']);d=np.r_[d[1:],d[0]]
        ax.semilogy(np.arange(len(d))-len(d)+1,d,'o-',color=col,lw=1.7,ms=4,label=f"{r['peaks']} peaks")
    ax.set(title='1  More returns near the same saddle',xlabel='Section returns before the excursion',ylabel='Distance to the saddle on the section')
    ax.legend(frameon=False,fontsize=10,loc='upper left');ax.grid(alpha=.16)
    ax=axs[0,1]
    with np.load(OUT/'P85-d85-s128.npz') as z:
        ax.semilogy(z['worst_times'],z['worst_history'],color=ORANGE,lw=2)
    ax.axhline(1,color=GREY,lw=1,ls='--')
    ax.text(.04,.95,'85-peak cycle\nPeak ≈ 620,000×',ha='left',va='top',transform=ax.transAxes,fontsize=11)
    ax.set(title='2  Large growth during an attracting cycle',xlabel='Elapsed time from a sensitive phase',ylabel='Maximum transverse linear gain')
    ax.grid(alpha=.16)
    ax=axs[1,0]
    widths=np.array([r['width'] for r in rows]);gains=np.array([r['gain'] for r in rows])
    xx=np.geomspace(widths[-1]*.65,widths[0]*1.5,200)
    ax.loglog(xx,np.sqrt(rows[7]['width']*rows[7]['gain']**2/xx),'--',color=GREY,label='85-peak calibration: G ∝ W⁻¹ᐟ²')
    ax.loglog(widths[:8],gains[:8],'o',color=BLUE,ms=6,label='15–85 peaks')
    ax.loglog(widths[8:],gains[8:],'D',color=ORANGE,ms=7,label='New tests: 95 and 105 peaks')
    for i,offset in [(0,(-32,10)),(7,(8,0)),(9,(10,-10))]:
        ax.annotate(str(rows[i]['peaks']),xy=(widths[i],gains[i]),xytext=offset,textcoords='offset points',fontsize=10)
    ax.set(title='3  Narrower windows, larger transient gain',xlabel='Symmetric-cycle stability width W in β',ylabel='Refined maximum linear gain G')
    ax.legend(frameon=False,fontsize=9,loc='lower left');ax.grid(alpha=.16)
    ax=axs[1,1]
    ss=np.linspace(-.5,1.5,150);ax.plot(ss,(ss+ss*ss)/2,color=GREY,lw=2,label='Prediction: (s + s²)/2')
    for id,col,mark in [('P55',BLUE,'o'),('P85',ORANGE,'x')]:
        r=read(OUT/f'{id}-normal-form-d75.json')['rows']
        ax.plot([x['s'] for x in r],[float(x['return_coordinate']) for x in r],mark,color=col,ms=6,label=id)
    ax.plot(ss,ss,':',color='#a3afb4',label='Fixed-point line')
    ax.annotate('Attracting cycle',(0,0),xytext=(-.48,.45),arrowprops={'arrowstyle':'-','color':GREY},fontsize=10)
    ax.annotate('Saddle cycle',(1,1),xytext=(.40,1.72),arrowprops={'arrowstyle':'-','color':GREY},fontsize=10)
    ax.set(title='4  The same return curve after rescaling',xlabel='Initial position between the two cycles, s',ylabel='Position after the return')
    ax.legend(frameon=False,fontsize=9,loc='lower right');ax.grid(alpha=.16)
    fig.text(.5,.025,'Numerical evidence for one period-adding sequence at α = 9, m₀ = −8/7, m₁ = −5/7.\nLinear gain is measured transverse to the flow; it is not a finite basin-escape threshold.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-return-explanation')


def phase():
    fig,axs=plt.subplots(2,2,figsize=(14.5,10))
    fig.subplots_adjust(left=.08,right=.92,bottom=.10,top=.88,hspace=.37,wspace=.27)
    fig.suptitle('The same cycle has very different sensitive phases',fontsize=21,y=.97)
    fig.text(.5,.92,'Color marks the largest linear growth available during the next period, starting at that point.',ha='center',fontsize=12)
    for col,(id,label) in enumerate([('S04','15 peaks'),('P85','85 peaks')]):
        with np.load(OUT/f'{id}-d85-s128.npz') as z:
            xyz=z['xyz'];g=z['gain'];t=z['t'];period=float(z['period'])
        ax=axs[0,col]
        points=xyz[:,[0,2]]
        segments=np.stack([points,np.roll(points,-1,axis=0)],axis=1)
        lc=LineCollection(segments,cmap='viridis',norm=LogNorm(10,1e6),linewidths=1.6)
        lc.set_array(g);ax.add_collection(lc);ax.autoscale();ax.set_aspect('equal',adjustable='datalim')
        ax.set(title=f'{label}: x–z projection',xlabel='x',ylabel='z')
        ax=axs[1,col];ax.semilogy(t/period,g,color=BLUE,lw=1.5)
        ax.set(title='Sensitivity depends on when the perturbation arrives',xlabel='Phase t/T',ylabel='Maximum transverse gain over the next T')
        ax.set_ylim(10,1e6);ax.grid(alpha=.18)
    cax=fig.add_axes([.94,.55,.015,.30]);fig.colorbar(lc,cax=cax,label='Linear gain (same color scale)')
    fig.text(.5,.025,'Both references have full-period Floquet spectral radius 0.25. This is asymptotic attraction,\nnot a guarantee that every perturbation shrinks by four in one period. Projected crossings are not orbit intersections.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-return-phase-sensitivity')


def forecasts(rows):
    forecast=read(OUT/'forecasts-frozen.json')['targets']
    fig,axs=plt.subplots(1,3,figsize=(15,5.4))
    fig.subplots_adjust(left=.06,right=.98,bottom=.23,top=.75,wspace=.34)
    fig.suptitle('Predictions published before the 95- and 105-peak searches',fontsize=20,y=.97)
    fig.text(.5,.89,'The same 85-peak calibration predicts three different properties of the new cycles.',ha='center',fontsize=12)
    specs=[('width','width','Stability width in β'),('gain','maximum_gain','Maximum transverse gain'),('anchor_pair_separation','companion_section_distance','Saddle separation on the inherited section')]
    for ax,(actual,predicted,title) in zip(axs,specs):
        vals=[]
        for target in forecast:
            r=next(x for x in rows if x['id']==target['id'])
            vals.append(r[actual]/float(target[predicted]))
        ax.bar([0,1],np.array(vals)-1,bottom=1,color=[BLUE,ORANGE],width=.55)
        ax.axhline(1,color=GREY,lw=1,ls='--')
        for x,v in enumerate(vals):ax.text(x,v+.0003,f'{100*(v-1):+.3f}%',ha='center',fontsize=12)
        ax.set_xticks([0,1],['95 peaks','105 peaks']);ax.set_title(title,fontsize=12)
        ax.set_ylim(.997,1.012);ax.set_ylabel('Measured / predicted');ax.grid(axis='y',alpha=.15)
    fig.text(.5,.07,'Forecast commit: 8df8da8. Locations were corrected by shooting and were not part of this forecast.\nBars use an expanded scale around 1 so the errors remain visible; no forecast is replaced after the search.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-return-forecast-check')


def catalogue():
    rows=[]
    for p in sorted(OUT.glob('*-optimized.json')):
        r=read(p)
        if r['id'] in ['P95','P105']:continue
        rows.append(r)
    rows.sort(key=lambda r:r['maximum']['gain'])
    fig,ax=plt.subplots(figsize=(13.8,10))
    fig.subplots_adjust(left=.10,right=.96,bottom=.14,top=.87)
    colors=[ORANGE if r['id'] in IDS else BLUE for r in rows]
    y=np.arange(len(rows));ax.barh(y,[r['maximum']['gain'] for r in rows],color=colors,height=.68)
    ax.set_yticks(y,[r['id'] for r in rows]);ax.set_xscale('log');ax.set_xlim(30,1e6)
    ax.set_xlabel('Refined maximum transverse linear gain over at most one period')
    ax.grid(axis='x',alpha=.15)
    fig.suptitle('All 34 displayed reference cycles can amplify transverse perturbations',fontsize=19,y=.97)
    fig.text(.5,.92,'Orange: the eight known members of the tested sequence. Blue: the other reference families.',ha='center',fontsize=12)
    fig.text(.5,.025,'Each cycle uses its own reference β. These are local refinements of sampled maxima in dimensionless Euclidean norm.\nThe sequence law does not classify the mechanisms organizing the other families.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-return-catalogue')


def main():
    setup();rows=sequence();overview(rows);phase();forecasts(rows);catalogue()
    print('Rendered four PNG/PDF figures',flush=True)


if __name__=='__main__':main()
