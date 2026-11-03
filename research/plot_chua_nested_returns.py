"""Figures distinguishing nested returns from ordinary period adding."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.publication/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from test_chua_nested_returns import OUT,read

BG='#fafaf7';INK='#182d39';GREY='#74858b';BLUE='#167c9d';PURPLE='#9a4f98';ORANGE='#cb6e30'
COLORS=[BLUE,PURPLE,ORANGE]


def save(fig,name):
    fig.savefig(ROOT/f'docs/assets/{name}.png',dpi=165,facecolor=BG)
    fig.savefig(OUT/f'{name}.pdf',facecolor=BG);plt.close(fig)


def profile(ax,u,blocks=None):
    u=np.array(u);x=np.arange(1,len(u)+1)
    ax.axhline(0,color=GREY,lw=.7)
    ax.plot(x,u,color=GREY,lw=.85,zorder=1)
    for j in range(len(u)):
        col=COLORS[(j//7)%3] if blocks else BLUE
        ax.vlines(x[j],0,u[j],color=col,lw=2)
        ax.scatter(x[j],u[j],color=col,s=22,zorder=3)
    if blocks:
        for j in range(blocks):
            ax.axvspan(j*7+.5,(j+1)*7+.5,color=COLORS[j%3],alpha=.055)
    ax.set_xlim(.5,len(u)+.5);ax.set_xticks(np.arange(1,len(u)+1,2 if len(u)>10 else 1))
    ax.set_xlabel('Successive normalized section crossings')
    ax.set_ylabel('Signed distance from the short saddle')
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))


def overview():
    fig=plt.figure(figsize=(14,10.5));gs=fig.add_gridspec(2,2,left=.085,right=.97,bottom=.10,top=.70,hspace=.58,wspace=.28)
    fig.suptitle('A long motion can organize still longer motions',fontsize=23,y=.975)
    fig.text(.5,.93,'One additional level tested in two Chua families; arbitrary depth remains unproved.',ha='center',fontsize=12)
    for x,title,body in [(.19,'5-peak short cycle','Earlier first-level sequence'),(.50,'35-peak family','Continue through its fold to a saddle'),(.81,'105 / 175 / 245 peaks','New cycles revisit that long saddle')]:
        fig.text(x,.845,title,ha='center',fontsize=17,color=INK,bbox=dict(boxstyle='round,pad=.65',fc='#e7eef0',ec='none'))
        fig.text(x,.784,body,ha='center',fontsize=10,color=GREY)
    fig.text(.345,.842,'→',fontsize=29,ha='center');fig.text(.66,.842,'→',fontsize=29,ha='center')
    a=read(OUT/'S06-k3-geometry.json');control=read(OUT/'negative-controls.json')['ordinary_105']
    ax=fig.add_subplot(gs[0,0]);profile(ax,a['short_saddle_coordinate'],3)
    ax.set_title('New 105-peak cycle: three seven-crossing motifs',loc='left',fontsize=12,pad=14)
    ax=fig.add_subplot(gs[0,1]);profile(ax,control['short_saddle_coordinate'])
    ax.set_title('Earlier 105-peak cycle: one longer sequence',loc='left',fontsize=12,pad=14)
    for ix,parent in enumerate(['S06','S11-n3']):
        ax=fig.add_subplot(gs[1,ix])
        for k,col in zip([3,5,7],COLORS):
            r=read(OUT/f'{parent}-k{k}-geometry.json');d=r['block_rms_distance_to_parent']
            ax.semilogy(np.arange(1,k+1),d,'o-',color=col,ms=5,label=f'{35*k} peaks ({k} blocks)')
        ax.axhline(.1,color=GREY,lw=.8,ls=':');ax.set_xticks(range(1,8))
        ax.set(title=f'{parent}: longer cycles follow the long saddle more closely',xlabel='Successive parent half-return blocks',ylabel='Distance / separation of the two parent branches')
        ax.legend(frameon=False,fontsize=9,loc='lower right');ax.grid(alpha=.15)
    fig.text(.5,.025,'Upper panels: each orbit uses its own same-parameter short saddle. Lower panels: each block compares all seven section points.\nThe dotted 0.1 level is a display guide, not an acceptance threshold. The figures show measured orbits, not an infinite construction.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-nested-returns-explanation')


def motifs():
    fig,axs=plt.subplots(3,1,figsize=(13.5,10),gridspec_kw=dict(height_ratios=[1,1.25,1.25]))
    fig.subplots_adjust(left=.09,right=.96,bottom=.11,top=.85,hspace=.70)
    fig.suptitle('The same peak count does not identify the same motion',fontsize=23,y=.97)
    fig.text(.5,.92,'S06 comparison · the repeated seven-crossing pattern resolves what an orbit outline hides',ha='center',fontsize=12)
    g=read(OUT/'S06-k3-geometry.json');c=read(OUT/'negative-controls.json')['ordinary_105']
    rows=[(g['parent_short_saddle_coordinate'],'35-peak organizing saddle: one motif',1),
          (g['short_saddle_coordinate'],'New 105-peak child: the motif returns three times',3),
          (c['short_saddle_coordinate'],'Earlier 105-peak cycle: a single 21-crossing word',None)]
    for ax,(u,title,b) in zip(axs,rows):
        profile(ax,u,b);ax.set_title(title,loc='left',pad=15)
    fig.text(.5,.025,'Signs record the side of the short saddle on a specified section, not left/right scroll occupancy.\nThe repeated binary word is a coarse diagnostic; distinct section points and independent orbit checks exclude multiple traversal.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-nested-returns-motifs')


def predictions():
    scores=read(OUT/'scores.json')['rows'];fig,axs=plt.subplots(1,3,figsize=(14,6.2))
    fig.subplots_adjust(left=.075,right=.975,bottom=.24,top=.73,wspace=.36)
    fig.suptitle('Frozen forecasts: widths and periods pass; four locations miss',fontsize=21,y=.97)
    fig.text(.5,.89,'Six targets · no coefficients refitted after the searches',ha='center',fontsize=12)
    for ix,(key,title,label,limit) in enumerate([
        ('location_in_widths','Location: 2 / 6 pass','Absolute error / measured window width',1),
        ('width_relative','Width: 6 / 6 pass','Absolute relative error (%)',30),
        ('period_absolute','Period: 6 / 6 pass','Absolute period error',.5)]):
        ax=axs[ix]
        for j,r in enumerate(scores):
            val=abs(r['errors'][key])*(100 if key=='width_relative' else 1)
            ax.bar(j,val,color=BLUE if j<3 else PURPLE,width=.65)
            ax.text(j,val*1.2,f'{val:.3g}',ha='center',va='bottom',fontsize=9)
        ax.set_yscale('log');ax.axhline(limit,color=ORANGE,ls='--',lw=1)
        ax.set(title=title,ylabel=label,xticks=range(6),xticklabels=['105','175','245']*2)
        ax.set_xlabel('Full-period peaks\nS06                         S11-n3')
        ax.set_ylim([.08,100000] if ix==0 else [.5,80] if ix==1 else [.001,1.3]);ax.grid(axis='y',alpha=.15)
    fig.text(.5,.055,'Dashed lines: the tolerances published before searching (one width, 30%, 0.5 time units).\nA small absolute β bias becomes tens of thousands of widths when the stability interval is extremely thin.',ha='center',fontsize=11,color=GREY)
    save(fig,'chua-nested-returns-predictions')


def maps():
    data=read(OUT/'parent-block-map.json')['rows'];moving=read(OUT/'moving-coordinate-map.json')['rows'];fig,axs=plt.subplots(2,2,figsize=(13.5,10))
    fig.subplots_adjust(left=.085,right=.97,bottom=.10,top=.84,hspace=.45,wspace=.28)
    fig.suptitle('The parent return resembles a quadratic map, with finite errors',fontsize=21,y=.97)
    fig.text(.5,.918,'Original calibration held fixed · a projected two-dimensional return, not a proved one-dimensional model',ha='center',fontsize=11)
    for col,row in enumerate(data):
        ax=axs[0,col]
        for c,color in zip(row['curves'],COLORS):
            x=np.array([r['x'] for r in c['samples']]);actual=np.array([r['actual'] for r in c['samples']])
            ax.plot(x,actual,color=color,lw=1.5,label=f"Actual return, M = {c['M']}")
        c=row['curves'][-1];ax.plot(x,2-x*x,color=INK,lw=1,ls='--',label='Quadratic prediction at M = 2')
        ax.plot(x,x,color=GREY,lw=.7,alpha=.5);ax.set(title=row['parent'],xlabel='Input X',ylabel='Projected return X');ax.legend(frameon=False,fontsize=9)
        ax=axs[1,col];r=moving[col];x=np.array([s['x'] for s in r['samples']])
        ax.plot(x,[s['actual'] for s in r['samples']],color=BLUE,lw=1.7,label='Actual return in moving coordinates')
        ax.plot(x,[s['quadratic'] for s in r['samples']],color=INK,lw=1,ls='--',label='Quadratic, calibrated after the search')
        ax.plot(x,x,color=GREY,lw=.7,alpha=.5)
        ax.text(.04,.08,f"Maximum sampled error: {r['maximum_scalar_error']:.3f}",transform=ax.transAxes,fontsize=10)
        ax.set(title='Diagnostic at the 245-peak child’s parameters',xlabel='X from the two same-parameter parent roots',ylabel='Projected return X');ax.legend(frameon=False,fontsize=8,loc='upper left')
    fig.text(.5,.025,'Top: original fixed secant and parameter scale. Bottom: new coordinates from the two parent roots and saddle multiplier.\nThe lower comparison is a post-search diagnostic. It does not replace the frozen forecasts or prove an invariant curve.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-nested-returns-maps')


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,'axes.labelcolor':INK,
        'axes.titleweight':'medium','axes.spines.top':False,'axes.spines.right':False,'axes.facecolor':BG,
        'xtick.color':GREY,'ytick.color':GREY,'axes.edgecolor':'#bdc8cb'})
    overview();motifs();predictions();maps()


if __name__=='__main__':main()
