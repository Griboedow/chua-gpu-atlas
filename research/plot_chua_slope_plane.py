"""Draw an unsmoothed m0-m1 parameter map and a magnified classical neighborhood."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch,Rectangle
from matplotlib.lines import Line2D

COLORS=['#d7dce1','#54b795','#e69b3f','#eddc6b','#365f9e','#161616','#a66ba8']
LABELS=['Escape cutoff reached','Equilibrium candidate','Regular candidate','Uncertain','Chaos candidate','Numerical failure','Different labels across starts']


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--assets',type=Path,required=True)
    args=p.parse_args();args.assets.mkdir(parents=True,exist_ok=True)
    data={k:np.load(args.source/(k+'.npz')) for k in ['wide','zoom']}
    cmap=ListedColormap(COLORS);norm=BoundaryNorm(np.arange(-.5,7),7)
    plt.rcParams.update({'font.size':11,'text.color':'#193449','axes.labelcolor':'#193449','axes.titleweight':'normal'})
    fig,axes=plt.subplots(1,2,figsize=(14,7.6));fig.subplots_adjust(left=.065,right=.98,top=.81,bottom=.22,wspace=.23)
    for ax,(name,d) in zip(axes,data.items()):
        c=d['labels'].reshape(len(d['m1']),len(d['m0']))
        ax.pcolormesh(d['m0'],d['m1'],c,cmap=cmap,norm=norm,shading='nearest',rasterized=True)
        ax.plot(-8/7,-5/7,marker='*',ms=14,mfc='white',mec='#13283b',mew=1.2)
        ax.set(xlabel='Inner slope m0',ylabel='Outer slope m1',xlim=(d['m0'][0],d['m0'][-1]),ylim=(d['m1'][0],d['m1'][-1]),
            title='Wide section' if name=='wide' else 'Zoom around the classical slopes')
        ax.ticklabel_format(useOffset=False,style='plain');ax.spines[['top','right']].set_visible(False)
    z=data['zoom'];axes[0].add_patch(Rectangle((z['m0'][0],z['m1'][0]),np.ptp(z['m0']),np.ptp(z['m1']),fill=False,edgecolor='white',lw=1.4))
    axes[0].annotate('Zoom',xy=(-8/7,-5/7),xytext=(-1.27,-.57),color='white',fontsize=11,
        arrowprops=dict(arrowstyle='->',color='white',lw=1.1))
    if (args.source/'closed-cycles.json').exists():
        for r in json.loads((args.source/'closed-cycles.json').read_text())['cycles']:
            c=r['cycle'];ax=axes[1]
            ax.plot(c['m0'],c['m1'],'o',ms=5,mfc='white',mec='#13283b')
            ax.annotate(f'{r["peak_count"]} peaks',xy=(c['m0'],c['m1']),xytext=(c['m0']+.0007,c['m1']+.0015),
                fontsize=10,bbox=dict(boxstyle='round,pad=.2',fc='white',ec='none',alpha=.95),
                arrowprops=dict(arrowstyle='->',color='white',lw=1.1))
    fig.suptitle('Chua system in the plane of the two diode slopes',fontsize=21,weight='bold',y=.97)
    fig.text(.5,.90,'Fixed alpha = 9, beta = 100/7. Both m0 and m1 vary independently.',ha='center',fontsize=13)
    present=sorted(set(np.unique(data['wide']['labels']))|set(np.unique(data['zoom']['labels'])))
    handles=[Patch(facecolor=COLORS[k],label=LABELS[k]) for k in present]
    handles.append(Line2D([],[],marker='*',ms=12,mfc='white',mec='#13283b',ls='',label='Classical point'))
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.085),ncol=4,frameon=False,fontsize=10)
    fig.text(.5,.035,'241 x 241 nodes per panel; three initial states at each node. No smoothing or interpolation.\n'
        'Regular and chaotic labels are finite-time candidates. Differing labels do not by themselves prove coexistence.',ha='center',fontsize=10,color='#526673')
    fig.savefig(args.assets/'chua-slope-plane.png',dpi=180,facecolor='white');plt.close(fig)
    if (args.source/'zoom-peaks.npz').exists():
        peaks=np.load(args.source/'zoom-peaks.npz');counts=np.full(len(z['labels']),np.nan)
        period=peaks['periods'];groups=np.full(len(period),4,dtype=float);groups[period==0]=np.nan
        for group,value in enumerate([5,7,10,14]):groups[period==value]=group
        counts[peaks['indices']]=groups
        fig,ax=plt.subplots(figsize=(8.2,7.4),layout='constrained')
        ax.set_facecolor('#e8edf2')
        peak_colors=['#2685ad','#e29432','#784fa3','#4e9b77','#34485d']
        ax.pcolormesh(z['m0'],z['m1'],counts.reshape(len(z['m1']),len(z['m0'])),cmap=ListedColormap(peak_colors),
            norm=BoundaryNorm(np.arange(-.5,5),5),shading='nearest')
        ax.plot(-8/7,-5/7,'*',ms=13,mfc='white',mec='black')
        ax.set(xlabel='Inner slope m0',ylabel='Outer slope m1',title='Repeating peak patterns inside the regular bands\nalpha = 9, beta = 100/7; seed zero')
        ax.ticklabel_format(useOffset=False,style='plain')
        ax.legend(handles=[Patch(facecolor=color,label=label) for color,label in zip(peak_colors,['5 peaks','7 peaks','10 peaks','14 peaks','Other resolved counts'])],
            loc='upper center',bbox_to_anchor=(.48,.99),framealpha=.95,fontsize=9)
        fig.supxlabel('Gray: not tested here or no repeated pattern resolved.\nPeak counts are recurrence diagnostics, not distinct family identities.',fontsize=10)
        fig.savefig(args.assets/'chua-slope-peak-patterns.png',dpi=180,facecolor='white');plt.close(fig)
    summary={name:dict(nodes=len(d['labels']),shape=[len(d['m1']),len(d['m0'])],m0_range=[float(d['m0'][0]),float(d['m0'][-1])],
        m1_range=[float(d['m1'][0]),float(d['m1'][-1])],counts={LABELS[k]:int((d['labels']==k).sum()) for k in range(7)}) for name,d in data.items()}
    (args.source/'plot-summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
