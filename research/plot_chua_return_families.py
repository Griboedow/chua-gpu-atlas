"""Publication figures for the comparison of three return sequences."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.publication/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from test_chua_return_families import OUT,read
from plot_chua_return_geometry import sequence

BG='#fafaf7';INK='#182d39';GREY='#75878c';COLORS=['#167c9d','#9a4f98','#cb6e30']


def save(fig,name):
    fig.savefig(ROOT/f'docs/assets/{name}.png',dpi=160,facecolor=BG)
    fig.savefig(OUT/f'{name}.pdf',facecolor=BG);plt.close(fig)


def sequences():
    old=sequence()[:5];calibration=read(OUT/'calibration.json');result=[old]
    for seq in calibration['sequences']:
        rows=[dict(id=r['id'],n=n+1,width=float(r['width']),gain=r['gain'],period=float(r['representative']['period'])) for n,r in enumerate(seq['calibration'])]
        for n in [3,4,5]:
            id=f"{seq['id']}-n{n}";s=next(r for r in read(OUT/'scores.json')['rows'] if r['id']==id)
            rows.append(dict(id=id,n=n,**s['measured']))
        result.append(rows)
    return result


def overview():
    fig,axs=plt.subplots(2,2,figsize=(14,9.5))
    fig.subplots_adjust(left=.085,right=.965,bottom=.12,top=.85,hspace=.49,wspace=.29)
    fig.suptitle('Layers record different repeating motions',fontsize=23,y=.975)
    fig.text(.5,.922,'In three tested sequences: extra turns near a short unstable cycle, followed by a global return.',ha='center',fontsize=12)
    ax=axs[0,0]
    with np.load(OUT/'S11-n5-gain-s128.npz') as z:u=z['xyz']
    with np.load(OUT/'S11-n5-parent.npz') as z:p=z['xyz']
    ax.plot(u[:,0],u[:,2],color=COLORS[1],lw=1,alpha=.7,label='Long attracting cycle: 55 peaks')
    ax.plot(p[:,0],p[:,2],color=INK,lw=1.65,label='Short unstable cycle: 5 peaks')
    ax.set(title='1  A long motion reuses a short orbit',xlabel='x',ylabel='z');ax.set_aspect('equal',adjustable='datalim')
    ax.legend(frameon=False,fontsize=9,loc='lower left')
    ax=axs[0,1]
    for id,col,peaks in [('S09',GREY,15),('S11-n3',COLORS[0],35),('S11-n5',COLORS[1],55)]:
        r=read(OUT/f'{id}-geometry.json');d=np.array(r['aligned_distance_to_parent']);d=np.r_[d,d[0]]
        ax.semilogy(np.arange(len(d)),d,'o-',color=col,ms=4,label=f'{peaks} peaks')
    ax.set(title='2  More section crossings near that orbit',xlabel='Successive crossings during one symmetric half',ylabel='Section distance to the short cycle')
    ax.legend(frameon=False,fontsize=10);ax.grid(alpha=.16)
    ax=axs[1,0];r=read(OUT/'S11-n5-directions.json');t=np.array(r['floquet_time']);a=np.array(r['floquet_amplitude'])
    ax.semilogy(t/t[-1],a,color=COLORS[1],lw=1.5);ax.axhline(1,color=GREY,ls='--',lw=.8)
    ax.scatter([.5,1],[.5,.25],color=INK,s=30,zorder=5)
    ax.annotate('End: 0.25 ×',xy=(1,.25),xytext=(.62,.07),textcoords='axes fraction',fontsize=10,arrowprops=dict(arrowstyle='-',color=GREY))
    ax.text(.04,.95,f'One fixed Floquet direction\nTemporary growth ≈ {a.max():,.0f} ×',transform=ax.transAxes,va='top',fontsize=10,bbox=dict(facecolor=BG,edgecolor='none',alpha=.85,pad=3))
    ax.set(title='3  Growth is undone during the return',xlabel='Fraction of the full period',ylabel='Transverse amplitude / initial amplitude');ax.grid(alpha=.16)
    ax=axs[1,1]
    for rows,col,label in zip(sequences(),COLORS,['S08 sequence','S11 sequence','E05 sequence']):
        ax.semilogy([r['n'] for r in rows],[r['width'] for r in rows],'o-',color=col,label=label,ms=5)
    ax.set(title='4  Additional turns need thinner windows',xlabel='Sequence member n',ylabel='Symmetric stability width in β',xticks=[1,2,3,4,5]);ax.grid(alpha=.16);ax.legend(frameon=False,fontsize=9)
    fig.text(.5,.035,'A layer is a range of parameter values supporting a repeating attracting motion.\nThese computations explain three sequences; they do not classify the whole atlas or give finite escape thresholds.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-return-families-explanation')


def comparison():
    fig,axs=plt.subplots(2,3,figsize=(15,9))
    fig.subplots_adjust(left=.065,right=.975,bottom=.13,top=.84,hspace=.55,wspace=.3)
    fig.suptitle('Different short cycles, similar local return geometry',fontsize=22,y=.975)
    fig.text(.5,.922,'Top: actual short saddle cycles at the long-cycle parameters. Bottom: independently evaluated long-cycle returns.',ha='center',fontsize=12)
    for j,(id,label,col) in enumerate(zip(['P55','S11-n5','E05-n5'],['S08 branch · 5 peaks','S11 branch · 5 peaks','E05 branch · 7 peaks'],COLORS)):
        with np.load(OUT/f'{id}-parent.npz') as z:xyz=z['xyz'];period=float(z['period'])
        ax=axs[0,j];ax.plot(xyz[:,0],xyz[:,2],color=col,lw=1.8);ax.set_aspect('equal',adjustable='datalim')
        ax.set(title=f'{label}\nShort period ≈ {period:.4f}',xlabel='x',ylabel='z')
        ax=axs[1,j];s=np.linspace(-.5,1.5,151);ax.plot(s,(s+s*s)/2,color=INK,lw=1.8,label='(s + s²)/2');ax.plot(s,s,':',color=GREY,label='Fixed-point line')
        if id=='P55':
            r=read(ROOT/'results/chua-return-geometry/P55-normal-form-d75.json');rows=r['rows'];xx=[v['s'] for v in rows];yy=[float(v['return_coordinate']) for v in rows];error=r['maximum_absolute_error']
        else:
            r=read(OUT/f'{id}-normal-form-d90.json');rows=r['samples'];xx=[v['s'] for v in rows];yy=[float(v['coordinate']) for v in rows];error=r['maximum_error']
        ax.plot(xx,yy,'o',mfc=BG,mec=col,ms=6,mew=1.7,label='Actual return samples')
        ax.set(title=f'{id}: largest discrepancy {error:.2g}',xlabel='Initial position between cycle and saddle, s',ylabel='Position after the signed half-return')
        ax.grid(alpha=.16)
        if j==0:ax.legend(frameon=False,fontsize=9,loc='lower right')
    fig.text(.5,.035,'s = 0: attracting cycle; s = 1: companion saddle at exactly the same parameters.\nNine secant samples per cycle; numerical agreement does not establish a global one-dimensional reduction.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-return-families-comparison')


def forecasts():
    rows=read(OUT/'scores.json')['rows'];fig,axs=plt.subplots(1,3,figsize=(15,6.5))
    fig.subplots_adjust(left=.08,right=.98,bottom=.21,top=.74,wspace=.38)
    fig.suptitle('Frozen forecasts: periods transfer; widths need better calibration',fontsize=21,y=.96)
    fig.text(.5,.875,'All six cycles were found. The failed width and gain forecasts are retained without refitting.',ha='center',fontsize=12)
    labels=[f"{r['id']}\n{r['peaks']} peaks" for r in rows];y=np.arange(6);cs=[COLORS[1]]*3+[COLORS[2]]*3
    for ax,key,title in zip(axs,['width','gain','period'],['Width: 0 / 6 within ±25%','Gain: 4 / 6 within ±20%','Period: 6 / 6 within ±0.5']):
        values=np.array([r['corrected_error'][key] for r in rows]);values*=1 if key=='period' else 100
        ax.barh(y,values,color=cs,height=.62);ax.axvline(0,color=GREY,lw=.8)
        if key!='period':
            gate=25 if key=='width' else 20;ax.axvline(-gate,ls='--',lw=1,color=GREY);ax.axvline(gate,ls='--',lw=1,color=GREY)
            ax.set_xlabel('Measured / forecast − 1 (%)')
        else:ax.set_xlabel('Measured − forecast (time units)')
        ax.set_yticks(y,labels if key=='width' else ['']*6);ax.invert_yaxis();ax.set_title(title,fontsize=12,pad=15);ax.grid(axis='x',alpha=.15)
    fig.text(.5,.065,'Forecasts published in commit a2444a0 before target searches; two old cycles calibrated each sequence.\nAcceptance gates are prescribed numerical tolerances, not confidence intervals. Window location was corrected by shooting.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-return-families-forecasts')


def contraction():
    fig,axs=plt.subplots(1,2,figsize=(13,5.7));fig.subplots_adjust(left=.08,right=.97,bottom=.18,top=.76,wspace=.27)
    fig.suptitle('Where a long attracting motion contracts',fontsize=21,y=.97)
    fig.text(.5,.88,'Follow one fixed perturbation direction through consecutive section returns during a symmetric half-period.',ha='center',fontsize=11)
    for id,col in [('S11-n5',COLORS[1]),('E05-n5',COLORS[2])]:
        r=read(OUT/f'{id}-audit.json')['section_variation'];f=np.array(r['floquet_direction_factors'])
        f=np.roll(f,-int(np.flatnonzero(f<1)[0]));a=np.r_[1.,np.cumprod(f)]
        axs[0].semilogy(np.arange(len(f))+1,f,'o-',color=col,label=id)
        axs[1].semilogy(np.arange(len(a)),a,'o-',color=col,label=id)
    axs[0].axhline(1,color=GREY,ls='--');axs[0].set(title='Two returns contract; nine expand',xlabel='Return number, aligned at the contraction',ylabel='Amplitude factor for this return')
    axs[1].axhline(.5,color=GREY,ls='--');axs[1].set(title='The complete half-period multiplies by 0.5',xlabel='Completed returns',ylabel='Amplitude / initial amplitude')
    for ax in axs:ax.legend(frameon=False,fontsize=10);ax.grid(alpha=.16)
    fig.text(.5,.045,'Dominant Floquet direction; Euclidean norm in section coordinates (y, z); cyclic phase alignment.\nThe number of contracting steps depends on the section; it is not a universal count of physical events.',ha='center',fontsize=10,color=GREY)
    save(fig,'chua-return-families-contraction')


if __name__=='__main__':
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':INK,'axes.labelcolor':INK,
        'axes.edgecolor':GREY,'axes.facecolor':BG,'figure.facecolor':BG,'axes.spines.top':False,
        'axes.spines.right':False,'axes.titleweight':'bold','xtick.color':INK,'ytick.color':INK})
    overview();comparison();forecasts();contraction()
