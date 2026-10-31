"""Publication figures comparing preregistered layer forecasts and new evidence."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from predict_chua_layers import ROOT,OUT

ASSETS=ROOT/'docs/assets'
INK='#163449';GRAY='#647681';BLUE='#2477bd';ORANGE='#dc6d35';TEAL='#168476';PALE='#edf3f7'


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
        'text.color':INK,'axes.labelcolor':INK,'xtick.color':GRAY,'ytick.color':GRAY,
        'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold',
        'axes.edgecolor':'#bfcbd2','grid.color':'#dbe4e9','grid.alpha':.65})


def save(fig,name):
    fig.savefig(ASSETS/(name+'.png'),dpi=170,facecolor='white',bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),facecolor='white',bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
    plt.close(fig)


def forecasts(data):
    fig=plt.figure(figsize=(13,8),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1,.44],width_ratios=[1.6,1])
    fig.suptitle('Predicted first. Searched second. Compared at the true scale.',fontsize=20,weight='bold')
    for i,r in enumerate(data['comparison']):
        ax=fig.add_subplot(gs[i,0]);left,right=r['predicted_interval_in_measured_width_units']
        ax.barh(.68,right-left,left=left,height=.19,color=BLUE,alpha=.8,label='Frozen forecast')
        ax.barh(.29,1,left=-1,height=.19,color=ORANGE,label='Measured stable core')
        ax.axvline(0,color=INK,lw=.8,ls=':');ax.set(xlim=(-1.5,4.1),ylim=(0,1),yticks=[.29,.68],yticklabels=['Found','Forecast'])
        ax.set_title(f'{r["peaks"]} peaks per complete cycle',loc='left',fontsize=15)
        ax.set_xlabel('Position relative to the measured fold / measured core width')
        ax.grid(axis='x');ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
        note='Forecast overlaps the core' if r['predicted_midpoint_inside_measured_core'] else 'Forecast misses the core; correction finds it'
        ax.text(.02,.96,note,transform=ax.transAxes,va='top',fontsize=10,color=TEAL if i==0 else ORANGE)
        ax=fig.add_subplot(gs[i,1]);ax.axis('off')
        txt=(f'WIDTH\nPredicted   {r["predicted_width"]:.4e}\nMeasured  {float(r["measured_width"]):.4e}\n'
             f'Difference  {100*r["width_error_relative_to_prediction"]:.2f}% of forecast\n\n'
             f'FOLD POSITION\nError: {r["signed_fold_error_in_measured_widths"]:.3f} measured widths\n'
             f'Beta = {float(r["measured_fold_beta"]):.14f}\n\n'
             f'PERIOD\nPredicted {r["predicted_period"]:.6f}  |  Measured {r["measured_period"]:.6f}')
        ax.text(0,.94,txt,va='top',fontsize=11,linespacing=1.45)
    ax=fig.add_subplot(gs[2,:]);ax.axis('off')
    ax.text(.02,.84,'Forecast locked in Git: 0fc47ea   |   Calibration: the previously known 15–65-peak layers',weight='bold',fontsize=12)
    ax.text(.02,.43,'Search: prescribed itinerary + orbit closing, followed by 55/75-digit comparison and independent full-field integrations.',fontsize=10)
    ax.text(.02,.05,'The blue intervals are point forecasts, not confidence bands. Their engineering location gate was wider than either new core.',fontsize=10,color=GRAY)
    save(fig,'chua-layer-prediction-proof')


def scaling(data):
    rows=data['rows'];n=np.array([r['n'] for r in rows]);w=np.array([r['width'] for r in rows]);period=np.array([r['fold_period'] for r in rows])
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    fig.suptitle('One short saddle cycle predicts three scales of the hierarchy',fontsize=19,weight='bold')
    ax=axes[0,0];ax.semilogy(n[:6],w[:6],'o-',color=GRAY,label='Known before forecast');ax.semilogy(n[6:],w[6:],'o',color=ORANGE,ms=8,label='Newly closed cycles')
    pred=[r['predicted_width'] for r in data['comparison']];ax.semilogy([7,8],pred,'x',ms=10,mew=2,color=BLUE,label='Frozen width forecasts')
    ax.set(title='A | Symmetric cores keep shrinking',xlabel='Sequence index n  (peaks = 10n + 5)',ylabel='Width in beta');ax.legend(fontsize=9);ax.grid()
    ax=axes[0,1];ax.plot(n[2:],[r['spacing_ratio'] for r in rows[2:]],'o-',color=BLUE,label='Successive fold-spacing ratio')
    ax.plot(n[1:],[r['width_ratio'] for r in rows[1:]],'s-',color=ORANGE,label='Successive core-width ratio')
    ax.axhline(data['expected_spacing_ratio'],color=BLUE,ls='--',label=f'1 / Lambda = {data["expected_spacing_ratio"]:.4f}')
    ax.axhline(data['expected_width_ratio'],color=ORANGE,ls='--',label=f'1 / Lambda² = {data["expected_width_ratio"]:.4f}')
    ax.set(title='B | Position and width have different scales',xlabel='Sequence index n',ylabel='Ratio to the previous member',ylim=(.10,.47));ax.legend(fontsize=9);ax.grid()
    ax=axes[1,0];delta=np.diff(period);ax.plot(n[1:],delta,'o-',color=ORANGE)
    ax.axhline(data['expected_period_increment'],color=INK,ls='--',label=f'2 Tp = {data["expected_period_increment"]:.6f}')
    ax.set(title='C | Two extra saddle turns per new member',xlabel='Sequence index n',ylabel='Increase in full period');ax.legend(fontsize=9);ax.grid()
    ax=axes[1,1];coord=np.array([r['quadratic_core_coordinate'] for r in rows]);ax.semilogy(n,abs(coord-1/16),'o-',color=TEAL)
    ax.set(title='D | An extra check of the quadratic return shape',xlabel='Sequence index n',ylabel='Error from predicted normalized position 1/16');ax.grid()
    for ax in axes.flat:ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.supxlabel('The quadratic-return mechanism is supported numerically. Its global tangency and an infinite hierarchy are not proved.',fontsize=10,color=GRAY)
    save(fig,'chua-layer-prediction-scales')


def descriptions():
    orbits=np.load(ROOT/'results/chua-layer-atlas/family-orbits.npz')
    records=json.loads((OUT/'motion-descriptions.json').read_text())['families']
    fig,axes=plt.subplots(3,2,figsize=(13,9),width_ratios=[1,2.4],layout='constrained')
    fig.suptitle('The same region route can carry different rhythms',fontsize=20,weight='bold')
    for i,id in enumerate(['S08','S11','E05']):
        r=next(r for r in records if r['id']==id);u=orbits[id+'_xyz'];t=orbits[id+'_t']
        ax=axes[i,0];ax.plot(u[:,0],u[:,2],color=INK,lw=1);ax.set(xlim=(-2.5,2.5),ylim=(-3.3,3.3),xlabel='x',ylabel='z',aspect='equal')
        ax.set_title(f'{id}: {r["peaks"]} peaks, T = {r["period"]:.3f}',loc='left',fontsize=12)
        ax=axes[i,1]
        ax.fill_between(t,-2.6,2.6,where=u[:,0]>1,color=ORANGE,alpha=.20,linewidth=0)
        ax.fill_between(t,-2.6,2.6,where=u[:,0]<-1,color=BLUE,alpha=.20,linewidth=0)
        ax.plot(t,u[:,0],color=INK,lw=1.2);ax.axhline(1,color=GRAY,ls=':',lw=.7);ax.axhline(-1,color=GRAY,ls=':',lw=.7)
        ax.set(xlim=(0,19),ylim=(-2.6,2.6),xlabel='Physical time within one complete cycle',ylabel='x')
        ax.set_title('Outer visits: R → R → L → L   |   Residence times distinguish the motion',loc='left',fontsize=10)
        ax.grid(axis='x')
    fig.supxlabel('Blue shading: x < -1. Orange: x > 1. Each reference uses its own beta; these are not three coexisting cycles.\n32 catalogue entries collapse to 12 full outer-region words, or just 4 primitive words. That is information loss, not fewer families.',fontsize=10,color=GRAY)
    save(fig,'chua-layer-motion-description')


def new_motions(data):
    arrays=np.load(OUT/'predicted-orbits.npz')
    fig,axes=plt.subplots(2,2,figsize=(14,7.5),width_ratios=[1,2.7],layout='constrained')
    fig.suptitle('The two motions found after the forecast was frozen',fontsize=20,weight='bold')
    for i,r in enumerate(data['comparison']):
        p=str(r['peaks']);u=arrays[p+'_xyz'];t=arrays[p+'_t'];half=(len(t)-1)//2
        ax=axes[i,0];ax.plot(u[:half+1,0],u[:half+1,2],color=BLUE,lw=.8);ax.plot(u[half:,0],u[half:,2],color=ORANGE,lw=.8)
        ax.set(xlim=(-2.5,2.5),ylim=(-3.3,3.3),xlabel='x',ylabel='z',aspect='equal')
        ax.set_title(f'{p} peaks | {r["peaks"]//5} section returns',loc='left',fontsize=12)
        ax=axes[i,1];ax.plot(t[:half+1],u[:half+1,0],color=BLUE,lw=.8);ax.plot(t[half:],u[half:,0],color=ORANGE,lw=.8)
        ax.set(xlim=(0,190),ylim=(-2.5,2.5),xlabel='Physical time',ylabel='x');ax.grid(axis='x')
        ax.set_title(f'One full cycle, T = {t[-1]:.6f}   |   Repeated-looking loops do not imply a shorter period',loc='left',fontsize=10)
    fig.supxlabel('Blue/orange mark inversion-related halves. The sampled paths illustrate independently checked cycles at their own parameters.\nThe complete return sequence is primitive; the overlapping orbit outlines alone conceal the difference.',fontsize=10,color=GRAY)
    save(fig,'chua-layer-predicted-motions')


def main():
    style();data=json.loads((OUT/'analysis.json').read_text());forecasts(data);scaling(data);descriptions();new_motions(data);organizer(data)
    print('Rendered five PNG/PDF figures.')


def organizer(data):
    report=json.loads((OUT/'organizer-approximants.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
    fig.suptitle('A candidate mechanism: linger near a saddle, leave, and return',fontsize=18,weight='bold')
    ax=axes[0];a=report['records'][0]['starting_coordinate']
    for curve,color in zip(report['curves'],[BLUE,INK,ORANGE]):
        # Magnify the vicinity where the moving parabola meets the stable line.
        x=(np.array(curve['initial_coordinate'])-a)*1e6
        y=np.array(curve['normal_distance'])*1e7
        ax.plot(x,y,color=color,lw=1.8,label=f'Beta offset {curve["beta_offset"]:+.0e}')
    ax.axhline(0,color=GRAY,ls=':',label='Tangent stable line')
    ax.set(xlim=(-4,4),ylim=(-16,4),xlabel='Displacement on initial unstable line (millionths)',
        ylabel='Final normal distance (units of 1e-7)',title='A | The returning branch develops a quadratic fold')
    ax.legend(fontsize=9,loc='lower center');ax.grid()
    ax=axes[1]
    for row,color in zip(data['rows'][-2:],[BLUE,ORANGE]):
        d=row['motion']['distances_to_parent'];ax.semilogy(range(len(d)),d,'o-',color=color,label=f'{row["peaks"]} peaks')
    ax.set(xlabel='Successive inverted half-section returns',ylabel='Distance to the short saddle on the section',
        title='B | The long cycle spends extra returns near that saddle')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True));ax.legend();ax.grid()
    fig.supxlabel(f'Finite-iterate tangency estimates cluster near beta = {report["approximate_beta"]:.11f}.\nSix approximations span {report["beta_spread"]:.2e}; this spread is not an error bound or a proof of an exact manifold tangency.',fontsize=10,color=GRAY)
    save(fig,'chua-layer-prediction-organizer')


if __name__=='__main__':main()
