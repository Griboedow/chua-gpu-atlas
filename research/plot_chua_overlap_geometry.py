"""Publication figures and a self-contained visual guide to coexistence geometry."""
import argparse
import base64
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

ROOT=Path(__file__).resolve().parents[1]
INK='#183142';MUTED='#526875';BLUE='#147e96';ROSE='#bd436c'
PURPLE='#7860b5';GOLD='#e4a23b';GREEN='#158a76'
COLORS=['#f5f7f9','#b6dce5','#edc3d0',PURPLE,GOLD,GREEN]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelcolor':INK,
    'text.color':INK,'axes.edgecolor':'#b1bec8','xtick.color':MUTED,'ytick.color':MUTED,
    'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white',
    'savefig.facecolor':'white','pdf.fonttype':42,'svg.fonttype':'none'})


def read(out,name):return json.loads((out/name).read_text())


def save(fig,out,assets,name):
    fig.savefig(assets/(name+'.png'),dpi=170)
    fig.savefig(out/(name+'.pdf'))
    plt.close(fig)


def curve(out,kind,code):
    rows=read(out,kind+'-'+code+'.json')['records']
    return np.array([[r['m1'],r['control'],r['m0']] for r in rows])


def projection(ax,out,kind='alpha',legend=True):
    curves={code:curve(out,kind,code) for code in ['FS','SF','FD','DF','SD','DS']}
    x=curves['FS'][:,0]
    ax.fill_between(x,curves['FD'][:,1],curves['DF'][:,1],color=GOLD,alpha=.25)
    ax.fill_between(x,curves['FS'][:,1],curves['SF'][:,1],color=PURPLE,alpha=.5)
    for code in ['FD','DF']:ax.plot(x,curves[code][:,1],color='#a96a1f',lw=1.6)
    for code in ['FS','SF']:ax.plot(x,curves[code][:,1],color=PURPLE,lw=1.6)
    for code in ['SD','DS']:ax.plot(x,curves[code][:,1],color=GREEN,lw=1.4,ls='--')
    ax.set(xlabel=r'outer diode slope $m_1$',ylabel=kind)
    ax.grid(alpha=.14,zorder=0)
    if legend:
        ax.legend(handles=[Patch(facecolor=GOLD,alpha=.4,label='Both families; either symmetry'),
            Patch(facecolor=PURPLE,alpha=.6,label='Both cycles symmetric'),
            Line2D([0],[0],color=GREEN,ls='--',label='Two asymmetric pairs possible')],
            loc='upper center',fontsize=8.5,frameon=True,facecolor='white',framealpha=.97)


def slice_arrays(data):
    x=np.array([r['m1'] for r in data['records']])
    edges={name:{role:np.array([r['edges'][name][role]['m0'] for r in data['records']])
                 for role in ['fold','symmetry','doubling']} for name in ['P5','P7']}
    center=(edges['P5']['fold']+edges['P5']['symmetry'])/2
    return x,edges,center


def phase(ax,data,labels=True):
    x,edges,center=slice_arrays(data)
    fine=np.linspace(x[0],x[-1],801)
    edges={name:{role:np.interp(fine,x,values) for role,values in e.items()} for name,e in edges.items()}
    center=np.interp(fine,x,center);x=fine
    height=np.linspace(-1.65,1.05,601)
    m0=center[None,:]+height[:,None]/1000
    count={}
    for name,e in edges.items():
        count[name]=((m0>e['symmetry'])&(m0<e['fold'])).astype(int)+2*((m0>e['doubling'])&(m0<e['symmetry'])).astype(int)
    a,b=count['P5'],count['P7']
    categories=np.select([(a>0)&(b==0),(b>0)&(a==0),(a==1)&(b==1),(a+b==3)&(a>0)&(b>0),(a==2)&(b==2)],
                         [1,2,3,4,5],default=0)
    ax.pcolormesh(x,height,categories,cmap=ListedColormap(COLORS),norm=BoundaryNorm(np.arange(-.5,6),6),shading='nearest',rasterized=True)
    for name,color in [('P5',BLUE),('P7',ROSE)]:
        for role,ls in [('fold','-'),('symmetry','--'),('doubling',':')]:
            ax.plot(x,1000*(edges[name][role]-center),color=color,ls=ls,lw=1)
    ax.set(xlim=(x[0],x[-1]),ylim=(height[0],height[-1]),xlabel=r'$m_1$',
        ylabel=r'$10^3(m_0-C_5)$' if labels else '')
    ax.set_title(r'$\alpha = $'+f'{data["alpha"]:g}',loc='left')
    ax.xaxis.set_major_locator(MaxNLocator(5));ax.yaxis.set_major_locator(MaxNLocator(5))
    return categories


def phase_legend(fig,y=.1):
    labels=['P5 only','P7 only','2: one symmetric cycle each','3: one pair + one cycle','4: two asymmetric pairs']
    fig.legend(handles=[Patch(facecolor=c,label=t) for c,t in zip(COLORS[1:],labels)],
        loc='lower center',bbox_to_anchor=(.5,y),ncol=3,frameon=False,fontsize=9)


def story(out,assets):
    fig=plt.figure(figsize=(14.5,9.3))
    gs=fig.add_gridspec(2,2,left=.075,right=.97,bottom=.20,top=.83,wspace=.26,hspace=.36,width_ratios=[1,1.05])
    ax=fig.add_subplot(gs[:,0]);projection(ax,out)
    ax.set(xlim=(-.785,-.718),ylim=(8.839,8.881),title='A | Both families attract at some common m0')
    for value in [8.85,8.86]:
        ax.axhline(value,color=INK,lw=1,ls=(0,(3,4)))
        ax.text(-.784,value+.00065,f'alpha = {value:g}',fontsize=9,weight='bold')
    tips=read(out,'turning-points.json')['records']
    for code in ['FS','SF']:
        t=next(r for r in tips if r['kind']=='alpha' and r['code']==code)
        ax.scatter(t['m1'],t['control'],c=PURPLE,s=30,edgecolors='white',zorder=6)
    ax=fig.add_subplot(gs[0,1]);phase(ax,read(out,'slice-8.85.json'))
    ax.set_title('B | alpha = 8.85: one symmetric-overlap island',loc='left',fontsize=10)
    ax=fig.add_subplot(gs[1,1]);phase(ax,read(out,'slice-8.86.json'))
    ax.set_title('C | alpha = 8.86: two purple islands, an orange bridge',loc='left',fontsize=10)
    ax.annotate('Three attracting cycles\nconnect the symmetric islands',xy=(-.7521,-.52),xytext=(-.765,-1.36),
        fontsize=9,color='#8c591a',arrowprops=dict(arrowstyle='->',color='#8c591a',lw=1.3),
        bbox=dict(boxstyle='round,pad=.3',facecolor='white',edgecolor='none',alpha=.95))
    fig.suptitle('Symmetry changes the apparent layer count',fontsize=23,weight='bold',y=.96)
    fig.text(.5,.899,'A numerical atlas of two Chua families: five and seven peaks per full period.  beta = 100/7.',ha='center',fontsize=11)
    phase_legend(fig,.095)
    fig.text(.5,.035,'Right: subtracting C5 = (fold5 + symmetry5)/2 straightens P5; vertical interval widths are preserved.\n'
        'Counts include these primary cycles before their first period doubling. White does not mean chaos or absence of other attractors.',
        ha='center',fontsize=9,color=MUTED)
    save(fig,out,assets,'chua-overlap-story')


def all_slices(out,assets):
    fig,axes=plt.subplots(2,2,figsize=(13.5,9.5));fig.subplots_adjust(left=.085,right=.97,bottom=.20,top=.83,hspace=.33,wspace=.23)
    subtitles=['No primary overlap','One connected overlap interval',
               'Symmetric islands, connected primary overlap','The primary overlap also splits']
    for ax,alpha,title in zip(axes.ravel(),[8.84,8.85,8.86,8.865],subtitles):
        phase(ax,read(out,f'slice-{alpha:g}.json'))
        ax.set_title(f'alpha = {alpha:g} | {title}',loc='left',fontsize=10)
    fig.suptitle('Four slices through the same family geometry',fontsize=22,weight='bold',y=.96)
    fig.text(.5,.90,'Continuous boundary tracking separates a change of section geometry from a change of orbit symmetry.',ha='center',fontsize=11)
    phase_legend(fig,.095)
    fig.text(.5,.035,'Vertical coordinate: 1000(m0 - C5), with C5 the midpoint of the P5 symmetric interval.  beta = 100/7.\n'
        'Solid = fold; dashed = symmetry loss; dotted = first period doubling. Filled cells interpolate the 101 corrected cross-sections.',ha='center',fontsize=9,color=MUTED)
    save(fig,out,assets,'chua-overlap-slices')
    # Individual views for the offline guide use the same checked data.
    for alpha in [8.84,8.85,8.86,8.865]:
        fig,ax=plt.subplots(figsize=(10,5.3));fig.subplots_adjust(left=.11,right=.97,bottom=.19,top=.88)
        phase(ax,read(out,f'slice-{alpha:g}.json'))
        fig.text(.5,.035,'P5: blue outlines. P7: rose outlines. Purple = 2, orange = 3, teal = 4 cycles from both families.',ha='center',fontsize=9)
        svg=out/f'slice-{alpha:g}.svg'
        fig.savefig(svg);plt.close(fig)
        svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')


def orbits(out,assets):
    witnesses=read(out,'witnesses.json')['records'];data=np.load(out/'witness-traces.npz')
    fig=plt.figure(figsize=(15,9.4));gs=fig.add_gridspec(2,3,left=.055,right=.98,bottom=.16,top=.81,hspace=.25,wspace=.22,height_ratios=[1.4,1])
    titles={2:'Two symmetric cycles',3:'One symmetric cycle + one pair',4:'Two asymmetric pairs'}
    for col,row in enumerate(witnesses):
        ax=fig.add_subplot(gs[0,col],projection='3d');wave=fig.add_subplot(gs[1,col])
        for cycle in row['cycles']:
            key=cycle['trace_key'];v=data[key+'_state'];t=data[key+'_time'];mirror=cycle['branch']=='asymmetric B'
            color=BLUE if cycle['id']=='P5' else ROSE
            ls='--' if mirror else '-'
            ax.plot(*v,color=color,lw=1.1,ls=ls,alpha=.85)
            total=np.linspace(0,40,3201)
            wave.plot(total,np.interp(total%cycle['period'],t,v[0]),color=color,ls=ls,lw=1,
                label=cycle['id']+(' B' if mirror else (' A' if not cycle['symmetric'] else '')))
        ax.set(xlabel='x',ylabel='y',zlabel='z');ax.set_box_aspect((1,1,1));ax.set_proj_type('ortho');ax.view_init(22,-58)
        ax.set_title(f'{row["count"]} attracting cycles\n'+titles[row['count']],fontsize=12,pad=15)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.set_major_locator(MaxNLocator(4));axis.set_pane_color((.97,.98,.99,1))
        ax.tick_params(labelsize=8,pad=0)
        wave.set(xlabel='time',ylabel='x',xlim=(0,40),ylim=(-3.4,3.4));wave.grid(alpha=.12)
        wave.legend(ncol=2,fontsize=8,frameon=False,loc='upper right')
        p=row['parameters'];wave.set_title(f'alpha = {p["alpha"]:g}; m0 = {p["m0"]:.10f}',fontsize=10)
    fig.suptitle('The layers represent different attracting motions',fontsize=22,weight='bold',y=.965)
    fig.text(.5,.89,'Blue: five peaks. Rose: seven peaks. Dashed traces are the inversion partners of solid asymmetric traces.',ha='center',fontsize=11)
    fig.text(.5,.065,'Within each column all orbits use identical parameters.  beta = 100/7; m1 = -0.7521.\n'
        'Every displayed orbit passes independent closure and stability checks. Counts are lower bounds; other families are not excluded.',ha='center',fontsize=10,color=MUTED)
    save(fig,out,assets,'chua-overlap-orbits')


def sections(out,assets):
    witnesses=read(out,'witnesses.json')['records']
    fig,axes=plt.subplots(2,3,figsize=(14,8));fig.subplots_adjust(left=.095,right=.97,bottom=.19,top=.80,hspace=.40,wspace=.31)
    for col,witness in enumerate(witnesses):
        for row,(family,color) in enumerate([('P5',BLUE),('P7',ROSE)]):
            ax=axes[row,col];cycles=[c for c in witness['cycles'] if c['id']==family]
            points=np.array([c['q'][:2] for c in cycles]);center=points.mean(axis=0)
            offsets=1000*(points-center)
            if len(points)==2:ax.plot(offsets[:,0],offsets[:,1],c='#aebbc4',lw=1,zorder=1)
            for i,(point,cycle) in enumerate(zip(offsets,cycles)):
                ax.scatter(*point,s=85,c=color,marker='o' if i==0 else '^',edgecolors='white',lw=.7,zorder=3)
                ax.annotate('symmetric' if cycle['symmetric'] else ('A' if i==0 else 'B'),point,xytext=(8,5),textcoords='offset points',fontsize=10)
            ax.axvline(0,c='#e4e9ee',lw=.7,zorder=0);ax.axhline(0,c='#e4e9ee',lw=.7,zorder=0)
            ax.set(xlim=(-1.7,1.7),ylim=(-2.8,2.8),xlabel=r'$10^3(y-y_c)$',ylabel=r'$10^3(z-z_c)$',
                title=f'{family}: {len(cycles)} '+('cycle' if len(cycles)==1 else 'distinct cycles'))
            ax.text(.5,-.29,f'center (yc, zc) = ({center[0]:.7f}, {center[1]:.7f})',ha='center',transform=ax.transAxes,fontsize=8,color=MUTED)
            if row==0:ax.text(.5,1.24,f'{witness["count"]} checked cycles at this point',ha='center',transform=ax.transAxes,fontsize=12,weight='bold')
    fig.suptitle('A section makes the distinct cycles visible',fontsize=22,weight='bold',y=.97)
    fig.text(.5,.902,'Recorded upward crossings of x = 0. Each panel is centered on its own family; all panels use the same scale.',ha='center',fontsize=11)
    fig.text(.5,.055,'The nearly overlapping 3D curves separate into distinct points on the section. A and B are related by state inversion and a phase shift.\n'
        'Centers are printed so the original section coordinates can be recovered. These are orbit points, not basin samples.',ha='center',fontsize=10,color=MUTED)
    save(fig,out,assets,'chua-overlap-sections')


def prediction(out,assets):
    tips=read(out,'turning-points.json')['records'];pred=read(out,'predictions.json')['records']
    tip=next(t for t in tips if t['kind']=='alpha' and t['code']=='FS')
    rows=[r for r in pred if r['kind']=='alpha' and r['code']=='FS']
    curvature=tip['curvatures'][-1]['second_derivative']
    fig,axes=plt.subplots(1,3,figsize=(15,5.8));fig.subplots_adjust(left=.065,right=.97,bottom=.23,top=.76,wspace=.32)
    ax=axes[0];x=np.linspace(-.006,.006,301)
    ax.plot(x,curvature*x*x/2,color=PURPLE,lw=2,label='Tip-curvature prediction')
    for row in rows:
        xx=[r['m1']-tip['m1'] for r in row['roots']]
        ax.plot(xx,[row['delta_control']]*2,color='#bac4cc',lw=.8)
        ax.scatter(xx,[row['delta_control']]*2,color=INK,s=22,zorder=3)
    ax.set(xlabel=r'$m_1-m_1^*$',ylabel=r'$\alpha-\alpha^*$',title='A | A horizontal cut meets the bend twice')
    ax.legend(loc='upper center',fontsize=8,frameon=False);ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax=axes[1];delta=np.array([r['delta_control'] for r in rows]);sep=np.array([r['corrected_separation'] for r in rows])
    xx=np.geomspace(delta.min(),delta.max(),100)
    ax.loglog(xx,2*np.sqrt(2*xx/curvature),color=PURPLE,lw=2,label='Prediction: square root')
    ax.scatter(delta,sep,color=INK,s=30,zorder=3,label='Corrected cut roots')
    ax.set(xlabel='distance from tip in alpha',ylabel='separation in m1',title='B | Prediction checked away from the tip')
    ax.legend(fontsize=8,frameon=False)
    ax=axes[2]
    for kind,color in [('alpha',BLUE),('beta',ROSE)]:
        for code,ls in zip(['FS','FF','SS','SF'],['-','--',':','-.']):
            rs=[r for r in pred if r['kind']==kind and r['code']==code]
            ax.semilogx([r['delta_control'] for r in rs],[1e6*r['relative_separation_error'] for r in rs],
                color=color,ls=ls,marker='o',ms=3,label=kind+' '+code)
    ax.set(xlabel='distance from tip in alpha or beta',ylabel='relative separation error (parts per million)',title='C | All eight symmetric-edge contact curves')
    ax.legend(ncol=2,fontsize=7.5,frameon=False,loc='upper left')
    for ax in axes:ax.grid(alpha=.15)
    fig.suptitle('The geometry predicts where the two islands will appear',fontsize=21,weight='bold',y=.96)
    fig.text(.5,.855,r'$\Delta m_1 \simeq 2\sqrt{2\delta/|a^{\prime\prime}(m_1^*)|}$'+'     |     curvature is measured at the tip, not fitted to the cut roots',ha='center',fontsize=12)
    fig.text(.5,.07,'32 held-out horizontal cuts; maximum measured separation error below 0.004%. F = fold; S = symmetry loss.\n'
        'This is a standard local geometric square-root law, verified here for these contact curves; the estimates are not interval-certified.',ha='center',fontsize=10,color=MUTED)
    save(fig,out,assets,'chua-overlap-prediction')


def volumes(out,assets):
    fig=plt.figure(figsize=(14,7.6));axes=[]
    for k,kind in enumerate(['alpha','beta']):
        ax=fig.add_subplot(1,2,k+1,projection='3d');axes.append(ax)
        for codes,color in [(['FD','DF'],GOLD),(['FS','SF'],PURPLE)]:
            a,b=[curve(out,kind,code) for code in codes]
            # A ruled ribbon joins contact loci; it is a geometric guide,
            # not an orbit boundary or a certified coexistence volume.
            u=np.linspace(0,1,13)
            m0=a[:,2,None]*(1-u)+b[:,2,None]*u
            m1=np.repeat(a[:,0,None],len(u),axis=1)
            control=a[:,1,None]*(1-u)+b[:,1,None]*u
            ax.plot_surface(m0,m1,control,color=color,alpha=.22,linewidth=0,shade=False)
            for code,r in zip(codes,[a,b]):ax.plot(r[:,2],r[:,0],r[:,1],color=color,lw=1.8,label=code)
        ax.set(xlabel='m0',ylabel='m1',zlabel=kind,title='beta = 100/7' if kind=='alpha' else 'alpha = 9')
        ax.set_box_aspect((1,1,1));ax.set_proj_type('ortho');ax.view_init(24,-52)
        ax.legend(fontsize=9,frameon=False,loc='upper left')
        ax.tick_params(labelsize=8,pad=1)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.set_major_locator(MaxNLocator(4));axis.set_pane_color((.97,.98,.99,1))
    fig.suptitle('Contact curves in the original three-dimensional coordinates',fontsize=21,weight='bold',y=.96)
    fig.text(.5,.875,'Purple: edges of simultaneous symmetric stability. Gold: edges including asymmetric primary branches.',ha='center',fontsize=11)
    fig.subplots_adjust(left=.02,right=.96,bottom=.16,top=.81,wspace=.1)
    fig.text(.5,.045,'The plotted lines are solved simultaneous events on two distinct cycles. Translucent ruled ribbons only connect those lines visually.\n'
        'F = fold; S = symmetry loss; D = first period doubling of an asymmetric cycle. The ribbons are not measured stability volumes.',ha='center',fontsize=10,color=MUTED)
    save(fig,out,assets,'chua-overlap-volumes')


def guide(out,assets):
    slices={str(alpha):(out/f'slice-{alpha:g}.svg').read_text(encoding='utf-8') for alpha in [8.84,8.85,8.86,8.865]}
    slices={key:value[value.index('<svg'):] for key,value in slices.items()}
    descriptions={
        '8.84':'No overlap of these two primary families in this slice. Other attractors have not been excluded.',
        '8.85':'One purple region supports both symmetric cycles. Orange margins include an asymmetric pair.',
        '8.86':'The two purple islands are joined by an orange region: one symmetric cycle and an asymmetric pair, for at least three attractors.',
        '8.865':'Even after including the asymmetric primary branches, this slice has two coexistence intervals. Period-doubled descendants are outside the count.'}
    images={name:base64.b64encode((assets/(name+'.png')).read_bytes()).decode() for name in ['chua-overlap-story','chua-overlap-orbits','chua-overlap-sections','chua-overlap-prediction','chua-overlap-volumes']}
    data=json.dumps(dict(slices=slices,descriptions=descriptions)).replace('</',r'<\/')
    html='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Chua: how layers overlap</title>
<style>body{font:17px/1.55 system-ui,sans-serif;color:#183142;background:#f5f7fa;margin:0}main{max-width:1120px;margin:auto;padding:38px 24px}h1{font-size:clamp(30px,5vw,48px);line-height:1.1}h2{margin-top:48px}p{max-width:900px}figure,.panel{margin:28px 0;background:white;padding:20px;border-radius:14px;box-shadow:0 2px 12px #1831420c}img,svg{max-width:100%;height:auto}select{font:inherit;padding:8px 14px;margin:0 14px;background:white;border:1px solid #9aaeba;border-radius:6px}.note{color:#526875;font-size:15px}.key{display:flex;gap:15px;flex-wrap:wrap;font-size:14px}.chip{display:inline-block;width:15px;height:15px;margin-right:5px;vertical-align:-2px}a{color:#147e96}#view svg{width:100%}</style>
<main><p class="note">REPRODUCIBLE NUMERICAL RESEARCH • CHUA CIRCUIT</p><h1>One family is not always one attracting motion.</h1>
<p>Two tracked families, with five and seven peaks, overlap. Each family can contain either a symmetric cycle or a pair of asymmetric cycles. This changes both the number of attractors and what a two-dimensional slice looks like.</p>
<figure><img alt="Overview of coexistence geometry and its symmetry-dependent slices" src="data:image/png;base64,STORY"></figure>
<h2>Move the section</h2><p>Choose a directly calculated alpha section. The other fixed parameter is beta = 100/7.</p>
<div class="panel"><label for="alpha">alpha</label><select id="alpha"><option>8.84</option><option>8.85</option><option selected>8.86</option><option>8.865</option></select><p id="description" aria-live="polite"></p><div id="view"></div>
<div class="key"><span><i class="chip" style="background:#b6dce5"></i>P5 only</span><span><i class="chip" style="background:#edc3d0"></i>P7 only</span><span><i class="chip" style="background:#7860b5"></i>2 cycles</span><span><i class="chip" style="background:#e4a23b"></i>3 cycles</span><span><i class="chip" style="background:#158a76"></i>4 cycles</span></div>
<p class="note">The vertical coordinate is 1000(m0 - C5), where C5 is the midpoint of the P5 symmetric interval at each m1. Subtraction straightens the band without changing its thickness. Colored cells interpolate corrected boundaries; white does not imply chaos. “P5 only” or “P7 only” can include an asymmetric pair.</p></div>
<h2>See the actual motions</h2><p>Within each column, all displayed trajectories use the same parameters. Every orbit has an independent closure and stability check.</p><figure><img alt="Independently verified examples with two, three, and four attracting cycles" src="data:image/png;base64,ORBITS"></figure><p>The mirrored cycles almost coincide in the broad orbit view. Their recorded section crossings separate them clearly:</p><figure><img alt="Magnified section points distinguish all nine checked cycles" src="data:image/png;base64,SECTIONS"></figure>
<h2>Make a prediction, then check it</h2><p>A smooth bend predicts a square-root separation of the two cut intersections. We used the curvature at the tip to predict 32 held-out cuts, then solved the intersections separately.</p><figure><img alt="Tip-curvature predictions compared with separately corrected intersections" src="data:image/png;base64,PREDICTION"></figure>
<h2>The original three-dimensional coordinates</h2><figure><img alt="Contact curves in m0, m1, alpha and m0, m1, beta" src="data:image/png;base64,VOLUMES"></figure>
<p class="note">Scope: two selected primary families through their first period doubling. These are numerical results, not an exhaustive attractor count or a proof of scientific priority. The square-root law and the underlying bifurcation mechanisms are established mathematics.</p>
<p><a href="https://github.com/Griboedow/chua-gpu-atlas/blob/main/docs/chua-overlap-geometry.md">Read the research report, equations, literature comparison and reproduction instructions</a></p>
<script>const data=DATA;const chooser=document.getElementById('alpha');function render(){document.getElementById('view').innerHTML=data.slices[chooser.value];document.getElementById('description').textContent=data.descriptions[chooser.value];}chooser.addEventListener('change',render);render();</script></main></html>'''
    for key,name in [('STORY','chua-overlap-story'),('ORBITS','chua-overlap-orbits'),('SECTIONS','chua-overlap-sections'),('PREDICTION','chua-overlap-prediction'),('VOLUMES','chua-overlap-volumes')]:html=html.replace(key,images[name])
    html=html.replace('DATA',data)
    (out/'guide.html').write_text(html,encoding='utf-8',newline='\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=ROOT/'results/chua-overlap-geometry')
    p.add_argument('--assets',type=Path,default=ROOT/'docs/assets')
    args=p.parse_args();args.assets.mkdir(parents=True,exist_ok=True)
    for action in [story,all_slices,orbits,sections,prediction,volumes,guide]:
        action(args.source,args.assets);print(action.__name__,flush=True)


if __name__=='__main__':main()
