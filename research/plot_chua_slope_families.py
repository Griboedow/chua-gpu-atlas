"""Show continued family identities, coexistence, and capture-map sensitivity."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
from continue_chua_slope_families import dump

COLORS = {'P5':'#197ba5', 'P7':'#b94d81'}


def read(source, name): return json.loads((source/name).read_text())


def sheet_figure(source, dpi=130):
    fig = plt.figure(figsize=(13.8,7.2), dpi=dpi, facecolor='white')
    fig.subplots_adjust(left=.025,right=.96,bottom=.17,top=.82,wspace=.02)
    axes=[]; overlaps=read(source,'overlaps.json')['records']
    for col, kind in enumerate(['alpha','beta'],1):
        ax=fig.add_subplot(1,2,col,projection='3d',computed_zorder=False); axes.append(ax); all_values=[]
        for name in ['P5','P7']:
            data=read(source,name+'-'+kind+'.json'); rows={tuple(r['index']):r for r in data['nodes']}
            polygons=[]
            for i in range(len(data['controls'])-1):
                for j in range(len(data['m1'])-1):
                    corners=[(i,j),(i+1,j),(i+1,j+1),(i,j+1)]
                    if not all(k in rows for k in corners):continue
                    v=[[rows[k]['representative'][p] for p in ['m0','m1',kind]] for k in corners]
                    polygons.extend([[v[0],v[1],v[2]],[v[0],v[2],v[3]]])
            ax.add_collection3d(Poly3DCollection(polygons,facecolors=COLORS[name],alpha=.72,linewidths=.12,edgecolors='white',zorder=1))
            values=np.array([[r['representative'][p] for p in ['m0','m1',kind]] for r in data['nodes']])
            all_values.append(values)
        values=np.vstack(all_values)
        ax.auto_scale_xyz(values[:,0],values[:,1],values[:,2])
        points=np.array([[r['parameters'][p] for p in ['m0','m1',kind]] for r in overlaps if r['kind']==kind])
        ax.scatter(points[:,0],points[:,1],points[:,2],s=24,c='#132b39',edgecolors='white',linewidths=.45,marker='D',depthshade=False,zorder=10)
        ax.set(xlabel='m0',ylabel='m1',zlabel=kind,
               title='Vary alpha; beta = 100/7' if kind=='alpha' else 'Vary beta; alpha = 9')
        ax.set_box_aspect((1,1,1));ax.set_proj_type('ortho');ax.view_init(elev=24,azim=-55,roll=0)
        ax.tick_params(labelsize=8,pad=1)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.set_major_locator(MaxNLocator(5))
            axis.set_pane_color((.96,.975,.985,1));axis._axinfo['grid']['color']=(.65,.7,.75,.25)
    fig.suptitle('Two continued rhythms can coexist',fontsize=22,weight='bold',y=.97)
    fig.text(.5,.89,'Blue: the five-peak family. Rose: the seven-peak family. Diamonds: independently checked overlap witnesses.',ha='center',fontsize=11)
    fig.text(.5,.075,'Surfaces connect midpoints of checked symmetric stability intervals; their tiny thickness is not exaggerated.\n'
             'Missing cells are left open. A continuation failure does not establish the end of a physical family.',ha='center',fontsize=10,color='#526673')
    return fig, axes


def coexistence(source, assets):
    w=read(source,'coexistence-witness.json'); traces=np.load(source/'coexisting-traces.npz')
    fig=plt.figure(figsize=(13.5,8.2)); grid=fig.add_gridspec(2,2,left=.07,right=.96,bottom=.12,top=.82,hspace=.38,wspace=.24)
    ax=fig.add_subplot(grid[:,0],projection='3d')
    for name in ['P5','P7']:
        v=traces[name+'_state']; ax.plot(*v,color=COLORS[name],lw=1.25,label=name+' cycle')
    ax.set(xlabel='x',ylabel='y',zlabel='z',title='Both closed orbits at identical parameters')
    ax.view_init(elev=24,azim=-55);ax.legend(loc='upper left',frameon=False)
    ax=fig.add_subplot(grid[0,1])
    for r in w['cycles']:
        name=r['id'];t=traces[name+'_time'];x=traces[name+'_state'][0]
        whole=np.linspace(0,40,4001); signal=np.interp(whole%r['period'],t,x)
        ax.plot(whole,signal,color=COLORS[name],lw=1.1,label=f"{name}: T = {r['period']:.4f}")
    ax.set(xlabel='time',ylabel='x',title='Different repeating rhythms');ax.legend(frameon=True,facecolor='white',framealpha=.95,fontsize=9,loc='upper right')
    ax.spines[['top','right']].set_visible(False)
    ax=fig.add_subplot(grid[1,1])
    edges=[]
    for n,name in enumerate(['P5','P7']):
        d=read(source,name+'-alpha.json');r=next(r for r in d['nodes'] if r['index']==w['index'])
        lo,hi=r['symmetry']['m0'],r['fold']['m0'];edges.extend([lo,hi])
        ax.plot([lo,hi],[1-n]*2,lw=14,color=COLORS[name],solid_capstyle='butt')
    ax.axvspan(*w['interval_m0'],color='#294356',alpha=.12,label='overlap')
    ax.axvline(w['parameters']['m0'],color='#132b39',ls='--',lw=1.2,label='two checked cycles')
    ax.set(yticks=[0,1],yticklabels=['7 peaks','5 peaks'],ylim=(-.65,1.65),xlabel='m0',title='Their stable symmetric intervals overlap')
    margin=(max(edges)-min(edges))*.12;ax.set_xlim(min(edges)-margin,max(edges)+margin)
    ax.ticklabel_format(axis='x',style='plain',useOffset=False);ax.tick_params(axis='x',labelsize=8)
    ax.spines[['top','right','left']].set_visible(False);ax.legend(frameon=False,fontsize=9,loc='lower left')
    p=w['parameters'];fig.suptitle('One parameter set, two attracting cycles',fontsize=22,weight='bold',y=.97)
    fig.text(.5,.90,f"alpha = {p['alpha']:g}   beta = 100/7   m0 = {p['m0']:.10f}   m1 = -5/7",ha='center',fontsize=12)
    fig.text(.5,.035,'Independent full-field closure and transverse multipliers confirm both numerical witnesses.\n'
        'The three original atlas starts all select P7 here, at both tested integration steps.',ha='center',fontsize=11,color='#526673')
    fig.savefig(assets/'chua-slope-coexistence.png',dpi=160);plt.close(fig)


def basins(source, assets):
    d=np.load(source/'basin.npz');r=np.load(source/'basin-refinement.npz');t=np.load(source/'basin-time-only.npz')
    labels=d['labels'][:-3].reshape(len(d['z']),len(d['y'])); values=np.select([labels==5,labels==7],[1,2],default=0)
    fig,axes=plt.subplots(1,2,figsize=(13,5.9));fig.subplots_adjust(left=.08,right=.96,bottom=.28,top=.79,wspace=.28)
    axes[0].pcolormesh(d['y'],d['z'],values,cmap=ListedColormap(['#d7dce1',COLORS['P5'],COLORS['P7']]),norm=BoundaryNorm([-.5,.5,1.5,2.5],3),shading='nearest')
    axes[0].set(xlabel='initial y (x = 0)',ylabel='initial z',title='Finite-time capture map: sensitive to step size')
    old=d['labels'][r['indices']];resolved=old!=0
    fractions=[np.mean(old[resolved]==t['labels'][resolved]),np.mean(old[resolved]==r['labels'][resolved])]
    axes[1].barh([1,0],fractions,color=['#519a7a','#c68849'],height=.5)
    axes[1].set(yticks=[1,0],yticklabels=['Longer run,\nsame step','Longer run,\nhalf step'],xlim=(0,1.13),xlabel='Fraction retaining the same cycle label',title=f'{resolved.sum()} previously resolved sampled starts')
    for y,v in zip([1,0],fractions):axes[1].text(v+.02,y,f'{100*v:.1f}%',va='center')
    axes[1].spines[['top','right']].set_visible(False)
    fig.suptitle('Coexistence is secure; detailed capture boundaries are not',fontsize=19,weight='bold')
    fig.legend(handles=[Patch(color=COLORS['P5'],label='Matched P5'),Patch(color=COLORS['P7'],label='Matched P7'),Patch(color='#d7dce1',label='Unresolved / other')],loc='lower center',bbox_to_anchor=(.5,.10),ncol=3,frameon=False)
    fig.text(.5,.035,'Stratified sensitivity sample, not whole-map accuracy. The speckled image does not establish a fractal basin boundary.',ha='center',fontsize=10,color='#526673')
    fig.savefig(assets/'chua-slope-basin-sensitivity.png',dpi=160);plt.close(fig)


def itinerary(source, assets):
    rows=read(source,'itinerary-change.json')['records']
    fig,axes=plt.subplots(1,2,figsize=(12,4.8));fig.subplots_adjust(left=.08,right=.96,bottom=.23,top=.77,wspace=.3)
    x=[r['cycle']['m1'] for r in rows]
    axes[0].plot(x,[r['spectral_radius'] for r in rows],'o-',color=COLORS['P7'],ms=4)
    axes[0].set(xlabel='m1',ylabel='Dominant transverse multiplier magnitude',title='The seven-peak cycle remains attracting',ylim=(0,1.05))
    axes[0].axhline(1,color='#9b6470',ls='--',lw=.8)
    counts=[r['cycle']['breakpoint_crossings'] for r in rows];change=int(np.flatnonzero(np.diff(counts))[0])
    axes[1].plot(x[:change+1],counts[:change+1],color=COLORS['P7'])
    axes[1].plot(x[change+1:],counts[change+1:],color=COLORS['P7'])
    axes[1].axvspan(x[change],x[change+1],color='#d7dce1',alpha=.6,label='Transition bracket')
    axes[1].scatter(x,counts,color=COLORS['P7'],s=22);axes[1].legend(frameon=False,fontsize=9,loc='upper left')
    axes[1].set(xlabel='m1',ylabel='Crossings of x = -1 or x = +1 per period',yticks=[8,12],title='Its route through diode regions changes')
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.ticklabel_format(axis='x',useOffset=False)
    fig.suptitle('A missing sheet edge need not mean a cycle disappears',fontsize=19,weight='bold')
    fig.text(.5,.07,'All 13 witnesses are closed with the actual piecewise vector field and independently checked.\n'
        'Alpha = 9 and beta = 100/7; m0 is adjusted along the path. The precise transition is not located.',ha='center',fontsize=10,color='#526673')
    fig.savefig(assets/'chua-slope-itinerary-change.png',dpi=160);plt.close(fig)


def animate(source, assets):
    fig,axes=sheet_figure(source,dpi=85);frames=[];palette=None
    for i in range(192):
        for ax in axes:ax.view_init(elev=24,azim=-55+20*np.sin(2*np.pi*i/192),roll=0)
        fig.canvas.draw();rgb=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy())
        if palette is None:palette=rgb.quantize(colors=256,method=Image.Quantize.MEDIANCUT,dither=Image.Dither.NONE)
        frames.append(rgb.quantize(palette=palette,dither=Image.Dither.NONE))
    durations=[10*(round((i+1)*12.5)-round(i*12.5)) for i in range(192)]
    path=assets/'chua-slope-family-sheets.gif'
    frames[0].save(path,save_all=True,append_images=frames[1:],duration=durations,loop=0,disposal=1,optimize=True);plt.close(fig)
    with Image.open(path) as gif:
        total=0
        for i in range(gif.n_frames):gif.seek(i);total+=gif.info['duration']
        assert total==24000
        dump(source/'animation.json',dict(frames=gif.n_frames,duration_seconds=total/1000,dimensions=list(gif.size),sweep_degrees=40,elevation=24,roll=0,changes_data=False))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--assets',type=Path,required=True)
    p.add_argument('--stage',choices=['static','animation','all'],default='all');a=p.parse_args();a.assets.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'text.color':'#193449','axes.labelcolor':'#193449'})
    if a.stage in ['static','all']:
        fig,_=sheet_figure(a.source);fig.savefig(a.assets/'chua-slope-family-sheets.png',dpi=160);plt.close(fig)
        coexistence(a.source,a.assets);basins(a.source,a.assets);itinerary(a.source,a.assets)
    if a.stage in ['animation','all']:animate(a.source,a.assets)


if __name__=='__main__':main()
