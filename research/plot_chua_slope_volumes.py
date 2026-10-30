"""Render measured 3D slope sections, interior crops, and a gentle camera sweep."""
import argparse,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes
from PIL import Image

KINDS=['alpha','beta']
COLORS=['#d7dce1','#54b795','#e69b3f','#eddc6b','#365f9e','#161616','#a66ba8']


def load(source):
    data={}
    for kind in KINDS:
        with np.load(source/(kind+'.npz')) as f:
            d={k:f[k] for k in ['m0','m1','control','labels']}
        d['labels']=d['labels'].reshape(len(d['control']),len(d['m1']),len(d['m0']))
        data[kind]=d
    return data


def crop(d,interior):
    if not interior:return d
    x=np.flatnonzero(d['m0']<=-1.06);y=np.flatnonzero((d['m1']>=-.90)&(d['m1']<=-.20))
    return dict(m0=d['m0'][x],m1=d['m1'][y],control=d['control'],labels=d['labels'][:,y][:,:,x])


def surface(d,label):
    mask=(d['labels']==label).astype(np.float32)
    if mask.min()==mask.max():return np.empty((0,3)),np.empty((0,3),dtype=int)
    axes=[d['control'],d['m1'],d['m0']]
    verts,faces,_,_=marching_cubes(mask,level=.5,spacing=tuple(float(a[1]-a[0]) for a in axes),allow_degenerate=False)
    verts+=np.array([a[0] for a in axes]);return verts[:,[2,1,0]],faces


def scene(fig,position,d,kind,interior):
    ax=fig.add_subplot(1,2,position,projection='3d');d=crop(d,interior)
    z,y,x=np.where(d['labels']==4);stride=max(1,math.ceil(len(x)/1800))
    ax.scatter(d['m0'][x[::stride]],d['m1'][y[::stride]],d['control'][z[::stride]],s=1.2,c=COLORS[4],alpha=.07,depthshade=False)
    vertices,faces=surface(d,2)
    ax.add_collection3d(Poly3DCollection(vertices[faces],facecolors=COLORS[2],linewidths=0,alpha=.96,shade=True))
    center=9. if kind=='alpha' else 100/7
    xx=[d['m0'][0],d['m0'][-1],d['m0'][-1],d['m0'][0],d['m0'][0]]
    yy=[d['m1'][0],d['m1'][0],d['m1'][-1],d['m1'][-1],d['m1'][0]]
    ax.plot(xx,yy,[center]*5,color='#7e909c',lw=.65,ls='--')
    ax.scatter([-8/7],[-5/7],[center],s=36,marker='*',c='white',edgecolors='#193449',depthshade=False)
    ax.set(xlim=(d['m0'][0],d['m0'][-1]),ylim=(d['m1'][0],d['m1'][-1]),zlim=(d['control'][0],d['control'][-1]),
        xlabel='m0',ylabel='m1',zlabel=kind)
    ax.set_title('Vary alpha; fixed beta = 100/7' if kind=='alpha' else 'Vary beta; fixed alpha = 9',fontsize=13,pad=14)
    ax.set_box_aspect((1,1,1.05));ax.set_proj_type('ortho');ax.view_init(elev=24,azim=-55,roll=0)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
        axis.set_major_locator(MaxNLocator(4));axis.set_pane_color((.965,.975,.985,1));axis._axinfo['grid']['color']=(.65,.7,.75,.25)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'));ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(labelsize=8,pad=1)
    return ax,dict(vertices=len(vertices),faces=len(faces),regular_nodes=int((d['labels']==2).sum()),
        ranges={key:[float(d[key][0]),float(d[key][-1])] for key in ['m0','m1','control']},context_points=len(x[::stride]))


def figure(data,interior,dpi=110):
    fig=plt.figure(figsize=(13.8,7.6),dpi=dpi,facecolor='white');fig.subplots_adjust(left=.01,right=.97,bottom=.16,top=.80,wspace=.03)
    axes=[];stats={}
    for k,kind in enumerate(KINDS,1):
        ax,s=scene(fig,k,data[kind],kind,interior);axes.append(ax);stats[kind]=s
    fig.suptitle('Inside the two slope volumes' if interior else 'Two three-dimensional slope sections',fontsize=22,weight='bold',y=.97)
    fig.text(.5,.90,'Interior crop exposes regular bands hidden by the outer regions.' if interior else 'Every volume contains 152,561 parameter points, with three initial states per point.',ha='center',fontsize=12)
    fig.text(.5,.07,'Orange: boundary of the sampled regular-candidate mask. Faint blue dots: sampled chaotic-candidate context.\n'
        'Dashed frame: the original 2D section. Binary-mask interpolation only; no smoothing or artificial end caps.',ha='center',fontsize=10,color='#526673')
    return fig,axes,stats


def static(data,source,assets):
    summary={}
    for interior in [False,True]:
        fig,_,stats=figure(data,interior);name='interior' if interior else 'overview';summary[name]=stats
        fig.savefig(assets/('chua-slope-volumes'+('-interior' if interior else '')+'.png'),dpi=160,facecolor='white');plt.close(fig)
    from matplotlib.colors import ListedColormap,BoundaryNorm
    from matplotlib.patches import Patch
    fig,axes=plt.subplots(2,3,figsize=(13.5,8.7));fig.subplots_adjust(left=.06,right=.98,top=.88,bottom=.16,hspace=.38,wspace=.25)
    for row,kind in enumerate(KINDS):
        d=data[kind]
        for col,index in enumerate([0,20,40]):
            ax=axes[row,col];ax.pcolormesh(d['m0'],d['m1'],d['labels'][index],cmap=ListedColormap(COLORS),norm=BoundaryNorm(np.arange(-.5,7),7),shading='nearest')
            ax.set(xlabel='m0',ylabel='m1',title=f'{kind} = {d["control"][index]:.5g}');ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Sections through the measured volumes',fontsize=20,weight='bold')
    handles=[Patch(facecolor=COLORS[k],label=name) for k,name in [(4,'Chaos candidate'),(2,'Regular candidate'),(1,'Equilibrium candidate'),(6,'Different labels across starts'),(3,'Uncertain'),(0,'Escape cutoff reached')]]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=3,frameon=False,fontsize=10)
    fig.text(.5,.012,'Top row: beta = 100/7. Bottom row: alpha = 9. The center panels reproduce the earlier slope section on this coarser mesh.',ha='center',fontsize=10)
    fig.savefig(assets/'chua-slope-volume-sections.png',dpi=160,facecolor='white');plt.close(fig)
    (source/'render-summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')


def animate(data,source,assets):
    fig,axes,_=figure(data,True,dpi=90);frames=[];palette=None;count=192;fps=8
    for i in range(count):
        angle=-55+20*np.sin(2*np.pi*i/count)
        for ax in axes:ax.view_init(elev=24,azim=angle,roll=0)
        fig.canvas.draw();rgb=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy())
        if palette is None:palette=rgb.quantize(colors=256,method=Image.Quantize.MEDIANCUT,dither=Image.Dither.NONE)
        frames.append(rgb.quantize(palette=palette,dither=Image.Dither.NONE))
        if (i+1)%24==0:print('animation',i+1,'/',count,flush=True)
    durations=[10*(round((i+1)*100/fps)-round(i*100/fps)) for i in range(count)]
    path=assets/'chua-slope-volumes.gif';frames[0].save(path,save_all=True,append_images=frames[1:],duration=durations,loop=0,disposal=1,optimize=True);plt.close(fig)
    with Image.open(path) as check:
        total=0
        for i in range(check.n_frames):check.seek(i);total+=check.info['duration']
        assert total==24000
        summary=dict(frames=check.n_frames,rendered_frames=count,duration_seconds=total/1000,dimensions=list(check.size),
            sweep_degrees=40,elevation_degrees=24,roll_degrees=0,bytes=path.stat().st_size,changes_data=False)
    (source/'animation.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--assets',type=Path,required=True)
    p.add_argument('--stage',choices=['static','animation','all'],default='all');a=p.parse_args();a.assets.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'text.color':'#193449','axes.labelcolor':'#193449'})
    data=load(a.source)
    if a.stage in ['static','all']:static(data,a.source,a.assets)
    if a.stage in ['animation','all']:animate(data,a.source,a.assets)


if __name__=='__main__':main()
