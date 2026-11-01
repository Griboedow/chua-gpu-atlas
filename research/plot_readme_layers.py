"""Make a family-by-family README tour and stability-window views at multiple scales."""
import argparse,colorsys,io,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
LOW,HIGH=14.24,14.4


def clip(poly,level,keep_above):
    result=[]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        ia=a[2]>=level if keep_above else a[2]<=level
        ib=b[2]>=level if keep_above else b[2]<=level
        if ia:result.append(a)
        if ia!=ib:result.append(a+(b-a)*(level-a[2])/(b[2]-a[2]))
    return np.array(result)


def geometry(sheet):
    grid=np.full((17,17,3),np.nan);beta=np.full((17,17),np.nan)
    for n in sheet['nodes']:
        i,j=n['index'];r=n['core'];grid[i,j]=[r['alpha'],(r['m0']+8/7)*1000,r['beta']];beta[i,j]=r['beta']
    faces=[]
    for i in range(16):
        for j in range(16):
            for corners in [[(i,j),(i+1,j),(i+1,j+1)],[(i,j),(i+1,j+1),(i,j+1)]]:
                poly=np.array([grid[a,b] for a,b in corners])
                if not np.isfinite(poly).all():continue
                poly=clip(poly,LOW,True)
                if len(poly)<3:continue
                poly=clip(poly,HIGH,False)
                if len(poly)>=3:faces.append(poly)
    points=grid.reshape(-1,3);points=points[np.isfinite(points).all(axis=1)&(points[:,2]>=LOW)&(points[:,2]<=HIGH)]
    return dict(faces=faces,points=points,beta=beta,offsets=(np.array(sheet['m0'])+8/7)*1000)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--assets',type=Path,required=True)
    args=p.parse_args();source=args.source;assets=args.assets;assets.mkdir(parents=True,exist_ok=True)
    registry=json.loads((ROOT/'results/chua-layer-atlas/family-registry.json').read_text())['families']
    ordered=sorted(registry,key=lambda r:r['core']['beta']);ids=sorted(r['id'] for r in registry)
    colors={id:colorsys.hsv_to_rgb((k*.61803398875)%1,.7,.72) for k,id in enumerate(ids)}
    sheets={id:json.loads((source/(id+'.json')).read_text()) for id in ids};geo={id:geometry(sheets[id]) for id in ids}
    allpoints=np.vstack([geo[id]['points'][::3] for id in ids]);audit=json.loads((source/'audit.json').read_text())
    assert audit['status']=='passed'
    assert all(g['faces'] and len(g['points']) for g in geo.values())
    plt.rcParams.update({'font.size':10,'text.color':'#193449','axes.labelcolor':'#193449','axes.titleweight':'normal'})

    def frame(family=None,index=0):
        fig=plt.figure(figsize=(12,7.7),dpi=100,facecolor='white')
        ax=fig.add_axes([.01,.13,.61,.70],projection='3d');right=fig.add_axes([.71,.22,.255,.56])
        if family is None:
            for id in ids:
                ax.add_collection3d(Poly3DCollection(geo[id]['faces'],facecolors=[(*colors[id],.30)],edgecolors=[(*colors[id],.55)],linewidths=.16))
            title='Overview: the 32 catalogued families';sub='The animation next isolates every family at the same scale.'
        else:
            id=family['id'];g=geo[id];color=colors[id]
            ax.scatter(*allpoints.T,s=.8,c='#cbd2d8',alpha=.12,depthshade=False)
            ax.add_collection3d(Poly3DCollection(g['faces'],facecolors=[(*color,.84)],edgecolors=[(*color,.9)],linewidths=.16))
            ax.scatter(*g['points'].T,s=5,c=[color],depthshade=False)
            title=f'{index:02d} / 32   {id}   |   {family["peak_count"]}-peak reference cycle'
            sub='One family highlighted; other families remain as faint context.'
            fig.text(.71,.12,f'Reference period T = {family["core"]["period"]:.4f}\n{len(g["points"])} checked samples inside the box',fontsize=10,color='#193449')
        ax.set(xlim=(8.94,9.06),ylim=(-6,6),zlim=(LOW,HIGH),xlabel='alpha',ylabel='1000 (m0 + 8/7)',zlabel='beta')
        ax.set_xticks([8.94,9,9.06]);ax.set_yticks([-4,0,4]);ax.set_zticks([14.24,14.28,14.32,14.36,14.4]);ax.tick_params(labelsize=9,pad=0)
        ax.view_init(elev=23,azim=-55);ax.set_box_aspect([1,1,1.15]);ax.grid(True)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_pane_color((.985,.99,.995,1))
        for id in ids:
            active=family is None or id==family['id'];g=geo[id]
            right.plot(g['offsets'],g['beta'][8],color=colors[id] if active else '#dce1e5',lw=1.1 if family is None else 2.8 if active else .65,
                marker='.' if active else None,ms=3.5,zorder=3 if active else 1)
        right.set(xlim=(-6,6),ylim=(LOW,HIGH),xlabel='1000 (m0 + 8/7)',ylabel='beta',title='Slice at alpha = 9')
        right.grid(alpha=.15);right.spines[['top','right']].set_visible(False)
        fig.text(.045,.94,'The layers, resolved by family',fontsize=21,weight='bold')
        fig.text(.045,.885,title,fontsize=16,color='#193449');fig.text(.045,.85,sub,fontsize=11,color='#526673')
        fig.text(.045,.045,'Dots are closed attracting cycles. Facets guide the eye between samples; their visible thickness is not window width.\n'
            'Fixed coordinates throughout. The separate zoom shows the widths of the narrowest stability intervals.',fontsize=10,color='#526673')
        buf=io.BytesIO();fig.savefig(buf,format='png',dpi=100);plt.close(fig);buf.seek(0)
        return Image.open(buf).convert('RGB')

    overview=frame();overview.save(assets/'chua-readme-layers-poster.png')
    frames=[overview.quantize(colors=192)]
    for index,family in enumerate(ordered,1):
        frames.append(frame(family,index).quantize(colors=192));print('frame',index,family['id'],flush=True)
    frames.append(frames[0].copy())
    frames[0].save(assets/'chua-readme-layers.gif',save_all=True,append_images=frames[1:],duration=[2200]+[1100]*32+[1600],loop=0,optimize=False,disposal=2)

    fig=plt.figure(figsize=(14.4,8.2),layout='constrained');grid=fig.add_gridspec(2,2,width_ratios=[1.18,1])
    ax=fig.add_subplot(grid[:,0]);old=np.load(ROOT/'results/chua-layer-census/atlas-comparison.npz');bi,mi=np.where(old['regular'][20])
    ax.scatter((old['m0'][mi]+8/7)*1000,old['beta'][bi],c='#d7dde2',s=18,marker='s',label='Old 41-point grid',zorder=0)
    for id in ids:
        ax.plot(geo[id]['offsets'],geo[id]['beta'][8],color=colors[id],lw=1,marker='.',ms=3.5,
            label=id if id in ['S08','S11','E06','E09','P55','P65'] else None)
    ax.set(xlim=(-6,6),ylim=(LOW,HIGH),xlabel='1000 (m0 + 8/7)',ylabel='beta',title='Family branches over the old sampled mask\nSlice at alpha = 9')
    ax.legend(loc='lower left',ncol=3,fontsize=8,framealpha=.92);ax.grid(alpha=.12)
    adding=json.loads((ROOT/'results/chua-layer-census/census-analysis.json').read_text())['adding_sequence']
    thin_report=json.loads((source/'thin-cores.json').read_text());assert thin_report['status']=='passed'
    thin=thin_report['families'];thin_counts=[len(r['nodes']) for r in thin]
    bands=[dict(id=id,maxima=r['maxima'],low=r['symmetry']['beta'],high=r['fold']['beta']) for id,r in zip(['S04','S05','S06','S07'],adding)]
    bands += [dict(id=r['id'],maxima=r['maxima'],low=r['symmetry_beta'],high=r['fold_beta']) for r in thin]
    ax=fig.add_subplot(grid[0,1])
    for j,b in enumerate(bands):
        ax.plot([b['low'],b['high']],[j,j],lw=9,solid_capstyle='butt',color=colors[b['id']])
        ax.plot((b['low']+b['high'])/2,j,'|',color=colors[b['id']],ms=13,mew=1.5)
    ax.set(yticks=range(6),yticklabels=[f'{b["id"]}: {b["maxima"]} peaks' for b in bands],ylim=(5.7,-.7),xlim=(14.30345,14.3065),
        xlabel='beta (same physical scale for all six rows)',title='Six windows crowd together\nOne categorical row per family')
    ax.ticklabel_format(axis='x',style='plain',useOffset=False);ax.grid(axis='x',alpha=.16)
    sub=grid[1,1].subgridspec(1,2,wspace=.13)
    for k,r in enumerate(thin):
        ax=fig.add_subplot(sub[0,k]);scale=1e7 if r['id']=='P55' else 1e8
        x=[(n['core']['beta']-r['fold_beta'])*scale for n in r['nodes']];rho=[n['core']['spectral_radius'] for n in r['nodes']]
        ax.axvspan(-r['width']*scale,0,color=colors[r['id']],alpha=.10)
        ax.plot(x,rho,'.-',color=colors[r['id']],lw=1.1,ms=4)
        ax.axhline(1,color='#7b8790',lw=.8,ls='--');ax.plot([-r['width']*scale,0],[1,1],'o',mfc='white',mec='#7b8790',ms=4)
        ax.text(.5,.95,'Open circles: estimated boundaries',transform=ax.transAxes,ha='center',va='bottom',fontsize=7,color='#526673')
        ax.set(xlim=(-r['width']*scale*1.08,r['width']*scale*.08),ylim=(-.03,1.1),xlabel=f'(beta - fold beta) x 10^{7 if k==0 else 8}',
            title=f'{r["id"]}: a true width zoom\nwidth = {r["width"]:.3e}')
        if k==0:ax.set_ylabel('Return multiplier magnitude\nBelow 1: attracting')
        else:ax.set_yticklabels([])
        ax.grid(alpha=.12)
    for ax in fig.axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('More points reveal the branches; magnification reveals their width',fontsize=17,weight='bold')
    fig.supxlabel(f'{audit["total_nodes"]:,} closed surface samples, including {audit["new_parameter_locations"]:,} new parameter locations. '
        f'The width zooms contain {thin_counts[0]} and {thin_counts[1]} newly closed cycles.\nThe old beta spacing was 0.004: about 126,000 times the P65 stability width. Connecting lines are visual guides.',fontsize=10)
    fig.savefig(assets/'chua-readme-layer-zoom.png',dpi=170,facecolor='white');fig.savefig(source/'readme-layer-zoom.pdf',facecolor='white');plt.close(fig)
    report=dict(families=32,gif_frames=len(frames),frame_order=[r['id'] for r in ordered],
        surface_samples=audit['total_nodes'],new_parameter_locations=audit['new_parameter_locations'],thin_core_samples=sum(thin_counts),
        face_counts={id:len(geo[id]['faces']) for id in ids},
        scope='Displayed facets interpolate closed representative points; they are not physical window thicknesses or validated interiors. The six-row chart uses a categorical vertical axis. The two detailed zooms retain physical beta widths.')
    (source/'render-summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
