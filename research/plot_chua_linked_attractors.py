"""Joint-topology figures and an offline viewer; components keep shared coordinates."""
import os,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.publication/matplotlib'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from plotly.offline import get_plotlyjs
from chua_link_diagrams import clean,crossings

OUT=ROOT/'results/chua-linked-attractors';ASSETS=ROOT/'docs/assets'
IDS=['F1a','F1b','F2a','F2b']
COLORS=['#197b9c','#d77536','#8b559b','#3a8d70']
BG='#fafaf7';INK='#1d3440';MUTED='#637a85'


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def load():
    curves=[np.load(OUT/f'{id}-d90-n384.npz')['xyz'] for id in IDS]
    controls=np.load(OUT/'standard-link-curves.npz')
    borromean=[controls[f'borromean_{i}'] for i in range(3)]
    roots=[read(OUT/f'{id}-d90-n384.json') for id in IDS]
    return curves,borromean,roots


def save(fig,stem):
    fig.savefig(ASSETS/f'{stem}.png',dpi=180,facecolor=BG)
    fig.savefig(ASSETS/f'{stem}.pdf',facecolor=BG)
    plt.close(fig)


def scene(ax,curves,labels):
    curves=[clean(c) for c in curves]
    p=np.vstack(curves);lo=p.min(0);hi=p.max(0);span=hi-lo
    for setter,l,h,s in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],lo,hi,span):setter(l-.06*s,h+.06*s)
    ax.set_box_aspect((1,1,1));ax.view_init(elev=23,azim=-58)
    # Determine every undercrossing in the same fixed-camera projection.
    # Whole-curve artist order would give incorrect over/under information.
    projected=np.array(proj3d.proj_transform(*p.T,ax.get_proj())).T
    projected[:,2]*=-1  # Larger depth coordinate now means nearer the viewer.
    sizes=np.array([len(c) for c in curves]);starts=np.r_[0,np.cumsum(sizes)[:-1]]
    components=np.repeat(np.arange(len(curves)),sizes)
    local=np.concatenate([np.arange(n) for n in sizes])
    ends=np.concatenate([s+(np.arange(n)+1)%n for s,n in zip(starts,sizes)])
    rr=crossings(projected,ends,components,local,sizes)
    lengths=np.linalg.norm(projected[ends,:2]-projected[:,:2],axis=1)
    intervals=[];midpoints=[]
    for start,n in zip(starts,sizes):
        arc=np.r_[0,np.cumsum(lengths[start:start+n])]
        intervals.append(arc);midpoints.append((arc[:-1]+arc[1:])/2)
    visible=np.ones(len(p),bool)
    gap=.00105 if len(curves)==4 else .00145
    for r in rr:
        side=1 if r[4]>0 else 0;edge=int(r[side]);fraction=r[side+2]
        component=int(components[edge]);start=starts[component];n=sizes[component]
        arc=intervals[component];position=arc[edge-start]+fraction*lengths[edge]
        delta=abs(midpoints[component]-position);delta=np.minimum(delta,arc[-1]-delta)
        visible[start:start+n]&=delta>gap+lengths[start:start+n]/2
    segments=np.stack([p,p[ends]],axis=1)
    depth=(projected[:,2]+projected[ends,2])/2
    selected=np.flatnonzero(visible);selected=selected[np.argsort(depth[selected])]
    colors=np.array(COLORS,dtype=object)[components[selected]]
    ax.add_collection3d(Line3DCollection(segments[selected],colors=colors,linewidths=1.15 if len(curves)==4 else 2.25,capstyle='round'))
    ax.set(xlabel='x',ylabel='y',zlabel='z');ax.tick_params(labelsize=8,pad=0)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
        axis.pane.fill=False;axis._axinfo['grid']['color']=(.72,.78,.80,.3)
    ax.legend([Line2D([],[],color=c,lw=2) for c in COLORS[:len(curves)]],labels,frameon=False,fontsize=9,loc='upper left',bbox_to_anchor=(.02,.94))


def overview(curves,borromean):
    fig=plt.figure(figsize=(14,9.5));gs=fig.add_gridspec(2,2,height_ratios=[3.7,1],left=.045,right=.96,top=.85,bottom=.12,wspace=.13,hspace=.12)
    fig.text(.05,.955,'Four coexisting Chua cycles: no joint link found',fontsize=23,weight='bold',color=INK)
    fig.text(.05,.906,'Actual cycles at one parameter point, compared with a known collective-link control.',fontsize=12,color=MUTED)
    ax=fig.add_subplot(gs[0,0],projection='3d');scene(ax,curves,IDS)
    ax.set_title('Chua: four numerically unlinked components',fontsize=13,pad=15,color=INK)
    ax=fig.add_subplot(gs[0,1],projection='3d');scene(ax,borromean,['Ring A','Ring B','Ring C'])
    ax.set_title('Control: Borromean rings',fontsize=13,pad=15,color=INK)
    for column,rows in enumerate([
        [('Individual cycles','4 / 4 unknots'),('Two-cycle sublinks','6 / 6 unlinks'),('Three-cycle sublinks','4 / 4 unlinks'),('All four together','Unlink')],
        [('Individual rings','3 / 3 unknots'),('Two-ring sublinks','3 / 3 unlinks'),('All three together','Nontrivial link'),('Detection','Exact diagram polynomial')]]):
        ax=fig.add_subplot(gs[1,column]);ax.set_axis_off()
        for j,(label,value) in enumerate(rows):
            y=.95-j*.25;ax.text(.10,y,label,color=MUTED,fontsize=11,va='center')
            ax.text(.60,y,value,color=INK if column==0 or j!=2 else COLORS[1],fontsize=11,weight='bold',va='center')
    fig.text(.05,.065,'Chua: α = 453/50, β = 71729/5000, m₀ = −8/7 + 3/500, m₁ = −5/7, γ = 0.',fontsize=10,color=MUTED)
    fig.text(.05,.030,'Each panel uses shared coordinates; small gaps mark undercrossings in this projection. Projection overlap is not linking.\nClassification concerns reconstructed polygonal curves; ambient isotopy to the exact ODE orbits is not certified.',fontsize=10,color=MUTED)
    save(fig,'chua-linked-attractors-overview')


def evidence(roots):
    files=sorted(OUT.glob('analysis-d*-n*.json'));summaries=[]
    for path in files:
        d=read(path);records=d['records'];digits=records[0]['digits'];n=records[0]['subdivisions']
        summaries.append((digits,n,records,d['gauss_integrals']))
    fig,axs=plt.subplots(1,2,figsize=(14,6.8));fig.subplots_adjust(left=.085,right=.97,bottom=.22,top=.76,wspace=.35)
    fig.text(.06,.94,'The negative result survives independent numerical checks',fontsize=22,weight='bold',color=INK)
    fig.text(.06,.875,'Unlink classifications use diagram reductions; a zero linking number alone would not be sufficient.',fontsize=11,color=MUTED)
    ax=axs[0];x=np.arange(len(summaries));labels=[]
    for j,(digits,n,records,gauss) in enumerate(summaries):
        total=len(records);passed=sum(r['status']=='unlink' for r in records)
        ax.bar(j,passed,color=COLORS[j%3],width=.57)
        ax.text(j,passed+.6,f'{passed} / {total}',ha='center',fontsize=12,weight='bold',color=INK)
        labels.append(f'{digits} digits\n{n} samples / arc\n{len(set(r["seed"] for r in records))} projections')
    ax.set(xticks=x,xticklabels=labels,ylabel='Diagrams reducing to crossing-free unlinks',ylim=(0,max(len(s[2]) for s in summaries)*1.22),title='All 15 nonempty subsets at every projection')
    ax.grid(axis='y',alpha=.17);ax.set_axisbelow(True)
    ax=axs[1]
    y=np.array([float(r['q'][0]) for r in roots]);z=np.array([float(r['q'][1]) for r in roots])
    for j,(yy,zz) in enumerate(zip(y,z)):
        ax.scatter(yy,zz,s=75,color=COLORS[j],zorder=3)
        offset=(10,-14) if j in [1,2] else (10,7)
        ax.annotate(IDS[j],(yy,zz),xytext=offset,textcoords='offset points',color=COLORS[j],weight='bold',fontsize=11)
    ax.set(xlabel='y on the common section x = 0',ylabel='z on the common section x = 0',title='Distinct reference points, shared coordinates')
    ax.set_xlim(y.min()-.005,y.max()+.01);ax.set_ylim(z.min()-.01,z.max()+.01);ax.grid(alpha=.17)
    maxgauss=max(abs(g['value']) for s in summaries for g in s[3])
    fig.text(.06,.10,f'Independent solid-angle Gauss integrals: all six pairs agree with zero; maximum absolute residual {maxgauss:.2g}.',fontsize=11,color=INK)
    fig.text(.06,.045,'Left: repeats across precision, polygon resolution and camera direction. Right: one chosen starting point per orbit.\nThese checks support the numerical classification; they do not provide a certified separation tube around every ODE orbit.',fontsize=10,color=MUTED)
    save(fig,'chua-linked-attractors-evidence')


def viewer(curves,borromean,roots):
    datasets={
        'chua':dict(names=IDS,curves=[np.round(c,12).tolist() for c in curves],title='Four coexisting Chua cycles',
            note='Every tested pair, triple and the four-component union reduces to an unlink. These are numerical polygonal classifications.',
            details='α = 453/50 · β = 71729/5000 · m₀ = −8/7 + 3/500 · m₁ = −5/7 · γ = 0',
            periods=[float(r['period']) for r in roots]),
        'borromean':dict(names=['Ring A','Ring B','Ring C'],curves=[np.vstack([c,c[:1]]).tolist() for c in borromean],
            title='Borromean rings — geometric control',
            note='Every pair is an unlink. All three together form a nontrivial link, detected by an exact diagram polynomial. These curves are not Chua trajectories.',
            details='Known ellipse construction · positive control for collective linking',periods=[])}
    source='''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Joint topology of Chua cycles</title>
<style>
:root{color-scheme:light;font-family:system-ui,-apple-system,Segoe UI,sans-serif;color:#1d3440;background:#f5f7f7}*{box-sizing:border-box}body{margin:0}main{max-width:1450px;margin:0 auto;padding:28px 28px 20px}h1{font-size:clamp(24px,3vw,34px);margin:0 0 10px;letter-spacing:-.025em}p{line-height:1.55}header p{margin:0;color:#637a85;max-width:950px}.layout{display:grid;grid-template-columns:285px minmax(0,1fr);gap:20px;margin-top:25px}.panel,.chart{background:#fff;border:1px solid #dbe3e6;border-radius:12px;overflow:hidden}.panel{padding:22px}.panel h2{font-size:15px;margin:22px 0 11px}.panel h2:first-child{margin-top:0}label{display:block;margin:9px 0;line-height:1.4}select,button{font:inherit;color:inherit;border:1px solid #bacbd2;border-radius:6px;background:#fff;padding:8px 10px}select{width:100%}button{cursor:pointer}button:hover{background:#edf3f5}button:focus-visible,select:focus-visible,input:focus-visible{outline:3px solid #197b9c;outline-offset:2px}.buttons{display:flex;gap:8px;flex-wrap:wrap}.dot{display:inline-block;width:11px;height:11px;border-radius:50%;margin:0 7px}.meta{font-size:12px;color:#637a85;line-height:1.5}#plot{height:650px;width:100%}.charthead{padding:20px 24px 0}.charthead h2{font-size:20px;margin:0 0 5px}.charthead p{font-size:13px;color:#637a85;margin:0}.note{border-top:1px solid #e2e9ec;padding:15px 24px;margin:0;min-height:80px;font-size:14px}.status{display:block;font-size:12px;color:#637a85;margin-top:12px}footer{max-width:1100px;color:#637a85;font-size:13px;padding-top:18px}a{color:#197b9c}details{margin-top:16px}summary{cursor:pointer;font-size:13px}.badges{display:flex;gap:8px;flex-wrap:wrap;margin-top:17px}.badge{font-size:12px;background:#eaf2f3;border-radius:20px;padding:5px 10px}.sr{position:absolute;left:-9999px}@media(max-width:900px){main{padding:20px 12px}.layout{grid-template-columns:1fr}.panel{display:grid;grid-template-columns:1fr 1fr;gap:8px 22px}.panel h2{margin-top:0}#plot{height:570px}}@media(max-width:540px){.panel{display:block}.panel h2{margin-top:20px}#plot{height:460px}.charthead{padding:16px 14px 0}.note{padding:14px}}
</style></head><body><main>
<header><h1>Joint topology of coexisting cycles</h1><p>Inspect all components in their shared coordinates. Hide a component to inspect a sublink. The Borromean control illustrates why individually unknotted curves and zero pairwise linking numbers do not determine the whole link.</p><div class="badges"><span class="badge">Offline viewer</span><span class="badge">Fixed camera until dragged</span><span class="badge">Common coordinates</span></div></header>
<div class="layout"><aside class="panel"><div><h2>Dataset</h2><label class="sr" for="dataset">Choose a dataset</label><select id="dataset"><option value="chua">Chua: four coexisting cycles</option><option value="borromean">Control: Borromean rings</option></select><p class="meta">Switching datasets changes the coordinate range for the entire scene. Components are never individually moved or rescaled.</p></div>
<div><h2>Visible components</h2><div id="components"></div><span id="visibleCount" class="status" aria-live="polite"></span></div>
<div><h2>Sublink presets</h2><label class="sr" for="subset">Choose a component subset</label><select id="subset"></select><p class="meta">The legend also shows component colors. Use the checkboxes to choose any pair, triple or individual cycle.</p></div>
<div><h2>View</h2><div class="buttons"><button id="reset">Reset camera</button><button id="top">View x–y</button><button id="front">View x–z</button></div><details><summary>Interaction and evidence</summary><p class="meta">Drag to rotate; scroll to zoom. No automatic rotation. Plotly's toolbar exports the current view as PNG.<br><br>Chua curves use 90-digit roots and 384 samples per affine section-to-section arc. Diagram simplification and sampling checks support an unlink result. An ambient-isotopy certificate for the exact ODE orbits is still absent.</p></details></div></aside>
<section class="chart" aria-label="Interactive three-dimensional link viewer"><div class="charthead"><h2 id="sceneTitle"></h2><p id="parameters"></p></div><div id="plot"></div><p id="note" class="note"></p></section></div>
<footer>Lines show actual reconstructed Chua coordinates or the explicitly labeled geometric control. Axes use one common scene scale per coordinate, shared by every component. A crossing in a two-dimensional view is not itself a link.<br><a href="https://github.com/Griboedow/chua-gpu-atlas/blob/main/docs/chua-linked-attractors.md">Study report</a> · <a href="https://github.com/Griboedow/chua-gpu-atlas/tree/main/results/chua-linked-attractors">Data and reproduction</a></footer>
</main><script>__PLOTLY__</script><script>
const DATA=__DATA__;const COLORS=__COLORS__;const defaults={eye:{x:1.55,y:1.7,z:1.05},up:{x:0,y:0,z:1}};let current='chua';let selected=[];
const plot=document.getElementById('plot');const componentBox=document.getElementById('components');const subset=document.getElementById('subset');
function updateVisibility(){Plotly.restyle(plot,{visible:selected});document.querySelectorAll('#components input').forEach((e,i)=>e.checked=selected[i]);document.getElementById('visibleCount').textContent=selected.filter(Boolean).length+' of '+selected.length+' components visible';}
function setSelection(indices){selected=selected.map((_,i)=>indices.includes(i));updateVisibility();}
function render(key){current=key;const d=DATA[key];selected=d.names.map(()=>true);document.getElementById('sceneTitle').textContent=d.title;document.getElementById('parameters').textContent=d.details;document.getElementById('note').textContent=d.note;componentBox.replaceChildren();
 d.names.forEach((name,i)=>{const label=document.createElement('label');const input=document.createElement('input');input.type='checkbox';input.checked=true;input.addEventListener('change',()=>{selected[i]=input.checked;subset.value='custom';updateVisibility()});const dot=document.createElement('span');dot.className='dot';dot.style.background=COLORS[i];label.append(input,dot,document.createTextNode(name+(d.periods.length?' · T = '+d.periods[i].toFixed(5):'')));componentBox.append(label)});
 subset.replaceChildren();const options=[['all','All '+d.names.length+' components']];d.names.forEach((name,i)=>options.push(['omit-'+i,'Omit '+name+' ('+(d.names.length-1)+' components)']));options.push(['custom','Custom selection']);options.forEach(([value,text])=>{const o=document.createElement('option');o.value=value;o.textContent=text;subset.append(o)});
 const traces=d.curves.map((c,i)=>({type:'scatter3d',mode:'lines',name:d.names[i],x:c.map(p=>p[0]),y:c.map(p=>p[1]),z:c.map(p=>p[2]),line:{color:COLORS[i],width:5},hovertemplate:d.names[i]+'<br>x %{x:.6f}<br>y %{y:.6f}<br>z %{z:.6f}<extra></extra>'}));
 const axis={showbackground:true,backgroundcolor:'#fafcfb',gridcolor:'#e2e9ec',zerolinecolor:'#cad8de'};Plotly.react(plot,traces,{margin:{l:8,r:8,t:5,b:5},paper_bgcolor:'#fff',font:{family:'system-ui, sans-serif',color:'#1d3440'},showlegend:true,legend:{x:.02,y:.98,bgcolor:'rgba(255,255,255,.8)'},scene:{xaxis:{...axis,title:'x'},yaxis:{...axis,title:'y'},zaxis:{...axis,title:'z'},aspectmode:'cube',camera:defaults},uirevision:key},{responsive:true,displaylogo:false,toImageButtonOptions:{format:'png',filename:'chua-joint-link-view',scale:2}}).then(()=>{plot.removeAllListeners('plotly_legendclick');plot.on('plotly_legendclick',e=>{selected[e.curveNumber]=!selected[e.curveNumber];subset.value='custom';updateVisibility();return false});updateVisibility()});}
document.getElementById('dataset').addEventListener('change',e=>render(e.target.value));subset.addEventListener('change',e=>{if(e.target.value==='all')setSelection(selected.map((_,i)=>i));else if(e.target.value.startsWith('omit-')){const omit=Number(e.target.value.slice(5));setSelection(selected.map((_,i)=>i).filter(i=>i!==omit))}});
document.getElementById('reset').addEventListener('click',()=>Plotly.relayout(plot,{'scene.camera':defaults}));document.getElementById('top').addEventListener('click',()=>Plotly.relayout(plot,{'scene.camera':{eye:{x:0,y:0,z:2.5},up:{x:0,y:1,z:0}}}));document.getElementById('front').addEventListener('click',()=>Plotly.relayout(plot,{'scene.camera':{eye:{x:0,y:-2.5,z:0},up:{x:0,y:0,z:1}}}));render('chua');
</script></body></html>'''
    source=source.replace('__PLOTLY__',get_plotlyjs()).replace('__DATA__',json.dumps(datasets,separators=(',',':'))).replace('__COLORS__',json.dumps(COLORS))
    source='\n'.join(line.rstrip() for line in source.splitlines())+'\n'
    (OUT/'linked-attractors.html').write_text(source,encoding='utf-8',newline='\n')


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.facecolor':BG,'axes.facecolor':BG,
        'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':INK,'text.color':INK,'xtick.color':MUTED,'ytick.color':MUTED})
    curves,borromean,roots=load();overview(curves,borromean);evidence(roots);viewer(curves,borromean,roots)
    print('Wrote two shared-coordinate figures and offline viewer',flush=True)


if __name__=='__main__':main()
