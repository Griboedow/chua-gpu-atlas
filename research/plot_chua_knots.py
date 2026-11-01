"""Scientific figures and an offline viewer for the numerical knot atlas."""
import json,os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'.publication/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial.transform import Rotation
from plotly.offline import get_plotlyjs
from chua_knot_diagrams import intersections
from reconstruct_chua_knot_orbits import rows

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/chua-knot-atlas';ASSET=ROOT/'docs/assets'
TYPES=['unknot','T(2,5)','T(3,5) mirror','T(2,3) # T(2,3) mirror','T(2,5) # T(2,5) mirror']
NAMES=['Unknot','T(2,5)','Mirror of T(3,5)','Square knot','T(2,5) # mirror']
COLORS=['#287f9c','#c66526','#7748a1','#bd466d','#398466']
SHORT=['Unknot','T(2,5)','T(3,5) mirror','Square knot','T(2,5) + mirror']


def load():
    records=[]
    for source in rows():
        id=source['id'];root=json.loads((OUT/f'{id}-d70-n192.json').read_text())
        d=json.loads((OUT/'diagrams'/f'{id}-d60-n96-s0.json').read_text())
        records.append(dict(id=id,peaks=source.get('peak_count',int(id[1:]) if id.startswith('P') else None),
                            beta=float(root['beta']),period=float(root['period']),symmetric=root['symmetric'],
                            type=d['knot_type'],group=TYPES.index(d['knot_type']),diagram=d))
    return sorted(records,key=lambda x:x['id'])


def portrait(ax,xyz,color):
    x=(xyz-xyz.mean(0))/np.ptp(xyz,axis=0)
    ax.plot(*x.T,color=color,lw=.65,alpha=.85)
    ax.set_box_aspect((1,1,1));ax.view_init(22,-62);ax.set_axis_off()
    ax.set(xlim=(-.57,.57),ylim=(-.57,.57),zlim=(-.57,.57))


def knot_drawing(ax,xyz,color):
    # Select a generic near-planar view. This changes only the camera.
    x=(xyz-xyz.mean(0))/np.ptp(xyz,axis=0)
    candidates=[]
    for a,b in ((.05,.07),(.11,.09),(.17,-.12),(-.07,.06),(.03,-.04)):
        p=x@Rotation.from_euler('xy',[a,b]).as_matrix()
        crossings=intersections(p)
        candidates.append((len(crossings),p,crossings))
    _,p,rr=min(candidates,key=lambda v:v[0])
    edges=np.roll(p,-1,axis=0)-p;length=np.linalg.norm(edges[:,:2],axis=1)
    arc=np.r_[0,np.cumsum(length)];mid=(arc[:-1]+arc[1:])/2;total=arc[-1]
    visible=np.ones(len(p),bool);gap=.035
    for r in rr:
        j=int(r[1] if r[4]>0 else r[0]);f=r[3] if r[4]>0 else r[2]
        position=arc[j]+f*length[j]
        delta=np.abs(mid-position);delta=np.minimum(delta,total-delta)
        visible &= delta>gap
    segments=np.stack((p[:,:2],np.roll(p,-1,axis=0)[:,:2]),axis=1)
    ax.add_collection(LineCollection(segments[visible],colors=color,linewidths=2.4,capstyle='round'))
    ax.set_xlim(p[:,0].min()-.08,p[:,0].max()+.08);ax.set_ylim(p[:,1].min()-.08,p[:,1].max()+.08)
    ax.set_aspect('equal');ax.axis('off')


def save(fig,name):
    fig.savefig(ASSET/f'chua-knots-{name}.png',dpi=190,facecolor=fig.get_facecolor())
    fig.savefig(ASSET/f'chua-knots-{name}.pdf',facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'figure.facecolor':'#fafaf7','axes.facecolor':'#fafaf7'})
    records=load();lookup={r['id']:r for r in records}
    standard=np.load(OUT/'standard-knot-curves.npz')
    controls=json.loads((OUT/'controls.json').read_text())
    shapes={name:standard[f'k{i}'] for i,name in enumerate(controls['types'])}
    examples=['P85','A04','A05','S01','E01']
    fig=plt.figure(figsize=(15,7.6));grid=fig.add_gridspec(2,5,left=.035,right=.985,bottom=.12,top=.81,hspace=.28,wspace=.12)
    fig.text(.04,.955,'Five knot types among 34 reference cycles',fontsize=24,weight='bold',color='#182d39')
    fig.text(.04,.905,'Numerical classification of attracting periodic orbits · α = 9, m₀ = −8/7, m₁ = −5/7 · β varies by family',fontsize=11)
    for j,id in enumerate(examples):
        count=sum(r['group']==j for r in records);r=lookup[id]
        ax=fig.add_subplot(grid[0,j],projection='3d');portrait(ax,np.load(OUT/f'{id}-d60-n96.npz')['xyz'],COLORS[j])
        ax.set_title(f'{id} · {r["peaks"]} peaks\n{count} reference'+(' cycles' if count!=1 else ' cycle'),fontsize=11,pad=0)
        ax=fig.add_subplot(grid[1,j]);knot_drawing(ax,shapes[TYPES[j]],COLORS[j]);ax.set_title(NAMES[j],fontsize=12,color=COLORS[j],weight='bold',pad=10)
    fig.text(.04,.075,'TOP  Actual Chua orbits; coordinates rescaled independently.    BOTTOM  Standard embeddings of the identified knots.',fontsize=10)
    fig.text(.04,.035,'The 15 → 25 → … → 85 peak sequence is unknotted in all eight tests. More peaks do not imply a more complicated knot.',fontsize=11,weight='bold')
    save(fig,'overview')

    fig=plt.figure(figsize=(15,17));grid=fig.add_gridspec(6,6,left=.02,right=.98,bottom=.055,top=.925,wspace=.06,hspace=.25)
    fig.text(.035,.975,'The reference atlas, coloured by knot type',fontsize=23,weight='bold',color='#182d39')
    fig.text(.035,.947,'Each panel shows one full periodic orbit. β differs between panels; these are not 34 simultaneous attractors.',fontsize=11)
    for index,r in enumerate(records):
        ax=fig.add_subplot(grid[index//6,index%6],projection='3d')
        portrait(ax,np.load(OUT/f'{r["id"]}-d60-n96.npz')['xyz'],COLORS[r['group']])
        ax.set_title(f'{r["id"]} · {r["peaks"]} peaks\n{SHORT[r["group"]]}',fontsize=10,pad=0,color=COLORS[r['group']])
    fig.text(.035,.025,'All 34 entries: reconstructed roots, projection checks and exact diagram invariants. Underlying ODE knot types remain numerical.',fontsize=10)
    save(fig,'gallery')

    fig,axes=plt.subplots(1,2,figsize=(14,6.7),gridspec_kw={'width_ratios':[1.2,1]})
    fig.subplots_adjust(left=.065,right=.975,top=.79,bottom=.15,wspace=.24)
    fig.suptitle('Peak count and knot type distinguish different features',fontsize=21,weight='bold',y=.97)
    for ax in axes:
        for j in range(5):
            subset=[r for r in records if r['group']==j]
            ax.scatter([r['beta'] for r in subset],[r['peaks'] for r in subset],s=60,c=COLORS[j],label=SHORT[j],edgecolors='white',linewidths=.5,zorder=3)
        ax.grid(alpha=.17);ax.set_xlabel('β at the reference orbit');ax.set_ylabel('Local maxima of x per full cycle')
    axes[0].set_title('All reference parameters',loc='left');axes[1].set_title('Enlargement near the accumulation',loc='left')
    axes[1].set_xlim(14.302,14.3087);axes[1].ticklabel_format(axis='x',style='plain',useOffset=False)
    axes[1].tick_params(axis='x',rotation=25)
    for id in ['S03','A01','A04','A05','S04','S05','S06','S07','P55','P65','P75','P85']:
        r=lookup[id];dy={'A01':7,'A04':-1,'A05':7,'S05':8}.get(id,4)
        axes[1].annotate(id,(r['beta'],r['peaks']),xytext=(5,dy),textcoords='offset points',fontsize=9,color=COLORS[r['group']])
    axes[0].legend(loc='upper left',bbox_to_anchor=(0,1.2),ncol=3,fontsize=9,frameon=False)
    fig.text(.065,.035,'At 25 peaks: A04 is T(2,5), whereas A01, S03, S05 and S10 are unknots. Equal knot type also does not identify a family.',fontsize=10)
    save(fig,'parameters')

    counter=json.loads((OUT/'sampling-counterexample.json').read_text())
    fig,axes=plt.subplots(1,3,figsize=(13,5.6));fig.subplots_adjust(left=.06,right=.98,bottom=.26,top=.71,wspace=.32)
    fig.suptitle('A coarse trajectory can change the apparent knot',fontsize=22,weight='bold',y=.96)
    fig.text(.06,.85,'A05 · the same 70-digit periodic root · only the placement and number of polygon vertices change',fontsize=11)
    checks=[('Uniform time\n4,000 vertices',counter['uniform'][0]),('Uniform time\n16,000 vertices',counter['uniform'][1]),
            ('Section-aligned\n11,520 vertices',json.loads((OUT/'diagrams/A05-d70-n192-s0.json').read_text()))]
    for ax,(title,r) in zip(axes,checks):
        coef=r['alexander_coefficients'];col='#ac4941' if len(coef)==13 else COLORS[2]
        ax.bar(range(len(coef)),coef,color=col,width=.7);ax.axhline(0,color='#333333',lw=.7)
        ax.set_title(title,fontsize=12);ax.set_xlabel('Power of t');ax.set_ylim(-1.4,1.4);ax.set_yticks([-1,0,1]);ax.set_ylabel('Alexander coefficient')
        ax.text(.5,-.36,'Different polynomial' if len(coef)==13 else 'Same T(3,5) polynomial',transform=ax.transAxes,ha='center',color=col,fontsize=11,weight='bold')
    fig.text(.06,.025,'Agreement of finer grids and projections is a numerical control, not an interval certificate of the ODE knot type.',fontsize=10)
    save(fig,'sampling')

    # Full precision coordinates are retained. Display arrays are also among the
    # classified polygons; an on-screen line width is not a separation bound.
    viewer=[]
    for r in records:
        xyz=np.load(OUT/f'{r["id"]}-d60-n96.npz')['xyz']
        viewer.append({k:v for k,v in r.items() if k!='diagram'}|dict(xyz=xyz.tolist()))
    body='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Chua periodic-orbit knot atlas</title>
<style>body{margin:0;background:#fafaf7;color:#182d39;font:16px system-ui}main{max-width:1400px;margin:auto;padding:26px}h1{font-size:32px;margin:0 0 10px}p{line-height:1.55}select{font:inherit;padding:8px;margin:0 14px 0 8px}.controls{padding:16px 0;border-block:1px solid #ccd5d9}.plots{display:grid;grid-template-columns:1fr 1fr}h2{font-size:20px;margin:18px 0 0} .plot{height:520px}#facts{font-variant-numeric:tabular-nums}footer{border-top:1px solid #ccd5d9;font-size:14px;margin-top:20px}@media(max-width:800px){.plots{grid-template-columns:1fr}.plot{height:420px}}</style>
<main><h1>Knots of attracting periodic orbits</h1><p>34 reference cycles, five numerically identified knot types. Drag either view to rotate it. The right-hand curve is a standard embedding of the knot, not a circuit trajectory.</p>
<div class="controls"><label for="family">Reference family</label><select id="family"></select><button id="reset">Reset views</button></div>
<div class="plots"><section><h2>Chua orbit</h2><div id="orbit" class="plot"></div></section><section><h2 id="kind"></h2><div id="standard" class="plot"></div></section></div>
<p id="facts"></p><p>All eight tested members of the 15 → 25 → … → 85 peak sequence are unknots. A04 and A01 both have 25 peaks but have different knot types.</p>
<footer><p>α = 9, m₀ = −8/7, m₁ = −5/7, γ = 0. Each family uses its own β. Asymmetric representatives have an inversion-related partner of mirror knot type. Both views rescale coordinate axes to expose the shape.</p><p>The classification uses high-precision shooting, section-aligned sampling, multiple projections, Reidemeister simplification and exact Alexander polynomials. It is numerical evidence for the orbit topology; no interval isotopy certificate is claimed. See <a href="https://github.com/Griboedow/chua-gpu-atlas/blob/main/docs/chua-knot-atlas.md">the report and checks</a>.</p></footer></main>
<script>__PLOTLY__</script><script>const data=__DATA__,shapes=__SHAPES__,colors=__COLORS__,names=__NAMES__;
const sel=document.getElementById('family');data.forEach((r,i)=>{let o=document.createElement('option');o.value=i;o.textContent=r.id+' · '+r.peaks+' peaks · '+names[r.group];sel.appendChild(o)});
function trace(xyz,c){return {type:'scatter3d',mode:'lines',x:xyz.map(p=>p[0]),y:xyz.map(p=>p[1]),z:xyz.map(p=>p[2]),line:{color:c,width:4},hoverinfo:'skip'}}
function layout(){return {paper_bgcolor:'#fafaf7',margin:{l:0,r:0,t:15,b:0},scene:{aspectmode:'cube',camera:{eye:{x:1.6,y:1.4,z:.9}},xaxis:{title:'x'},yaxis:{title:'y'},zaxis:{title:'z'}}}}
function draw(){let r=data[+sel.value];Plotly.react('orbit',[trace(r.xyz,colors[r.group])],layout(),{responsive:true,displaylogo:false});let l=layout();l.scene.xaxis.title='';l.scene.yaxis.title='';l.scene.zaxis.title='';Plotly.react('standard',[trace(shapes[r.type],colors[r.group])],l,{responsive:true,displaylogo:false});document.getElementById('kind').textContent=names[r.group];document.getElementById('facts').textContent=r.id+' · β = '+r.beta.toPrecision(16)+' · period = '+r.period.toFixed(8)+' · '+r.peaks+' peaks · '+(r.symmetric?'half-period inversion symmetry':'asymmetric reference');}
sel.value=data.findIndex(r=>r.id==='A04');sel.addEventListener('change',draw);document.getElementById('reset').addEventListener('click',draw);draw();</script></html>'''
    body=body.replace('__PLOTLY__',get_plotlyjs()).replace('__DATA__',json.dumps(viewer,separators=(',',':')))
    body=body.replace('__SHAPES__',json.dumps({k:np.vstack((v,v[0])).tolist() for k,v in shapes.items()},separators=(',',':')))
    body=body.replace('__COLORS__',json.dumps(COLORS)).replace('__NAMES__',json.dumps(NAMES))
    (OUT/'knot-atlas.html').write_text(body,encoding='utf-8')
    print('Wrote four PNG/PDF figures and the offline viewer')


if __name__=='__main__':main()
