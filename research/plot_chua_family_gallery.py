"""Make 32 individual orbit cards and an offline two-family visual comparison.

Uses the published reference orbits without rerunning discovery or changing
catalogue entries. The time series and peak order distinguish nearly identical
projected outlines. Every card uses the same state-coordinate scales.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.signal import find_peaks
from plotly.offline import get_plotlyjs

ROOT=Path(__file__).resolve().parents[1]
BLUE='#176f93'; ORANGE='#d56a36'; INK='#193449'; GRAY='#647580'


def portrait(ax, xyz, ticks=True):
    half=(len(xyz)-1)//2
    ax.plot(xyz[:half+1,0],xyz[:half+1,2],color=BLUE,lw=.85)
    ax.plot(xyz[half:,0],xyz[half:,2],color=ORANGE,lw=.85)
    ax.scatter(xyz[0,0],xyz[0,2],s=15,color=INK,zorder=5)
    ax.axhline(0,color='#dbe2e6',lw=.6,zorder=0);ax.axvline(0,color='#dbe2e6',lw=.6,zorder=0)
    ax.set(xlim=(-2.5,2.5),ylim=(-3.45,3.45),aspect='equal')
    if ticks:ax.set(xlabel='x',ylabel='z')
    else:ax.set(xticks=[],yticks=[])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/chua-layer-atlas')
    p.add_argument('--assets',type=Path,default=ROOT/'docs/assets/chua-families')
    args=p.parse_args();out=args.output;assets=args.assets;assets.mkdir(parents=True,exist_ok=True)
    records=sorted(json.loads((out/'family-behavior.json').read_text()),key=lambda r:r['reference_beta'])
    orbits=np.load(out/'family-orbits.npz');payload={};checks=[]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
        'text.color':INK,'axes.labelcolor':INK,'xtick.color':GRAY,'ytick.color':GRAY})
    for r in records:
        id=r['id'];xyz=orbits[id+'_xyz'];t=orbits[id+'_t'];peaks,_=find_peaks(xyz[:,0]);half=(len(t)-1)//2
        assert np.isfinite(xyz).all() and np.isfinite(t).all()
        assert np.all(np.diff(t)>0) and abs(t[-1]-r['period'])<1e-8
        assert len(peaks)==r['maxima'],id+' sampled peaks disagree with the independently counted maxima'
        closure=float(max(abs(xyz[-1]-xyz[0])));assert closure<1e-4
        checks.append(dict(id=id,samples=len(t),sampled_peaks=len(peaks),expected_peaks=r['maxima'],endpoint_gap=closure))
        symmetry='half-period inversion symmetry' if r['symmetric'] else 'asymmetric representative'
        visit='both outer regions' if r['both_outer_regions'] else 'one outer region'
        fig=plt.figure(figsize=(14,4.8),layout='constrained');grid=fig.add_gridspec(1,3,width_ratios=[.85,1.65,1.15])
        ax=fig.add_subplot(grid[0]);portrait(ax,xyz);ax.set_title('Orbit outline (x-z projection)',fontsize=11)
        ax=fig.add_subplot(grid[1]);ax.plot(t[:half+1],xyz[:half+1,0],color=BLUE,lw=1.05,label='First half')
        ax.plot(t[half:],xyz[half:,0],color=ORANGE,lw=1.05,label='Second half')
        ax.set(xlim=(0,t[-1]),ylim=(-2.5,2.5),xlabel='Time within one complete cycle',ylabel='x',title='The complete repeating signal')
        ax.grid(alpha=.13)
        ax=fig.add_subplot(grid[2]);index=np.arange(1,len(peaks)+1)
        ax.plot(index,xyz[peaks,0],color='#a6b2ba',lw=.6,zorder=1)
        ax.scatter(index,xyz[peaks,0],c=np.where(t[peaks]<t[-1]/2,BLUE,ORANGE),s=21,zorder=2)
        ax.set(xlim=(.5,len(peaks)+.5),ylim=(-2.5,2.5),xlabel='Peak number in this cycle',ylabel='Height of x at a local maximum',title=f'{len(peaks)} peaks before repetition')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True,nbins=7))
        ax.grid(alpha=.13)
        fig.suptitle(f'{id}  |  {r["maxima"]} peaks  |  T = {r["period"]:.6f}  |  {symmetry}',fontsize=16,weight='bold')
        fig.supxlabel(f'alpha = 9, beta = {r["reference_beta"]:.11f}, m0 = -8/7, m1 = -5/7  |  {visit}\n'
            'One reference orbit per family, at different beta values. Blue/orange mark the two halves; projected lines may overlap.',fontsize=9)
        fig.savefig(assets/(id+'.png'),dpi=150,facecolor='white');plt.close(fig)
        # Two-to-one decimation is only for browser display. Full arrays and
        # independently counted peaks remain available in the source evidence.
        payload[id]=dict(**r,xyz=np.round(xyz[::2],7).tolist(),peak_values=np.round(xyz[peaks,0],7).tolist())
    for page in range(2):
        fig,axes=plt.subplots(4,4,figsize=(11.8,13.6),layout='constrained')
        for ax,r in zip(axes.flat,records[page*16:(page+1)*16]):
            portrait(ax,orbits[r['id']+'_xyz'],False)
            label='symmetric' if r['symmetric'] else 'asymmetric'
            ax.set_title(f'{r["id"]}  |  {r["maxima"]} peaks\nT = {r["period"]:.3f}  |  {label}',fontsize=10)
        fig.suptitle(f'Chua family portraits  {page+1}/2',fontsize=21,weight='bold')
        fig.supxlabel('Same x-z scales in every panel. Blue = first half of the cycle; orange = second half.\n'
            'An outline can hide repeated passes. Individual cards show the full signal and ordered peaks.',fontsize=10)
        fig.savefig(assets/f'overview-{page+1}.png',dpi=165,facecolor='white')
        fig.savefig(out/f'family-portraits-{page+1}.pdf',facecolor='white');plt.close(fig)
    html=HTML.replace('__PLOTLY__',get_plotlyjs()).replace('__DATA__',json.dumps(payload,separators=(',',':')))
    html='\n'.join(line.rstrip() for line in html.splitlines())+'\n'
    (out/'family-gallery.html').write_text(html,encoding='utf-8',newline='\n')
    report=dict(families=len(records),checks=checks,passed=True,
        maximum_plot_endpoint_gap=max(r['endpoint_gap'] for r in checks),
        note='Display-array checks, not new existence or stability certificates. Small endpoint gaps are integration errors below the plotting scale.',
        source_sha256={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['family-behavior.json','family-orbits.npz']})
    (out/'gallery-checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    lines=['# Individual Chua family portraits','',
        'The figures compare one numerically checked reference orbit from each of 32 catalogued families. A family is defined by periodic-orbit continuation and recorded branch relationships, rather than by peak count or orbit outline alone.',
        '', '## How to read the differences','',
        '- **Where it moves:** E06 stays near one outer scroll; the other reference representatives visit both outer regions.',
        '- **How long until it repeats:** the full period T ranges from about 9.57 to 140.53. A longer cycle can retrace almost the same visible loops many times.',
        '- **What repeats:** the full signal and the ordered heights of its local maxima reveal patterns hidden by the orbit outline. Two five-peak families, S08 and S11, have different geometry and periods of about 10.89 and 14.50.',
        '- **Whether the halves mirror:** for a symmetric representative, the second half equals the first with all state coordinates reversed. Twenty-three reference representatives have this property; nine are asymmetric. Symmetry can change along a family\'s descendant branches.',
        '', 'Continuation, branch relationships and parameter location distinguish families whose orbit outlines are visually similar.',
        '', 'The [knot atlas](chua-knot-atlas.md) adds spatial topology for these 32 representatives and the later P75/P85 cycles. It compares each orbit with a standard knot embedding and explains why dense, section-aligned samples are needed to distinguish close strands. Knot type supplements the recurrence patterns shown here; it does not identify a family by itself.',
        '', '[Offline interactive comparison](../results/chua-layer-atlas/family-gallery.html): download and open the HTML in a browser. Select two families, rotate their 3D paths and compare signals. A phase slider shifts the second signal; the saved starting crossing need not be the same position within a repeated multi-loop pattern.',
        '', 'The interactive time axis is the fraction of a full cycle, so equal panel width does not mean equal physical period. The peak-order panel uses the saved starting phase and does not shift with the signal slider. The orange/blue on the static cards mark cycle halves; in the interactive comparison they identify the selected families.',
        '', 'All cards fix alpha = 9, m0 = -8/7, m1 = -5/7 and use the family-specific reference beta shown. The E-family reference roots can lie outside the original beta box; their continued families enter that box elsewhere. These are comparisons at different parameters, not a claim of 32 simultaneous attractors. Asymmetric cards show one partner. The normalized transverse multiplier is chosen by convention, so this gallery does not rank stability or basin size.',
        '', '![First sixteen family outlines](assets/chua-families/overview-1.png)',
        '', '![Second sixteen family outlines](assets/chua-families/overview-2.png)',
        '', '## All 32 cards','', '| Family | Peaks | Full period T | Reference symmetry |', '|---|---:|---:|---|']
    for r in records:lines.append(f'| [{r["id"]}](#{r["id"].lower()}) | {r["maxima"]} | {r["period"]:.6f} | {"Symmetric" if r["symmetric"] else "Asymmetric"} |')
    for r in records:
        lines.extend(['',f'## {r["id"]}','',f'![{r["id"]}: orbit, complete signal and peak sequence](assets/chua-families/{r["id"]}.png)'])
    lines.extend(['','Source: [family behaviour data](../results/chua-layer-atlas/family-behavior.json), [sampled orbits](../results/chua-layer-atlas/family-orbits.npz), and [display checks](../results/chua-layer-atlas/gallery-checks.json). Scientific interpretation and numerical limitations are in the [refined atlas report](chua-refined-atlas.md).','',
        'Regenerate with `python research/plot_chua_family_gallery.py`. The PNGs and PDFs use the full stored arrays; the interactive paths use every second sample for display.',''])
    (ROOT/'docs/chua-family-gallery.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
    print('Rendered 32 cards, two overview PNG/PDF sheets and the offline comparison; all sampled peak counts agree.')


HTML='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Compare Chua window families</title><style>
body{font:16px system-ui;margin:24px;color:#193449;background:#fff}h1{margin:0 0 10px}p{max-width:1150px;line-height:1.5}.controls{display:flex;flex-wrap:wrap;gap:24px;align-items:center;padding:16px;background:#f0f5f7;border-radius:9px}select{font:inherit;padding:5px}#facts{display:flex;gap:24px;margin:14px 0;flex-wrap:wrap}.fact{padding:12px 16px;border-left:4px solid;line-height:1.6;background:#f8fafb;min-width:280px}.note{font-size:14px;color:#526674}#plot{height:880px}label{white-space:nowrap}
</style><script>__PLOTLY__</script></head><body>
<h1>Compare the 32 Chua window families</h1>
<p>A family has a repeating motion pattern. Its outline can look almost identical to another family's, while the sequence or number of loops before repetition differs. Compare the <b>3D path</b>, the <b>complete signal</b>, and its <b>ordered peaks</b>.</p>
<div class="controls"><label>Blue family <select id="first"></select></label><label>Orange family <select id="second"></select></label>
<label>Shift orange signal <input id="phase" type="range" min="0" max="1" step="0.001" value="0"><output id="shift">0%</output></label></div>
<div id="facts"></div><div id="plot"></div>
<p class="note">One reference representative per family, at alpha = 9, m0 = -8/7 and m1 = -5/7. Beta differs between families. These cycles are not all coexisting at one parameter point. E-family reference betas may lie outside the original atlas box.</p>
<p class="note">The time axis is a fraction of a full period: longer cycles are compressed to the same width. Saved starting phases can differ; the slider shifts only the orange signal. The peak sequence retains its saved starting phase. A 2D projection can overlap itself; rotate the 3D view. Stability multipliers were normalized for identity comparison and do not rank basin size.</p>
<script>
const D=__DATA__, first=document.querySelector('#first'), second=document.querySelector('#second'), slider=document.querySelector('#phase');
const ids=Object.keys(D).sort((a,b)=>D[a].reference_beta-D[b].reference_beta);
for(const id of ids){for(const el of [first,second]){const o=document.createElement('option');o.value=id;o.textContent=id+' | '+D[id].maxima+' peaks';el.appendChild(o)}}
const none=document.createElement('option');none.value='';none.textContent='None';second.prepend(none);first.value='S08';second.value='S11';
const colors=['#176f93','#d56a36'];let initialized=false;
function facts(r,i){return '<div class="fact" style="border-color:'+colors[i]+'"><b>'+r.id+'</b> · '+r.maxima+' peaks · T = '+r.period.toFixed(6)+'<br>'+(r.symmetric?'Half-period inversion symmetry':'Asymmetric representative')+' · '+(r.both_outer_regions?'Both outer regions':'One outer region')+'<br>beta = '+r.reference_beta.toFixed(11)+'</div>'}
function signal(r,shift){const a=r.xyz,n=a.length-1;return a.map((v,i)=>{let s=((i/n+shift)%1)*n,k=Math.floor(s),f=s-k;return a[k][0]*(1-f)+a[Math.min(k+1,n)][0]*f})}
function render(){const selected=[D[first.value],second.value?D[second.value]:null],traces=[];document.querySelector('#facts').innerHTML=selected.map((r,i)=>r?facts(r,i):'').join('');document.querySelector('#shift').textContent=Math.round(slider.value*100)+'%';
 selected.forEach((r,i)=>{if(!r)return;const a=r.xyz,phase=a.map((_,j)=>j/(a.length-1));
 traces.push({type:'scatter3d',x:a.map(v=>v[0]),y:a.map(v=>v[1]),z:a.map(v=>v[2]),mode:'lines',name:r.id,line:{color:colors[i],width:4},scene:'scene',hovertemplate:r.id+'<br>x=%{x:.3f}<br>y=%{y:.3f}<br>z=%{z:.3f}<extra></extra>'});
 traces.push({type:'scatter',x:phase,y:signal(r,i?Number(slider.value):0),mode:'lines',name:r.id,xaxis:'x',yaxis:'y',showlegend:false,line:{color:colors[i],width:1.5},hovertemplate:r.id+'<br>fraction of period=%{x:.3f}<br>x=%{y:.4f}<extra></extra>'});
 traces.push({type:'scatter',x:r.peak_values.map((_,j)=>j+1),y:r.peak_values,mode:'lines+markers',name:r.id,xaxis:'x2',yaxis:'y2',showlegend:false,line:{color:colors[i],width:1},marker:{size:6},hovertemplate:r.id+'<br>peak %{x}<br>height=%{y:.5f}<extra></extra>'});});
 const layout={height:880,margin:{l:65,r:30,t:45,b:60},font:{family:'system-ui',color:'#193449'},uirevision:'keep-camera',legend:{x:0,y:1.03,orientation:'h'},
 scene:{domain:{x:[0,.48],y:[.40,1]},xaxis:{title:{text:'x'},range:[-2.5,2.5]},yaxis:{title:{text:'y'},range:[-.45,.45]},zaxis:{title:{text:'z'},range:[-3.45,3.45]},aspectmode:'cube'},
 xaxis:{domain:[.58,1],anchor:'y',title:{text:'Fraction of one complete period'},range:[0,1]},yaxis:{domain:[.48,.96],anchor:'x',title:{text:'x'},range:[-2.5,2.5]},
 xaxis2:{domain:[0,1],anchor:'y2',title:{text:'Peak number (saved starting phase)'},dtick:Math.max(...selected.filter(Boolean).map(r=>r.maxima))<=15?1:10},yaxis2:{domain:[0,.25],anchor:'x2',title:{text:'Peak height'},range:[-2.5,2.5]},
 annotations:[{text:'Rotate the 3D orbit',xref:'paper',yref:'paper',x:.20,y:1.02,showarrow:false},{text:'Complete repeating signal',xref:'paper',yref:'paper',x:.82,y:1.02,showarrow:false},{text:'Ordered peaks reveal the repetition pattern',xref:'paper',yref:'paper',x:.5,y:.32,showarrow:false}]};
 const action=initialized?Plotly.react:Plotly.newPlot;action('plot',traces,layout,{responsive:true,displaylogo:false});initialized=true;}
first.addEventListener('change',()=>{slider.value=0;render()});second.addEventListener('change',()=>{slider.value=0;render()});slider.addEventListener('input',render);render();
</script></body></html>'''


if __name__=='__main__':main()
