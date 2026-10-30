"""Build a self-contained interactive viewer of both measured slope volumes."""
import argparse,json
from pathlib import Path
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plot_chua_slope_volumes import load,crop,surface,COLORS,KINDS


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);a=p.parse_args();data=load(a.source)
    fig=make_subplots(rows=2,cols=2,specs=[[{'type':'scene'},{'type':'scene'}],[{'type':'xy'},{'type':'xy'}]],
        row_heights=[.66,.34],vertical_spacing=.11,horizontal_spacing=.075,
        subplot_titles=['Vary alpha; beta = 100/7','Vary beta; alpha = 9','Full m0-m1 section through the alpha volume','Full m0-m1 section through the beta volume'])
    scale=[pair for k,color in enumerate(COLORS) for pair in [(k/7,color),((k+1)/7,color)]];indices={};payload={}
    for col,kind in enumerate(KINDS,1):
        d=data[kind];indices[kind]={}
        for interior in [False,True]:
            scope='interior' if interior else 'full';cut=crop(d,interior)
            for label in [2,4]:
                v,f=surface(cut,label);indices[kind][scope+str(label)]=len(fig.data)
                fig.add_trace(go.Mesh3d(x=v[:,0],y=v[:,1],z=v[:,2],i=f[:,0],j=f[:,1],k=f[:,2],color=COLORS[label],
                    visible=interior and label==2,opacity=.90 if label==2 else .24,
                    name=('Regular' if label==2 else 'Chaos')+' candidate boundary',showlegend=False,
                    hovertemplate=f'm0: %{{x:.5f}}<br>m1: %{{y:.5f}}<br>{kind}: %{{z:.5f}}<extra></extra>',
                    lighting=dict(ambient=.72,diffuse=.65,specular=.12,roughness=.8),flatshading=False),row=1,col=col)
        indices[kind]['plane']=len(fig.data)
        fig.add_trace(go.Surface(x=d['m0'],y=d['m1'],z=np.full((len(d['m1']),len(d['m0'])),d['control'][20]),surfacecolor=d['labels'][20],
            cmin=-.5,cmax=6.5,colorscale=scale,showscale=False,visible=False,opacity=.96,hoverinfo='skip',showlegend=False),row=1,col=col)
        indices[kind]['map']=len(fig.data)
        fig.add_trace(go.Heatmap(x=d['m0'],y=d['m1'],z=d['labels'][20],zmin=-.5,zmax=6.5,colorscale=scale,showscale=False,
            hovertemplate='m0: %{x:.5f}<br>m1: %{y:.5f}<br>Class: %{z}<extra></extra>',name=kind+' section'),row=2,col=col)
        fig.update_xaxes(title_text='m0',row=2,col=col);fig.update_yaxes(title_text='m1',row=2,col=col)
        cut=crop(d,True);scene='scene' if col==1 else 'scene2'
        fig.update_layout(**{scene:dict(xaxis=dict(title='m0',range=[cut['m0'][0],cut['m0'][-1]],nticks=5),
            yaxis=dict(title='m1',range=[cut['m1'][0],cut['m1'][-1]],nticks=5),zaxis=dict(title=kind,range=[d['control'][0],d['control'][-1]],nticks=5),
            aspectmode='cube',camera=dict(eye=dict(x=1.6,y=-1.8,z=1.2),up=dict(x=0,y=0,z=1),projection=dict(type='orthographic')),
            uirevision='keep-view')})
        payload[kind]={k:d[k].tolist() for k in ['m0','m1','control','labels']}
    fig.update_layout(height=1030,margin=dict(l=48,r=25,t=45,b=45),paper_bgcolor='white',plot_bgcolor='white',font=dict(family='Arial, sans-serif',color='#193449'),uirevision='keep-view')
    script=r'''
const model=__MODEL__, traces=__TRACES__, graph=document.getElementById('volumes');
let updating=false, requested=false;
async function refresh(){
 requested=true;if(updating)return;updating=true;
 try{do{requested=false;await renderView();}while(requested);}finally{updating=false;}
}
async function renderView(){
 const scope=document.getElementById('scope').value, mode=document.getElementById('mode').value;
 const showPlane=document.getElementById('show-plane').checked, visibility=new Array(graph.data.length).fill(false), layout={};
 const levels={alpha:Number(document.getElementById('alpha-slice').value),beta:Number(document.getElementById('beta-slice').value)};
 for(const [n,kind] of ['alpha','beta'].entries()){
  const d=model[kind], t=traces[kind], i=levels[kind], level=d.control[i];
  document.getElementById(kind+'-value').textContent=level.toFixed(5);
  visibility[t[scope+'2']]=mode!=='chaos';visibility[t[scope+'4']]=mode!=='regular';visibility[t.plane]=showPlane;visibility[t.map]=true;
  const xi=d.m0.map((v,k)=>k).filter(k=>scope==='full'||d.m0[k]<=-1.06);
  const yi=d.m1.map((v,k)=>k).filter(k=>scope==='full'||(d.m1[k]>=-.90&&d.m1[k]<=-.20));
  const x=xi.map(k=>d.m0[k]),y=yi.map(k=>d.m1[k]),z=yi.map(()=>xi.map(()=>level)),colors=yi.map(j=>xi.map(k=>d.labels[i][j][k]));
  await Plotly.restyle(graph,{x:[x],y:[y],z:[z],surfacecolor:[colors]},[t.plane]);
  await Plotly.restyle(graph,{z:[d.labels[i]]},[t.map]);
  const scene=n===0?'scene':'scene2';layout[scene+'.xaxis.range']=[x[0],x.at(-1)];layout[scene+'.yaxis.range']=[y[0],y.at(-1)];
 }
 await Plotly.restyle(graph,{visible:visibility});await Plotly.relayout(graph,layout);
 document.getElementById('status').textContent=(scope==='full'?'Full volume':'Interior crop')+' | '+mode+' candidates | alpha = '+document.getElementById('alpha-value').textContent+' | beta = '+document.getElementById('beta-value').textContent;
 document.getElementById('status').dataset.ready='true';
}
for(const id of ['scope','mode','show-plane'])document.getElementById(id).addEventListener('change',refresh);
for(const kind of ['alpha','beta'])document.getElementById(kind+'-slice').addEventListener('input',refresh);
document.getElementById('reset').addEventListener('click',()=>{
 const camera={eye:{x:1.6,y:-1.8,z:1.2},up:{x:0,y:0,z:1},projection:{type:'orthographic'}};
 Plotly.relayout(graph,{'scene.camera':camera,'scene2.camera':camera});
});
refresh();
'''.replace('__MODEL__',json.dumps(payload,separators=(',',':'))).replace('__TRACES__',json.dumps(indices))
    chart=fig.to_html(full_html=False,include_plotlyjs=True,div_id='volumes',post_script=script,config=dict(displaylogo=False,responsive=True,scrollZoom=False))
    page='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Chua: two three-dimensional slope sections</title><style>
body{margin:0;background:#f1f5f8;color:#193449;font:16px/1.45 Arial,sans-serif}main{max-width:1500px;margin:auto;padding:24px}
h1{font-size:28px;margin:0 0 8px}p{margin:8px 0 16px}.controls{display:flex;gap:18px;align-items:center;flex-wrap:wrap;background:white;padding:16px;border-radius:10px}
label{display:flex;gap:8px;align-items:center}select,button{font:inherit;padding:5px 9px}input[type=range]{width:180px}output{min-width:78px;font-variant-numeric:tabular-nums}
.legend{display:flex;gap:16px;flex-wrap:wrap;font-size:13px;margin-top:14px}.legend span:before{content:"";display:inline-block;width:12px;height:12px;margin-right:5px;background:var(--color)}
#status{font-size:13px;margin:10px 0;color:#526673}.note{font-size:13px;color:#526673}.plotly-graph-div{border-radius:10px}details{background:white;padding:12px;margin-top:15px}
</style></head><body><main><h1>Chua: two three-dimensional slope sections</h1>
<p>Drag either 3D view to rotate it. Compare the full volume with an interior crop, then move the sections through alpha or beta. Everything works offline.</p>
<div class="controls"><label>View <select id="scope"><option value="interior">Interior crop</option><option value="full">Full volume</option></select></label>
<label>Show <select id="mode"><option value="regular">Regular candidates</option><option value="chaos">Chaos candidates</option><option value="both">Both</option></select></label>
<label><input id="show-plane" type="checkbox"> Show section planes in 3D</label><button id="reset">Reset cameras</button>
<label>alpha <input id="alpha-slice" type="range" min="0" max="40" value="20"><output id="alpha-value">9.00000</output></label>
<label>beta <input id="beta-slice" type="range" min="0" max="40" value="20"><output id="beta-value">14.28571</output></label></div>
<div class="legend"><span style="--color:#d7dce1">0 Escape cutoff</span><span style="--color:#54b795">1 Equilibrium candidate</span>
<span style="--color:#e69b3f">2 Regular candidate</span><span style="--color:#eddc6b">3 Uncertain</span><span style="--color:#365f9e">4 Chaos candidate</span>
<span style="--color:#161616">5 Numerical failure</span><span style="--color:#a66ba8">6 Different labels across starts</span></div><div id="status">Loading measured volumes...</div>
'''+chart+'''
<details><summary>What these surfaces do and do not mean</summary><p>The alpha volume fixes beta = 100/7; the beta volume fixes alpha = 9. Each contains 61 x 61 x 41 sampled parameter points with three starts per point. The lower panels show actual sampled sections, with m0 along the horizontal axis.</p>
<p>Surfaces interpolate a binary candidate mask at level 0.5, without smoothing, convex hulls or artificial end caps. They are not proven bifurcation surfaces. Thin layers can fall between grid nodes. The interior crop keeps m0 &le; -1.06 and -0.90 &le; m1 &le; -0.20; nothing outside it is displayed in that mode.</p>
<p>Classes describe finite observation intervals. Different labels across starts do not by themselves establish coexistence. See the accompanying report for step/time sensitivity and the low-alpha control whose original chaos label was not retained in the longer checks.</p></details>
</main></body></html>'''
    # The bundled Plotly shaders contain trailing spaces; keep generated diffs clean.
    page='\n'.join(line.rstrip(' \t') for line in page.split('\n'))
    target=a.source/'volumes.html';target.write_text(page,encoding='utf-8',newline='\n')
    (a.source/'viewer-summary.json').write_text(json.dumps(dict(traces=len(fig.data),indices=indices,section_levels=41,bytes=target.stat().st_size,
        self_contained=True,scope='All volume nodes contribute to the mask meshes; the section sliders show saved grid planes.'),indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Wrote interactive viewer:',len(fig.data),'traces,',target.stat().st_size,'bytes',flush=True)


if __name__=='__main__':main()
