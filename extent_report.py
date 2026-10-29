"""Render the coarse multiscale point cloud and explain its finite-domain limits."""
from pathlib import Path
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from report import style
from explore_extent import OUT,ROOT
from run import save_json

def original_box(logarithmic=False):
    limits=[(7,11),(100/7-4,100/7+4),(-8/7-.16,-8/7+.16)]
    verts=np.array([[limits[0][i],limits[1][j],limits[2][k]] for i in [0,1] for j in [0,1] for k in [0,1]])
    if logarithmic:verts[:,:2]=np.log2(verts[:,:2])
    lines=[]
    for i in range(8):
        for bit in [1,2,4]:
            j=i^bit
            if i<j:lines.extend([verts[i],verts[j],[np.nan]*3])
    return np.asarray(lines)

def cloud_figure(scout,data,logarithmic=False):
    keep=(data['classes']==4).any(axis=1)
    p=data['parameters'][keep]
    estimates=np.max(np.where(data['classes'][keep]==4,data['lle'][keep],-np.inf),axis=1)
    colors=p[:,2]
    selected=np.linspace(0,len(scout['parameters'])-1,min(12000,len(scout['parameters'])),dtype=int)
    other=scout['parameters'][selected]
    transformed=p[:,:3].copy();background=other[:,:3].copy()
    if logarithmic:
        transformed[:,:2]=np.log2(transformed[:,:2]);background[:,:2]=np.log2(background[:,:2])
    fig=go.Figure()
    fig.add_trace(go.Scatter3d(x=background[:,0],y=background[:,1],z=background[:,2],mode='markers',
        marker=dict(size=1,color='#70839b',opacity=.07),name='Scouted points',visible='legendonly',hoverinfo='skip'))
    fig.add_trace(go.Scatter3d(x=transformed[:,0],y=transformed[:,1],z=transformed[:,2],mode='markers',
        marker=dict(size=3,color=colors,colorscale='Viridis',opacity=.85,colorbar=dict(title='m₀',thickness=12)),
        customdata=np.c_[p[:,:3],estimates],name='Repeated GPU chaos candidates',
        hovertemplate='α=%{customdata[0]:.5g}<br>β=%{customdata[1]:.6g}<br>m₀=%{customdata[2]:.5f}<br>Rough GPU λ₁=%{customdata[3]:.4f}<extra></extra>'))
    check_file=OUT/('full_field_checks.json' if (OUT/'full_field_checks.json').exists() else 'event_resolved_checks.json')
    if check_file.exists():
        checks=json.loads(check_file.read_text())['checks']
        good=[r for r in checks if min(r['windows'])>.005 and r['lle']>.01]
        if good:
            actual=np.array([r['parameters'][:3] for r in good]);shown=actual.copy()
            if logarithmic:shown[:,:2]=np.log2(shown[:,:2])
            fig.add_trace(go.Scatter3d(x=shown[:,0],y=shown[:,1],z=shown[:,2],mode='markers',
                marker=dict(size=7,color='#f6be57',symbol='diamond'),name='Full-field positive-LLE checks' if check_file.name=='full_field_checks.json' else 'Historical v1.0 CPU check locations',
                customdata=np.c_[actual,[r['lle'] for r in good]],
                hovertemplate='α=%{customdata[0]:.4f}<br>β=%{customdata[1]:.3f}<br>m₀=%{customdata[2]:.5f}<br>Event-resolved λ₁=%{customdata[3]:.3f}<extra></extra>'))
    box=original_box(logarithmic)
    fig.add_trace(go.Scatter3d(x=box[:,0],y=box[:,1],z=box[:,2],mode='lines',
                              line=dict(color='#f5b75a',width=5),name='Original classic box',hoverinfo='skip'))
    axes=dict(xaxis_title='α',yaxis_title='β',zaxis_title='m₀',aspectmode='cube',
              camera=dict(eye=dict(x=1.8,y=-1.8,z=1.2),projection=dict(type='orthographic')))
    if logarithmic:
        axes.update(xaxis=dict(title='α · log scale',tickvals=[1,3,5,7,8],ticktext=['2','8','32','128','256']),
                    yaxis=dict(title='β · log scale',tickvals=[0,4,8,12,14],ticktext=['1','16','256','4096','16384']))
    style(fig,'Wide point cloud · '+('logarithmic α, β' if logarithmic else 'linear coordinates'))
    fig.update_layout(scene=axes,height=680,margin=dict(l=5,r=45,t=65,b=70))
    return fig

def envelope_plot(data):
    keep=(data['classes']==4).any(axis=1);p=data['parameters'][keep]
    alpha=np.geomspace(2,256,500);critical=5*alpha/7+10*alpha**2/49
    fig=go.Figure(go.Scattergl(x=p[:,0],y=p[:,1],mode='markers',marker=dict(size=4,color=p[:,2],colorscale='Plasma',opacity=.7,
        colorbar=dict(title='m₀',thickness=12)),name='GPU candidates',customdata=p[:,2],
        hovertemplate='α=%{x:.4g}<br>β=%{y:.5g}<br>m₀=%{customdata:.4f}<extra></extra>'))
    fig.add_trace(go.Scatter(x=alpha,y=critical,mode='lines',line=dict(color='#f6be57',dash='dash',width=2),
        name='Outer equilibria change stability'))
    style(fig,'The bend follows a growing parameter scale')
    fig.update_layout(xaxis=dict(title='α',type='log'),yaxis=dict(title='β',type='log'),height=500,
                      margin=dict(l=70,r=70,t=65,b=65),legend=dict(orientation='h',y=-.17))
    return fig

def static_figure(data):
    keep=(data['classes']==4).any(axis=1);p=data['parameters'][keep]
    plt.rcParams.update({'figure.facecolor':'#f5f9fd','axes.facecolor':'#f5f9fd','font.size':10})
    fig=plt.figure(figsize=(16,6.5),layout='constrained')
    for column,logarithmic in enumerate([False,True]):
        ax=fig.add_subplot(1,2,column+1,projection='3d',computed_zorder=False)
        xyz=p[:,:3].copy()
        if logarithmic:xyz[:,:2]=np.log2(xyz[:,:2])
        dots=ax.scatter(*xyz.T,c=p[:,2],s=7,cmap='viridis',alpha=.75,rasterized=True)
        box=original_box(logarithmic);ax.plot(*box.T,color='#bf7e1a',lw=1.6,zorder=10)
        if logarithmic:
            ax.set_xticks([1,3,5,7,8],labels=['2','8','32','128','256'])
            ax.set_yticks([0,3,6,9,12,14],labels=['1','8','64','512','4096','16384'])
        ax.set(xlabel='α' if not logarithmic else 'α (log spacing)',ylabel='β' if not logarithmic else 'β (log spacing)',
               zlabel='m₀',title='Linear parameter coordinates' if not logarithmic else 'Logarithmic α and β reveal the small-scale part')
        ax.view_init(23,-55);ax.set_box_aspect((1,1,1))
    fig.colorbar(dots,ax=fig.axes,shrink=.65,label='Inner diode slope m₀')
    fig.suptitle('CHUA / THE REGION EXTENDS FAR BEYOND THE ORIGINAL BOX\n32,768 scouting points; refined candidates; fixed m₁ = −5/7',fontsize=17)
    fig.savefig(OUT/'wide_extent.png',dpi=180)
    fig.savefig(OUT/'wide_extent.pdf')
    plt.close(fig)

def build():
    scout=dict(np.load(OUT/'scout.npz'));data=dict(np.load(OUT/'verified.npz'))
    summary=json.loads((OUT/'summary.json').read_text())
    check_file=OUT/('full_field_checks.json' if (OUT/'full_field_checks.json').exists() else 'event_resolved_checks.json')
    evidence=json.loads(check_file.read_text()) if check_file.exists() else {'checks':[]}
    historical=check_file.name=='event_resolved_checks.json'
    correction=('These are historical v1.0 CPU values. The affine-branch event routine later failed at a distant point after a missed crossing; these values are not definitive reference measurements. See the corrected full-field controls in the expanded study.' if historical else 'These checks use the corrected DOP853 reference, evaluating the full piecewise field at every stage. Crossing events are diagnostics and do not select a branch.')
    keep=(data['classes']==4).any(axis=1);p=data['parameters'][keep]
    crit=5*p[:,0]/7+10*p[:,0]**2/49
    above=int((p[:,1]>crit).sum())
    coarse=np.max(np.where(scout['classes'][data['indices']]==4,scout['lle'][data['indices']],-np.inf),axis=1)
    fine=np.max(np.where(data['classes']==4,data['lle'],-np.inf),axis=1)
    consistent=keep&(np.abs(fine-coarse)<np.maximum(.05,.4*np.maximum(abs(fine),abs(coarse))))
    analysis={'candidate_nodes':int(keep.sum()),'lle_consistent_between_runs':int(consistent.sum()),
              'candidates_above_outer_stability_boundary':above,
              'alpha_above_192':int((p[:,0]>192).sum()),'m0_below_minus_6':int((p[:,2]<-6).sum()),
              'boundary_equation':'beta = 5*alpha/7 + 10*alpha**2/49',
              'scope':'The stability boundary is not itself a proven chaos boundary. No global closure/unboundedness claim.'}
    save_json(OUT/'interpretation.json',analysis)
    static_figure(data)
    def div(fig,ident):return fig.to_html(full_html=False,include_plotlyjs=False,div_id=ident,config={'responsive':True,'displaylogo':False})
    rows=''.join(f"<tr><td>{r['parameters'][0]:.4f}</td><td>{r['parameters'][1]:.3f}</td><td>{r['parameters'][2]:.5f}</td><td>{r['gpu_lle']:.3f}</td><td>{r['lle']:.3f}</td></tr>" for r in evidence['checks'])
    body=f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>Chua — the outer extent of the parameter region</title><script src="plotly.min.js"></script><style>
    *{{box-sizing:border-box}}body{{margin:0;background:#091321;color:#e3edf7;font:16px/1.7 system-ui,sans-serif}}
    main{{max-width:1550px;margin:auto;padding:42px 30px}}h1{{font-size:46px;line-height:1.12}}h2{{font-size:27px}}p,li{{color:#b8cadd}}a{{color:#64d5ca}}
    .pair{{display:grid;grid-template-columns:1fr 1fr;gap:20px}}.pair>*{{min-width:0}}.box{{padding:20px 26px;border-left:4px solid #f5ba5b;background:#172a42;border-radius:10px}}
    section{{margin-top:38px;border-top:1px solid #293c53;padding-top:24px}}table{{border-collapse:collapse;width:100%}}td,th{{text-align:left;padding:10px;border-bottom:1px solid #2b4056}}
    .note{{font-size:13px}}@media(max-width:1100px){{.pair{{grid-template-columns:1fr}}main{{padding:24px 14px}}h1{{font-size:34px}}}}
    </style></head><body><main><p>CHUA / WIDE PARAMETER EXPLORATION / GPU</p><h1>Where does this region end?</h1>
    <div class="box"><b>No outer cutoff has yet been located toward larger α or more negative m₀ within the explored range.</b>
    <p>Candidates occur close to α=256 and m₀=−8, the new search limits. Observed β values reach {p[:,1].max():,.0f}.
    The original model was a small fragment. A finite point cloud cannot prove global boundedness or unboundedness.</p></div>
    <p>Every point represents parameters (α,β,m₀), not a trajectory point (x,y,z).
    This entire map fixes m₁=−5/7, matching the classic example. The second original example uses the slightly rounded m₁=−0.714;
    its section is kept separate.</p>
    <p>The scout sampled {len(scout['parameters']):,} parameter tuples with three starts each.
    After repeating candidates with half the time step and longer observation, {int(keep.sum()):,} candidates remain.
    No surface is interpolated between them; neither connections nor absence of other regimes in the gaps are established.</p>
    <div class="pair">{div(cloud_figure(scout,data),'extent-linear')}{div(cloud_figure(scout,data,True),'extent-log')}</div>
    <p class="note">The gold wireframe marks the original classic box. The right plot uses logarithmic α and β, changing the apparent geometry.
    Toggle “Scouted points” in the legend to see other tested points. Color denotes m₀;
    gold diamonds identify the recorded CPU check locations; see the reference-solver correction below.</p>
    <section><h2>Why the region bends</h2>
    <p>Outer equilibria change linear stability on the surface
    <b>β = 5α/7 + 10α²/49 ≈ 0.714α + 0.204α²</b>, introducing a β scale that grows faster than α.
    This helps interpret the observed bend. All {int(keep.sum()):,} repeated candidates lie below that surface;
    the surface itself is not a proven global boundary of chaos.</p>
    {div(envelope_plot(data),'stability-envelope')}
    <details><summary>Derivation of the stability surface</summary><p>On an outer branch the diode derivative is m₁.
    The Jacobian characteristic polynomial is λ³ + [1+α(1+m₁)]λ² + (β+αm₁)λ + αβ(1+m₁).
    Substitution of m₁=−5/7 into the cubic stability criterion gives β=5α/7+10α²/49.
    This is local linear equilibrium stability, not a calculation of the complete global attractor.</p></details></section>
    <section><h2>Independent checks of distant points</h2>
    <div class="box"><b>Reference-solver correction</b><p>{correction}</p>
    <p><a href="../expanded/index.html">Expanded study and corrected controls</a></p></div>
    <table><tr><th>α</th><th>β</th><th>m₀</th><th>GPU λ₁</th><th>Recorded CPU λ₁</th></tr>{rows}</table>
    <p class="note">Quantitative convergence of the original fixed-step cloud is not established.
    Read this as a coarse candidate map, not a precise map of chaos intensity.
    <a href="{check_file.name}">Recorded checks</a> · <a href="summary.json">GPU settings</a>.</p></section>
    <section><h2>Scope of the search</h2>
    <p>Sobol points sample logarithmic coordinates: 2≤α≤256, 0.25≤β/α≤64, 0.003≤−1−m₀≤7.
    Only the m₀&lt;−1 continuation is explored. Low-discrepancy sampling reduces large gaps but can still miss thin windows and separate regions.</p>
    <p>The scout uses transient=100 and observation=200; the repeated pass uses 300 and 800.
    Time steps are grouped by the system rates and then halved. The coordinate escape cutoff is 1000 because distant attractors can be much larger.
    Initial conditions include neighborhoods of the parameter-dependent outer equilibrium.</p>
    <p>The sample contains {analysis['alpha_above_192']} candidates with α&gt;192 and {analysis['m0_below_minus_6']} with m₀&lt;−6.
    Truncation in these directions therefore does not look like a located natural edge.
    Continuous connection of all distant points to the original region, global closure, and unboundedness remain unresolved.</p>
    <p>Boundedness of an individual attractor in state space does not require boundedness of the parameter set supporting attractors.</p>
    <p><a href="../detailed/index.html">Refined atlas: thin windows and coexisting attractors</a> ·
    <a href="wide_extent.png">Overview PNG</a> · <a href="wide_extent.pdf">PDF</a> · <a href="verified.npz">Repeated candidate calculations</a></p></section>
    </main><script>window.addEventListener('load',()=>requestAnimationFrame(()=>document.querySelectorAll('.js-plotly-plot').forEach(p=>Plotly.Plots.resize(p))));</script></body></html>"""
    shutil.copyfile(ROOT/'outputs'/'plotly.min.js',OUT/'plotly.min.js')
    (OUT/'index.html').write_text(body,encoding='utf-8')
    print('Wide exploration report saved',analysis,flush=True)

if __name__=='__main__':build()
