"""Compare the original wide cloud with the fourfold expansion, without a hull."""
import json
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from report import style
from explore_expanded import ROOT,OUT
from run import save_json

def candidates(data):
    return data['parameters'][(data['classes']==4).any(axis=1)]

def old_domain():
    vertices=np.array([[a,a*r,m] for a in [2,256] for r in [.25,64] for m in [-8,-1.003]])
    lines=[]
    for i in range(8):
        for bit in [1,2,4]:
            if i<(j:=i^bit):
                # Subdivide alpha edges before the nonlinear beta/beta_crit view.
                # Transforming only their endpoints would draw false straight borders.
                first,last=vertices[i],vertices[j]
                if first[0]!=last[0]:
                    a=np.geomspace(first[0],last[0],100)
                    fraction=(a-first[0])/(last[0]-first[0])
                else:fraction=np.linspace(0,1,2)
                lines.extend(first[None,:]+fraction[:,None]*(last-first)[None,:])
                lines.append([np.nan]*3)
    return np.array(lines)

def transform(p,mode):
    x=p[:,:3].copy()
    if mode=='log':x[:,:2]=np.log10(x[:,:2])
    if mode=='scaled':
        x[:,1]=x[:,1]/(5*x[:,0]/7+10*x[:,0]**2/49)
        x[:,0]=np.log10(x[:,0])
    return x

def cloud(p,old,scout,checks,mode):
    fig=go.Figure();xyz=transform(p,mode);older=transform(old,mode)
    fig.add_trace(go.Scatter3d(x=older[:,0],y=older[:,1],z=older[:,2],mode='markers',
        marker=dict(size=2,color='#9caac0',opacity=.25),name='Previous fixed-step candidates',hoverinfo='skip'))
    fig.add_trace(go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='markers',
        marker=dict(size=3,color=p[:,2],cmin=-29,cmax=-1,colorscale='Viridis',opacity=.85,
                    colorbar=dict(title='m₀',thickness=12)),customdata=p,
        name='Expanded adaptive candidates',hovertemplate='α=%{customdata[0]:.6g}<br>β=%{customdata[1]:.6g}<br>m₀=%{customdata[2]:.6g}<extra></extra>'))
    shown=transform(old_domain(),mode)
    fig.add_trace(go.Scatter3d(x=shown[:,0],y=shown[:,1],z=shown[:,2],mode='lines',
        line=dict(color='#f6be57',width=4),name='Previous search limits',hoverinfo='skip'))
    good=[c for c in checks if c['positive_in_all_independent_windows']]
    if good:
        actual=np.array([c['parameters'] for c in good]);shown=transform(actual,mode)
        fig.add_trace(go.Scatter3d(x=shown[:,0],y=shown[:,1],z=shown[:,2],mode='markers',
            marker=dict(size=7,color='white',symbol='diamond'),customdata=actual,
            name='Independent positive-LLE examples',hovertemplate='α=%{customdata[0]:.6g}<br>β=%{customdata[1]:.6g}<br>m₀=%{customdata[2]:.6g}<extra></extra>'))
    shown=transform(scout,mode)
    fig.add_trace(go.Scatter3d(x=shown[:,0],y=shown[:,1],z=shown[:,2],mode='markers',
        marker=dict(size=1,color='#90a0b4',opacity=.06),name='All new scouting points',visible='legendonly',hoverinfo='skip'))
    xa=dict(title='α');ya=dict(title='β')
    if mode in ['log','scaled']:
        ticks=np.array([2,8,32,128,512,1024]);xa=dict(title='α · log scale',tickvals=np.log10(ticks),ticktext=ticks.astype(str))
    if mode=='log':
        ticks=np.array([.01,1,100,10000,300000]);ya=dict(title='β · log scale',tickvals=np.log10(ticks),ticktext=['0.01','1','100','10,000','300,000'])
    if mode=='scaled':ya=dict(title='β / βcrit',range=[0,1.35])
    title={'linear':'Linear coordinates: the expanded outline','log':'Logarithmic α and β: old and new scales','scaled':'Divide β by the equilibrium-stability scale'}[mode]
    style(fig,title)
    fig.update_layout(height=700,margin=dict(l=10,r=40,t=70,b=100),scene=dict(xaxis=xa,yaxis=ya,zaxis_title='m₀',
        aspectmode='cube',camera=dict(eye=dict(x=1.8,y=-1.8,z=1.2),projection=dict(type='orthographic'))),
        legend=dict(orientation='h',y=-.08))
    return fig

def static(p,old):
    plt.rcParams.update({'figure.facecolor':'#f5f9fd','axes.facecolor':'#f5f9fd','font.size':10})
    fig=plt.figure(figsize=(16,7),layout='constrained')
    for i,mode in enumerate(['linear','log']):
        ax=fig.add_subplot(1,2,i+1,projection='3d');xyz=transform(p,mode);o=transform(old,mode)
        ax.scatter(*o.T,c='#a0adbd',s=3,alpha=.2,rasterized=True)
        points=ax.scatter(*xyz.T,c=p[:,2],cmap='viridis',s=9,alpha=.8,vmin=-29,vmax=-1,rasterized=True)
        box=transform(old_domain(),mode);ax.plot(*box.T,color='#bc851c',lw=1.5)
        ax.set(xlabel='α',ylabel='β',zlabel='m₀',title='Linear parameter coordinates' if i==0 else 'Logarithmic α and β')
        if mode=='log':
            ticks=[2,8,32,128,512,1024];ax.set_xticks(np.log10(ticks),labels=list(map(str,ticks)))
            ticks=[1,100,10000,300000];ax.set_yticks(np.log10(ticks),labels=['1','100','10,000','300,000'])
        ax.view_init(23,-55);ax.set_box_aspect((1,1,1))
    fig.colorbar(points,ax=fig.axes,shrink=.6,label='Inner diode slope m₀')
    fig.suptitle(f'CHUA / FOUR TIMES FARTHER\n{len(p):,} repeated adaptive GPU candidates; gold = previous search limits',fontsize=18)
    fig.savefig(OUT/'expanded_extent.png',dpi=170);fig.savefig(OUT/'expanded_extent.pdf');plt.close(fig)

def build():
    summary=json.loads((OUT/'summary.json').read_text());data=dict(np.load(OUT/'verified.npz'))
    p=candidates(data);old=candidates(np.load(ROOT/'outputs/extent/verified.npz'));scout=np.load(OUT/'scout.npz')['parameters']
    evidence=json.loads((OUT/'independent_checks.json').read_text()) if (OUT/'independent_checks.json').exists() else {'checks':[]}
    checks=evidence['checks'];static(p,old)
    def div(fig,ident):return fig.to_html(full_html=False,include_plotlyjs=False,div_id=ident,config={'responsive':True,'displaylogo':False})
    rows=''.join(f"<tr><td>{c['parameters'][0]:.4f}</td><td>{c['parameters'][1]:,.2f}</td><td>{c['parameters'][2]:.5f}</td><td>{c['matched_adaptive_lle']:.4f}</td><td>{c['independent']['lle']:.4f}</td><td>{'Positive in all windows' if c['positive_in_all_independent_windows'] else 'Not confirmed in all windows'}</td></tr>" for c in checks)
    shells=[]
    for lo,hi in [(2,256),(256,512),(512,768),(768,1024)]:
        selected=p[(p[:,0]>=lo)&(p[:,0]<hi)]
        shells.append({'alpha_interval':[lo,hi],'count':len(selected),
                       'minimum_m0':float(selected[:,2].min()) if len(selected) else None,
                       'maximum_beta':float(selected[:,1].max()) if len(selected) else None})
    save_json(OUT/'shells.json',shells)
    shell_rows=''.join(f"<tr><td>{s['alpha_interval'][0]}–{s['alpha_interval'][1]}</td><td>{s['count']}</td><td>{s['minimum_m0']}</td><td>{s['maximum_beta']}</td></tr>" for s in shells)
    times=summary['scout']['timing']['wall_seconds']+summary['verification']['timing']['wall_seconds']
    body=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>Chua — four times farther</title><script src="plotly.min.js"></script><style>
    *{{box-sizing:border-box}}body{{margin:0;background:#091321;color:#e3edf7;font:16px/1.7 system-ui,sans-serif}}
    main{{max-width:1550px;margin:auto;padding:42px 30px}}h1{{font-size:48px;line-height:1.1}}h2{{font-size:27px}}p,li{{color:#b8cadd}}a{{color:#64d5ca}}
    .pair{{display:grid;grid-template-columns:1fr 1fr;gap:20px}}.pair>*{{min-width:0}}.box{{padding:20px 26px;border-left:4px solid #f5ba5b;background:#172a42;border-radius:10px}}
    section{{margin-top:38px;border-top:1px solid #293c53;padding-top:24px}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%}}td,th{{text-align:left;padding:10px;border-bottom:1px solid #2b4056}}
    @media(max-width:1100px){{.pair{{grid-template-columns:1fr}}main{{padding:24px 14px}}h1{{font-size:34px}}}}
    </style></head><body><main><p>CHUA / EXPANDED GPU SURVEY</p><h1>Four times farther</h1>
    <div class="box"><b>Search limits: α ≤ 1024, m₀ ≥ −29, β up to about {summary['theoretical_beta_max']:,.0f}.</b>
    <p>Compared with the previous search, the α limit and depth −1−m₀ increase fourfold.
    The upper β search limit increases about 17-fold because its natural reference scale grows quadratically with α.</p>
    <p>{len(scout):,} new parameter tuples × three initial conditions; {len(p):,} candidates survive tighter tolerances and longer observation.
    Bulk scout and verification take {times:.1f} s on the GTX 1080, excluding independent CPU checks.</p></div>
    <p>The repeated cloud reaches α={p[:,0].max():.3f}, β={p[:,1].max():,.1f} and m₀={p[:,2].min():.5f} (separate coordinate extrema).
    There are {summary['alpha_above_256']:,} candidates beyond α=256 and {summary['m0_below_minus_8']:,} beyond m₀=−8.
    Near the new limits, {summary['alpha_above_768']} candidates have α&gt;768 and {summary['m0_below_minus_22']} have m₀&lt;−22.</p>
    <div class="pair">{div(cloud(p,old,scout,checks,'linear'),'expanded-linear')}{div(cloud(p,old,scout,checks,'log'),'expanded-log')}</div>
    <p>Color represents m₀. Gray points are the previous 2,249 candidates; the gold wireframe is the previous search domain.
    White diamonds mark selected independently positive examples. No surface fills the gaps. The logarithmic view changes the apparent geometry.
    Toggle “All new scouting points” in the legend to see tested locations outside the retained cloud.</p>
    <section><h2>Does the bend continue?</h2>
    <p>The third view uses q=β/βcrit, where βcrit=5α/7+10α²/49 is the outer-equilibrium linear-stability threshold.
    This removes the leading quadratic scale from β and makes changes in the cloud easier to compare across α.
    It is a coordinate transformation, not evidence that the full family is exactly self-similar.</p>
    {div(cloud(p,old,scout,checks,'scaled'),'expanded-scaled')}
    <div class="scroll"><table><tr><th>α interval</th><th>Repeated candidates</th><th>Lowest sampled m₀</th><th>Highest sampled β</th></tr>{shell_rows}</table></div>
    <p>Counts depend on the logarithmic sampling measure. They are not volumes or probabilities of chaos.</p></section>
    <section><h2>Independent checks at distant points</h2>
    <p>Selected extremes start from the end of the verified GPU trajectory, then discard 20 time units and measure 80 more.
    An independent DOP853 calculation evaluates the full piecewise field at every stage; a further GPU run uses rtol=10⁻⁹, atol=10⁻¹² with matched initial state and intervals.</p>
    <div class="scroll"><table><tr><th>α</th><th>β</th><th>m₀</th><th>Matched GPU λ₁</th><th>Full-field DOP853 λ₁</th><th>Independent sign check</th></tr>{rows}</table></div>
    <p>Agreement of signs is evidence for the selected finite-time regimes. Exponent magnitudes can differ, and neither a sign check nor this sparse cloud proves asymptotic strange-attractor existence.
    <a href="independent_checks.json">Full independent records</a> · <a href="adaptive_validation.json">Adaptive solver controls</a>.</p></section>
    <section><h2>What changed in the computation?</h2>
    <p>The expanded survey uses adaptive Dormand–Prince 5(4) on the GPU. Error control includes the state and tangent vector; normalization occurs after every accepted step.
    The scout uses rtol=10⁻⁶, atol=10⁻⁹, transient=100 and observation=200. Candidate verification uses rtol=10⁻⁸, atol=10⁻¹¹, transient=200 and observation=600.</p>
    <p>The coordinate escape cutoff increases from 1000 to 10⁶ so larger oscillations are not automatically excluded.
    Three starts include parameter-dependent neighborhoods of the outer equilibrium. Sampling is uniform in log α, log q and log(−1−m₀),
    with 2≤α≤1024, 0.01≤q≤1.3 and 0.003≤−1−m₀≤28. The fixed outer slope remains m₁=−5/7.</p>
    <p>The previous and expanded surveys differ in integrator, observation, escape cutoff and sampling coordinates. Their overlay compares geometric extent; it is not a controlled estimate of changed chaos volume.</p>
    <p>The original v1.0 CPU routine extended a single affine branch until an event. It failed at the new maximum-α example: a missed crossing could leave it on the wrong branch.
    The corrected reference always evaluates the actual piecewise field. Crossing events are recorded as diagnostics and do not control branch selection.
    Earlier v1.0 event-reference values are retained as historical records, but should not be treated as definitive reference measurements.</p>
    <p>A finite search still cannot settle whether the complete parameter family terminates or extends indefinitely, or whether all distant candidates connect.
    Bounded motion of one attractor in state space does not imply a bounded set of parameters supporting attractors.</p>
    <p><a href="../extent/index.html">Previous wide survey</a> · <a href="expanded_extent.png">Comparison PNG</a> · <a href="expanded_extent.pdf">PDF</a> ·
    <a href="verified.npz">Repeated data</a> · <a href="summary.json">Exact settings and counts</a></p></section>
    </main><script>window.addEventListener('load',()=>requestAnimationFrame(()=>document.querySelectorAll('.js-plotly-plot').forEach(p=>Plotly.Plots.resize(p))));</script></body></html>'''
    shutil.copy2(ROOT/'outputs/plotly.min.js',OUT/'plotly.min.js');(OUT/'index.html').write_text(body,encoding='utf-8')
    print(json.dumps({'candidates':len(p),'alpha_max':float(p[:,0].max()),'beta_max':float(p[:,1].max()),'m0_min':float(p[:,2].min()),'shells':shells},indent=2),flush=True)

if __name__=='__main__':build()
