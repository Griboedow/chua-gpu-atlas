"""Interpret the refined GPU survey and build an interactive parameter atlas."""
from __future__ import annotations
import html
import json
from pathlib import Path
import shutil
import numpy as np
from scipy.ndimage import label
from skimage.measure import marching_cubes
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from report import style, COLORS
from run import save_json

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"/"detailed"

def load(name):
    return dict(np.load(OUT/f"{name}.npz"))

def meta(name):
    return json.loads((OUT/f"{name}.json").read_text(encoding="utf-8"))

def aggregate(classes):
    result=np.full(len(classes),3,dtype=np.uint8)
    result[(classes==0).all(axis=1)]=0
    result[(classes==1).all(axis=1)]=1
    result[np.isin(classes,[0,1,2]).all(axis=1)&(classes==2).any(axis=1)]=2
    result[(classes==4).any(axis=1)]=4
    result[(classes==5).any(axis=1)&~(classes==4).any(axis=1)]=5
    return result

def mesh(mask,data,color,name,visible=True,export=None,opacity=.65):
    axes=[data[k] for k in ["alpha","beta","m0"]]
    volume=np.asarray(mask).reshape(tuple(len(a) for a in axes))
    if volume.min()==volume.max():
        return go.Mesh3d(x=[],y=[],z=[],name=name,visible=visible)
    verts,faces,_,_=marching_cubes(volume.astype(np.float32),level=.5,
                    spacing=tuple(float(a[1]-a[0]) for a in axes),allow_degenerate=False)
    verts+=np.array([a[0] for a in axes])
    if export:
        np.savez_compressed(OUT/f"{export}_mesh.npz",vertices=verts,faces=faces,
                            coordinate_names=np.array(["alpha","beta","m0"]))
    return go.Mesh3d(x=verts[:,0],y=verts[:,1],z=verts[:,2],i=faces[:,0],j=faces[:,1],k=faces[:,2],
        color=color,opacity=opacity,name=name,visible=visible,showlegend=True,hoverinfo="skip",
        lighting=dict(ambient=.65,diffuse=.7,specular=.15,roughness=.7),flatshading=False)

def scene(data):
    return dict(xaxis=dict(title="α",range=[float(data['alpha'][0]),float(data['alpha'][-1])]),
                yaxis=dict(title="β",range=[float(data['beta'][0]),float(data['beta'][-1])]),
                zaxis=dict(title="m₀",range=[float(data['m0'][0]),float(data['m0'][-1])]),
                aspectmode="cube",camera=dict(eye=dict(x=1.8,y=1.8,z=1.2),projection=dict(type="orthographic")))

def region(data,metadata):
    c=data["classes"];chaos=c==4;both=data["visits_both_outer_regions"]
    masks=[chaos.any(axis=1),chaos.all(axis=1),(chaos&~both).any(axis=1),(chaos&both).any(axis=1)]
    titles=["Any seed","All seeds","One outer region","Both outer regions"]
    colors=["#27bbba","#f6be57","#8f96ff","#f475a0"]
    fig=go.Figure()
    for idx,(mask,color,title) in enumerate(zip(masks,colors,titles)):
        fig.add_trace(mesh(mask,data,color,title,idx==0,metadata['name'] if idx==0 else None))
    idx=np.flatnonzero(masks[0]&~masks[1])
    if len(idx)>15000:idx=idx[np.linspace(0,len(idx)-1,15000,dtype=int)]
    p=data["parameters"][idx]
    fig.add_trace(go.Scatter3d(x=p[:,0],y=p[:,1],z=p[:,2],mode="markers",visible=False,
        marker=dict(size=2,color="#f6be57",opacity=.55),name="Seed disagreement",
        hovertemplate="α=%{x:.4f}<br>β=%{y:.4f}<br>m₀=%{z:.5f}<extra></extra>"))
    center=metadata["configuration"]["parameters"]
    fig.add_trace(go.Scatter3d(x=[center[0]],y=[center[1]],z=[center[2]],mode="markers",
        marker=dict(color="white",size=6,symbol="diamond"),name="Reference"))
    fig.add_trace(go.Scatter3d(x=[center[0]]*2,y=[data['beta'][0],data['beta'][-1]],z=[center[2]]*2,
                              mode="lines",line=dict(color="white",width=4),name="β transect"))
    xx,yy=np.meshgrid([data['alpha'][0],data['alpha'][-1]],[data['beta'][0],data['beta'][-1]])
    fig.add_trace(go.Surface(x=xx,y=yy,z=np.full((2,2),-1.),visible="legendonly",showlegend=True,
                            name="m₀ = −1",showscale=False,colorscale=[[0,"#b5c5da"],[1,"#b5c5da"]],opacity=.3))
    buttons=[]
    for k,title in enumerate(titles+["Seed disagreement"]):
        vis=[False]*5+[True,True,"legendonly"]
        vis[k]=True
        buttons.append(dict(label=title,method="update",args=[{"visible":vis}]))
    style(fig,f"{metadata['grid']}³ parameter grid · m₁ = {center[3]:.6f}")
    fig.update_layout(scene=scene(data),height=650,margin=dict(l=10,r=10,t=110,b=70),
        updatemenus=[dict(type="dropdown",buttons=buttons,x=0,y=1.01,xanchor="left",yanchor="bottom",
                         bgcolor="#d4deeb",font=dict(color="#13273b"))])
    return fig

def slice_slider(data,metadata):
    n=metadata["grid"];volume=aggregate(data['classes']).reshape(n,n,n)
    scale=[]
    for k,color in enumerate(COLORS):scale.extend([[k/6,color],[(k+1)/6,color]])
    common=dict(x=data['alpha'],y=data['beta'],zmin=-.5,zmax=5.5,colorscale=scale,
                colorbar=dict(thickness=12,tickvals=list(range(6)),ticktext=["Escaped","Equilibrium","Near zero","Uncertain","Chaos","Failure"]),
                hovertemplate="α=%{x:.5f}<br>β=%{y:.5f}<br>Class=%{z}<extra></extra>")
    fig=go.Figure(go.Heatmap(z=volume[:,:,n//2].T,**common))
    steps=[dict(label=f"{data['m0'][k]:.5f}",method="restyle",
                args=[{"z":[volume[:,:,k].T.tolist()]},[0]]) for k in range(n)]
    style(fig,"Move the slice through m₀")
    fig.update_layout(height=520,margin=dict(l=65,r=100,t=65,b=105),xaxis_title="α",yaxis_title="β",
        sliders=[dict(active=n//2,currentvalue=dict(prefix="m₀ = "),pad=dict(t=45),steps=steps)])
    return fig

def beta_plot(data,name):
    fig=make_subplots(rows=3,cols=1,shared_xaxes=True,vertical_spacing=.08,row_heights=[.32,.5,.18],
                      subplot_titles=["Largest Lyapunov exponent · three starts","Last local maxima of x(t)","Repeated peak pattern length"])
    beta=data['beta']
    for j,color in enumerate(["#f6be57","#90a9c1","#657a92"]):
        fig.add_trace(go.Scatter(x=beta,y=data['lle'][:,j],mode="lines",line=dict(width=1,color=color),
                                name=f"Start {j+1}"),row=1,col=1)
    peaks=data['peaks'][:,-64:]
    fig.add_trace(go.Scattergl(x=np.repeat(beta,64),y=peaks.ravel(),mode="markers",
        marker=dict(size=1.4,color="#52c9bf",opacity=.5),name="x maxima",hoverinfo="skip"),row=2,col=1)
    period=data['peak_period'].astype(float);period[period==0]=np.nan
    fig.add_trace(go.Scatter(x=beta,y=period,mode="markers",marker=dict(size=3,color="#f475a0"),
                            name="Peak recurrence"),row=3,col=1)
    fig.add_hline(y=0,row=1,col=1,line=dict(color="#8ba1b8",width=1,dash="dot"))
    windows=[(14.31,14.33)] if name=="classic" else [(25.165625,25.371875),(28.86875,30.040625)]
    for lo,hi in windows:
        fig.add_vrect(x0=lo,x1=hi,fillcolor="#f475a0",opacity=.12,line_width=0,row="all",col=1)
    style(fig,"β transect · α and m₀ fixed at the reference")
    fig.update_layout(height=700,margin=dict(l=60,r=20,t=80,b=60),legend=dict(orientation="h",y=-.11))
    fig.update_xaxes(title_text="β",row=3,col=1)
    fig.update_yaxes(title_text="λ₁",row=1,col=1)
    fig.update_yaxes(title_text="x maxima",row=2,col=1)
    return fig

def zoom_figure(data,metadata):
    c=data['classes']
    regular=np.isin(c,[1,2]).all(axis=1)&(c==2).any(axis=1)
    fig=go.Figure(mesh(regular,data,"#f475a0","Near-zero LLE · all starts",export="window_regular",opacity=.85))
    fig.add_trace(mesh(data['candidate_all'],data,"#27bbba","Chaos · all starts",opacity=.12))
    fig.add_trace(go.Scatter3d(x=[9,9],y=[14.24,14.40],z=[-8/7]*2,mode="lines",line=dict(color="white",width=5),name="β transect"))
    fig.add_trace(go.Scatter3d(x=[9],y=[100/7],z=[-8/7],mode="markers",marker=dict(size=7,color="#f6be57"),name="Classic reference"))
    fig.add_trace(go.Scatter3d(x=[9],y=[14.32],z=[-8/7],mode="markers",marker=dict(size=7,color="white"),name="Five-peak example"))
    style(fig,"Local 3D zoom · thin regular window inside chaos")
    fig.update_layout(scene=scene(data),height=650,margin=dict(l=10,r=10,t=65,b=70))
    return fig

def basin_plot(data):
    scale=[[0,"#25374b"],[.249,"#25374b"],[.25,"#43c6bf"],[.499,"#43c6bf"],
           [.5,"#ee8caa"],[.749,"#ee8caa"],[.75,"#f6be57"],[1,"#f6be57"]]
    # The delivered basin contains labels 1,2 and the unresolved exact origin (3).
    fig=go.Figure(go.Heatmap(x=data['x'],y=data['y'],z=data['basin'].T,zmin=0,zmax=3,
        colorscale=[[0,"#25374b"],[.1666,"#25374b"],[.1667,"#43c6bf"],[.4999,"#43c6bf"],
                    [.5,"#ee8caa"],[.8333,"#ee8caa"],[.8334,"#f6be57"],[1,"#f6be57"]],
        colorbar=dict(tickvals=[0,1,2,3],ticktext=["Escaped","Left attractor","Right attractor","Unresolved"],thickness=15),
        hovertemplate="x₀=%{x:.4f}<br>y₀=%{y:.4f}<br>Basin=%{z}<extra></extra>"))
    style(fig,"Two basins at the same parameters · α=15.6, β=32")
    fig.update_layout(xaxis_title="Initial x₀",yaxis_title="Initial y₀ · z₀=0",margin=dict(l=60,r=120,t=65,b=60),height=550)
    return fig

def example_plot(items):
    titles=["β=28 · double-scroll chaos","β=29 · recurrent window","β=32 · two single-scroll attractors","β=36 · two regular cycles"]
    fig=make_subplots(rows=2,cols=2,subplot_titles=titles,horizontal_spacing=.1,vertical_spacing=.15)
    for idx,item in enumerate(items):
        data=load("example_"+item['name'])
        for j,color in enumerate(["#43c6bf","#ee8caa"]):
            xyz=data['xyz'][j,::4]
            fig.add_trace(go.Scatter(x=xyz[:,0],y=xyz[:,2],mode="lines",line=dict(width=.8,color=color),
                name="Start +0.7" if j==0 else "Start −0.7",showlegend=idx==0,hoverinfo="skip"),row=idx//2+1,col=idx%2+1)
        fig.update_xaxes(title_text="x",row=idx//2+1,col=idx%2+1)
        fig.update_yaxes(title_text="z",row=idx//2+1,col=idx%2+1)
    style(fig,"Same equations, four different long-time regimes")
    fig.update_layout(height=750,margin=dict(l=60,r=20,t=90,b=60))
    return fig

def summarize(name,data,metadata):
    mask=data['candidate_any'];c=data['classes']; n=metadata['grid']
    components,ncomp=label(mask.reshape(n,n,n))
    sizes=np.bincount(components.ravel())[1:]
    return {"name":name,"grid":n,"nodes":len(mask),"candidate_nodes":int(mask.sum()),
        "candidate_fraction":float(mask.mean()),"all_seeds_fraction_of_candidates":float(data['candidate_all'].sum()/mask.sum()),
        "largest_6_connected_component_fraction_of_candidates":float(sizes.max()/mask.sum()),
        "connected_components":ncomp,"single_node_components":int((sizes==1).sum()),
        "maximum_candidate_m0":float(data['parameters'][mask,2].max()),
        "chaotic_seed_runs_visiting_both_outer_regions":float(data['visits_both_outer_regions'][c==4].mean()),
        "grid_step":[float(data[k][1]-data[k][0]) for k in ['alpha','beta','m0']],
        "scope":"Grid connectivity and fractions, not invariant topology, basin probabilities, or exact volumes."}

def static_findings(classic,reference,basin,items):
    plt.rcParams.update({"figure.facecolor":"#f6f9fd","axes.facecolor":"#f6f9fd","font.size":11})
    fig,axes=plt.subplots(2,2,figsize=(15,10),layout="constrained")
    ax=axes[0,0]
    ax.plot(classic['beta'],classic['lle'][:,0],color="#147c85",lw=1.5)
    ax.axvspan(14.31,14.33,color="#ef89ab",alpha=.35,label="Repeated five-peak orbit near β=14.32")
    ax.axvline(100/7,color="#a56b08",ls="--",label="Original reference β=100/7")
    ax.axhline(0,color="#697b8c",lw=.7)
    ax.set(xlim=(14.24,14.40),ylim=(-.03,.42),xlabel="β",ylabel="Largest Lyapunov exponent",
           title="A 0.24% change of β can remove chaos")
    ax.legend(fontsize=8,loc="upper right")
    ax=axes[0,1]
    peaks=reference['peaks'][:,-48:]
    ax.scatter(np.repeat(reference['beta'],48),peaks.ravel(),s=.3,alpha=.35,c="#147c85",rasterized=True)
    ax.axvspan(28.86875,30.040625,color="#ef89ab",alpha=.25)
    ax.set(xlabel="β",ylabel="Local maxima of x(t)",title="Regular windows interrupt the chaotic bands",xlim=(23,38))
    ax=axes[1,0]
    xyz=load('example_single_scroll_pair')['xyz']
    for j,color in enumerate(["#147c85","#c4547d"]):
        ax.plot(xyz[j,::2,0],xyz[j,::2,2],lw=.3,color=color,label="Start +0.7" if j==0 else "Start −0.7")
    ax.set(xlabel="x",ylabel="z",title="β=32: two chaotic attractors coexist")
    ax.legend(loc="upper right",fontsize=9)
    ax=axes[1,1]
    cmap=ListedColormap(["#25374b","#43bfb8","#e77da0","#f1be55","#8661bc"])
    ax.pcolormesh(basin['x'],basin['y'],basin['basin'].T,cmap=cmap,vmin=0,vmax=4,shading="nearest",rasterized=True)
    ax.set(xlabel="Initial x₀",ylabel="Initial y₀ (z₀=0)",title="Initial conditions choose the left or right attractor")
    fig.suptitle("CHUA / WHAT THE REFINED MAP REVEALS",fontsize=20)
    fig.savefig(OUT/"findings.png",dpi=180)
    fig.savefig(OUT/"findings.pdf")
    plt.close(fig)

def build():
    records={n:(load(n+'_scan'),meta(n+'_scan')) for n in ['classic','reference','window_zoom']}
    lines={n:load(n+'_beta') for n in ['classic','reference']}
    basin=load('basin');items=meta('examples')
    summaries=[summarize(n,*records[n]) for n in ['classic','reference']]
    save_json(OUT/'analysis_summary.json',summaries)
    static_findings(lines['classic'],lines['reference'],basin,items)
    def div(fig,ident):
        return fig.to_html(full_html=False,include_plotlyjs=False,div_id=ident,config={"responsive":True,"displaylogo":False})
    global_parts=[]
    for name,summary in zip(['classic','reference'],summaries):
        data,metadata=records[name]
        global_parts.append(f"""<article><h2>{html.escape(metadata['configuration']['label'])}</h2>
        <p>{summary['candidate_nodes']:,} chaos candidate nodes ({summary['candidate_fraction']:.2%} of the grid).
        All three starts pass at {summary['all_seeds_fraction_of_candidates']:.1%} of these nodes.</p>
        {div(region(data,metadata),'region-'+name)}{div(slice_slider(data,metadata),'slices-'+name)}
        <p class="note">Clipped at the scan limits. The white line is the β transect shown below.
        One/Both outer regions records visits to x&lt;−1 and x&gt;1 during observation.
        Seed disagreement means different finite-time classifications; this alone does not prove coexisting attractors.</p>
        <p><a href="{name}_mesh.npz">Triangular surface mesh</a> · <a href="{name}_scan.json">Computation settings</a></p></article>""")
    body=f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>Chua — refined atlas and findings</title><script src="plotly.min.js"></script><style>
    *{{box-sizing:border-box}}body{{margin:0;background:#091321;color:#e1edf8;font:16px/1.65 system-ui,sans-serif}}
    main{{max-width:1550px;margin:auto;padding:40px 30px}}h1{{font-size:48px;line-height:1.1}}h2{{font-size:26px}}h3{{font-size:21px}}
    p,li{{color:#b7c9dc}}a{{color:#63d3cc}}.eyebrow{{color:#63d3cc;letter-spacing:.15em;font-size:12px;text-transform:uppercase}}
    .pair{{display:grid;grid-template-columns:1fr 1fr;gap:22px}}.pair>*{{min-width:0}}section{{margin-top:45px;border-top:1px solid #2c3e54;padding-top:26px}}
    .callout{{background:#192b43;border-left:4px solid #f475a0;border-radius:8px;padding:20px 26px;margin:25px 0}}
    .note{{font-size:13px}}.plotly-graph-div{{border-radius:12px;overflow:hidden}}table{{border-collapse:collapse;width:100%}}td,th{{padding:10px 14px;border-bottom:1px solid #2b3d52;text-align:left}}
    @media(max-width:1100px){{.pair{{grid-template-columns:1fr}}main{{padding:24px 14px}}h1{{font-size:36px}}}}
    </style></head><body><main><div class="eyebrow">Python / OpenCL / refined parameter atlas</div><h1>Structure inside chaos</h1>
    <p>Each point of the 3D parameter map represents (α, β, m₀).
    A separate trajectory in state coordinates (x,y,z) shows the motion produced by one parameter tuple.</p>
    <div class="callout"><b>Nearby parameters can produce very different long-time behavior.</b>
    <p>The classic β=100/7≈14.2857 example is chaotic. At β=14.32, a change of just 0.24%,
    the calculation finds a repeating five-peak pattern and λ₁ close to zero. The old β grid spacing of 0.16 missed this window.</p></div>
    <p>The general grids increase from 51³ to 81³: 3,188,646 runs with three initial conditions.
    A separate 41³ zoom uses dt=0.005, transient=600 and observation=1600.
    Each β transect contains 1,601 points; retained local maxima provide recurrence diagnostics.</p>
    <section><h2>1. Refined 3D parameter regions</h2><p>Switch views in the menu and move the m₀ slice.
    The surface uses every calculated grid node, with no smoothing to erase thin windows.</p>
    <div class="pair">{''.join(global_parts)}</div></section>
    <section><h2>2. A thin regular window inside chaos</h2>
    <p>The pink surface encloses nodes with near-zero/nonpositive λ₁ for all three starts.
    Teal marks the chaos criterion for all starts. At β=14.32 on the central transect, peak analysis additionally finds a five-peak pattern.
    Exact periodicity was not checked at every pink node: this is a numerical model of a regular window, not a proven bifurcation surface.</p>
    <div class="pair">{div(zoom_figure(*records['window_zoom']),'window-zoom')}{div(slice_slider(*records['window_zoom']),'window-slices')}</div>
    <p class="note">Local spacings: Δα=0.003, Δβ=0.004, Δm₀=0.0003. Boundary precision is limited by these spacings and the observation duration.</p></section>
    <section><h2>3. Moving along β</h2>
    <p>Top: λ₁. Middle: local maxima of x(t); separate branches indicate repeated peaks, while dense bands contain many distinct maxima.
    Bottom: the detected peak-pattern length. Missing points mean no repeating pattern was found at the chosen tolerance; they do not independently prove chaos.</p>
    <div class="pair">{div(beta_plot(lines['classic'],'classic'),'beta-classic')}{div(beta_plot(lines['reference'],'reference'),'beta-reference')}</div>
    <p>At α=15.6, m₀=−1.143 and m₁=−0.714, a broad regular window spans approximately β=28.87…30.04 on this transect.
    At β=32, chaotic motion returns around one outer center. Increasing β does not monotonically increase or decrease chaos.</p></section>
    <section><h2>4. Same parameters, different attractors</h2>
    <p>At α=15.6, β=32, m₀=−1.143, m₁=−0.714, symmetric starts produce two distinct single-scroll attractors.
    Both have λ₁≈0.418 after transient=1000 and observation=2400: strong numerical evidence of coexisting chaotic attractors.</p>
    {div(example_plot(items),'examples')}
    {div(basin_plot(basin),'basin')}
    <p>The lower map varies initial (x₀,y₀) with z₀=0: teal and pink select different scrolls.
    It contains 151²=22,801 initial states. The only unassigned state is the exact equilibrium (0,0,0).
    Equal colored-node counts reflect the symmetric grid and equations, not universal basin probabilities.</p></section>
    <section><h2>5. What the calculations establish, and what remains open</h2>
    <p>Over 99.9% of the candidate nodes in each general grid belong to one component under six-neighbor grid connectivity.
    At this resolution the dominant structure is one large component interrupted by thin windows.
    Isolated single nodes do not establish separate physical islands.</p>
    <p>The equations explain the special role of m₀=−1. For m₀&lt;−1 and fixed m₁&gt;−1, the origin is accompanied by
    outer equilibria (±x*,0,∓x*), with x*=(m₁−m₀)/(1+m₁)&gt;1.
    At m₀=−1, the inner segment contains a continuum of equilibria. This changes the equilibrium structure,
    but does not make the plane a proven boundary of all chaotic behavior.</p>
    <p>A jagged mesh does not prove a fractal boundary. Three starts cannot cover every basin,
    and finite-time positive exponents cannot exclude long chaotic transients. Candidates reach beyond both general scan boxes.</p>
    <p>Independent adaptive SciPy DOP853 confirms the five-peak pattern at β=14.32 and six-peak pattern at β=29 in their respective sections.
    A repeated GPU check with dt=0.0025, transient=1200 and observation=3200 preserves the category at 252 of 256 selected zoom nodes.
    This checks a stratified sample, not the complete boundary. <a href="validation_details.json">Validation record</a>.</p>
    <p>Periodic windows and changes between single and double scrolls are established Chua phenomena;
    this project measures their locations in the selected sections.
    <a href="https://people.eecs.berkeley.edu/~chua/papers/Kennedy93b.pdf">Kennedy: Chua’s circuit primer</a> ·
    <a href="https://www.scholarpedia.org/article/Chua_circuit">Leon Chua’s survey</a>.</p>
    <p><a href="findings.png">Findings PNG</a> · <a href="findings.pdf">PDF</a> ·
    <a href="analysis_summary.json">Numerical summaries</a> · <a href="../index.html">Original atlas</a></p></section>
    </main><script>window.addEventListener('load',()=>requestAnimationFrame(()=>document.querySelectorAll('.js-plotly-plot').forEach(p=>Plotly.Plots.resize(p))));</script></body></html>"""
    shutil.copyfile(ROOT/'outputs'/'plotly.min.js',OUT/'plotly.min.js')
    (OUT/'index.html').write_text(body,encoding="utf-8")
    print(json.dumps(summaries,indent=2),flush=True)
    print("Detailed report saved",flush=True)

if __name__=='__main__':build()
