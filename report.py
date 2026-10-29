"""Offline Plotly report and publication-friendly Matplotlib PNG/PDF figures."""
from __future__ import annotations

import html
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs

COLORS=["#243448","#5c7793","#acbfd1","#f6be57","#27bbba","#ed6080"]


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def style(fig,title):
    fig.update_layout(template="plotly_dark",paper_bgcolor="#101d30",plot_bgcolor="#101d30",
                      font=dict(family="Arial, sans-serif",size=13,color="#dce7f3"),
                      title=dict(text=title,font=dict(size=19)),margin=dict(l=15,r=15,b=15,t=65),
                      height=580,legend=dict(orientation="h",y=-0.03))
    return fig


def region_figure(data,meta):
    p=data["parameters"]
    candidate=data["candidate_any"]
    everywhere=data["candidate_all"]
    classes=data["classes"]
    selected=np.flatnonzero(candidate)
    fig=go.Figure()
    # Iso-surfaces interpolate a binary mask over the computed grid, preserving holes.
    # They are clipped at the scanned box; do not extrapolate a closed envelope.
    for mask,color,name in [(candidate,"#27bbba","At least one seed"),
                            (everywhere,"#f6be57","All three seeds")]:
        fig.add_trace(go.Isosurface(x=p[:,0],y=p[:,1],z=p[:,2],value=mask.astype(float),
            isomin=.5,isomax=.5,surface_count=1,caps=dict(x_show=False,y_show=False,z_show=False),
            colorscale=[[0,color],[1,color]],showscale=False,opacity=.6,
            name=name,showlegend=True,visible=True if name=="At least one seed" else False,
            hoverinfo="skip"))
    # Deterministic display subsampling only; every point remains in the NPZ dataset.
    if len(selected)>35000:
        selected=selected[np.linspace(0,len(selected)-1,35000,dtype=int)]
    lle=np.max(np.where(classes==4,data["lle"],-np.inf),axis=1)
    fig.add_trace(go.Scatter3d(x=p[selected,0],y=p[selected,1],z=p[selected,2],mode="markers",
        marker=dict(size=2,opacity=.65,color=lle[selected],colorscale="Viridis",
                    colorbar=dict(title="λ₁",thickness=12),cmin=0),
        customdata=np.column_stack([lle[selected],np.sum(classes[selected]==4,axis=1)]),
        hovertemplate="α=%{x:.4f}<br>β=%{y:.4f}<br>m₀=%{z:.5f}<br>λ₁=%{customdata[0]:.4f}"
                      "<br>Candidate seeds: %{customdata[1]}/3<extra></extra>",
        name="Computed candidate nodes",visible=False))
    uncertain=np.flatnonzero((classes==3).any(axis=1)&~candidate)
    if len(uncertain)>15000:
        uncertain=uncertain[np.linspace(0,len(uncertain)-1,15000,dtype=int)]
    fig.add_trace(go.Scatter3d(x=p[uncertain,0],y=p[uncertain,1],z=p[uncertain,2],mode="markers",
        marker=dict(size=2,color="#f6be57",opacity=.4),name="Uncertain nodes",visible=False))
    c=meta["configuration"]["parameters"]
    fig.add_trace(go.Scatter3d(x=[c[0]],y=[c[1]],z=[c[2]],mode="markers",
        marker=dict(size=7,color="#ffffff",symbol="diamond",line=dict(color="#0a1423",width=2)),
        name="Reference point",hovertemplate="Reference<br>α=%{x}<br>β=%{y}<br>m₀=%{z}<extra></extra>"))
    buttons=[]
    for label,vis in [("Any seed",[True,False,False,False,True]),
                      ("All seeds",[False,True,False,False,True]),
                      ("LLE points",[False,False,True,False,True]),
                      ("Uncertain",[True,False,False,True,True])]:
        buttons.append(dict(label=label,method="update",args=[{"visible":vis}]))
    fig.update_layout(scene=dict(xaxis=dict(title="α",range=[float(p[:,0].min()),float(p[:,0].max())]),
                                 yaxis=dict(title="β",range=[float(p[:,1].min()),float(p[:,1].max())]),
                                 zaxis=dict(title="m₀",range=[float(p[:,2].min()),float(p[:,2].max())]),
                                 aspectmode="cube",camera=dict(eye=dict(x=1.8,y=1.8,z=1.2),
                                 center=dict(x=0,y=0,z=0),projection=dict(type="orthographic"))),
                      updatemenus=[dict(type="buttons",direction="right",buttons=buttons,x=0,y=1.01,
                                       xanchor="left",yanchor="bottom",
                                       bgcolor="#cbd5e1",font=dict(color="#14273c"))])
    style(fig,f"Parameter region · m₁ = {c[3]:.6f}")
    fig.update_layout(margin=dict(l=15,r=15,b=55,t=110))
    return fig


def trajectory_figure(data,meta):
    xyz=data["xyz"]
    fig=go.Figure(go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode="lines",
        line=dict(width=2,color=data["t"],colorscale="Tealgrn"),hoverinfo="skip",name="Trajectory"))
    fig.update_layout(scene=dict(xaxis_title="x",yaxis_title="y",zaxis_title="z",aspectmode="manual",
                                 aspectratio=dict(x=1,y=.65,z=1.25),
                                 camera=dict(eye=dict(x=1.8,y=1.8,z=1.2),projection=dict(type="orthographic"))))
    return style(fig,f"Phase space · λ₁ = {meta['lle']:.4f}")


def slice_figure(data,meta):
    n=meta["grid"]
    classes=data["classes"]
    # Existence screen = strongest evidence among tested seeds; all escaped = 0.
    aggregate=np.full(len(classes),3,dtype=int)
    aggregate[(classes==0).all(axis=1)]=0
    aggregate[(classes==1).all(axis=1)]=1
    aggregate[np.isin(classes,[0,1,2]).all(axis=1)&(classes==2).any(axis=1)]=2
    aggregate[(classes==4).any(axis=1)]=4
    aggregate[(classes==5).any(axis=1)&~(classes==4).any(axis=1)]=5
    volume=aggregate.reshape((n,n,n))
    indices=[n//4,n//2,3*n//4]
    fig=make_subplots(rows=1,cols=3,subplot_titles=[f"m₀ = {data['m0'][k]:.4f}" for k in indices],
                      horizontal_spacing=.06)
    scale=[]
    for k,color in enumerate(COLORS):
        scale.extend([[k/6,color],[(k+1)/6,color]])
    for j,k in enumerate(indices):
        fig.add_trace(go.Heatmap(x=data["alpha"],y=data["beta"],z=volume[:,:,k].T,
            zmin=-.5,zmax=5.5,colorscale=scale,showscale=j==2,
            colorbar=dict(thickness=12,tickvals=list(range(6)),
                          ticktext=["Escaped","Equilibrium","Near-zero / negative","Uncertain","Chaos candidate","Failure"]),
            hovertemplate="α=%{x:.3f}<br>β=%{y:.3f}<br>Class=%{z}<extra></extra>"),row=1,col=j+1)
        fig.update_xaxes(title_text="α",row=1,col=j+1)
        fig.update_yaxes(title_text="β",row=1,col=j+1)
    style(fig,"Slices through the computed grid")
    fig.update_layout(height=380,margin=dict(l=45,r=150,b=45,t=70))
    return fig


def static_overview(records,directory):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                         "figure.facecolor":"#f8fbfe","axes.facecolor":"#f8fbfe",
                         "axes.labelcolor":"#243448","text.color":"#152740"})
    fig=plt.figure(figsize=(15,6*len(records)),layout="constrained")
    for row,(name,scan,sm,traj,tm) in enumerate(records):
        ax=fig.add_subplot(len(records),2,2*row+1,projection="3d")
        xyz=traj["xyz"]
        ax.plot(*xyz.T,lw=.24,color="#007f84",alpha=.85)
        ax.set(xlabel="x",ylabel="y",zlabel="z",title=f"{sm['configuration']['label']}  |  λ₁ = {tm['lle']:.3f}")
        ax.view_init(22,50)
        ax.set_box_aspect((1,.65,1))
        ax=fig.add_subplot(len(records),2,2*row+2,projection="3d",computed_zorder=False)
        idx=np.flatnonzero(scan["candidate_any"])
        # Plot only a deterministic display subset for legibility.
        if len(idx)>12000:
            idx=idx[np.linspace(0,len(idx)-1,12000,dtype=int)]
        p=scan["parameters"][idx]
        lle=np.max(np.where(scan["classes"][idx]==4,scan["lle"][idx],-np.inf),axis=1)
        dots=ax.scatter(*p[:,:3].T,c=lle,cmap="viridis",s=2,alpha=.5,rasterized=True)
        c=sm["configuration"]["parameters"]
        ax.scatter(*c[:3],c="#d23e5c",s=65,marker="*",edgecolor="white",linewidth=.7,zorder=100)
        ax.set(xlabel="α",ylabel="β",zlabel="m₀",title=f"{sm['grid']}³ grid · m₁ = {c[3]:.6f} · {sm['candidate_fraction_of_nodes']:.1%} candidate nodes")
        ax.view_init(23,43)
        fig.colorbar(dots,ax=ax,shrink=.55,pad=.02,label="Largest finite-time Lyapunov exponent")
    fig.suptitle("CHUA  /  GPU PARAMETER ATLAS\nBounded positive-Lyapunov candidates from three initial conditions",fontsize=19)
    fig.savefig(directory/"overview.png",dpi=170)
    fig.savefig(directory/"overview.pdf")
    plt.close(fig)
    # Unsmoothed categorical slices show periodic windows and uncertain cells.
    for name,scan,sm,traj,tm in records:
        n=sm["grid"]
        a=scan["classes"]
        agg=np.full(len(a),3)
        agg[(a==0).all(axis=1)]=0
        agg[(a==1).all(axis=1)]=1
        agg[np.isin(a,[0,1,2]).all(axis=1)&(a==2).any(axis=1)]=2
        agg[(a==4).any(axis=1)]=4
        agg[(a==5).any(axis=1)&~(a==4).any(axis=1)]=5
        fig,axes=plt.subplots(1,3,figsize=(15,4.5),layout="constrained")
        for ax,k in zip(axes,[n//4,n//2,3*n//4]):
            im=ax.pcolormesh(scan["alpha"],scan["beta"],agg.reshape((n,)*3)[:,:,k].T,
                            cmap=ListedColormap(COLORS),norm=BoundaryNorm(np.arange(-.5,6),6),shading="nearest")
            ax.set(xlabel="α",ylabel="β",title=f"m₀ = {scan['m0'][k]:.4f}")
        cb=fig.colorbar(im,ax=axes,ticks=range(6),shrink=.9)
        cb.ax.set_yticklabels(["Escaped","Equilibrium","Near-zero / negative","Uncertain","Chaos candidate","Failure"])
        fig.suptitle(f"{sm['configuration']['label']} — finite-time classification")
        fig.savefig(directory/f"{name}_slices.png",dpi=170)
        plt.close(fig)


def build_report(directory):
    directory=Path(directory)
    records=[]
    for file in sorted(directory.glob("*_scan.npz")):
        name=file.stem[:-5]
        trajectory=directory/f"{name}_trajectory.npz"
        if not trajectory.exists():
            continue
        records.append((name,dict(np.load(file)),read_json(file.with_suffix(".json")),
                        dict(np.load(trajectory)),read_json(trajectory.with_suffix(".json"))))
    if not records:
        raise RuntimeError("No complete scan + trajectory results; run 'python run.py all' first")
    static_overview(records,directory)
    (directory/"plotly.min.js").write_text(get_plotlyjs(),encoding="utf-8")
    sections=[]
    for name,scan,sm,traj,tm in records:
        cfg=sm["configuration"]
        c=cfg["parameters"]
        settings=sm["settings"]
        plots=[trajectory_figure(traj,tm),region_figure(scan,sm),slice_figure(scan,sm)]
        divs=[f.to_html(full_html=False,include_plotlyjs=False,
                        config={"responsive":True,"displaylogo":False,"toImageButtonOptions":{"format":"png","scale":2}})
              for f in plots]
        faces=", ".join(k for k,v in sm["boundary_faces_with_candidates"].items() if v)
        sections.append(f"""<section><div class="eyebrow">{html.escape(name)} / {sm['grid']} × {sm['grid']} × {sm['grid']}</div>
        <h2>{html.escape(cfg['label'])}</h2>
        <p>Center: α = {c[0]:.6g}, β = {c[1]:.6g}, m₀ = {c[2]:.6g}; fixed m₁ = {c[3]:.6g}.
        <a href="{html.escape(cfg['source'],quote=True)}">Parameter source ↗</a></p>
        <div class="metrics"><div><b>{sm['candidate_nodes_any_seed']:,}</b><span>chaos candidate nodes</span></div>
        <div><b>{sm['candidate_nodes_all_seeds']:,}</b><span>candidates for all three starts</span></div>
        <div><b>{sm['timing']['wall_seconds']:.1f} s</b><span>GPU scan time</span></div>
        <div><b>{tm['lle']:.4f}</b><span>reference-point λ₁</span></div></div>
        <div class="pair">{divs[0]}{divs[1]}</div><div class="slices">{divs[2]}</div>
        <p class="note">Candidates touch these box faces: {html.escape(faces or 'none')}.
        The region is clipped at the search limits; its complete outer boundary is unknown.
        The surface interpolates a binary grid mask. Points and slices display the computed nodes.
        Point color is the maximum λ₁ among starts passing the candidate criterion.</p>
        <details><summary>Settings and data</summary><p>dt = {settings['dt']}; transient = {settings['transient']};
        observation = {settings['duration']}; four windows of {settings['duration']/4}; λ₁ &gt; {settings['threshold']},
        every window &gt; {settings['threshold']/2}. Initial states: {html.escape(str(settings['seeds']))}.</p>
        <p><a href="{name}_scan.npz">Complete NPZ data</a> · <a href="{name}_scan.json">JSON metadata</a> ·
        <a href="{name}_trajectory.npz">Trajectory NPZ</a> · <a href="{name}_slices.png">Slices PNG</a></p></details></section>""")
    total=sum(r[2]["grid_nodes"]*len(r[2]["settings"]["seeds"]) for r in records)
    device=records[0][2]["device"]
    validation=""
    vf=directory/"validation.json"
    if vf.exists():
        validation=f'<p><a href="validation.json">Independent checks and refinement →</a></p>'
    followups=[]
    for folder,label in [("detailed","Refined atlas: thin windows and coexisting attractors"),
                         ("extent","Outer extent: a wide point cloud"),
                         ("expanded","Four times farther: expanded adaptive GPU survey")]:
        if (directory/folder/"index.html").exists():
            followups.append(f'<a href="{folder}/index.html">{label} →</a>')
    navigation='<p>'+"<br>".join(followups)+'</p>' if followups else ''
    text=f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>Chua · GPU parameter atlas</title><script src="plotly.min.js"></script><style>
    *{{box-sizing:border-box}}body{{margin:0;background:#091321;color:#dce7f3;font:16px/1.65 system-ui,sans-serif}}
    main{{max-width:1550px;margin:auto;padding:48px 36px}}h1{{font-size:clamp(34px,4vw,64px);line-height:1.1;margin:14px 0 24px}}
    h2{{font-size:29px;margin:10px 0}}a{{color:#69d4cf}}p{{max-width:1150px;color:#aebfd2}}.eyebrow{{color:#65cec6;letter-spacing:.16em;text-transform:uppercase;font-size:12px}}
    .lead{{font-size:20px;max-width:1050px}}.formula{{padding:22px 28px;border:1px solid #304259;border-radius:12px;font:18px/1.85 ui-monospace,monospace;background:#101d30}}
    section{{border-top:1px solid #283b51;margin-top:52px;padding-top:35px}}.pair{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}
    .pair>div,.slices{{min-width:0;border-radius:14px;overflow:hidden;background:#101d30}}.slices{{margin-top:18px}}
    .metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:16px;margin:25px 0}}.metrics div{{background:#122238;border-radius:10px;padding:18px 22px}}
    .metrics b{{display:block;font-size:29px;color:#edf7ff}}.metrics span{{font-size:13px;color:#a9bed2}}
    .note{{font-size:13px}}details{{background:#101d30;border-radius:10px;padding:16px 22px}}summary{{cursor:pointer}}footer{{margin-top:40px;color:#8fa5bd;font-size:13px}}
    @media(max-width:1000px){{main{{padding:24px 14px}}.pair{{grid-template-columns:1fr}}.metrics{{grid-template-columns:1fr 1fr}}}}
    </style></head><body><main><div class="eyebrow">Nonlinear dynamics / Python + OpenCL / float64</div>
    <h1>Chua: mapping the chaotic regimes</h1>{navigation}<p class="lead">Double-scroll trajectories and a three-dimensional parameter map (α, β, m₀).
    {total:,} trajectories computed on {html.escape(device['name'])}. Rotate and zoom the plots.</p>
    <div class="formula">dx/dt = α [y − x − h(x)]<br>dy/dt = x − y + z<br>dz/dt = −βy<br>
    h(x) = m₁x + ½(m₀ − m₁)(|x + 1| − |x − 1|)</div>
    <p>Colored regions are numerical chaos candidates: the trajectory stays within |x|, |y|, |z| ≤ 100,
    keeps moving and has a positive largest finite-time Lyapunov exponent in four post-transient windows.
    A finite observation and a small set of initial conditions cannot prove attractor existence or absence.
    The screen includes both single-scroll and double-scroll chaotic regimes.</p>
    <p>“Any seed” requires one successful start; “All seeds” requires all three.
    “LLE points” colors nodes by λ₁; “Uncertain” adds unresolved nodes.
    A negative or near-zero exponent alone does not establish periodicity.</p>
    {''.join(sections)}<footer><a href="overview.png">Overview PNG</a> · <a href="overview.pdf">Overview PDF</a>
    {validation}<p>RK4 jointly advances the state and variational equation; the tangent vector is normalized every ten steps.
    No fast-math flags or silent CPU fallback. Boundary accuracy is limited by the grid, integrator and observation time.</p></footer>
    </main><script>window.addEventListener('load',()=>{{
      requestAnimationFrame(()=>document.querySelectorAll('.js-plotly-plot').forEach(p=>Plotly.Plots.resize(p)));
    }});</script></body></html>"""
    (directory/"index.html").write_text(text,encoding="utf-8")
    print("Report saved: index.html in the selected output directory",flush=True)


if __name__=="__main__":
    import sys
    build_report(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/"outputs")
