"""Postprocess one completed full-state N=3 CPU audit; no integration."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT=Path(__file__).resolve().parents[1]


def main(gif=False,summary_only=False):
    stem=ROOT/'results/motor-rows/cpu-candidate-N3-783';report=json.loads(stem.with_suffix('.json').read_text())
    if report['status']!='complete':raise RuntimeError('Completed CPU audit required')
    data=np.load(stem.with_suffix('.npz'));time=data['time'];state=data['states'];mask=time>report['settings']['warm']
    t=time[mask];z=state[mask];c=z[:,:3];s=z[:,3:6];X=np.c_[z[:,6:],-z[:,6:].sum(axis=1)]
    sec=data['section_states'];st=data['section_time'];pairs=[(0,1),(0,2),(1,2)]
    distances=np.stack([np.sqrt((c[:,i]-c[:,j])**2+(s[:,i]-s[:,j])**2+(X[:,i]-X[:,j])**2) for i,j in pairs],axis=1)
    closest=np.argmin(distances,axis=1);changes=np.r_[0,np.flatnonzero(np.diff(closest)!=0)+1,len(closest)]
    windows=report['windows'];rates=np.array([w['spectrum'] for w in windows]);durations=np.array([w['duration'] for w in windows]);elapsed=np.cumsum(durations)
    cumulative=np.cumsum(rates*durations[:,None],axis=0)/elapsed[:,None]
    blocks=[]
    for start in range(0,len(windows),5):
        chosen=windows[start:start+5];length=sum(w['duration'] for w in chosen)
        blocks.append(dict(duration=length,spectrum=np.average([w['spectrum'] for w in chosen],weights=[w['duration'] for w in chosen],axis=0).tolist()))
    # A fixed physical-coordinate normal, not a PCA-remapped visualization.
    covariance=np.cov(s,rowvar=False,ddof=1);normal=np.ones(3)/np.sqrt(3.)
    mean_s=np.mean(s,axis=1);normal_variance=float(normal@covariance@normal)
    summary=dict(N=3,dimension=8,parameters=report['parameters'],source_row=783,
        observed_tau=report['completed_observation_tau'],spectrum=report['cumulative_spectrum'],five_hundred_tau_blocks=blocks,
        observed_upward_X1_returns=len(sec),trace_error=report['sum_trace_error'],flow_rate=report['flow_exponent'],
        telescopic_flow_rate=report['telescopic_flow_exponent'],sampled_domain=report['sampled_domain'],
        row_pair_distances=report['row_pair_distances'],
        harmonic_projection_geometry=dict(coordinates=['s1','s2','s3'],
            sample_selection='Saved CPU states with original tau > 200; no new integration',
            sample_count=len(s),covariance_and_std_normalization='Sample covariance and standard deviation, ddof=1',
            individual_std=np.sqrt(np.diag(covariance)).tolist(),
            mean_s_range=[float(np.min(mean_s)),float(np.max(mean_s))],mean_s_std=float(np.std(mean_s,ddof=1)),
            covariance=covariance.tolist(),covariance_eigenvalues_ascending=np.linalg.eigvalsh(covariance).tolist(),
            fixed_normal=normal.tolist(),variance_along_fixed_normal=normal_variance,
            total_variance=float(np.trace(covariance)),fixed_normal_variance_fraction=normal_variance/float(np.trace(covariance)),
            qualification='The small variance along (1,1,1)/sqrt(3) explains flattening of the physical harmonic projection. The common mode is not constant. This is not a dimensional reduction, a PCA-remapped figure, or evidence of new topology.'),
        clustering=dict(metric='Euclidean distance in normalized (c_i,s_i,X_i), sampled about every .2 tau',
            pair_labels=[[i+1,j+1] for i,j in pairs],nearest_pair_fraction=[float(np.mean(closest==j)) for j in range(3)],
            nearest_pair_switches=int(len(changes)-2),
            fractions_below_threshold={str(th):np.mean(distances<th,axis=0).tolist() for th in [1e-3,1e-2,.1]},
            qualification='Nearest pair changes are geometric diagnostics, not proof of an invariant synchronized cluster or a symbolic dynamics.'),
        source_npz_exact_sha256=hashlib.sha256(stem.with_suffix('.npz').read_bytes()).hexdigest(),
        source_metadata_canonical_lf_sha256=hashlib.sha256(stem.with_suffix('.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        caveat='Eight-dimensional exact existing model. Three-dimensional pictures are projections. No claim of scientific priority, topological novelty or a proved asymptotic attractor.')
    stem.with_name(stem.name+'-summary').with_suffix('.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    if summary_only:return
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    fig=plt.figure(figsize=(14,9),layout='constrained');fig.suptitle('Three motor rows: numerical candidate in an existing 8D physical model\nAll 3D views are projections; normalized variables, no apparatus-specific units',fontsize=14)
    for panel,(values,labels) in enumerate([(np.c_[X[:,:2],s[:,0]],['X1','X2','s1']),
            (np.c_[c[:,0],s[:,0],X[:,0]],['c1','s1','X1']),(s,['s1','s2','s3'])],1):
        ax=fig.add_subplot(2,3,panel,projection='3d');ax.scatter(*values.T,s=.25,c='#296e87',alpha=.4,depthshade=False)
        ax.set(xlabel=labels[0],ylabel=labels[1],zlabel=labels[2],title='Projection of the full 8D trajectory');ax.view_init(24,35)
    ax=fig.add_subplot(234);pts=ax.scatter(sec[:,7],sec[:,0],c=st,s=4,cmap='viridis');fig.colorbar(pts,ax=ax,label='tau')
    ax.set(xlabel='X2 on X1=0',ylabel='c1 on X1=0',title=f'Upward section: {len(sec)} returns');ax.grid(alpha=.15)
    ax=fig.add_subplot(235)
    for j in range(8):ax.plot(elapsed,cumulative[:,j],lw=1,label=f'{j+1}')
    ax.axhline(0,color='black',lw=.6);ax.set(xlabel='Observed tau after 200 discarded',ylabel='Full QR finite-time exponent',title='One positive, one near zero, six negative')
    ax.legend(ncol=4,fontsize=8);ax.grid(alpha=.15)
    ax=fig.add_subplot(236);late=t>=t[-1]-100
    for j in range(3):ax.plot(t[late]-t[late][0],X[late,j],lw=.7,label=f'X{j+1}')
    ax.set(xlabel='Final 100 tau',ylabel='Normalized row shift',title='Full-state row motions');ax.legend(fontsize=8);ax.grid(alpha=.15)
    fig.savefig(stem.with_suffix('.png'),dpi=170);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained');fig.suptitle('Row similarity and sampled Poincare diagnostics; no fixed cluster reduction imposed',fontsize=13)
    ax=axs[0,0];zoom=t>=t[-1]-50
    for j,(i,k) in enumerate(pairs):ax.plot(t[zoom],distances[zoom,j],lw=.7,label=f'Rows {i+1},{k+1}')
    ax.set(xlabel='tau',ylabel='Normalized row-state distance',title='All pairwise distances: final 50 tau');ax.legend(fontsize=8)
    ax=axs[0,1];ax.step(t[zoom],closest[zoom],where='post',lw=.7);ax.set(xlabel='tau',ylabel='Currently closest pair',yticks=[0,1,2],yticklabels=['1,2','1,3','2,3'],title='Closest-pair identity: final 50 tau')
    ax=axs[1,0];ax.scatter(sec[:-1,7],sec[1:,7],s=5,c='#b44b2a');ax.set(xlabel='X2 at return n',ylabel='X2 at return n+1',title='Return pairs: projection, not a scalar map proof')
    ax=axs[1,1]
    for j in range(3):ax.plot(t,c[:,j]**2+s[:,j]**2,lw=.35,label=f'Row{j+1}')
    ax.axhline(2,c='black',ls='--',lw=.8);ax.set(xlabel='tau',ylabel='c_i squared + s_i squared',title='Invariant density disks (bound 2)');ax.legend(fontsize=8)
    for ax in axs.ravel():ax.grid(alpha=.15)
    fig.savefig(stem.with_name(stem.name+'-diagnostics').with_suffix('.png'),dpi=170);plt.close(fig)
    if gif:
        fig=plt.figure(figsize=(5.5,5.5));ax=fig.add_subplot(111,projection='3d')
        ax.scatter(X[:,0],X[:,1],s[:,0],s=.5,c='#296e87',alpha=.45,depthshade=False);ax.set(xlabel='X1',ylabel='X2',zlabel='s1')
        fig.suptitle('Three motor rows: projection of an 8D candidate\nNormalized row shifts and motor harmonic',fontsize=10)
        def update(frame):ax.view_init(22+16*np.sin(frame*2*np.pi/60),frame*6);return []
        FuncAnimation(fig,update,frames=60,interval=80).save(stem.with_suffix('.gif'),writer=PillowWriter(fps=12),dpi=95);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--gif',action='store_true');parser.add_argument('--summary-only',action='store_true')
    args=parser.parse_args();main(args.gif,args.summary_only)
