"""Render the saved exact-motor survey; performs no new integration."""
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


def main(gif=False):
    folder=ROOT/'results/motor-pair';assets=ROOT/'docs/assets';assets.mkdir(exist_ok=True)
    scout=np.load(folder/'gpu-scout.npz');tight=np.load(folder/'gpu-refined-complete.npz')
    cpu=np.load(folder/'gpu-cpu-444-extended.npz');report=json.loads((folder/'gpu-cpu-444-extended.json').read_text())
    if report['status']!='complete':raise RuntimeError('CPU continuation must finish before final rendering')
    p=cpu['parameters'];time=cpu['time'];x=cpu['balanced_states'];physical=np.c_[p[1]*(x[:,0]+x[:,1]),x[:,1:]]
    selected=time>=max(100.,time[-1]-2000.);trajectory=physical[selected];times=time[selected]
    stride=max(1,len(trajectory)//18000);trajectory=trajectory[::stride];times=times[::stride]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    fig=plt.figure(figsize=(13,10),layout='constrained')
    ax=fig.add_subplot(221);valid=scout['status']==0;pre=valid&(scout['lle']>.01)&np.all(scout['lle_windows']>.005,axis=1)
    ax.scatter(scout['parameters'][valid,0],scout['parameters'][valid,1],s=8,c='#c7d6de',label='Completed scout')
    ax.scatter(scout['parameters'][~valid,0],scout['parameters'][~valid,1],s=13,c='#8b8b8b',marker='x',label='Step budget: unresolved')
    ax.scatter(scout['parameters'][pre,0],scout['parameters'][pre,1],s=45,c='#c04b24',label='Retained after tighter continuation')
    ax.set(xscale='log',yscale='log',xlabel=r'$\rho$',ylabel=r'$r$',title='2,880 trajectories; 252 scout lanes censored')
    ax.legend(fontsize=8,loc='upper left');ax.grid(alpha=.15)
    ax=fig.add_subplot(222)
    for ident,pa,w in zip(tight['source_ids'],tight['parameters'],tight['lle_windows']):
        ax.plot(np.arange(1,5),w,'o-',ms=3,lw=1,label=f'{ident}: r={pa[1]:.1f}')
    ax.axhline(0,c='black',lw=.7);ax.set(xlabel='Observation window (each 200 / rho in tau)',ylabel='Leading FTLE per original tau',
        title='All eight tightened trajectories remain positive',xticks=[1,2,3,4]);ax.legend(ncol=2,fontsize=7);ax.grid(alpha=.15)
    ax=fig.add_subplot(223,projection='3d');ax.scatter(*trajectory.T,s=.25,c='#286379',alpha=.4,depthshade=False)
    ax.set(xlabel='h',ylabel='Y',zlabel='Z',title='Exact physical state coordinates: representative row 444')
    ax.view_init(elev=24,azim=32)
    ax=fig.add_subplot(224);windows=report['windows'];durations=np.array([v['duration'] for v in windows]);t=np.cumsum(durations)
    values=np.array([v['spectrum'] for v in windows]);cumulative=np.cumsum(values*durations[:,None],axis=0)/t[:,None]
    flow=np.cumsum(np.array([v['flow_exponent'] for v in windows])*durations)/t
    for j in range(3):ax.plot(t,cumulative[:,j],label=f'QR exponent {j+1}')
    ax.plot(t,flow,'--',c='black',lw=1,label='Exact flow-direction finite-time rate');ax.axhline(0,c='gray',lw=.7)
    ax.set(xlabel='Observed original tau (after 100 discarded)',ylabel='Finite-time rate',title='Independent CPU full-spectrum convergence check')
    ax.legend(fontsize=8);ax.grid(alpha=.15)
    fig.suptitle('Collective molecular motors: a numerical chaotic candidate in an existing exact model\nFinite-time evidence; scientific priority and asymptotic attraction remain unproved',fontsize=14)
    fig.savefig(assets/'motor_pair_gpu_survey.png',dpi=170);plt.close(fig)
    fig=plt.figure(figsize=(14,8),layout='constrained')
    ax=fig.add_subplot(231,projection='3d');ax.scatter(*trajectory.T,s=.25,c='#286379',alpha=.4,depthshade=False);ax.view_init(24,32)
    ax.set(xlabel='h',ylabel='Y',zlabel='Z',title='Full 3D state')
    for number,(i,j) in enumerate([(0,1),(0,2),(1,2)],2):
        ax=fig.add_subplot(2,3,number);ax.scatter(trajectory[:,i],trajectory[:,j],s=.25,c='#286379',alpha=.4)
        ax.set(xlabel=['h','Y','Z'][i],ylabel=['h','Y','Z'][j],title='Physical coordinate projection');ax.grid(alpha=.15)
    sec=cpu['section_states'];sec_time=cpu['section_time'];mask=sec_time>=500.;sec=sec[mask];sec_time=sec_time[mask]
    ax=fig.add_subplot(235);dots=ax.scatter(sec[:,1],sec[:,2],c=sec_time,s=30,cmap='viridis');fig.colorbar(dots,ax=ax,label='Continuation time tau')
    ax.set(xlabel='Y on h=0',ylabel='Z on h=0',title=f'Upward Poincare section: {len(sec)} new returns');ax.grid(alpha=.15)
    ax=fig.add_subplot(236)
    if len(sec)>1:ax.plot(sec[:-1,2],sec[1:,2],'o',ms=4,c='#c04b24')
    ax.set(xlabel='Z at return n',ylabel='Z at return n+1',title='Sampled return pairs');ax.grid(alpha=.15)
    fig.suptitle(f'Existing motor-pair model, rho={p[0]:.8f}, r={p[1]:.6f}\nDifferent views are geometric diagnostics, not evidence of a new attractor type',fontsize=14)
    fig.savefig(assets/'motor_pair_candidate_geometry.png',dpi=170);plt.close(fig)
    if gif:
        fig=plt.figure(figsize=(5.5,5.5));ax=fig.add_subplot(111,projection='3d')
        ax.scatter(*trajectory[::2].T,s=.4,c='#286379',alpha=.5,depthshade=False);ax.set(xlabel='h',ylabel='Y',zlabel='Z')
        fig.suptitle('Existing motor-pair model: numerical candidate\nPhysical coordinates; finite-time point samples',fontsize=10)
        def update(frame):ax.view_init(22+16*np.sin(frame*2*np.pi/60),frame*6);return []
        FuncAnimation(fig,update,frames=60,interval=80).save(assets/'motor_pair_candidate_rotation.gif',writer=PillowWriter(fps=12),dpi=95)
        plt.close(fig)
    data=dict(scout=dict(lanes=len(valid),completed=int(valid.sum()),step_budget_incomplete=int(np.sum(scout['status']==3)),
        preliminary_count=int(pre.sum())),tightened=dict(completed=int(np.sum(tight['status']==0)),
        minimum_window=float(np.min(tight['lle_windows'])),maximum_window=float(np.max(tight['lle_windows'])),
        rows=[dict(source_id=int(i),rho=float(pa[0]),r=float(pa[1]),mean_lle=float(le),windows=w.tolist(),
            mean_divergence=float(div)) for i,pa,le,w,div in zip(tight['source_ids'],tight['parameters'],tight['lle'],tight['lle_windows'],tight['mean_divergence_original'])]),
        cpu=dict(status=report['status'],cumulative_spectrum=report['cumulative_spectrum'],extension_spectrum=report['extension_spectrum'],
            observed_tau=report['cumulative_observed_tau'],section_returns=report['observed_section_hits'],
            extension_section_returns=int(len(sec)),sum_trace_error=report['cumulative_sum_trace_error'],
            flow_exponent=report['cumulative_flow_exponent'],telescopic_flow_exponent=report['cumulative_telescopic_flow_exponent']),
        source_hash_policy='Exact NPZ bytes; canonical LF text files',source_sha256={},
        limitations=['252 scout trajectories remain step-budget censored; they are not negative results.',
            'Only one representative has an independent long full-spectrum CPU check.',
            'Finite-time bounded persistence does not prove asymptotic attraction, scientific priority or new topology.',
            'The physical parameter rectangle is exploratory, not an apparatus calibration.',
            'Incomplete historical NPZ mean_divergence_original fields are full-request-normalized partial integrals, not valid means; use only status 0.'])
    for name in ['gpu-scout.npz','gpu-refined.npz','gpu-refined-complete.npz','gpu-cpu-444.json','gpu-cpu-444-extended.json','gpu-cpu-444-extended.npz']:
        b=(folder/name).read_bytes();data['source_sha256'][name]=hashlib.sha256(b if name.endswith('.npz') else b.replace(b'\r\n',b'\n')).hexdigest()
    (folder/'gpu-summary.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--gif',action='store_true');main(parser.parse_args().gif)
