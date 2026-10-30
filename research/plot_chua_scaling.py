"""Figures from saved data only; this script performs no integrations."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from chua_scaling_gpu import OUT, ROOT
from analyze_chua_scaling import verify_resume_source


def plot():
    verify_resume_source()
    scout=np.load(OUT/'scout.npz');ref=np.load(OUT/'refined.npz');p=scout['parameters']
    retained=np.zeros(len(p),bool);retained[ref['source_indices'][ref['retained_candidate']]]=True
    values=np.full(len(p),np.nan);values[ref['source_indices']]=ref['lle']
    resumed_path=OUT/'refined-clock-resume.npz'
    if resumed_path.exists():
        resumed=np.load(resumed_path)
        retained[resumed['source_indices']]=resumed['retained_candidate']
        values[resumed['source_indices']]=resumed['lle']
    fig,axes=plt.subplots(1,3,figsize=(14,4.4),sharex=True,sharey=True,layout='constrained')
    for ax,m in zip(axes,[-8/7,-2.,-30.]):
        select=abs(p[:,2]-m)<1e-8
        tuples=np.unique(p[select,:2],axis=0)
        for a,k in tuples:
            q=select&(p[:,0]==a)&(p[:,1]==k)
            survives=q&retained
            if survives.any():
                ax.scatter(a,k,c=[min(values[survives].max(),3)],vmin=0,vmax=3,cmap='viridis',s=32)
            elif np.any(q&scout['preliminary_candidate']):ax.scatter(a,k,c='#cc6a26',marker='x',s=22)
            elif np.any(q&(scout['status']==0)):ax.scatter(a,k,c='#c5cad0',s=7)
            else:ax.scatter(a,k,c='#d63a49',marker='s',s=13)
        aa=np.geomspace(9,2048,100);ax.plot(aa,10/49+5/(7*aa),'--',color='black',lw=1)
        ax.set(xscale='log',xlabel='alpha',title=f'Fixed inner slope m0 = {m:.6g}',ylim=(.02,.29))
        ax.grid(alpha=.18)
    axes[0].set_ylabel('kappa = beta / alpha²')
    sm=plt.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,3),cmap='viridis')
    fig.colorbar(sm,ax=axes,label='Largest retained finite-time LLE (original time; capped at 3)')
    fig.suptitle('Fixed-slope Chua sections: finite observations, not a chaos boundary',fontsize=14)
    fig.text(.5,-.06,'Colored: repeated positive screen. Orange ×: preliminary signal not retained. Gray: no preliminary signal.\nRed square: all three starts reached magnitude cutoff. Dashed: outer linear stability line.',ha='center',fontsize=10)
    fig.savefig(ROOT/'docs/assets/chua-scaling-sections.png',dpi=170,bbox_inches='tight');plt.close(fig)
    files=[OUT/'cpu-504.npz',OUT/('cpu-2198-aligned.npz' if (OUT/'cpu-2198-aligned.npz').exists() else 'cpu-2198.npz')]
    available=[f for f in files if f.exists()]
    if not available:return
    fig=plt.figure(figsize=(12,5*len(available)),layout='constrained')
    for row,file in enumerate(available):
        x=np.load(file);y=x['states'];t=x['time'];p=x['parameters'];section=x['section_states']
        ax=fig.add_subplot(len(available),2,2*row+1,projection='3d')
        skip=max(1,len(y)//30000);v=y[::skip]
        ax.plot(v[:,0],v[:,1],v[:,2],lw=.2,color='#3568af',alpha=.65)
        ax.set(xlabel='x',ylabel='y',zlabel='Z = z / alpha',title=f'alpha={p[0]:g}, kappa={p[1]:g}, m0={p[2]:.6g}')
        ax=fig.add_subplot(len(available),2,2*row+2)
        if len(section):ax.scatter(section[:,0],section[:,2]/np.sqrt(p[1]),s=3,c='#3568af',alpha=.5)
        ax.set(xlabel='x at y=0, upward crossing',ylabel='R = Z / sqrt(kappa)',title=f'CPU Poincare section, {len(section):,} returns')
        ax.grid(alpha=.2)
    fig.suptitle('Independent CPU trajectories in linearly scaled physical coordinates',fontsize=14)
    fig.savefig(ROOT/'docs/assets/chua-scaling-trajectories.png',dpi=170);plt.close(fig)


if __name__=='__main__':plot()
