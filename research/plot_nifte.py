"""Render saved NIFTE controls and survey diagnostics; no new integrations."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nifte import ROOT


def main():
    plt.rcParams.update({'font.size':10,'axes.titlesize':11})
    fig=plt.figure(figsize=(13,10),layout='constrained');gs=fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[0,0],projection='3d')
    for name,label,color,style in [('source-5d-cycle','5D stable','#2166ac','-'),('source-9d-cycle','9D stable','#1b9e77','-'),('source-5d-unstable-cycle','5D unstable','#d95f02','--')]:
        a=np.load(ROOT/(name+'.npz'));s=a['states'];p=a['parameters'];m=p[-2]/(621*9.81)*1000
        ax.plot(s[:,1]*m,s[:,2]*m,s[:,0]*p[-2]/1000,label=label,color=color,ls=style,lw=1.25)
    ax.set(xlabel='Displacer displacement (mm)',ylabel='Power displacement (mm)',zlabel='Pressure deviation (kPa)',title='Published periodic regimes, independently recovered\nThree physical coordinates; 5D/9D state projections')
    ax.legend(loc='upper left',fontsize=9);ax.view_init(elev=24,azim=40)
    a=np.load(ROOT/'focus-5d.npz');lf=a['Lf_grid'];g=a['gain_grid'];amp=(a['stats'][:,3]-a['stats'][:,2]).reshape(len(lf),len(g),2)
    diff=np.abs(amp[:,:,0]-amp[:,:,1])/(np.max(amp,axis=2)+1e-12)
    ax=fig.add_subplot(gs[0,1]);im=ax.pcolormesh(lf/1e6,g,diff.T,vmin=0,vmax=1,cmap='viridis',shading='nearest')
    ax.plot(.286255183178673,302.871676391795,'r*',ms=13,label='Two marginal linear modes')
    ax.set(xlabel='Feedback inertance Lf (10⁶ kg m⁻⁴)',ylabel='Thermal gain G = F K',title='Two starts can approach different oscillation sizes\nAmplitude difference is not a chaos criterion')
    ax.legend(loc='upper left',fontsize=9);fig.colorbar(im,ax=ax,label='Relative difference of displacer ranges')
    for index,n in enumerate([5,9]):
        a=np.load(ROOT/f'scout-{n}d.npz');lf=a['Lf_grid'];g=a['gain_grid'];value=a['lle'].reshape(len(lf),len(g),2).max(axis=2)
        ax=fig.add_subplot(gs[1,index]);im=ax.pcolormesh(lf/1e6,g,value.T,cmap='RdBu_r',vmin=-.01,vmax=.01,shading='nearest')
        ax.set(xscale='log',yscale='log',xlabel='Feedback inertance Lf (10⁶ kg m⁻⁴)',ylabel='Thermal gain G = F K',title=f'{n}D: larger finite-time estimate from two starts\n500 s discarded + 1,000 s observed')
        ax.plot(.474,452.1,'ko',ms=4);fig.colorbar(im,ax=ax,label='Finite-time largest exponent (s⁻¹)',extend='both')
    fig.suptitle('NIFTE: source-cycle controls and a search for persistent chaos',fontsize=15)
    path=Path(__file__).resolve().parents[1]/'docs/assets/nifte-controls.png';path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=170);plt.close(fig);print(path)
    scout=np.load(ROOT/'random-5d.npz');spec=np.load(ROOT/'cpu-row215-spectrum.npz');cycle=np.load(ROOT/'cpu-row215-cycle.npz')
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    x=np.arange(1,5);ax=axes[0]
    ax.plot(x,scout['lle_windows'][215],'o-',label='Scout: four 250 s windows',color='#d95f02')
    ax.plot(x,spec['spectrum_windows_per_s'][:,0],'s-',label='CPU: four later 1,000 s windows',color='#2166ac')
    ax.plot(x,spec['flow_neutral_window_rates_per_s'],'kx',ms=8,label='CPU neutral-flow control')
    ax.axhline(0,color='.5',lw=.7);ax.set(xticks=x,xlabel='Window number within each separate run',ylabel='Finite-time exponent (s⁻¹)',title='All four scout windows were positive')
    ax.legend(fontsize=8,loc='upper right')
    s=cycle['states'];p=cycle['parameters'];ax=axes[1]
    ax.plot(s[:,1]*p[-2]/(621*9.81)*1000,s[:,0]*p[-2]/1000,color='#2166ac')
    ax.set(xlabel='Displacer displacement (mm)',ylabel='Pressure deviation (kPa)',title='Independent closed stable cycle: T = 1.99955 s\nLargest transverse multiplier = 0.909785')
    fig.suptitle('A weak positive estimate can come from a regular periodic orbit',fontsize=13)
    path=path.with_name('nifte-neutral-control.png');fig.savefig(path,dpi=170);plt.close(fig);print(path)


if __name__=='__main__':main()
