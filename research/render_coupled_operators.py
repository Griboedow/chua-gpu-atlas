"""Physical projections and sections, not faithful 3D embeddings of a 6D flow."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from coupled_operators import initial_modes, segment


def dense_tail(times, states, p, duration=30):
    """Densify between CPU-verified event/QR knots using the exact branch flow.

    Coarse QR endpoints are sufficient for spectra but visibly alias fast
    mechanics. Every interval ends at or before the next switching event.
    """
    dt=min(.01,2*np.pi/np.sqrt(p[0]-p[1]**2/4)/64)
    first=max(0,int(np.searchsorted(times,times[-1]-duration))-1)
    dense_t=[];dense_y=[];error=0.
    for j in range(first,len(times)-1):
        h=times[j+1]-times[j];modes=initial_modes(states[j],p)
        mesh=np.linspace(0,h,max(2,int(np.ceil(h/dt))+1))
        ys=segment(states[j],p,modes,mesh).T
        error=max(error,float(np.max(abs(ys[-1]-states[j+1]))))
        dense_t.extend(times[j]+mesh[:-1]);dense_y.extend(ys[:-1])
    if error>1e-7:raise RuntimeError(f'CPU knot/analytic branch mismatch: {error}')
    print(f'Dense plot branch-end error: {error:.3g}',flush=True)
    return np.asarray(dense_t),np.asarray(dense_y)


def draw(prefix, output, gif=False):
    data=np.load(str(prefix)+'.npz');meta=json.loads(Path(str(prefix)+'.json').read_text())
    section=data['sections'];k,b,dm,dp=meta['parameters'];rootk=np.sqrt(k)
    t,s=dense_tail(data['times'],data['states'],np.array(meta['parameters']))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(13,10),layout='constrained')
    ax=fig.add_subplot(221,projection='3d')
    take=t>t[-1]-30
    ax.plot(s[take,0],s[take,3],s[take,2]-s[take,5],lw=.45,alpha=.8,color='#008d91')
    ax.set(xlabel=r'$x_-$',ylabel=r'$x_+$',zlabel=r'$y_- - y_+$',title='Physical projection of a six-state flow')
    ax.view_init(24,42)
    ax=fig.add_subplot(222)
    if len(section):
        im=ax.scatter(section[:,3],section[:,4]/rootk,c=section[:,2],s=5,alpha=.7,cmap='viridis')
        fig.colorbar(im,ax=ax,label=r'$y_-$')
    ax.set(xlabel=r'$x_+$',ylabel=r'$v_+/\sqrt{k}$',title=r'Projected section: $x_-=\delta_-$, $v_->0$')
    ax=fig.add_subplot(223);take=t>t[-1]-10
    ax.plot(t[take]-t[-1]+10,s[take,0],lw=.9,label=r'$x_-$');ax.plot(t[take]-t[-1]+10,s[take,3],lw=.9,label=r'$x_+$')
    ax.axhline(dm,c='C0',ls='--',alpha=.4);ax.axhline(dp,c='C1',ls='--',alpha=.4)
    ax.set(xlabel='Time / thermal relaxation time',ylabel='Normalized displacement',title='Two material displacements; constant input powers');ax.legend()
    ax=fig.add_subplot(224);spec=np.array(meta['spectrum']);win=np.array(meta['spectrum_windows'])
    ax.axhline(0,c='.6',lw=.8);ax.plot(np.arange(1,7),spec,'o-',color='#574696',label='Full observation')
    for row in win:ax.plot(np.arange(1,7),row,'.',alpha=.35,color='#574696')
    ax.set(xlabel='Lyapunov exponent index',ylabel='Exponent / thermal time',title='Independent CPU spectrum; dots show four windows');ax.legend()
    fig.suptitle(f'Coupled optical operators | k={k:.4g}, b={b:.4g}, delta=({dm:.4f}, {dp:.4f})\nNumerical chaotic candidate; novelty and experimental realization unconfirmed',fontsize=13)
    Path(output).parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(output)+'.png',dpi=155);plt.close(fig)
    if gif:
        fig=plt.figure(figsize=(7,6.2));ax=fig.add_subplot(111,projection='3d')
        take=t>t[-1]-25;xyz=np.column_stack([s[take,0],s[take,3],s[take,2]-s[take,5]])
        ax.plot(*xyz.T,lw=.5,color='#008d91',alpha=.85)
        ax.set(xlabel=r'$x_-$',ylabel=r'$x_+$',zlabel=r'$y_- - y_+$')
        ax.set_title('Coupled optical operators\n3D projection of a six-dimensional flow',pad=18)
        fig.text(.5,.04,'Camera rotation only | scientific novelty unconfirmed',ha='center',fontsize=9)
        fig.subplots_adjust(left=0,right=.94,bottom=.10,top=.86)
        from matplotlib.ticker import MaxNLocator
        for axis in (ax.xaxis,ax.yaxis,ax.zaxis):axis.set_major_locator(MaxNLocator(4))
        def update(i):ax.view_init(25+18*np.sin(4*np.pi*i/72),360*i/72)
        FuncAnimation(fig,update,frames=72,interval=65).save(str(output)+'.gif',writer=PillowWriter(fps=15),dpi=95)
        plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prefix',default='results/light-operators/coupled-cpu-1026');p.add_argument('--output',default='docs/assets/coupled_operators_1026');p.add_argument('--gif',action='store_true');a=p.parse_args();draw(a.prefix,a.output,a.gif)
