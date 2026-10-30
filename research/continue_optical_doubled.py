"""Follow the nontrivial 2T branch toward the first optical flip candidate."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment
from continue_optical_cycle import close,shoot,ROOT,OUT
from smooth_operators import rhs_and_jac


def run():
    attempts=json.loads((OUT/'doubled-probes.json').read_text())['attempts']
    original=next(r for r in attempts if r['classification']=='nontrivial doubled-orbit candidate')['orbit']
    p=np.array(original['parameters']);y=np.array(original['section_state']);T=original['period']
    scale=np.array([1,np.sqrt(p[0]),1,1,np.sqrt(p[0]),1]);rows=[]
    for w in [.0064,.00646,.00650,.00654,.00657,.00658]:
        p[4]=w;r=close(p,y,T,rtol=2e-12,max_step=.0025);y=np.array(r['section_state']);T=r['period']
        half,_=shoot(p,y,T/2,rtol=2e-12,max_step=.0025)
        nonclosure=float(max(abs((half-y)/scale)))
        if nonclosure<1e-5:raise RuntimeError('Collapsed to repeated parent orbit')
        end,M=shoot(p,y,T,method='Radau',rtol=1e-12,max_step=.0025)
        f=rhs_and_jac(end,p)[0];S=((np.eye(6)-np.outer(f,np.eye(6)[0])/f[0])@M)[1:,1:]
        ev=np.linalg.eigvals(S);expected=np.array(r['multipliers']);expected=expected[:,0]+1j*expected[:,1]
        ii,jj=linear_sum_assignment(abs(ev[:,None]-expected[None,:]))
        def event(t,q):return q[0]-p[2]
        event.direction=1
        sol=solve_ivp(lambda t,q:rhs_and_jac(q,p)[0],(0,T),y,rtol=1e-12,atol=1e-14,
                      method='DOP853',max_step=.0025,events=event)
        inside=sol.t_events[0][(sol.t_events[0]>1e-6)&(sol.t_events[0]<T-1e-6)]
        if not sol.success or len(inside)!=1:raise RuntimeError('Not a two-crossing cycle')
        r.update(half_period_nonclosure=nonclosure,upward_crossings_per_period=2,
                 independent_radau_closure=float(max(abs((end-y)/scale))),
                 independent_radau_multiplier_error=float(max(abs(ev[ii]-expected[jj]))))
        if r['independent_radau_closure']>1e-8 or r['independent_radau_multiplier_error']>1e-6:raise RuntimeError(r)
        rows.append(r);print(json.dumps(dict(w=w,T=T,nonclosure=nonclosure,rho=r['spectral_radius'])),flush=True)
        (OUT/'doubled-branch.json').write_text(json.dumps(dict(records=rows,interpretation='Unstable two-crossing periodic branch on the stable-parent side, numerically approaching the first flip. Supports a subcritical local picture; not a rigorously certified normal form or a chaos-crisis mechanism.'),indent=2)+'\n',encoding='utf8')
    plot(rows)


def plot(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator,FormatStrFormatter
    critical=json.loads((OUT/'analysis.json').read_text())['crossings'][0]['parameters'][4]
    widths=np.array([r['parameters'][4] for r in rows]);a=np.array([r['half_period_nonclosure'] for r in rows]);rho=np.array([r['spectral_radius'] for r in rows])
    fig,axes=plt.subplots(1,2,figsize=(10.5,3.9),layout='constrained')
    axes[0].plot(critical-widths,a*a,'o-',color='#a85444');axes[0].set(xlabel='Distance below first flip width',ylabel='Squared half-period nonclosure',title='Doubled orbit approaches the parent')
    axes[1].plot(widths,rho,'o-',color='#a85444');axes[1].axhline(1,color='#555',lw=1)
    axes[1].axvline(critical,color='#a87b27',ls='--',lw=1);axes[1].set(xlabel='Assumed optical width w',ylabel='Largest Poincare multiplier magnitude',title='Doubled branch remains unstable')
    axes[0].xaxis.set_major_locator(MaxNLocator(4))
    axes[0].ticklabel_format(axis='x',style='sci',scilimits=(-3,3),useMathText=True)
    axes[1].xaxis.set_major_locator(MaxNLocator(4))
    axes[1].xaxis.set_major_formatter(FormatStrFormatter('%.5f'))
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('A local subcritical-flip picture: numerical periodic orbits, not an attractor',fontsize=12)
    fig.savefig(ROOT/'docs/assets/optical_doubled_branch.png',dpi=170);plt.close(fig)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--plot',action='store_true');a=parser.parse_args()
    if a.plot:plot(json.loads((OUT/'doubled-branch.json').read_text())['records'])
    else:run()
