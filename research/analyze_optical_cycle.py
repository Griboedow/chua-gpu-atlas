"""Locate local cycle-stability crossings and check selected saved cycles."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq,linear_sum_assignment
from scipy.integrate import solve_ivp
from continue_optical_cycle import close,shoot,ROOT,OUT
from smooth_operators import rhs_and_jac


def eig(record):
    a=np.array(record['multipliers']);return a[:,0]+1j*a[:,1]


def run():
    narrow=json.loads((OUT/'from-narrow.json').read_text())['records']
    broad=json.loads((OUT/'from-broad.json').read_text())['records']
    scale=np.array([1,np.sqrt(narrow[0]['parameters'][0]),1,1,np.sqrt(narrow[0]['parameters'][0]),1])
    state_error=period_error=eig_error=0.
    for a,b in zip(narrow,broad[::-1]):
        assert abs(a['parameters'][4]-b['parameters'][4])<1e-12
        state_error=max(state_error,float(max(abs((np.array(a['section_state'])-b['section_state'])/scale))))
        period_error=max(period_error,abs(a['period']-b['period']))
        ea,eb=eig(a),eig(b);ii,jj=linear_sum_assignment(abs(ea[:,None]-eb[None,:]))
        eig_error=max(eig_error,float(max(abs(ea[ii]-eb[jj]))))
    crossings=[];controls=[]
    for lo,hi in [(.0065,.007),(.0135,.014)]:
        anchor=min(narrow,key=lambda r:abs(r['parameters'][4]-lo));cache={}
        def at(w):
            if w not in cache:
                p=np.array(anchor['parameters']);p[4]=w
                cache[w]=close(p,np.array(anchor['section_state']),anchor['period'],rtol=3e-12,max_step=.005)
            return cache[w]
        def target(w):
            e=eig(at(w));v=e[np.argmin(abs(e+1))]
            if abs(v.imag)>1e-6:raise RuntimeError('Crossing is not real')
            return float(v.real+1)
        w=brentq(target,lo,hi,xtol=1e-11);r=at(w)
        h=1e-6;slope=(target(w+h)-target(w-h))/(2*h)
        e=eig(r);critical=np.argmin(abs(e+1))
        r.update(width_bracket=[lo,hi],critical_multiplier=e[critical].real,
                 other_multiplier_max=float(max(abs(np.delete(e,critical)))),
                 multiplier_derivative_width=slope,
                 interpretation='Simple real multiplier crossing -1; flip candidate, nonlinear nondegeneracy not established.')
        crossings.append(r)
        print(json.dumps(dict(width=w,period=r['period'],slope=slope)),flush=True)
    selected=[min(narrow,key=lambda r:abs(r['parameters'][4]-w)) for w in [.003,.005,.01,.04]]+crossings
    for r in selected:
        p=np.array(r['parameters']);y=np.array(r['section_state']);T=r['period']
        end,M=shoot(p,y,T,method='Radau',rtol=2e-12,max_step=.0025)
        f=rhs_and_jac(end,p)[0]
        section=((np.eye(6)-np.outer(f,np.eye(6)[0])/f[0])@M)[1:,1:]
        ev=np.linalg.eigvals(section);expected=eig(r)
        ii,jj=linear_sum_assignment(abs(ev[:,None]-expected[None,:]))
        def event(t,z):return z[0]-p[2]
        event.direction=1
        sol=solve_ivp(lambda t,z:rhs_and_jac(z,p)[0],(0,T),y,rtol=2e-12,atol=2e-14,
                      method='DOP853',max_step=.0025,events=event,dense_output=True)
        interior=sol.t_events[0][(sol.t_events[0]>1e-6)&(sol.t_events[0]<T-1e-6)]
        row=dict(width=float(p[4]),method='Radau',rtol=2e-12,max_step=.0025,
                 closure_balanced_max=float(max(abs((end-y)/scale))),
                 multiplier_max_difference=float(max(abs(ev[ii]-expected[jj]))),
                 liouville_logdet_error=float(abs(np.linalg.slogdet(M)[1]+2*(p[1]+1)*T)),
                 upward_section_crossings_per_period=1+len(interior),
                 half_period_nonclosure=float(max(abs((sol.sol(T/2)-y)/scale))))
        controls.append(row)
        if row['closure_balanced_max']>1e-8 or row['multiplier_max_difference']>1e-6 or not sol.success:
            raise RuntimeError(row)
        print(json.dumps(row),flush=True)
    result=dict(common_widths=len(narrow),bidirectional_max_balanced_state_error=state_error,
                bidirectional_max_period_error=period_error,bidirectional_max_matched_multiplier_error=eig_error,
                crossings=crossings,independent_solver_controls=controls,
                limitations='Only this cycle branch was continued. No global basin result, no rigorous nonlinear flip proof, no full period-doubling route or crisis, and no interpretation of fixed-parameter Newton failure as a fold.')
    (OUT/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    plot(narrow,result)


def plot(rows,analysis):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    w=np.array([r['parameters'][4] for r in rows]);period=np.array([r['period'] for r in rows])
    rho=np.array([r['spectral_radius'] for r in rows]);stable=rho<1
    fig,axes=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
    for ax,data in zip(axes,[period,rho]):
        ax.plot(w,data,color='#777',lw=1)
        ax.scatter(w[stable],data[stable],s=13,color='#237e86',label='Stable cycle')
        ax.scatter(w[~stable],data[~stable],s=20,color='#c95d46',label='Unstable cycle')
        for r in analysis['crossings']:ax.axvline(r['parameters'][4],color='#aa7c28',ls='--',lw=1)
        ax.set_xlabel('Assumed optical width w');ax.grid(alpha=.2)
    axes[0].set(title='Same cycle continued from both ends',ylabel='Period / thermal relaxation time')
    axes[1].axhline(1,color='#222',lw=.8);axes[1].set(title='Two local stability crossings',ylabel='Largest Poincare multiplier magnitude')
    axes[0].legend(fontsize=8);axes[1].legend(fontsize=8)
    fig.suptitle('Optical capture cycle: stable, unstable, then stable again\nDashed lines mark real multipliers at -1; no global chaos boundary is claimed',fontsize=12)
    fig.savefig(ROOT/'docs/assets/optical_cycle_continuation.png',dpi=170);plt.close(fig)


if __name__=='__main__':run()
