"""Exact equal-cluster cycles and transverse Floquet blocks of motor rows.

No model modification and no GPU.  Five-state cluster closure is checked
against the independent source force solve; complete N=8 monodromy supplies
an independent decomposition check.  This is a local stability calculation,
not a proof of a global route to chaos.  Default CPU wall budget: 600 seconds.
"""
import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag
from scipy.optimize import brentq

from motor_rows import field_jac as full_field_jac
from motor_rows import rhs as full_rhs
from motor_rows import source_force_constraint_field, source_parameters


def field_jac(state, parameters):
    h, y, z, u, v = state
    rho, r, b = parameters
    a = -rho+r*(1-z)
    field = np.array([a*h-r*y, h*(1-z)-y, h*y-z,
                      -u+h*v, -v-h*u+b])
    jac = np.array([[a, -r, -r*h, 0, 0],
                    [1-z, -1, -h, 0, 0],
                    [y, h, -1, 0, 0],
                    [v, 0, 0, -1, h],
                    [-u, 0, 0, -h, -1]])
    return field, jac


def embedding(state, parameters, n=8):
    if n % 2:
        raise ValueError('Equal clusters require even N')
    h, y, z, u, v = state
    rho, r, _ = parameters
    sign = np.r_[np.ones(n//2), -np.ones(n//2)]
    c, s, x = 1-z+sign*u, v-sign*y, sign*(h-r*y)/rho
    return np.r_[c, s, x[:-1]]


def embedding_matrix(parameters, n=8):
    zero = np.zeros(5)
    return np.column_stack([embedding(zero+q, parameters, n)-embedding(zero, parameters, n)
                            for q in np.eye(5)])


def transverse_matrix(state, parameters, sign):
    h, y, z, u, v = state
    rho, r, _ = parameters
    c, s, velocity = 1-z+sign*u, v-sign*y, -sign*h
    return np.array([[-1, -velocity-r*s, rho*s],
                     [velocity, -1+r*c, -rho*c], [0, r, -rho]])


def transverse_embedding(n, first, second):
    # Embed (delta c,delta s,delta X) on two rows with opposite signs.
    out = np.zeros((3*n-1, 3))
    for row, sign in ((first, 1), (second, -1)):
        out[row, 0] = sign
        out[n+row, 1] = sign
        if row < n-1:
            out[2*n+row, 2] = sign
    return out


def audit(samples=64):
    rng = np.random.default_rng(2026103008)
    errors = dict(base_jacobian=0., source_embedding=0., base_tangent_embedding=0.,
                  transverse_A_embedding=0., transverse_B_embedding=0.)
    for _ in range(samples):
        p = (10**rng.uniform(-1, 1), 10**rng.uniform(-1, 1), rng.uniform(-2, 2))
        y = rng.normal(size=5)
        f, jac = field_jac(y, p)
        e = embedding_matrix(p)
        full = embedding(y, p)
        ffull = source_force_constraint_field(full, p, 8)
        jfull = full_field_jac(full, p, 8)[1]
        cj = np.column_stack([np.imag(field_jac(y.astype(complex)+1e-28j*q, p)[0])/1e-28
                              for q in np.eye(5)])
        comparisons = dict(base_jacobian=(jac, cj), source_embedding=(ffull, e@f),
                           base_tangent_embedding=(jfull@e, e@jac))
        for sign, first, second, name in ((1, 0, 3, 'A'), (-1, 4, 7, 'B')):
            t = transverse_embedding(8, first, second)
            comparisons['transverse_'+name+'_embedding'] = (jfull@t, t@transverse_matrix(y,p,sign))
        for name,(left,right) in comparisons.items():
            errors[name] = max(errors[name], float(np.max(np.abs(left-right)/(1+np.abs(right)))))
    assert max(errors.values()) < 2e-11, errors
    return dict(seed=2026103008, samples=samples, complex_step=1e-28,
                normalized_error='max(abs(left-right)/(1+abs(right)))', **errors)


class Integrator:
    def __init__(self, wall_seconds):
        self.started = time.monotonic()
        self.deadline = self.started+wall_seconds
        self.calls = []

    def integrate(self, state, duration, p, method='DOP853', augmented=False,
                  full=False, dense=False, events=None, max_step=.05):
        start = time.monotonic()
        if full:
            d = 23
            fun0 = lambda y: full_field_jac(y,p,8)
        else:
            d = 5
            fun0 = lambda y: field_jac(y,p)
        def fun(t, q):
            if time.monotonic() > self.deadline:
                raise TimeoutError('Declared CPU wall budget reached')
            f,j = fun0(q[:d])
            if not augmented:
                return f
            values = [f, (j@q[d:d+d*d].reshape(d,d)).ravel()]
            if not full:
                off = d+d*d
                for sign in (1,-1):
                    m = q[off:off+9].reshape(3,3)
                    values.append((transverse_matrix(q[:d],p,sign)@m).ravel())
                    off += 9
            return np.concatenate(values)
        if augmented:
            init = np.r_[state, np.eye(d).ravel()]
            if not full:
                init = np.r_[init, np.eye(3).ravel(), np.eye(3).ravel()]
        else:
            init = state
        kw = dict(jac=lambda t,y:fun0(y)[1]) if method=='Radau' and not augmented else {}
        sol = solve_ivp(fun,(0.,float(duration)),init,method=method,rtol=2e-11,
                        atol=2e-13,max_step=max_step,dense_output=dense,events=events,**kw)
        self.calls.append(dict(method=method,duration=float(duration),dimension=len(init),
                               success=bool(sol.success),message=sol.message,nfev=sol.nfev,
                               njev=sol.njev,wall_seconds=time.monotonic()-start))
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol


def serial_eigen(matrix, period):
    eigen = np.linalg.eigvals(matrix)
    order = np.argsort(-np.abs(eigen))
    eigen = eigen[order]
    return dict(multipliers=[dict(real=float(z.real),imag=float(z.imag),modulus=float(abs(z))) for z in eigen],
                exponents=[float(np.log(abs(z))/period) if abs(z)>0 else None for z in eigen])


def shoot(integrator, guess, period, p, max_iterations=6):
    state = np.asarray(guess,dtype=float).copy()
    state[0] = 0.
    history = []
    for it in range(max_iterations):
        sol = integrator.integrate(state,period,p,augmented=True)
        final = sol.y[:5,-1]
        residual = final-state
        norm = float(np.max(np.abs(residual)))
        history.append(dict(iteration=it,period=float(period),residual_max_abs=norm,state=state.tolist()))
        if norm < 5e-11:
            break
        matrix = sol.y[5:30,-1].reshape(5,5)
        shooting_jac = np.column_stack(((matrix-np.eye(5))[:,1:],field_jac(final,p)[0]))
        correction = np.linalg.solve(shooting_jac,-residual)
        if np.max(np.abs(correction[:4]))>1 or abs(correction[-1])>.25*period:
            raise RuntimeError('Shooting correction exceeded declared bounded step')
        state[1:] += correction[:4]
        period += correction[-1]
    if history[-1]['residual_max_abs']>5e-9:
        raise RuntimeError('Periodic closure failed declared 5e-9 gate: '+repr(history))
    return state,float(period),history


def run(output, wall_seconds=600):
    output.mkdir(parents=True,exist_ok=True)
    clock = Integrator(wall_seconds)
    record = dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  model='Exact equal-cluster restriction of published motor-row model',
                  original_time='tau',N=8,cluster_sizes=[4,4],
                  transverse_block_multiplicity_each=3,
                  tolerances=dict(rtol=2e-11,atol=2e-13,max_step=.05),
                  wall_budget_seconds=wall_seconds,algebra=audit(),cases=[],failed_attempts=[])
    path = output/'analysis.json'
    def save():
        record['solver_calls']=clock.calls
        record['wall_seconds']=time.monotonic()-clock.started
        record['script_sha256']=hashlib.sha256(Path(__file__).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()
        path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    source=source_parameters(.1)
    r,b=source['r'],1.
    record['source_positive_control']=source
    # Ordinary physical near-equilibrium start, declared here rather than
    # inferred from the paper's unspecified random draws.
    p=(source['rho'],r,b)
    initial=np.array([.02,0.,0.,0.,1.])
    def upward(t,y): return y[0]
    upward.direction=1
    try:
        convergence=clock.integrate(initial,300.,p,events=upward)
        times=convergence.t_events[0]
        states=convergence.y_events[0]
        if len(times)<5:
            raise RuntimeError('Insufficient source-control section returns')
        period=float(times[-1]-times[-2])
        state=states[-1]
        record['initial_control']=dict(initial=initial.tolist(),return_count=len(times),
            last_return_intervals=np.diff(times[-6:]).tolist(),last_section_state=state.tolist())
        grid=np.r_[np.linspace(source['rho'],2.934090909090909,18)]
        for index,rho in enumerate(grid):
            p=(float(rho),r,b)
            try:
                state,period,history=shoot(clock,state,period,p)
                sol=clock.integrate(state,period,p,augmented=True,dense=True)
                matrix=sol.y[5:30,-1].reshape(5,5)
                ma,mb=sol.y[30:39,-1].reshape(3,3),sol.y[39:48,-1].reshape(3,3)
                sample_t=np.linspace(0,period,501)
                sample=sol.sol(sample_t)[:5]
                physical=np.column_stack([embedding(q,p) for q in sample.T])
                radius=physical[:8]**2+physical[8:16]**2
                f0=field_jac(state,p)[0]
                case=dict(index=index,parameters=list(p),period=period,initial_state=state.tolist(),
                    shooting_history=history,closure_max_abs=float(np.max(np.abs(sol.y[:5,-1]-state))),
                    phase_condition='h(0)=0 with h_prime>0',h_prime_at_section=float(f0[0]),
                    flow_multiplier_residual=float(np.max(np.abs(matrix@f0-f0))),
                    base=serial_eigen(matrix,period),transverse_A=serial_eigen(ma,period),
                    transverse_B=serial_eigen(mb,period),max_row_disk_squared=float(radius.max()))
                np.savez_compressed(output/f'cycle-{index:02d}.npz',time=sample_t,base=sample,
                    full=physical,base_monodromy=matrix,transverse_A=ma,transverse_B=mb,
                    parameters=np.array(p))
                if index in (0,len(grid)-1):
                    independent=clock.integrate(state,period,p,method='Radau',dense=True)
                    case['Radau_base_check']=dict(closure_max_abs=float(np.max(np.abs(independent.y[:,-1]-state))),
                        dense_curve_max_abs=float(np.max(np.abs(independent.sol(sample_t)-sample))))
                    # Independent complete 23D variational dynamics, not a
                    # projected leading tangent, checks all repeated blocks.
                    full=clock.integrate(embedding(state,p),period,p,full=True,augmented=True)
                    fullmon=full.y[23:,-1].reshape(23,23)
                    basis=np.column_stack([embedding_matrix(p)]+[transverse_embedding(8,i,3) for i in range(3)]
                                           +[transverse_embedding(8,i,7) for i in range(4,7)])
                    predicted=basis@block_diag(matrix,ma,ma,ma,mb,mb,mb)@np.linalg.inv(basis)
                    case['full_23D_check']=dict(closure_max_abs=float(np.max(np.abs(full.y[:23,-1]-embedding(state,p)))),
                        monodromy_decomposition_max_abs=float(np.max(np.abs(fullmon-predicted))),
                        basis_condition=float(np.linalg.cond(basis)),spectrum=serial_eigen(fullmon,period))
                    np.savez_compressed(output/f'full-monodromy-{index:02d}.npz',monodromy=fullmon,
                                        decomposed=predicted,basis=basis)
                record['cases'].append(case)
                save()
                print(json.dumps(dict(index=index,rho=float(rho),period=period,
                    base_exponents=case['base']['exponents'],
                    transverse_exponents=case['transverse_A']['exponents'],
                    closure=case['closure_max_abs'])),flush=True)
            except Exception as e:
                record['failed_attempts'].append(dict(index=index,parameters=list(p),error=repr(e)))
                save()
                raise
        record['status']='completed_declared_continuation'
    except Exception as e:
        record['status']='partial_or_failed'
        record['terminal_error']=repr(e)
    finally:
        record['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        save()
    return record


def refine_and_plot(output, wall_seconds=600):
    """Independent variational checks and bounded local crossing refinement.

    Reads, but never replaces, the initial continuation artifact.  No branch
    beyond that interval is sought and no global transition is inferred.
    """
    source_path=output/'analysis.json'
    original=json.loads(source_path.read_text(encoding='utf-8'))
    clock=Integrator(wall_seconds)
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                source_analysis_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                source_driver='source-versions/initial-driver.py',
                independent_variational_checks=[],root_evaluations=[],failed_attempts=[])
    def save():
        record['wall_seconds']=time.monotonic()-clock.started
        record['solver_calls']=clock.calls
        record['script_sha256']=hashlib.sha256(Path(__file__).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()
        (output/'refinement.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    try:
        for case in (original['cases'][0],original['cases'][-1]):
            p=case['parameters']; period=case['period']; state=np.array(case['initial_state'])
            original_arrays=np.load(output/f"cycle-{case['index']:02d}.npz")
            sol=clock.integrate(state,period,p,method='Radau',augmented=True,dense=True)
            sample_t=np.linspace(0,period/2,301)
            parity=np.array([-1,-1,1,-1,1])[:,None]
            symmetry_error=float(np.max(np.abs(sol.sol(sample_t+period/2)[:5]-parity*sol.sol(sample_t)[:5])))
            ma=sol.y[30:39,-1].reshape(3,3)
            mb=sol.y[39:48,-1].reshape(3,3)
            row=dict(index=case['index'],parameters=p,period=period,
                     Radau_base_closure_max_abs=float(np.max(np.abs(sol.y[:5,-1]-state))),
                     base_monodromy_difference_max_abs=float(np.max(np.abs(sol.y[5:30,-1].reshape(5,5)-original_arrays['base_monodromy']))),
                     transverse_A_monodromy_difference_max_abs=float(np.max(np.abs(ma-original_arrays['transverse_A']))),
                     transverse_B_monodromy_difference_max_abs=float(np.max(np.abs(mb-original_arrays['transverse_B']))),
                     Radau_transverse_A=serial_eigen(ma,period),Radau_transverse_B=serial_eigen(mb,period),
                     half_period_cluster_exchange_max_abs=symmetry_error)
            record['independent_variational_checks'].append(row)
            save()
        left,right=original['cases'][6],original['cases'][5]
        r,b=left['parameters'][1:]
        def rootfun(rho):
            near=min(original['cases'],key=lambda q:abs(q['parameters'][0]-rho))
            p=(rho,r,b)
            state,period,history=shoot(clock,np.array(near['initial_state']),near['period'],p)
            sol=clock.integrate(state,period,p,augmented=True)
            ma=sol.y[30:39,-1].reshape(3,3)
            eigen=np.linalg.eigvals(ma)
            chosen=eigen[np.argmin(eigen.real)]
            if abs(chosen.imag)>1e-7:
                raise RuntimeError('Crossing multiplier is not real within declared tolerance')
            item=dict(rho=float(rho),period=period,multiplier=float(chosen.real),
                      root_function=float(chosen.real+1),initial_state=state.tolist(),
                      closure_max_abs=float(np.max(np.abs(sol.y[:5,-1]-state))),shooting_history=history)
            record['root_evaluations'].append(item)
            save()
            return float(chosen.real+1)
        rho_cross=brentq(rootfun,left['parameters'][0],right['parameters'][0],xtol=2e-10,rtol=2e-12,maxiter=16)
        nearest=min(record['root_evaluations'],key=lambda q:abs(q['rho']-rho_cross))
        # A bounded symmetric derivative supplies a transversality estimate,
        # not the nonlinear normal-form coefficient of an equivariant flip.
        epsilon=1e-4
        slope=(rootfun(rho_cross+epsilon)-rootfun(rho_cross-epsilon))/(2*epsilon)
        rho_source=original['source_positive_control']['rho']
        delta=1.1*(rho_source/rho_cross)**2-1
        record['crossing']=dict(rho=float(rho_cross),r=float(r),b=float(b),
            source_activity_delta=float(delta),period=nearest['period'],multiplier=nearest['multiplier'],
            d_multiplier_d_rho_centered=float(slope),difference_epsilon=epsilon,
            multiplicity_full_N8=6,base_cycle_stable_in_restriction=True,
            classification='Symmetry-degenerate transverse multiplier -1 crossing; local flip candidate only',
            not_established=['Nonlinear nondegeneracy','Existence/stability of doubled branches',
                'Global route to observed hyperchaos','Identity between periodic Floquet and chaotic Lyapunov directions'])
        record['status']='completed_bounded_crossing_and_independent_checks'
        save()
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        cases=original['cases']
        rho=np.array([q['parameters'][0] for q in cases])
        activity=1.1*(rho_source/rho)**2-1
        trans=np.array([q['transverse_A']['exponents'][0] for q in cases])
        base=np.array([q['base']['exponents'][1] for q in cases])
        fig,axes=plt.subplots(1,2,figsize=(11,4.4),constrained_layout=True)
        axes[0].plot(activity,trans,'o-',label='Leading transverse exponent (six copies)')
        axes[0].plot(activity,base,'s-',label='Leading non-neutral cluster exponent')
        axes[0].axhline(0,color='black',lw=.8)
        axes[0].axvline(delta,color='#ad3636',ls='--',label='Transverse multiplier crosses -1')
        axes[0].set(xlabel='Source activity increment delta',ylabel='Floquet exponent per normalized time',
                    title='Periodic 4+4 branch: local stability')
        axes[0].legend(fontsize=8)
        for index,label,color in ((0,'Source periodic control','#276a9f'),(len(cases)-1,'At hyperchaotic-candidate parameters','#ad3636')):
            trace=np.load(output/f'cycle-{index:02d}.npz')
            full=trace['full']
            axes[1].plot(full[8],full[12],color=color,label=label)
        axes[1].set(xlabel='Cluster A sine coefficient',ylabel='Cluster B sine coefficient',
                    title='Same periodic branch survives but becomes unstable')
        axes[1].legend(fontsize=8)
        for ax in axes: ax.grid(alpha=.2)
        fig.suptitle('Exact N=8 motor rows: transverse cluster instability, not a proven global chaos route')
        fig.savefig(output/'transverse-stability.png',dpi=180)
        plt.close(fig)
    except Exception as e:
        record['status']='partial_or_failed'
        record['failed_attempts'].append(dict(error=repr(e)))
    finally:
        record['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        save()
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('results/motor-rows/cluster-cycles'))
    parser.add_argument('--wall-seconds',type=float,default=600.)
    parser.add_argument('--refine-and-plot',action='store_true')
    args=parser.parse_args()
    report=(refine_and_plot if args.refine_and_plot else run)(args.output,args.wall_seconds)
    print(json.dumps(dict(status=report['status'],cases=len(report.get('cases',[])),wall_seconds=report['wall_seconds'])))
