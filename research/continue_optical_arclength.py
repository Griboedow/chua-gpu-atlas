"""Bounded pseudo-arclength continuation of the known optical capture cycle.

Only the explicitly assumed logistic width varies. This is orbit shooting, not
a search for chaos or a basin computation. Six balanced closure equations plus
one arclength hyperplane determine five section coordinates, period and width.
The width sensitivity is integrated analytically alongside the state and STM.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from smooth_operators import rhs_and_jac

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/light-operators/cycle-arclength'
SOURCE = ROOT / 'results/light-operators/cycle-continuation/from-narrow.json'
P = np.array([74.1567393414624, 1.784180421282182, .4139517339970917, .25683679452165964, .04])
SCALE = np.array([1., np.sqrt(P[0]), 1., 1., np.sqrt(P[0]), 1.])
T_SCALE = 2.
W_SCALE = .01
BEGAN = 0.
WALL_BUDGET = 780.


def budget():
    if time.perf_counter()-BEGAN > WALL_BUDGET:
        raise TimeoutError('Overall CPU wall budget reached before next orbit solve')


def width_derivative(y, p, jac):
    fw = np.zeros(6)
    fw[2] = -(y[3]-p[3])/p[4]*jac[2, 3]
    fw[5] = (p[2]-y[0])/p[4]*jac[5, 0]
    return fw


def pack(state, period, width):
    return np.r_[np.asarray(state)[1:]/SCALE[1:], period/T_SCALE, width/W_SCALE]


def unpack(u):
    y = np.r_[P[2], np.asarray(u[:5])*SCALE[1:]]
    return y, float(u[5]*T_SCALE), float(u[6]*W_SCALE)


def shoot(u, method='DOP853', rtol=2e-11, max_step=.01, dense=False):
    budget()
    y, period, width = unpack(u)
    if not .3 < period < 6 or not .0001 < width < .1:
        raise ValueError('Trial leaves the explicitly bounded period/width study')
    p = P.copy(); p[4] = width
    def augmented(t, z):
        f, jac = rhs_and_jac(z[:6], p)
        return np.r_[f, (jac@z[6:42].reshape(6, 6)).ravel(),
            jac@z[42:48]+width_derivative(z[:6], p, jac)]
    sol = solve_ivp(augmented, (0, period), np.r_[y, np.eye(6).ravel(), np.zeros(6)],
        method=method, rtol=rtol, atol=rtol*.01, max_step=max_step, dense_output=dense)
    if not sol.success:
        raise RuntimeError(sol.message)
    end = sol.y[:6, -1]
    matrix = sol.y[6:42, -1].reshape(6, 6)
    sensitivity = sol.y[42:48, -1]
    f = rhs_and_jac(end, p)[0]
    residual = (end-y)/SCALE
    state_jac = (matrix-np.eye(6))/SCALE[:, None]*SCALE[None, :]
    closure_jac = np.column_stack([state_jac[:, 1:], T_SCALE*f/SCALE, W_SCALE*sensitivity/SCALE])
    return dict(state=y, end=end, period=period, width=width, parameters=p,
        matrix=matrix, sensitivity=sensitivity, residual=residual,
        jacobian=closure_jac, final_field=f, solution=sol if dense else None)


def geometry(result, previous=None):
    _, singular, vh = np.linalg.svd(result['jacobian'], full_matrices=True)
    tangent = vh[-1]
    if previous is not None and np.dot(tangent, previous) < 0:
        tangent = -tangent
    return tangent, singular


def multipliers(result):
    f = result['final_field']
    if abs(f[0]) < 1e-7:
        raise RuntimeError('Section near grazing; change phase condition rather than declaring a fold')
    derivative = ((np.eye(6)-np.outer(f, np.eye(6)[0])/f[0])@result['matrix'])[1:, 1:]
    return np.linalg.eigvals(derivative)


def record(u, result, tangent, singular, **extra):
    eig = multipliers(result)
    fixed_singular = np.linalg.svd(result['jacobian'][:, :6], compute_uv=False)
    first_field = rhs_and_jac(result['state'], result['parameters'])[0]
    return dict(u=np.asarray(u).tolist(), state=result['state'].tolist(), period=result['period'], width=result['width'],
        closure_balanced_max=float(np.max(abs(result['residual']))),
        tangent=np.asarray(tangent).tolist(), width_tangent=float(tangent[-1]*W_SCALE),
        full_row_singular_values=singular.tolist(), fixed_width_singular_values=fixed_singular.tolist(),
        multipliers=[[float(z.real), float(z.imag)] for z in eig], spectral_radius=float(max(abs(eig))),
        crossing_velocity=float(result['state'][1]),
        neutral_flow_error=float(np.linalg.norm(result['matrix']@first_field-result['final_field'])/np.linalg.norm(result['final_field'])),
        liouville_logdet_error=float(abs(np.linalg.slogdet(result['matrix'])[1]+2*(P[1]+1)*result['period'])),
        **extra)


def correct(predictor, tangent):
    u = predictor.copy()
    history = []
    for iteration in range(14):
        result = shoot(u)
        residual = np.r_[result['residual'], np.dot(tangent, u-predictor)]
        norm = float(np.max(abs(residual)))
        history.append(norm)
        if norm < 3e-11:
            if result['state'][1] <= 1e-6:
                raise RuntimeError('Section orientation lost; a general phase condition is needed')
            return u, result, history
        jac = np.vstack([result['jacobian'], tangent])
        delta = np.linalg.solve(jac, -residual)
        accepted = False
        for power in range(8):
            trial = u+delta*2.**(-power)
            try:
                value = shoot(trial)
            except ValueError:
                continue
            rr = np.r_[value['residual'], np.dot(tangent, trial-predictor)]
            if max(abs(rr)) < norm:
                u = trial
                accepted = True
                break
        if not accepted:
            raise RuntimeError(f'Arclength Newton line search stopped at residual {norm:g}')
    raise RuntimeError(f'Arclength Newton iteration budget: residual {history[-1]:g}')


def derivative_controls(u):
    rng = np.random.default_rng(2026103081)
    errors = []
    for _ in range(64):
        y = rng.uniform(-.2, 1.2, 6)
        p = P.copy(); p[4] = 10**rng.uniform(-2.5, -1.)
        exact = width_derivative(y, p, rhs_and_jac(y, p)[1])
        h = p[4]*1e-5
        plus = p.copy(); plus[4] += h
        minus = p.copy(); minus[4] -= h
        finite = (rhs_and_jac(y, plus)[0]-rhs_and_jac(y, minus)[0])/(2*h)
        errors.append(float(np.max(abs(finite-exact))/(1+np.max(abs(exact)))))
    base = shoot(u, rtol=2e-12, max_step=.005)
    h = 2e-5
    plus = u.copy(); plus[-1] += h
    minus = u.copy(); minus[-1] -= h
    derivative = (shoot(plus, rtol=2e-12, max_step=.005)['residual']-
        shoot(minus, rtol=2e-12, max_step=.005)['residual'])/(2*h)
    endpoint_error = float(np.max(abs(derivative-base['jacobian'][:, -1])))
    return dict(pointwise_width_derivative_relative_max=max(errors), endpoint_width_column_max_error=endpoint_error,
        endpoint_scaled_width_difference_step=h,
        passed=bool(max(errors) < 1e-7 and endpoint_error < 2e-6))


def fold_refine(seed):
    cache = {}
    def evaluate(u):
        if 'u' not in cache or not np.array_equal(u, cache['u']):
            result = shoot(u, rtol=3e-12, max_step=.005)
            eig = multipliers(result)
            critical = eig[np.argmin(abs(eig-1))]
            if abs(critical.imag) > 1e-5:
                raise RuntimeError('Nearest unit multiplier is not a real isolated fold candidate')
            cache.update(u=np.array(u, copy=True), result=result,
                residual=np.r_[result['residual'], critical.real-1])
        return cache['result'], cache['residual']
    def fun(u):
        return evaluate(u)[1]
    def jac(u):
        result, residual = evaluate(u)
        top = result['jacobian'].copy()
        bottom = np.empty(7)
        for j in range(7):
            h = 1e-5*max(1., abs(u[j]))
            plus = u.copy(); plus[j] += h
            minus = u.copy(); minus[j] -= h
            bottom[j] = (fun(plus)[-1]-fun(minus)[-1])/(2*h)
        return np.vstack([top, bottom])
    solved = root(fun, seed, jac=jac, method='hybr', options={'xtol': 2e-10, 'maxfev': 100})
    result, residual = evaluate(solved.x)
    tangent, singular = geometry(result)
    # A finite-difference nondegeneracy check at fixed width, in the declared
    # scaled shooting coordinates. Signs of left singular vectors are arbitrary.
    left, fixed_sv, right = np.linalg.svd(result['jacobian'][:, :6])
    ell = left[:, -1]
    direction = right[-1]
    parameter_coefficient = float(ell@result['jacobian'][:, -1])
    quadratic = []
    for h in [.002, .001, .0005]:
        offset = np.r_[h*direction, 0.]
        second = (shoot(solved.x+offset, rtol=3e-12, max_step=.005)['residual']-
            2*result['residual']+shoot(solved.x-offset, rtol=3e-12, max_step=.005)['residual'])/h**2
        coefficient = float(.5*ell@second)
        quadratic.append(dict(h=h, coefficient=coefficient,
            predicted_physical_width_coefficient=float(-W_SCALE*coefficient/parameter_coefficient)))
    row = record(solved.x, result, tangent, singular,
        solver_success=bool(solved.success), solver_message=solved.message,
        augmented_residual_max=float(max(abs(residual))), nfev=int(solved.nfev),
        fixed_width_smallest_singular=float(fixed_sv[-1]),
        parameter_transversality_coefficient=parameter_coefficient,
        quadratic_directional_checks=quadratic,
        verified_numerically=bool(max(abs(residual)) < 2e-8 and singular[-1] > 1e-6),
        caveat='Numerical fold and nondegeneracy evidence, not a computer-assisted bifurcation theorem.')
    return row


def independent_control(u, label):
    dop = shoot(u, rtol=3e-12, max_step=.005)
    radau = shoot(u, method='Radau', rtol=3e-12, max_step=.003, dense=True)
    de = multipliers(dop); re = multipliers(radau)
    from scipy.optimize import linear_sum_assignment
    a, b = linear_sum_assignment(abs(de[:, None]-re[None, :]))
    tangent, singular = geometry(radau)
    times = np.linspace(0, radau['period'], 2001)
    trajectory = radau['solution'].sol(times)[:6].T
    half = radau['solution'].sol(radau['period']/2)[:6]
    quarter = radau['solution'].sol(radau['period']/4)[:6]
    # Count transverse upward roots on a dense mesh, excluding the initial root.
    level = trajectory[:, 0]-P[2]
    crossings = np.flatnonzero((level[:-1] < 0) & (level[1:] >= 0))
    row = record(u, radau, tangent, singular, label=label, method='Radau', rtol=3e-12, max_step=.003,
        independent_multiplier_max_error=float(np.max(abs(de[a]-re[b]))),
        half_period_nonclosure=float(np.max(abs((half-radau['state'])/SCALE))),
        quarter_period_nonclosure=float(np.max(abs((quarter-radau['state'])/SCALE))),
        sampled_upward_crossings_excluding_initial=len(crossings),
        thermal_range=[float(trajectory[:, [2, 5]].min()), float(trajectory[:, [2, 5]].max())])
    np.savez_compressed(OUT / f'control-{label}.npz', time=times, trajectory=trajectory,
        parameters=radau['parameters'], section_state=radau['state'], monodromy=radau['matrix'], width_sensitivity=radau['sensitivity'])
    return row


def save(data):
    data['wall_seconds'] = time.perf_counter()-BEGAN
    (OUT / 'continuation.json').write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def run(max_steps=100, wall=780.):
    global BEGAN, WALL_BUDGET
    BEGAN = time.perf_counter(); WALL_BUDGET = wall
    if (OUT / 'continuation.json').exists():
        raise RuntimeError('Output exists; choose a fresh --output directory')
    OUT.mkdir(parents=True, exist_ok=True)
    source = json.loads(SOURCE.read_text())['records']
    first = next(r for r in source if abs(r['parameters'][-1]-.0395) < 1e-10)
    last = next(r for r in source if abs(r['parameters'][-1]-.04) < 1e-10)
    previous_u = pack(first['section_state'], first['period'], .0395)
    u = pack(last['section_state'], last['period'], .04)
    tangent_hint = (u-previous_u)/np.linalg.norm(u-previous_u)
    controls = derivative_controls(u)
    if not controls['passed']:
        raise RuntimeError(f'Width sensitivity validation failed: {controls}')
    result = shoot(u)
    tangent, singular = geometry(result, tangent_hint)
    source_hash = hashlib.sha256(SOURCE.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    data = dict(model='Assumed logistic optical profile, fixed reference-A k,b,delta_minus,delta_plus',
        source_file=SOURCE.relative_to(ROOT).as_posix(), source_canonical_lf_sha256=source_hash,
        driver_canonical_lf_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
        settings=dict(state_scale=SCALE.tolist(), period_scale=T_SCALE, width_scale=W_SCALE,
            closure_tolerance=3e-11, method='DOP853', rtol=2e-11, atol=2e-13, max_step=.01,
            initial_arc_step=.005, maximum_arc_step=.012, minimum_arc_step=1e-5,
            max_steps=max_steps, wall_budget=wall),
        derivative_controls=controls, records=[record(u, result, tangent, singular, arclength=0., arc_step=0.)],
        rejected_attempts=[], turning_brackets=[], duplicates=[], stop_reason=None)
    step = .005; distance = 0.; after_turn = 0
    save(data)
    for index in range(max_steps):
        try:
            predictor = u+step*tangent
            new_u, new_result, history = correct(predictor, tangent)
            correction = float(np.linalg.norm(new_u-predictor))
            if correction > 2*step:
                raise RuntimeError('Corrector moved more than two predictor steps; possible branch jump')
            new_tangent, new_singular = geometry(new_result, tangent)
        except Exception as exc:
            data['rejected_attempts'].append(dict(step_index=index, arc_step=step,
                starting_u=u.tolist(), reason=f'{type(exc).__name__}: {exc}'))
            step *= .5
            save(data)
            if step < 1e-5 or isinstance(exc, TimeoutError):
                data['stop_reason'] = 'minimum arc step or wall budget after retained failure'
                break
            continue
        distance += float(np.linalg.norm(new_u-u))
        row = record(new_u, new_result, new_tangent, new_singular, arclength=distance,
            arc_step=step, newton_history=history, predictor_correction=correction)
        if tangent[-1]*new_tangent[-1] < 0:
            data['turning_brackets'].append([len(data['records'])-1, len(data['records'])])
        data['records'].append(row)
        if len(data['records']) > 4:
            old_u = np.array([r['u'] for r in data['records'][:-3]])
            matches = np.flatnonzero(np.linalg.norm(old_u-new_u, axis=1) < 1e-6)
            if len(matches):
                data['duplicates'].append(dict(record=len(data['records'])-1, previous=matches.tolist()))
        u, tangent = new_u, new_tangent
        if len(history) <= 4 and correction < .3*step:
            step = min(.012, step*1.2)
        elif len(history) > 7:
            step = max(1e-5, step*.65)
        save(data)
        print(f'arc={distance:.5f}, w={row["width"]:.10f}, T={row["period"]:.9f}, rho={row["spectral_radius"]:.7f}, dw/ds={row["width_tangent"]:.6g}, wall={data["wall_seconds"]:.1f}', flush=True)
        if data['turning_brackets']:
            after_turn += 1
            highest = max(r['width'] for r in data['records'])
            if after_turn >= 15 and row['width'] < highest-.002:
                data['stop_reason'] = 'resolved turning branch followed for more than .002 width'
                break
        if row['width'] > .052 or row['width'] < .025:
            data['stop_reason'] = 'declared local width study boundary reached'
            break
    if data['stop_reason'] is None:
        data['stop_reason'] = 'arc-step count budget reached'
    save(data)
    if data['turning_brackets']:
        a, b = data['turning_brackets'][0]
        try:
            fold = fold_refine((np.array(data['records'][a]['u'])+np.array(data['records'][b]['u']))/2)
            data['fold'] = fold
        except Exception as exc:
            data['fold_failure'] = f'{type(exc).__name__}: {exc}'
        save(data)
    selected = [('initial', np.array(data['records'][0]['u'])),
        ('last', np.array(data['records'][-1]['u']))]
    if data.get('fold', {}).get('verified_numerically'):
        selected.append(('fold', np.array(data['fold']['u'])))
    independent = []
    for label, coordinate in selected:
        try:
            independent.append(independent_control(coordinate, label))
        except Exception as exc:
            independent.append(dict(label=label, failure=f'{type(exc).__name__}: {exc}'))
        (OUT / 'independent-controls.json').write_text(json.dumps(independent, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    save(data)
    print(json.dumps(dict(records=len(data['records']), turns=data['turning_brackets'], stop=data['stop_reason'],
        fold=data.get('fold', {}).get('width'), wall=data['wall_seconds']), indent=2), flush=True)


def analyze():
    """Add one exact-width partner and robust section controls, then render.

    Existing independent results are retained. The original one-period dense
    samples can miss an endpoint crossing by rounding; integrate slightly past
    the period with an event solver instead of counting endpoint signs.
    """
    global BEGAN, WALL_BUDGET
    BEGAN = time.perf_counter(); WALL_BUDGET = 300.
    data = json.loads((OUT / 'continuation.json').read_text())
    rows = data['records']
    anchor_path = OUT / 'unstable-w004.json'
    if not anchor_path.exists():
        after = rows[data['turning_brackets'][0][1]:]
        pair = next((a, b) for a, b in zip(after[:-1], after[1:]) if a['width'] >= .04 >= b['width'])
        a, b = pair
        fraction = (.04-a['width'])/(b['width']-a['width'])
        seed = np.array(a['u'])+fraction*(np.array(b['u'])-np.array(a['u']))
        y, period, _ = unpack(seed)
        # Use the independently existing fixed-width corrector on the newly
        # located unstable side. It writes no parent continuation files.
        from continue_optical_cycle import close
        p = P.copy(); p[4] = .04
        closed = close(p, y, period, rtol=3e-12, max_step=.005)
        u = pack(closed['section_state'], closed['period'], .04)
        control = independent_control(u, 'unstable-w004')
        control['fixed_width_newton_history'] = closed['newton_history']
        control['source_arclength_width_bracket'] = [a['width'], b['width']]
        control['source_interpolated_u'] = seed.tolist()
        control['u'] = u.tolist()
        anchor_path.write_text(json.dumps(control, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    anchor = json.loads(anchor_path.read_text())
    section_path = OUT / 'section-audit.json'
    if not section_path.exists():
        targets = [('initial', rows[0]), ('last', rows[-1]), ('fold', data['fold']), ('unstable-w004', anchor)]
        checks = []
        for label, row in targets:
            y = np.array(row['state']); period = row['period']
            p = P.copy(); p[4] = row['width']
            def event(t, state):
                return state[0]-p[2]
            event.direction = 1
            sol = solve_ivp(lambda t, state: rhs_and_jac(state, p)[0], (0, 1.1*period), y,
                method='Radau', jac=lambda t, state: rhs_and_jac(state, p)[1],
                rtol=3e-12, atol=3e-14, max_step=.003, events=event)
            events = sol.t_events[0]
            noninitial = events[events > 1e-7]
            checks.append(dict(label=label, solver_success=bool(sol.success), method='Radau',
                duration=1.1*period, rtol=3e-12, atol=3e-14, max_step=.003,
                upward_crossings_excluding_initial=int(len(noninitial)), event_times=events.tolist(),
                expected_period=period, period_error=None if not len(noninitial) else float(abs(noninitial[0]-period)),
                crossing_velocities=sol.y_events[0][:, 1].tolist(),
                passed=bool(sol.success and len(noninitial)==1 and abs(noninitial[0]-period)<1e-8)))
            np.savez_compressed(OUT / f'control-{label}-section.npz', time=sol.t, trajectory=sol.y.T,
                event_times=events, event_states=sol.y_events[0], parameters=p)
        result = dict(explanation='The original dense one-period samples sometimes report zero upward crossings because a tiny negative closure error places the endpoint just before the root. This separate event audit integrates to 1.1 periods and does not overwrite those original sampled counts.',
            all_passed=all(r['passed'] for r in checks), rows=checks)
        section_path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        if not result['all_passed']:
            raise RuntimeError('A section control needs investigation')
    render()
    print(json.dumps(dict(unstable_at_w004=anchor['period'], rho=anchor['spectral_radius'],
        closure=anchor['closure_balanced_max'], analysis_wall=time.perf_counter()-BEGAN), indent=2), flush=True)


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data = json.loads((OUT / 'continuation.json').read_text())
    rows = data['records']; fold = data['fold']
    anchor = json.loads((OUT / 'unstable-w004.json').read_text())
    split = data['turning_brackets'][0][1]
    pieces = [(rows[:split]+[fold], '#2878b5', 'Stable side'),
        ([fold]+rows[split:], '#c44e35', 'Unstable side')]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), layout='constrained')
    for part, color, label in pieces:
        widths = np.array([r['width'] for r in part])
        periods = [r['period'] for r in part]
        positive = [max(complex(*z).real for z in r['multipliers']) for r in part]
        axes[0].plot(widths, periods, '.-', color=color, label=label, lw=1.2, ms=3)
        axes[1].plot(widths, positive, '.-', color=color, lw=1.2, ms=3)
    axes[0].scatter([fold['width']], [fold['period']], marker='*', s=100, color='#7c4da6', zorder=5, label='Refined fold')
    axes[1].scatter([fold['width']], [1], marker='*', s=100, color='#7c4da6', zorder=5)
    axes[0].scatter([.04, .04], [rows[0]['period'], anchor['period']], marker='s', s=30, facecolors='none', edgecolors='black', zorder=6)
    for label, color, title in [('initial', '#2878b5', 'Stable cycle, w=.04'),
        ('unstable-w004', '#c44e35', 'Saddle cycle, w=.04'), ('fold', '#7c4da6', 'Fold cycle')]:
        orbit = np.load(OUT / f'control-{label}.npz')['trajectory']
        axes[2].plot(orbit[:, 0], orbit[:, 3], color=color, lw=1.1, label=title)
    axes[0].set(xlabel='Assumed logistic width w', ylabel='Cycle period', title='The branch turns before w=.0405')
    axes[1].set(xlabel='Assumed logistic width w', ylabel='Positive real Poincare multiplier', title='Local attraction changes at +1')
    axes[1].axhline(1, color='#666', linestyle='--', lw=1)
    axes[2].set(xlabel='Displacement x_minus', ylabel='Displacement x_plus', title='Physical projection of six-state cycles')
    axes[0].legend(fontsize=8)
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('A numerical fold of the optical capture cycle\nLocal orbit structure does not establish a chaotic crisis or basin boundary', fontsize=13)
    target = (ROOT / 'docs/assets/optical_cycle_arclength.png' if
        OUT.resolve() == (ROOT / 'results/light-operators/cycle-arclength').resolve()
        else OUT / 'optical-cycle-arclength.png')
    fig.savefig(target, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--wall', type=float, default=780.)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--analyze', action='store_true')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.output:
        OUT = args.output.resolve()
    if args.analyze:
        analyze()
    elif args.plot:
        render()
    else:
        run(args.steps, args.wall)
