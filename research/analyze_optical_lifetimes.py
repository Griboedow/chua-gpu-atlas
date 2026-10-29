"""Independently check periodic captures and plot the finite-time cohort.

No extrapolated escape law or asymptotic attractor classification is fitted.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from smooth_operators import rhs_and_jac

ROOT = Path('results/light-operators/lifetimes')


def audit_captures():
    data = np.load(ROOT / 'ensemble.npz')
    controls = json.loads((ROOT / 'periodic-controls.json').read_text())
    rows = []
    for index in np.flatnonzero(data['status'] == 1):
        p = data['parameters'][index]
        control = next(c for c in controls if c['parameters'][4] == p[4])
        scale = np.array([1, np.sqrt(p[0]), 1, 1, np.sqrt(p[0]), 1])
        def event(t, y):
            return y[0] - p[2]
        event.direction = 1
        period = control['period']
        sol = solve_ivp(lambda t, y: rhs_and_jac(y, p)[0],
            (0, 4 * period), data['final_states'][index], method='DOP853',
            rtol=1e-12, atol=1e-14, max_step=.005, events=event, dense_output=True)
        if not sol.success or len(sol.t_events[0]) < 3:
            raise RuntimeError(f'CPU capture integration failed at lane {index}')
        reference = np.array(control['section_state'])
        section_error = float(np.max(abs((sol.y_events[0][-2:] - reference) / scale)))
        period_error = float(np.max(abs(np.diff(sol.t_events[0][-3:]) - period)))
        closure = float(np.max(abs((sol.sol(4 * period) - sol.sol(3 * period)) / scale)))
        rows.append(dict(index=int(index), width=float(p[4]),
            capture_checkpoint=float(data['capture_time'][index]),
            section_error=section_error, period_error=period_error,
            scaled_flow_closure=closure,
            agrees_with_checked_cycle=bool(max(section_error, period_error, closure) < 1e-7)))
    result = dict(method='Independent DOP853, rtol=1e-12, atol=1e-14, max_step=.005; four periods per captured endpoint',
        scope='Every classified endpoint is compared with the same upward Poincare section of an independently checked attracting cycle. This does not prove global convergence or a chaotic saddle.',
        captured_count=len(rows), all_passed=all(r['agrees_with_checked_cycle'] for r in rows), rows=rows)
    (ROOT / 'capture-cpu-audit.json').write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}), flush=True)
    if not result['all_passed']:
        raise RuntimeError('At least one periodic capture needs investigation')


def summarize_and_plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    data = np.load(ROOT / 'ensemble.npz')
    times = data['checkpoint_times']
    widths = np.unique(data['parameters'][:, 4])
    groups = []
    fig = plt.figure(figsize=(12, 9), layout='constrained')
    grid = fig.add_gridspec(3, 2, width_ratios=[1, 1.8])
    survival = fig.add_subplot(grid[:2, 0])
    evidence = fig.add_subplot(grid[2, 0])
    evidence.axis('off')
    colors = ['#2680b9', '#cf6945', '#238b75']
    for j, width in enumerate(widths):
        indices = np.flatnonzero(data['parameters'][:, 4] == width)
        statuses = data['status'][indices]
        active = indices[statuses == 0]
        captures = indices[statuses == 1]
        late = data['lle_windows'][-1, active]
        groups.append(dict(width=float(width), starts=len(indices), captured=len(captures),
            still_observed=len(active), unresolved=int(np.sum(statuses > 1)),
            last_block_all_windows_positive=int(np.sum(np.all(late > .01, axis=1))),
            last_block_lle_range=None if not len(active) else [float(np.min(data['exponents'][-1, active])), float(np.max(data['exponents'][-1, active]))],
            latest_capture_checkpoint=None if not len(captures) else float(np.max(data['capture_time'][captures]))))
        counts = np.sum(data['status_history'][:, indices] == 0, axis=1)
        survival.step(np.r_[0, times], np.r_[len(indices), counts], where='post', color=colors[j], label=f'w = {width:g}')
        ax = fig.add_subplot(grid[j, 1])
        order = indices[np.argsort(np.nan_to_num(data['capture_time'][indices], nan=np.inf))]
        array = np.ma.masked_invalid(data['exponents'][:, order].T)
        cmap = plt.get_cmap('RdBu_r').copy()
        cmap.set_bad('#d5d5d5')
        mesh = ax.imshow(array, origin='lower', aspect='auto', interpolation='nearest',
            extent=[0, times[-1], -.5, len(order)-.5], cmap=cmap,
            norm=TwoSlopeNorm(vmin=-.05, vcenter=0, vmax=.4))
        ax.set_title(f'w = {width:g}: {len(captures)}/{len(indices)} checked periodic captures')
        ax.set_ylabel('Initial condition\n(sorted by capture)')
        ax.set_xlabel('Time / thermal relaxation time')
    survival.set(xlabel='Time / thermal relaxation time', ylabel='Trajectories still under observation', ylim=(0, 33),
        title='Fixed cohort; no extrapolation')
    survival.legend()
    survival.grid(alpha=.2)
    evidence.text(0, .96, '96 independent initial conditions\n500-unit observation blocks\n\nGray: integration stopped after\ntwo periodic-capture confirmations.\nEach captured endpoint is checked\nwith an independent CPU solver.\n\nFinite survival does not prove\na chaotic attractor.', va='top', fontsize=10)
    fig.colorbar(mesh, ax=fig.axes[2:], shrink=.8, label='Largest finite-time exponent in each block')
    fig.suptitle('Long chaotic transients in coupled optical operators\nSix-state model with assumed logistic illumination', fontsize=15)
    target = Path('docs/assets/coupled_operators_lifetimes.png')
    fig.savefig(target, dpi=170)
    plt.close(fig)
    result = dict(observation_horizon=float(times[-1]), groups=groups,
        caveat='Cohort counts are descriptive, not basin probabilities. Survivors are censored at the observation horizon; no escape-rate fit or asymptotic classification is made.')
    (ROOT / 'summary.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--audit', action='store_true')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.audit:
        audit_captures()
    if args.plot:
        summarize_and_plot()
