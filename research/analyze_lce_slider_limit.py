"""Close a periodic control and render the normalized relay-limit surveys."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root

from lce_slider import return_map
from check_lce_slider import dop853_crossing


def control():
    data = np.load('results/lce-slider/relay-bias/scan.npz')
    index = np.flatnonzero(data['parameters'][:, 3] == .5)[0]
    p = data['parameters'][index]
    initial = data['final_sections'][index]
    def residual(s):
        end, _, jac = return_map(s, p, 512, True)
        return end-s, jac-np.eye(2)
    solved = root(residual, initial, jac=True, tol=1e-11)
    s = solved.x
    end, period, jac = return_map(s, p, 512, True)
    middle, t1, j1 = dop853_crossing(s, p, True, 1e-12, True)
    independent, t2, j2 = dop853_crossing(middle, p, False, 1e-12, True)
    values = np.linalg.eigvals(jac)
    result = dict(parameters=p.tolist(), section=s.tolist(), period=float(period),
        transverse_multipliers=[[float(x.real), float(x.imag)] for x in values],
        transverse_exponents=(np.log(abs(values))/period).tolist(),
        closure_error=float(max(abs(end-s))),
        dop853_closure_error=float(max(abs(independent-s))),
        dop853_period_error=float(abs(t1+t2-period)),
        dop853_jacobian_error=float(np.max(abs(j2@j1-jac))),
        half_period_symmetry_error=float(max(abs(middle-np.array([-s[0], 1-s[1]])))),
        interpretation='Numerically stable one-return symmetric periodic orbit of the normalized I=0 limit; not a validated existence proof.')
    assert result['closure_error'] < 1e-10
    assert result['dop853_closure_error'] < 1e-8
    assert result['dop853_jacobian_error'] < 1e-7
    assert max(abs(values)) < 1
    out = Path('results/lce-slider/relay-limit-control')
    out.mkdir(parents=True, exist_ok=True)
    (out/'summary.json').write_text(json.dumps(result, indent=2))
    return result


def render():
    intensity = json.loads(Path('results/lce-slider/intensity-limit/summary.json').read_text())
    bias = json.loads(Path('results/lce-slider/relay-bias/summary.json').read_text())
    trace = np.load('results/lce-slider/intensity-limit-cpu/trace.npz')['section']
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for ax in axes.flat:
        ax.grid(alpha=.2)
    ax = axes[0, 0]
    levels = np.array([g['intensity'] for g in intensity['groups']])
    values = np.array([g['lle'] for g in intensity['groups']])
    mean = values.mean(axis=1)
    ax.errorbar(levels, mean, yerr=[mean-values.min(axis=1), values.max(axis=1)-mean],
                fmt='o', capsize=3, color='#276eab', ms=4)
    ax.axhline(0, color='black', lw=.8)
    ax.set_xscale('symlog', linthresh=1e-6, linscale=1.2)
    ax.set_xlim(-.15e-6, 15)
    ax.set_xticks([0, 1e-5, .001, .1, 10])
    ax.set_xticklabels(['0 (limit)', '$10^{-5}$', '$10^{-3}$', '$10^{-1}$', '10'])
    ax.set(xlabel='Intensity I; threshold retuned to fixed delta', ylabel='Leading transverse exponent / time',
           title='A  Weak-light continuation: 16 starts per sample')
    ax = axes[0, 1]
    levels = np.array([g['delta'] for g in bias['groups']])
    values = np.array([g['lle'] for g in bias['groups']])
    mean = values.mean(axis=1)
    ax.errorbar(levels, mean, yerr=[mean-values.min(axis=1), values.max(axis=1)-mean],
                fmt='.', capsize=2, color='#276eab', ms=4)
    selected = np.array([g['retained'] > 0 for g in bias['groups']])
    ax.scatter(levels[selected], mean[selected], c='#c94242', s=18, zorder=5, label='All four windows positive')
    ax.axhline(0, color='black', lw=.8)
    ax.set(xlabel='Bias delta (other half related by exact reflection)', ylabel='Leading transverse exponent / time',
           title='B  Limit I=0: 102 offsets, 8 starts each', xlim=(0, .51))
    ax.legend(fontsize=8, loc='upper right')
    ax = axes[1, 0]
    ax.scatter(trace[:, 1], trace[:, 2], s=2, c='#276eab', alpha=.4, rasterized=True)
    ax.set(xlabel='Velocity V at upward crossing U=0', ylabel='Photo state Y at crossing',
           title='C  Independent DOP853: 4,000 late returns')
    ax = axes[1, 1]
    ax.scatter(trace[:-1, 2], trace[1:, 2], s=2, c='#276eab', alpha=.4, rasterized=True)
    low, high = trace[:, 2].min(), trace[:, 2].max()
    ax.plot([low, high], [low, high], color='gray', lw=.7)
    ax.set(xlabel='$Y_n$', ylabel='$Y_{n+1}$', title='D  Return-map projection (not a proven 1D map)')
    fig.suptitle('Light-switching slider: a biased relay limit\nNormalized coordinates; physical amplitude vanishes as I approaches zero', fontsize=14)
    path = Path('docs/assets/lce_slider_relay_limit.png')
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(path)


if __name__ == '__main__':
    print(json.dumps(control(), indent=2))
    render()
