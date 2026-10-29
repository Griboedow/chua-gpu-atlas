"""Render original evidence figures from saved numerical results; no GPU needed."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from synchronematic_observables import collective_order, pair_signal

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/synchronematic'
OUT = ROOT / 'docs/assets'


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    cpu = DATA / 'n4-cpu-long'
    if (DATA / 'n4-cpu-extended/spectrum.json').exists():
        cpu = DATA / 'n4-cpu-extended'
    a = np.load(cpu / 'trace.npz')
    meta = json.loads((cpu / 'spectrum.json').read_text())
    state, time, history, pos = a['state'], a['time'], a['history'], a['positions']
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.facecolor': '#f9fbfd', 'axes.facecolor': '#ffffff', 'savefig.facecolor': '#f9fbfd'})
    fig = plt.figure(figsize=(13, 9), layout='constrained')
    gs = fig.add_gridspec(2, 2)
    ax = fig.add_subplot(gs[0, 0])
    for i in range(4):
        for j in range(i):
            d = np.linalg.norm(pos[i] - pos[j])
            ax.plot(pos[[i, j], 0], pos[[i, j], 1], color='#bfd1db', lw=1 + 2 / d**3)
    ax.scatter(pos[:, 0], pos[:, 1], s=180, color=['#167d9a', '#e36b37', '#8b5ab8', '#2e9b72'], zorder=3)
    for i, xy in enumerate(pos):
        ax.annotate(str(i + 1), xy, ha='center', va='center', color='white', weight='bold')
    ax.set(aspect='equal', xlabel='Normalized horizontal position', ylabel='Normalized vertical position',
           title='Four fixed centers, eight angular variables')
    ax.set_aspect('equal', adjustable='datalim')
    ax.text(.02, .98, 'Coupling scales as distance$^{-3}$\n7 dimensions after removing\nthe common phase',
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(facecolor='white', edgecolor='none', alpha=.9))
    ax = fig.add_subplot(gs[0, 1])
    colors = ['#c34e32', '#227d9c', '#8176b5']
    for i in range(3):
        ax.plot(history[:, 0], history[:, i + 1], label=fr'$\lambda_{i+1}$', color=colors[i])
    ax.axhline(0, color='#465466', lw=.8)
    ax.set(xlabel='Observation time s (after discarded transient)', ylabel='Running Lyapunov estimate',
           title='Independent DOP853 + full tangent-space QR')
    ax.legend(loc='upper right', ncol=3)
    spectrum = np.array(meta['spectrum'])
    ax.text(.02, .03, f'Final largest exponent: {spectrum[0]:.4f}\nSecond exponent approaches the neutral flow direction',
            transform=ax.transAxes, fontsize=9)
    ax = fig.add_subplot(gs[1, 0], projection='3d')
    xyz = np.column_stack([pair_signal(state, 0, 1, odd=True), pair_signal(state, 0, 2, odd=True), np.sin(2 * state[:, 4])])
    stride = max(1, len(xyz) // 12000)
    ax.scatter(*xyz[::stride].T, c=time[::stride], cmap='viridis', s=.5, alpha=.45, rasterized=True)
    ax.set(xlabel=r'$P_{12}$', ylabel=r'$P_{13}$', zlabel=r'$\sin(2\chi_1)$',
           title='Gauge-invariant projection of the 7D trajectory')
    ax.text2D(.02, .98, r'$P_{ij}=\sin(\theta_j-\theta_i)\cos(\chi_j-\chi_i)$', transform=ax.transAxes, fontsize=9, va='top')
    ax.view_init(elev=24, azim=38)
    ax = fig.add_subplot(gs[1, 1])
    span = (time >= time[-1] - 60)
    r, nematic = collective_order(state)
    s = abs(nematic)
    ax.plot(time[span] - time[span][0], r[span], color='#167d9a', label='Coherence C', lw=1)
    ax.plot(time[span] - time[span][0], s[span], color='#e36b37', label='Nematic order S', lw=1)
    ax.set(xlabel='Time within the final observation segment', ylabel='Order parameter', ylim=(-.02, 1.02),
           title='Persistent fluctuations in collective order')
    ax.legend(loc='lower left', ncol=2)
    fig.suptitle('Synchronematic chaos candidate: literal printed-equation variant', fontsize=15)
    fig.savefig(OUT / 'synchronematic_evidence.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout='constrained')
    for ax, name, label in [(axes[0], 'geometry', 'Maximum coordinate displacement / original minimum spacing'),
                             (axes[1], 'coefficients', 'Maximum relative coefficient perturbation')]:
        groups = json.loads((DATA / f'n4-{name}/groups.json').read_text())
        x = np.array([v['level'] for v in groups])
        rate = np.array([v['retained'] / v['count'] for v in groups])
        ax.plot(x, rate, 'o-', color='#167d9a')
        for xx, rr, item in zip(x, rate, groups):
            ax.annotate(f"{item['retained']}/{item['count']}", (xx, rr), xytext=(0, 7), textcoords='offset points',
                        ha='center', fontsize=9)
        ax.set(xlabel=label, ylabel='Fraction retaining the finite-time criterion', ylim=(-.04, 1.14), title=name.title() + ' sensitivity')
        ax.grid(alpha=.15)
    ax = axes[2]
    for name, label, color in [('solenoidal', 'Solenoidal residual scale', '#167d9a'),
                                ('nonreciprocity', 'Even phase terms scale', '#e36b37')]:
        file = DATA / f'n4-{name}/groups.json'
        if file.exists():
            groups = json.loads(file.read_text())
            ax.plot([v['level'] for v in groups], [v['median_lle'] for v in groups], '.-', label=label, color=color, ms=3)
    ax.axhline(0, color='#333', lw=.6)
    ax.axvline(1, color='#777', ls=':', lw=.8)
    ax.set(xlabel='Control value (1 = published model)', ylabel='Median largest exponent', title='Controlled changes to the equations')
    ax.legend(fontsize=8)
    fig.suptitle('Robustness experiments: finite ensembles, not probabilities or rigorous stability bounds', fontsize=13)
    fig.savefig(OUT / 'synchronematic_robustness.png', dpi=180)
    plt.close(fig)

    a = np.load(ROOT / 'results/photofeedback/stability/curves.npz')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    axes[0].semilogx(a['tau'], a['eigenvalue_real'], color='#167d9a')
    axes[0].axhline(0, color='#333', lw=.7)
    axes[0].set(xlabel='Assumed photoresponse time constant', ylabel='Largest real part of equilibrium eigenvalues',
                title='Instability at intermediate response times')
    for tau, color in [('0.05', '#8176b5'), ('0.5', '#e36b37'), ('3.0', '#167d9a')]:
        mask = a['time'] > 1480
        axes[1].plot(a['time'][mask] - 1480, a['state_' + tau][mask, 0], label=fr'$\tau={tau}$', color=color, lw=1.1)
    axes[1].set(xlabel='Late time segment', ylabel='Dimensionless height', title='Oscillation between the stability boundaries')
    axes[1].legend()
    fig.suptitle('Optical-feedback mechanism test (uncalibrated model)', fontsize=13)
    fig.savefig(OUT / 'photofeedback_stability.png', dpi=180)
    plt.close(fig)
    transition = DATA / 'n4-transition/scan.npz'
    if transition.exists():
        a = np.load(transition)
        labels = a['experiment_labels'][:, 0]
        fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
        for ax in axes:
            ax.scatter(labels, a['lle'], s=6, alpha=.45, color='#167d9a')
            ax.axhline(0, color='#333', lw=.7)
            ax.set(xlabel='Residual scale η (1 = printed model)', ylabel='Largest finite-time exponent')
        axes[0].set(title='126 levels, eight identical starts per level')
        axes[0].axvline(1, color='#777', ls=':', lw=.8)
        axes[1].set(xlim=(1.052, 1.101), title='Regular windows inside the irregular region')
        axes[1].axvspan(1.074, 1.078, color='#269b72', alpha=.15)
        fig.suptitle('Controlled equation change: gradient part + η × divergence-free residual', fontsize=13)
        fig.savefig(OUT / 'synchronematic_transition.png', dpi=180)
        plt.close(fig)
    print('Saved evidence figures to', OUT)


if __name__ == '__main__':
    main()
