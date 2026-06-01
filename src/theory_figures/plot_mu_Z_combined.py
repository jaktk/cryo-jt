import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'data_analysis'))
from FluidProps import FluidProps

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')

PRESSURES_MPa = np.arange(0, 21, 4)
FLUIDS = [('helium', '-'), ('neon', '-.')]

YLIM_MU = (-4, 4)
YLIM_Z = (0, 2)


def _isobar(FP, p, T_grid):
    try:
        T_sat = FP.get_saturation_T(p)
    except Exception:
        T_sat = -np.inf

    Ts, mus, Zs = [], [], []
    for T in T_grid:
        if T <= T_sat:
            continue
        try:
            mus.append(FP.get_JT_coefficient(p=p, T=T))
            Zs.append(FP.get_compressibility_factor(p=p, T=T))
            Ts.append(T)
        except Exception:
            continue
    return np.asarray(Ts), np.asarray(mus), np.asarray(Zs)


def plot_mu_Z():
    fig = plt.figure(figsize=(5.4, 4.4))
    fig.set_layout_engine('none')
    gs = GridSpec(
        2, 2, figure=fig,
        width_ratios=[30, 1],
        hspace=0.08, wspace=0.04,
        left=0.12, right=0.90, bottom=0.11, top=0.97,
    )
    ax_mu = fig.add_subplot(gs[0, 0])
    ax_Z = fig.add_subplot(gs[1, 0], sharex=ax_mu)
    cax = fig.add_subplot(gs[:, 1])

    norm = mpl.colors.Normalize(vmin=0, vmax=PRESSURES_MPa.max())
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.viridis)

    T_grid = np.arange(4.0, 300.1, 0.5)
    for fluid_name, ls in FLUIDS:
        FP = FluidProps(fluid_name)
        for p in PRESSURES_MPa:
            Ts, mus, Zs = _isobar(FP, p, T_grid)
            if Ts.size == 0:
                continue
            c = sm.to_rgba(p)
            ax_mu.plot(Ts, mus, color=c, lw=1.0, ls=ls)
            ax_Z.plot(Ts, Zs, color=c, lw=1.0, ls=ls)

    ax_mu.axhline(0, color='0.4', lw=0.6, zorder=0)
    ax_Z.axhline(1, color='0.4', lw=0.6, zorder=0)
    ax_mu.set_xlim(0, 300)
    ax_mu.set_ylim(*YLIM_MU)
    ax_Z.set_ylim(*YLIM_Z)

    ax_mu.tick_params(labelbottom=False)
    ax_Z.set_xlabel(r'$T$/K')
    ax_mu.set_ylabel(r'$\mu_{\mathrm{JT}}$/(K$\cdot$MPa$^{-1}$)')
    ax_Z.set_ylabel(r'$Z$')

    legend_handles = [
        Line2D([0], [0], color='0.3', lw=1.2, ls=ls, label=label)
        for (label, ls) in FLUIDS
    ]
    ax_mu.legend(handles=legend_handles, loc='upper right', frameon=True,
                 edgecolor='black', facecolor='white', framealpha=1.0,
                 fontsize=9)

    cb = fig.colorbar(sm, cax=cax, ticks=PRESSURES_MPa)
    cb.set_label(r'$p$/MPa')

    fig.savefig(os.path.join(os.path.dirname(__file__), '..', '..', 'img', 'mu_Z_He_Ne.pdf'))


def main():
    plt.style.use(STYLE)
    plot_mu_Z()


if __name__ == '__main__':
    main()
    plt.show()
