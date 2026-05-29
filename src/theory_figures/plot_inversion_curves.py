import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'data_analysis'))
from FluidProps import FluidProps

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'img')

MEAS_P_RANGE = (1.0, 10.0)    # MPa
MEAS_T_RANGE = (65.0, 80.0)   # K

def inversion_pressure(FP, T, p_low=0.05, p_high=45.0):
    """Find p (MPa) where mu_JT = 0 at fixed T, or NaN if no sign change."""
    try:
        mu_lo = FP.get_JT_coefficient(p=p_low, T=T)
        mu_hi = FP.get_JT_coefficient(p=p_high, T=T)
    except Exception:
        return np.nan
    if not (np.isfinite(mu_lo) and np.isfinite(mu_hi)):
        return np.nan
    if mu_lo * mu_hi > 0:
        return np.nan
    try:
        return brentq(lambda p: FP.get_JT_coefficient(p=p, T=T),
                      p_low, p_high)
    except Exception:
        return np.nan


def plot_inversion_mix():
    FP = FluidProps(('Helium', 'Neon'))
    x_He = np.round(np.arange(0.0, 1.01, 0.1), 2)
    T_grid = np.arange(4.0, 223.5, 2.0)

    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.viridis)

    fig, ax = plt.subplots(figsize=(4.5, 4.0))

    # measurement region (shaded rectangle)
    rect = mpatches.Rectangle(
        (MEAS_P_RANGE[0], MEAS_T_RANGE[0]),
        MEAS_P_RANGE[1] - MEAS_P_RANGE[0],
        MEAS_T_RANGE[1] - MEAS_T_RANGE[0],
        facecolor='0.85', edgecolor='none', alpha=0.7, zorder=0)
    ax.add_patch(rect)

    for x in x_He:
        if 0.0 < x < 1.0:
            FP.set_composition_from_1st_fraction(x)
        elif x == 0.0:
            FP_pure = FluidProps('neon')
            p_inv = np.array([inversion_pressure(FP_pure, T) for T in T_grid])
            m = np.isfinite(p_inv)
            ax.plot(p_inv[m], T_grid[m], color=sm.to_rgba(x), lw=1.0)
            continue
        elif x == 1.0:
            FP_pure = FluidProps('helium')
            p_inv = np.array([inversion_pressure(FP_pure, T) for T in T_grid])
            m = np.isfinite(p_inv)
            ax.plot(p_inv[m], T_grid[m], color=sm.to_rgba(x), lw=1.0)
            continue
        p_inv = np.array([inversion_pressure(FP, T) for T in T_grid])
        m = np.isfinite(p_inv)
        ax.plot(p_inv[m], T_grid[m], color=sm.to_rgba(x), lw=1.0)

    cb = fig.colorbar(sm, ax=ax, pad=0.02, ticks=np.arange(0, 1.05, 0.2))
    cb.set_label(r'$x_{\mathrm{He}}$ / molar')

    ax.set_xlabel(r'$p$ / MPa')
    ax.set_ylabel(r'$T$ / K')
    ax.set_xlim(0, 35)
    ax.set_ylim(0, 250)
    out = os.path.join(OUT_DIR, 'Tinv_mix.pdf')
    fig.savefig(out)


def main():
    plt.style.use(STYLE)
    os.makedirs(OUT_DIR, exist_ok=True)
    plot_inversion_mix()


if __name__ == '__main__':
    main()
    plt.show()
