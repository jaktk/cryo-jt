import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import Chebyshev

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')
DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'derived_data', 'jt_coeffs', 'Helium-Neon_65K_5MPa_JT.csv')
OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'img', 'toc.pdf')


def main():
    plt.style.use(STYLE)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    mpl.rcParams.update({
        "font.size": 8,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })
    df = pd.read_csv(DATA).sort_values('p/MPa')
    p = df['p/MPa'].to_numpy()
    T = df['T/K'].to_numpy()

    # smooth Chebyshev fit of the isenthalp for the curve
    deg = min(3, len(p) - 1)
    fit = Chebyshev.fit(p, T, deg)
    pp = np.linspace(p.min(), p.max(), 200)
    TT = fit(pp)

    fig, ax = plt.subplots(figsize=(3.25, 1.75))

    ax.plot(pp, TT, '-', color='#0072B2', lw=1.6, zorder=2)
    ax.plot(p, T, 'o', color='#0072B2', ms=4.5, mec='black', mew=0.5,
            zorder=3)

    # expansion / cooling arrow, offset above the isenthalp (high p -> low p)
    off = 0.9  # K, vertical offset so the arrow does not overlap the data
    ax.annotate('', xy=(pp[20], TT[20] + off), xytext=(pp[-20], TT[-20] + off),
                arrowprops=dict(arrowstyle='-|>', color='#D55E00', lw=1.8))
    ax.text(0.46, 0.78, 'throttling: $p\\downarrow,\\ T\\downarrow$',
            transform=ax.transAxes, color='#D55E00', fontsize=7.5,
            ha='center', va='center', rotation=14)

    ax.text(0.97, 0.30,
            r'$\mu_{\rm JT}=\left(\dfrac{\partial T}{\partial p}\right)_{h}$',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9)
    ax.text(0.03, 0.92, r'$^4$He-Ne, 65 K', transform=ax.transAxes,
            fontsize=8, color='0.25', va='top')

    ax.set_xlabel(r'$p$ / MPa', fontsize=8.5, labelpad=1)
    ax.set_ylabel(r'$T$ / K', fontsize=8.5, labelpad=1)
    ax.tick_params(labelsize=7.5, pad=2)
    ax.margins(0.08)

    fig.savefig(OUT, bbox_inches='tight')


if __name__ == '__main__':
    main()
    plt.show()
