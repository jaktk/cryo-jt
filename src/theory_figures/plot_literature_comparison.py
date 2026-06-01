import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'data_analysis'))
from FluidProps import FluidProps

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')
DATA = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'literature', 'literature_measurements.csv')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'img')

AUTHOR_COLORS = {
    'Roebuck':              '#009E73',  # green
    'Tkaczuk':              '#CC79A7',  # magenta
    'Smith':                '#D55E00',  # red-orange
    '(King, Potter)':       '#0072B2',  # blue
    '(Bier,Ernst,Maurer)':  '#E69F00',  # orange
}

FLUID_MARKERS = {
    'propylene':       's',
    'nitrogen':        'o',
    'helium':          'D',
    'argon':           'p',
    'n-butane':        'v',
    'air':             '^',
    'helium-nitrogen': '<',
    'helium-neon':     '>',
}

FLUID_LABELS = {
    'propylene':       'propylene',
    'nitrogen':        'nitrogen',
    'helium':          'helium',
    'argon':           'argon',
    'n-butane':        r'$n$-butane',
    'air':             'air',
    'helium-nitrogen': r'helium-nitrogen',
    'helium-neon':     r'helium-neon',
}


def relative_error(group):
    fluid = group['fluid'].iloc[0]
    if fluid.lower() in ('helium-neon', 'helium-nitrogen'):
        components = tuple(c.capitalize() for c in fluid.split('-'))
        FP = FluidProps(components)
        x_mean = float(np.mean(group['x1']))
        FP.set_composition_from_1st_fraction(x_mean)
    else:
        FP = FluidProps(fluid)
    mu_calc = np.array([FP.get_JT_coefficient(p=p, T=T)
                        for p, T in zip(group['p/MPa'], group['T/K'])])
    rel = group['rel_err']
    if rel.isnull().all():
        return (mu_calc - group['uJT_meas/(K/MPa)']) / mu_calc * 100.0
    return rel * 100.0


def main():
    plt.style.use(STYLE)
    os.makedirs(OUT_DIR, exist_ok=True)
    df = pd.read_csv(DATA)

    fig, ax = plt.subplots(figsize=(6, 3))

    ax.add_patch(mpatches.Rectangle((0, -5), 600, 10, facecolor='0.88', edgecolor='none', zorder=-2))
    ax.axhline(0, color='k', lw=0.6, zorder=-1)

    seen_fluids = set()
    legend_handles = []

    for (author, fluid), group in df.groupby(['authors', 'fluid'], sort=False):
        try:
            err = relative_error(group)
        except Exception as e:
            print(f"  skip {author} / {fluid}: {e}")
            continue
        color = AUTHOR_COLORS.get(author, '0.4')
        marker = FLUID_MARKERS.get(fluid, 'x')
        ax.scatter(group['T/K'], err,
                   c=color, edgecolor='k', linewidths=0.7,
                   s=60, alpha=0.85, marker=marker,
                   zorder=0 if author == 'Roebuck' else 1)
        if fluid not in seen_fluids:
            seen_fluids.add(fluid)
            legend_handles.append(Line2D(
                [0], [0], marker=marker, color='k', linestyle='',
                markerfacecolor='w', markeredgecolor='k', markersize=6,
                label=FLUID_LABELS.get(fluid, fluid)))

    ax.set_yscale('symlog', linthresh=1)
    ax.set_xlim(0, 600)
    ax.set_ylim(-1000, 1000)
    ax.yaxis.set_major_locator(FixedLocator([-1000, -100, -10, -1, 0, 1, 10, 100, 1000]))
    ax.yaxis.set_major_formatter(FixedFormatter([r'$-10^3$', r'$-10^2$', r'$-10^1$', r'$-10^0$', '0',
                                                 r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$']))
    ax.set_xlabel(r'$T$/K')
    ax.set_ylabel(r'$(\mu_{\rm JT}^{\rm EOS} - \mu_{\rm JT}^{\rm meas}) / \mu_{\rm JT}^{\rm EOS} \times 100$')
    ax.legend(handles=legend_handles,
              loc='lower center', bbox_to_anchor=(0.455, 1.02), ncol=4,
              frameon=True, edgecolor='black', facecolor='white',
              framealpha=1.0, handletextpad=0.1,
              labelspacing=0.1, columnspacing=0.7)

    out = os.path.join(OUT_DIR, 'literature_review_plus_this_work.pdf')
    fig.savefig(out)


if __name__ == '__main__':
    main()
    plt.show()
