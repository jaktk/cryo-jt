import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')
MARKERS = ['o', 's', '^', 'D', 'v', '>', '<', 'p', '*', 'h']
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'derived_data', 'jt_coeffs')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'img')
SUMMARY = os.path.join(DATA_DIR, 'jt_coefficients_summary.csv')

K = 2.0  # coverage factor for the plotted error bars

PURE_FLUID_COLS = [
    ('Nitrogen', r'N$_2$'),
    ('Argon',    'Ar'),
    ('Helium',   r'$^4$He'),
]
PURE_YLIMS = {
    'Nitrogen': {'mu': (0.0, 20.0),  'dev': (-5.0, 5.0)},
    'Argon':    {'mu': (3.0, 12.0),  'dev': (-5.0, 5.0)},
    'Helium':   {'mu': (-0.6, 0),  'dev': (-40.0, 40.0)},
}

MIXTURE_LIMS = {
    'HeNe': {'mu': (0.0, 5.0),  'dev': (-20.0, 20.0), 'x': (0.0, 10.0)},
    'HeN2': {'mu': (0.0, 12.0), 'dev': (-5.0, 20.0), 'x': (0.0, 10.0)},
}

def parse_name(path):
    base = os.path.basename(path).replace('_JT.csv', '')
    m = re.match(r'^(.+?)_(\d+)K_(\d+)MPa$', base)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return base, None, None


def load_isenthalps(prefix, summary):
    """ Concatenate all *_JT.csv files for the given fluid prefix """
    files = sorted(glob.glob(os.path.join(DATA_DIR, prefix + '_*_JT.csv')))
    frames = []
    for f in files:
        fluid, T, p = parse_name(f)
        df = pd.read_csv(f)
        base = os.path.basename(f).replace('_JT.csv', '.csv')
        sm_row = summary[summary['filename'] == base]
        if sm_row.empty:
            continue
        sm_row = sm_row.iloc[0]
        if fluid == 'Helium-Neon':
            x_he = float(sm_row['x1'])
        elif fluid == 'Nitrogen-Helium':
            x_he = float(sm_row['x2'])
        else:
            x_he = float('nan')
        df['fluid'] = fluid
        df['T_nom'] = T
        df['p_nom'] = p
        df['x_he'] = x_he
        df['iso'] = f"{T} K, {p} MPa"
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def group_by_composition(df, bins):
    """ Partition the isenthalps into composition columns """
    groups = []
    for label, lo, hi in bins:
        sub = df[(df['x_he'] >= lo) & (df['x_he'] < hi)]
        if not sub.empty:
            groups.append((label, sub))
    return groups


def make_figure(groups, outname, reference='eos', limits=None):
    """Multi-column mixture figure.

    ``limits``: optional dict with keys ``'mu'`` (top-row y limits),
    ``'dev'`` (bottom-row y limits), and ``'x'`` (shared x limits). Each
    value is a ``(lo, hi)`` tuple. Limits are applied to every column.
    """
    n_cols = len(groups)
    if n_cols == 0:
        return

    if n_cols == 1:
        figsize = (4, 6)
    elif n_cols == 2:
        figsize = (4.6, 4.6)
    else:
        figsize = (6.5, 4.3)

    fig, axes = plt.subplots(2, n_cols, figsize=figsize, sharex='col', sharey='row',
                             gridspec_kw={'height_ratios': [2, 1]})
    if n_cols == 1:
        axes = np.array([[axes[0]], [axes[1]]])

    for j, (col_label, df) in enumerate(groups):
        ax_mu, ax_dev = axes[0, j], axes[1, j]
        iso_groups = list(df.groupby('iso', sort=False))
        n_iso = len(iso_groups)
        for i, (iso, g) in enumerate(iso_groups):
            g = g.sort_values('p/MPa')
            if 'used' in g.columns:
                g = g[g['used']]
            if g.empty:
                continue
            c = plt.cm.viridis(i / max(n_iso - 1, 1) * 0.85)
            mk = MARKERS[i % len(MARKERS)]
            ax_mu.plot(g['p/MPa'], g['JT_eos/(K/MPa)'], '-', color=c,
                       lw=1.0, zorder=2)
            ax_mu.errorbar(g['p/MPa'], g['JT_meas/(K/MPa)'],
                           yerr=K * g['JT_UNC/(K/MPa)'],
                           fmt=mk, color=c, ms=5.5, lw=0, elinewidth=0.7,
                           capsize=1.8, mec='black', mew=0.5,
                           label=iso, alpha=0.8, zorder=3)
            denom = g['JT_eos/(K/MPa)'] if reference == 'eos' \
                else g['JT_meas/(K/MPa)']
            dev = 100.0 * (g['JT_meas/(K/MPa)'] - g['JT_eos/(K/MPa)']) / denom
            ax_dev.plot(g['p/MPa'], dev, mk, color=c, ms=5.5,
                        mec='black', mew=0.5, alpha=0.8, zorder=3)

        ax_dev.axhline(0, color='0.2', lw=0.6, zorder=1)
        if col_label is not None and n_cols > 1:
            ax_mu.set_title(col_label, fontsize=9)
        ax_mu.legend(ncol=2 if n_cols == 1 else 1,
                     loc='best',
                     fontsize=9,
                     frameon=True, edgecolor='black',
                     facecolor='white', framealpha=1.0)
        ax_dev.set_xlabel(r'$p$ / MPa')

    axes[0, 0].set_ylabel(r'$\mu_{\mathrm{JT}}$ / K MPa$^{-1}$')
    _ref = (r'\mu_{\mathrm{JT}}^{\mathrm{EOS}}' if reference == 'eos'
            else r'\mu_{\mathrm{JT}}^{\mathrm{meas}}')
    axes[1, 0].set_ylabel(
        r'$(\mu_{\mathrm{JT}}^{\mathrm{meas}}-\mu_{\mathrm{JT}}^{\mathrm{EOS}})/'
        + _ref + r'\times 100$')

    if limits is not None:
        if 'mu' in limits:
            for ax in axes[0, :]:
                ax.set_ylim(*limits['mu'])
        if 'dev' in limits:
            for ax in axes[1, :]:
                ax.set_ylim(*limits['dev'])
        if 'x' in limits:
            for ax in axes[1, :]:
                ax.set_xlim(*limits['x'])

    fig.align_ylabels(axes[:, 0])
    out = os.path.join(OUT_DIR, outname)
    fig.savefig(out)


def make_pure_fluids_figure(summary, outname, ylims=None):
    if ylims is None:
        ylims = PURE_YLIMS

    fig, axes = plt.subplots(
        2, 3, figsize=(7.6, 4.6), sharex='col',
        gridspec_kw={'height_ratios': [2, 1]},
    )

    for j, (fluid_key, title) in enumerate(PURE_FLUID_COLS):
        ax_mu, ax_dev = axes[0, j], axes[1, j]
        ax_mu.set_title(title, fontsize=10)
        ax_dev.set_xlabel(r'$p$ / MPa')
        ax_dev.axhline(0, color='0.2', lw=0.6, zorder=1)

        df = load_isenthalps(fluid_key, summary)
        if not df.empty:
            iso_groups = list(df.groupby('iso', sort=False))
            n_iso = len(iso_groups)
            for i, (iso, g) in enumerate(iso_groups):
                g = g.sort_values('p/MPa')
                if 'used' in g.columns:
                    g = g[g['used']]
                if g.empty:
                    continue
                c = plt.cm.viridis(i / max(n_iso - 1, 1) * 0.85)
                mk = MARKERS[i % len(MARKERS)]
                ax_mu.plot(g['p/MPa'], g['JT_eos/(K/MPa)'], '-', color=c,
                           lw=1.0, zorder=2)
                ax_mu.errorbar(g['p/MPa'], g['JT_meas/(K/MPa)'],
                               yerr=K * g['JT_UNC/(K/MPa)'],
                               fmt=mk, color=c, ms=5.5, lw=0, elinewidth=0.7,
                               capsize=1.8, mec='black', mew=0.5,
                               label=iso, alpha=0.8, zorder=3)
                dev = 100.0 * (g['JT_meas/(K/MPa)'] - g['JT_eos/(K/MPa)']) \
                      / g['JT_eos/(K/MPa)']
                ax_dev.plot(g['p/MPa'], dev, mk, color=c, ms=5.5,
                            mec='black', mew=0.5, alpha=0.8, zorder=3)
            ax_mu.legend(loc='best', fontsize=9, frameon=True,
                         edgecolor='black', facecolor='white', framealpha=1.0)

        if fluid_key in ylims:
            ax_mu.set_ylim(*ylims[fluid_key]['mu'])
            ax_dev.set_ylim(*ylims[fluid_key]['dev'])
            ax_mu.set_xlim(0, 10)

    axes[0, 0].set_ylabel(r'$\mu_{\mathrm{JT}}$ / K MPa$^{-1}$')
    axes[1, 0].set_ylabel(
        r'$(\mu_{\mathrm{JT}}^{\mathrm{meas}}-\mu_{\mathrm{JT}}^{\mathrm{EOS}})'
        r'/\mu_{\mathrm{JT}}^{\mathrm{EOS}} \times 100$'
    )
    fig.align_ylabels(axes[:, 0])
    fig.savefig(os.path.join(OUT_DIR, outname))


def main():
    plt.style.use(STYLE)
    os.makedirs(OUT_DIR, exist_ok=True)
    summary = pd.read_csv(SUMMARY)

    make_pure_fluids_figure(summary, 'measurements_pures_errorbars.pdf')

    # He-Ne
    df_hene = load_isenthalps('Helium-Neon', summary)
    if not df_hene.empty:
        groups = group_by_composition(df_hene, [
            (r'(a) $x_{\mathrm{He}} \approx 0.22$', 0.00, 0.25),
            (r'(b) $x_{\mathrm{He}} \approx 0.30$', 0.25, 0.38),
            (r'(c) $x_{\mathrm{He}} \approx 0.47$', 0.38, 1.01),
        ])
        make_figure(groups, 'measurements_HeNe_errorbars.pdf',
                    reference='meas', limits=MIXTURE_LIMS['HeNe'])

    # He-N2
    df_hen2 = load_isenthalps('Nitrogen-Helium', summary)
    if not df_hen2.empty:
        groups = group_by_composition(df_hen2, [
            (r'(a) $x_{\mathrm{He}} \approx 0.15$', 0.00, 0.30),
            (r'(b) $x_{\mathrm{He}} \approx 0.50$', 0.30, 1.01),
        ])
        make_figure(groups, 'measurements_HeN2_errorbars.pdf',
                    reference='meas', limits=MIXTURE_LIMS['HeN2'])


if __name__ == '__main__':
    main()
    plt.show()
