import sys
import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(__file__))
from theoretical_jt_uncertainty import (
    conventional_relative_uncertainty,
    _conv_uncert_for_slope,
)

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'derived_data', 'jt_coeffs')
PT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'derived_data', 'p_T_pairs')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'img')

FLUID_STYLE = {
    'Nitrogen':       {'label': r'N$_2$',       'color': '#D62728', 'marker': 'o'},
    'Argon':          {'label': 'Ar',           'color': '#FF7F0E', 'marker': 's'},
    'Helium':         {'label': 'He',           'color': '#1F77B4', 'marker': '^'},
    'Helium-Neon':    {'label': 'He-Ne',        'color': '#2CA02C', 'marker': 'D'},
    'Nitrogen-Helium':{'label': r'He-N$_2$',    'color': '#9467BD', 'marker': 'v'},
}


def parse_filename(fname):
    """ Extract fluid name, temperature, and pressure from filename """
    base = os.path.basename(fname).replace('_JT.csv', '')
    m = re.match(r'^(.+?)_(\d+)K_(\d+)MPa$', base) # Pattern: Fluid_TK_PMPa  e.g. Nitrogen_160K_6MPa
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return base, None, None


def load_all_measurements():
    """ Load all JT coefficient CSV files and return combined DataFrame """
    pattern = os.path.join(DATA_DIR, '*_JT.csv')
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No *_JT.csv files found in {DATA_DIR}")

    rows = []
    for f in files:
        fluid, T_nom, p_nom = parse_filename(f)
        df = pd.read_csv(f)
        # inlet state of each point, from the p-T pairs (same row order)
        pt = pd.read_csv(os.path.join(PT_DIR, os.path.basename(f).replace('_JT.csv', '.csv')))
        df['p_in/MPa'] = pt['PT101/MPa'].values
        df['T_in/K'] = pt['TT101/K'].values
        # keep only the interior points used for the mu_JT evaluation
        # (the isenthalp extremities are excluded from the analysis)
        if 'used' in df.columns:
            df = df[df['used']]
        df['fluid'] = fluid
        df['T_nom'] = T_nom
        df['p_nom'] = p_nom
        df['label'] = f"{fluid} {T_nom}K {p_nom}MPa"
        rows.append(df)

    data = pd.concat(rows, ignore_index=True)
    # Compute |mu_JT| from measured values
    data['abs_mu_meas'] = data['JT_meas/(K/MPa)'].abs()
    data['abs_mu_eos'] = data['JT_eos/(K/MPa)'].abs()
    return data


def compute_theoretical_uncertainty_per_point(data):
    """
    Conventional (sensor-level) relative uncertainty of each measured point,
    from its own pressure drop p_in - p_out and its inlet and outlet
    temperatures.
    """
    theory_pct = []
    for _, row in data.iterrows():
        abs_mu = abs(row['JT_meas/(K/MPa)'])
        T_in, T_out = row['T_in/K'], row['T/K']
        dp = row['p_in/MPa'] - row['p/MPa']
        if abs_mu > 1e-6 and dp > 0 and T_in != T_out:
            u_conv = conventional_relative_uncertainty(
                mu_jt=abs_mu, delta_T=abs(T_in - T_out), delta_p=dp,
                T_in=T_in, T_out=T_out) * 100
        else:
            u_conv = np.nan
        theory_pct.append(u_conv)

    data['theory_unc_pct'] = theory_pct
    return data


def make_figure(data):
    fig, ax = plt.subplots(figsize=(3.1, 3.0))

    # Theoretical uncertainty curve
    slopes = np.concatenate([
        np.linspace(0.02, 0.5, 30),
        np.linspace(0.6, 2.0, 20),
        np.linspace(2.5, 8.0, 20),
        np.linspace(9.0, 22.0, 15),
    ])

    # Band for T_in = 65..160 K
    dp_ref = 8.5
    conv_65 = np.array([_conv_uncert_for_slope(s, dp_ref, 65) for s in slopes])
    conv_160 = np.array([_conv_uncert_for_slope(s, dp_ref, 160) for s in slopes])
    conv_80 = np.array([_conv_uncert_for_slope(s, dp_ref, 80) for s in slopes])

    ax.fill_between(slopes, conv_65, conv_160, color='0.80', alpha=0.4, zorder=1, label='Theoretical sensor uncertainty')
    ax.plot(slopes, conv_80, 'k-', lw=1.2, zorder=2, label=r'Conventional ($T_{\rm in}$ = 80 K)')

    # Measured deviations
    for fluid, style in FLUID_STYLE.items():
        mask = data['fluid'] == fluid
        if not mask.any():
            continue
        sub = data[mask]
        ax.scatter(sub['abs_mu_eos'], sub['rel_abs_err_perc'] / 100,
                   c=style['color'], marker=style['marker'],
                   s=18, alpha=0.8, zorder=5,
                   edgecolor='k', linewidths=0.5,
                   label=style['label'])

    # Annotation
    n2_mask = (data['fluid'] == 'Nitrogen') & (data['T_nom'] == 160)
    if n2_mask.any():
        sub = data[n2_mask]
        idx = sub['abs_mu_eos'].sub(7.0).abs().idxmin()
        x_pt = sub.loc[idx, 'abs_mu_eos']
        y_meas = sub.loc[idx, 'rel_abs_err_perc'] / 100
        y_theo = sub.loc[idx, 'theory_unc_pct'] / 100
        if y_meas > y_theo and not np.isnan(y_theo):
            ax.annotate('',
                        xy=(x_pt, y_theo), xytext=(x_pt, y_meas),
                        arrowprops=dict(arrowstyle='<->', color='0.3',
                                        lw=1.2, shrinkA=2, shrinkB=2))
            ax.text(x_pt - 0.3, 0.002, 'systematic\nerror gap', color='0.3', ha='right', va='center')

    # Formatting
    ax.set_xlabel(r'$|\mu_{\rm JT}^{\rm EOS}|$/(K$\cdot$MPa$^{-1}$)')
    ax.set_ylabel(r'$|\mu_{\rm JT}^{\rm meas} - \mu_{\rm JT}^{\rm EOS}| / \mu_{\rm JT}^{\rm EOS}$, $U_{\rm r}(\mu_{\rm JT})$')
    ax.set_yscale('log')
    ax.set_xlim(0, 20)
    ax.set_ylim(1e-5, 1)

    # Fluid regions
    ax.axvspan(0, 0.5, alpha=0.08, color='C0')
    ax.axvspan(4, 8, alpha=0.06, color='C3')
    ax.text(0.6, 0.55, 'He', color='C0', ha='left', weight='bold')
    ax.text(6.0, 0.55, r'N$_2$, Ar', color='C3', ha='center', weight='bold')

    # Reference lines
    ax.axhline(0.05, color='0.5', ls=':', lw=0.7)
    ax.axhline(0.02, color='0.5', ls=':', lw=0.7)
    ax.text(19.7, 0.053, '0.05', color='0.4', ha='right', va='bottom')
    ax.text(19.7, 0.0215, '0.02', color='0.4', ha='right', va='bottom')

    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=True,
              edgecolor='black', facecolor='white', framealpha=1.0)

    return fig


def print_summary_table(data):
    """ Print a summary table comparing theory and experiment """
    print(f"\n{'Measurement':<30s} {'|μ_JT|':>10s} {'Theory':>10s} "
          f"{'Actual':>10s} {'Gap':>8s}")
    print(f"{'':30s} {'(K/MPa)':>10s} {'Ur (%)':>10s} "
          f"{'err (%)':>10s} {'factor':>8s}")
    print('-' * 70)

    for label in sorted(data['label'].unique()):
        sub = data[data['label'] == label]
        mu_range = f"{sub['abs_mu_eos'].min():.1f}-{sub['abs_mu_eos'].max():.1f}"
        theory = sub['theory_unc_pct'].mean()
        actual = sub['rel_abs_err_perc'].mean()
        gap = actual / theory if theory > 0 else np.inf
        print(f"{label:<30s} {mu_range:>10s} {theory:10.2f} "
              f"{actual:10.1f} {gap:8.1f}x")


def main():
    plt.style.use(STYLE)
    data = load_all_measurements()
    data = compute_theoretical_uncertainty_per_point(data)
    print_summary_table(data)
    fig = make_figure(data)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    outpath = os.path.join(OUTPUT_DIR, 'theory_vs_experiment_uncertainty.pdf')
    fig.savefig(outpath)


if __name__ == '__main__':
    main()

