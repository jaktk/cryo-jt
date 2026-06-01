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
    combined_temperature_uncertainty,
    _conv_uncert_for_slope,
    U_P,
)

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'derived_data', 'jt_coeffs')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'img')

FLUID_STYLE = {
    'Nitrogen':       {'label': r'N$_2$',       'color': '#D62728', 'marker': 'o'},
    'Argon':          {'label': 'Ar',           'color': '#FF7F0E', 'marker': 's'},
    'Helium':         {'label': 'He',           'color': '#1F77B4', 'marker': '^'},
    'Helium-Neon':    {'label': 'He-Ne',        'color': '#2CA02C', 'marker': 'D'},
    'Nitrogen-Helium':{'label': r'N$_2$-He',    'color': '#9467BD', 'marker': 'v'},
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
    Compute theoretical (conventional) relative uncertainty at each
    measured point using the actual pressure drop and temperatures.
    """
    theory_pct = []
    for _, row in data.iterrows():
        # Use the actual measured T and p for that point
        T_out = row['T/K']
        p_out = row['p/MPa']
        mu_meas = row['JT_meas/(K/MPa)']
        abs_mu = abs(mu_meas)
        T_unc = row['T_UNC/K']
        p_unc = row['p_UNC/MPa']
        U_T = combined_temperature_uncertainty(T_out) # temperature measuremnets uncertainty

        if abs_mu > 1e-6:
            # Use equipment uncertainties with a representative dp
            # matching the theoretical analysis (p_in=10, p_out=1.5)
            dp_ref = 8.5  # same as theoretical analysis
            u_conv = _conv_uncert_for_slope(abs_mu, dp_ref, T_out) * 100
        else:
            u_conv = np.nan

        theory_pct.append(u_conv)

    data['theory_unc_pct'] = theory_pct
    return data


def make_figure(data):
    fig, ax = plt.subplots(figsize=(5, 4))

    # Theoretical uncertainty curve
    slopes = np.concatenate([
        np.linspace(0.02, 0.5, 30),
        np.linspace(0.6, 2.0, 20),
        np.linspace(2.5, 8.0, 20),
        np.linspace(9.0, 22.0, 15),
    ])

    # Band for T_in = 65..160 K
    dp_ref = 8.5
    conv_65 = np.array([_conv_uncert_for_slope(s, dp_ref, 65) * 100 for s in slopes])
    conv_160 = np.array([_conv_uncert_for_slope(s, dp_ref, 160) * 100 for s in slopes])
    conv_80 = np.array([_conv_uncert_for_slope(s, dp_ref, 80) * 100 for s in slopes])

    ax.fill_between(slopes, conv_65, conv_160, color='0.80', alpha=0.4, zorder=1, label='Theoretical sensor uncertainty')
    ax.plot(slopes, conv_80, 'k-', lw=1.2, zorder=2, label=r'Conventional ($T_{\rm in}$ = 80 K)')

    # Measured deviations
    for fluid, style in FLUID_STYLE.items():
        mask = data['fluid'] == fluid
        if not mask.any():
            continue
        sub = data[mask]
        ax.scatter(sub['abs_mu_eos'], sub['rel_abs_err_perc'],
                   c=style['color'], marker=style['marker'],
                   s=40, alpha=0.8, zorder=5,
                   edgecolor='k', linewidths=0.7,
                   label=style['label'])

    # Annotation
    n2_mask = (data['fluid'] == 'Nitrogen') & (data['T_nom'] == 160)
    if n2_mask.any():
        sub = data[n2_mask]
        idx = sub['abs_mu_eos'].sub(7.0).abs().idxmin()
        x_pt = sub.loc[idx, 'abs_mu_eos']
        y_meas = sub.loc[idx, 'rel_abs_err_perc']
        y_theo = sub.loc[idx, 'theory_unc_pct']
        if y_meas > y_theo and not np.isnan(y_theo):
            ax.annotate('',
                        xy=(x_pt, y_theo), xytext=(x_pt, y_meas),
                        arrowprops=dict(arrowstyle='<->', color='0.3',
                                        lw=1.2, shrinkA=2, shrinkB=2))
            ax.text(x_pt + 0.3, np.sqrt(y_meas * y_theo), 'systematic\nerror gap', color='0.3', va='center')

    # Formatting
    ax.set_xlabel(r'$|\mu_{\rm JT}^{\rm EOS}|$/(K$\cdot$MPa$^{-1}$)')
    ax.set_ylabel(r'$(\mu_{\rm JT}^{\rm EOS} - \mu_{\rm JT}^{\rm meas}) / \mu_{\rm JT}^{\rm EOS} \times 100$ at $k=2$')
    ax.set_yscale('log')
    ax.set_xlim(0, 20)
    ax.set_ylim(0.05, 100)

    # Fluid regions
    ax.axvspan(0, 0.5, alpha=0.08, color='C0')
    ax.axvspan(4, 8, alpha=0.06, color='C3')
    ax.text(0.8, 0.06, 'He', color='C0', ha='center', weight='bold')
    ax.text(6.0, 0.06, r'N$_2$, Ar', color='C3', ha='center', weight='bold')

    # Reference lines
    ax.axhline(5, color='0.5', ls=':', lw=0.7)
    ax.axhline(2, color='0.5', ls=':', lw=0.7)
    ax.text(18, 5.3, r'5 \%', color='0.4', va='bottom')
    ax.text(18, 2.15, r'2 \%', color='0.4', va='bottom')

    ax.legend(loc='upper right', ncol=2, frameon=True, fontsize=11,
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
    fig.savefig(outpath, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
    plt.show()

