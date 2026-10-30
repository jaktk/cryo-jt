"""
Heat-loss and kinetic-energy corrections to the measured path coefficient,

    mu_JT = (dT/dp)_path - Q / (m c_p dp) + de_k / (c_p dp),

evaluated for every point used in the mu_JT evaluation.

The mass flow rate m is not measured; it is interpolated at the measured outlet
pressure from the capillary design calculations of each isenthalp
(data/derived_data/capillary_mass_flow, D = 0.365 mm, L = 6 m). The heat
exchange is a constant Q = 27 mW (positive when the fluid loses heat), and
de_k = (u2^2 - u1^2) / 2 with u = m / (rho A). c_p is the mean of the inlet and
outlet values.

Outputs (data/derived_data/jt_coeffs):
  * corrections_per_point.csv  per-point corrections in K/MPa and % of |mu_JT|
  * corrections_summary.csv    mean |corrections| per fluid group, the values of
                               CORRECTIONS_PERC in calculate_jt_coefficients.py
and prints the deviation from the EOS with and without the corrections, and the
heat exchange that best fits each fluid group.

Run after calculate_jt_coefficients.py.
"""
import os
import numpy as np
import pandas as pd
from scipy import optimize
from FluidProps import FluidProps
from repo_root import REPO_ROOT

Q_LOSS = 27e-3                      # W
AREA = np.pi * (0.365e-3)**2 / 4    # m^2, equivalent internal diameter

ROOT = REPO_ROOT
JT_DIR = os.path.join(ROOT, 'data', 'derived_data', 'jt_coeffs')
PT_DIR = os.path.join(ROOT, 'data', 'derived_data', 'p_T_pairs')
MF_DIR = os.path.join(ROOT, 'data', 'derived_data', 'capillary_mass_flow')

# Capillary design calculation of each measured isenthalp
GRIDS = {
    'Nitrogen_150K_6MPa': 'N2_D0.365_L6_p6.0_T149.8.csv',
    'Nitrogen_160K_5MPa': 'N2_D0.365_L6_p5.0_T160.0.csv',
    'Nitrogen_160K_6MPa': 'N2_D0.365_L6_p6.0_T160.0.csv',
    'Nitrogen_160K_7MPa': 'N2_D0.365_L6_p7.0_T160.0.csv',
    'Nitrogen_160K_9MPa': 'N2_D0.365_L6_p9.0_T160.0.csv',
    'Nitrogen_160K_12MPa': 'N2_D0.365_L6_p12.0_T160.0.csv',
    'Argon_180K_5MPa': 'Ar_D0.365_L6_p5.0_T180.0.csv',
    'Argon_180K_12MPa': 'Ar_D0.365_L6_p12.0_T180.0.csv',
    'Helium_65K_5MPa': 'He_D0.365_L6_p5.0_T65.0.csv',
    'Helium_65K_7MPa': 'He_D0.365_L6_p7.0_T65.0.csv',
    'Helium_141K_10MPa': 'He_D0.365_L6_p10.0_T140.6.csv',
    'Helium-Neon_65K_5MPa': 'HeNe_D0.365_L6_p5.0_T65.0_0.215He.csv',
    'Helium-Neon_65K_7MPa': 'HeNe_D0.365_L6_p7.0_T65.0_0.213He.csv',
    'Helium-Neon_65K_8MPa': 'HeNe_D0.365_L6_p8.0_T65.0_0.467He.csv',
    'Helium-Neon_65K_10MPa': 'HeNe_D0.365_L6_p10.0_T65.1_0.311He.csv',
    'Helium-Neon_66K_10MPa': 'HeNe_D0.365_L6_p10.0_T65.7_0.461He.csv',
    'Helium-Neon_80K_7MPa': 'HeNe_D0.365_L6_p7.0_T80.0_0.467He.csv',
    'Helium-Neon_80K_8MPa': 'HeNe_D0.365_L6_p8.0_T80.4_0.327He.csv',
    'Helium-Neon_80K_10MPa': 'HeNe_D0.365_L6_p10.0_T79.9_0.467He.csv',
    'Nitrogen-Helium_140K_8MPa': 'HeN2_D0.365_L6_p8.0_T140.0_0.503He.csv',
    'Nitrogen-Helium_160K_5MPa': 'HeN2_D0.365_L6_p5.0_T160.0_0.145He.csv',
    'Nitrogen-Helium_160K_8MPa': 'HeN2_D0.365_L6_p8.0_T160.0_0.145He.csv',
}

FLUIDS = {
    'Nitrogen': 'Nitrogen', 'Argon': 'Argon', 'Helium': 'Helium',
    'Helium-Neon': ('Helium', 'Neon'), 'Nitrogen-Helium': ('Nitrogen', 'Helium'),
}


def group_name(fluid, x_he=None):
    """ Fluid group of the corrections table; He-N2 is split by composition. """
    if fluid == 'Nitrogen-Helium':
        return 'Nitrogen-Helium-0.14' if x_he < 0.3 else 'Nitrogen-Helium-0.50'
    return fluid


def per_point_corrections():
    rows = []
    for stem, grid_file in GRIDS.items():
        fluid = stem.split('_')[0]
        jt = pd.read_csv(os.path.join(JT_DIR, f'{stem}_JT.csv'))
        pt = pd.read_csv(os.path.join(PT_DIR, f'{stem}.csv'))
        grid = pd.read_csv(os.path.join(MF_DIR, grid_file)).sort_values('p2/MPa')
        fp = FluidProps(FLUIDS[fluid])
        x_he = None
        if fluid in ('Helium-Neon', 'Nitrogen-Helium'):
            x1 = jt['x1'].mean()
            fp.set_composition_from_1st_fraction(x1)
            x_he = x1 if fluid == 'Helium-Neon' else 1.0 - x1
        for _, r in jt[jt['used'] == True].iterrows():
            q = pt.loc[(pt['PT102/MPa'] - r['p/MPa']).abs().idxmin()]
            p1, t1, p2, t2 = q['PT101/MPa'], q['TT101/K'], q['PT102/MPa'], q['TT102/K']
            m = np.interp(p2, grid['p2/MPa'], grid['m/g_per_s']) / 1e3
            cp = np.mean([fp.get_cp(p1, t1), fp.get_cp(p2, t2)])
            dp = (p1 - p2) * 1e6
            u1 = m / fp.get_rhomass(p1, t1) / AREA
            u2 = m / fp.get_rhomass(p2, t2) / AREA
            rows.append({
                'isenthalp': stem, 'group': group_name(fluid, x_he),
                'p/MPa': p2, 'm/(g/s)': m * 1e3,
                'JT_meas/(K/MPa)': r['JT_meas/(K/MPa)'], 'JT_eos/(K/MPa)': r['JT_eos/(K/MPa)'],
                # heat-loss correction per watt of Q, K/MPa per W
                'dQ_per_W/(K/MPa/W)': -1e6 / (m * cp * dp),
                'dk/(K/MPa)': (u2**2 - u1**2) / 2 / (cp * dp) * 1e6,
            })
    df = pd.DataFrame(rows)
    df['dQ/(K/MPa)'] = Q_LOSS * df['dQ_per_W/(K/MPa/W)']
    df['dQ/%'] = 100 * df['dQ/(K/MPa)'].abs() / df['JT_meas/(K/MPa)'].abs()
    df['dk/%'] = 100 * df['dk/(K/MPa)'].abs() / df['JT_meas/(K/MPa)'].abs()
    return df


def deviation(df, q):
    """ Relative deviation in % from the EOS, referenced like the AAD of the paper. """
    mu = df['JT_meas/(K/MPa)'] + q * df['dQ_per_W/(K/MPa/W)'] + (df['dk/(K/MPa)'] if q else 0)
    mixture = df['group'].str.contains('-')
    ref = np.where(mixture, mu, df['JT_eos/(K/MPa)'])
    return 100 * (mu - df['JT_eos/(K/MPa)']) / ref


def main():
    df = per_point_corrections()
    df.to_csv(os.path.join(JT_DIR, 'corrections_per_point.csv'), index=False)

    summary = df.groupby('group', sort=False).agg(**{
        'dQ/(K/MPa)': ('dQ/(K/MPa)', lambda s: s.abs().mean()),
        'dQ/%': ('dQ/%', 'mean'),
        'dk/(K/MPa)': ('dk/(K/MPa)', lambda s: s.abs().mean()),
        'dk/%': ('dk/%', 'mean'),
    })
    df['dev_raw'] = deviation(df, 0.0).abs()
    df['dev_corr'] = deviation(df, Q_LOSS).abs()
    summary['AAD_raw/%'] = df.groupby('group', sort=False)['dev_raw'].mean()
    summary['AAD_corrected/%'] = df.groupby('group', sort=False)['dev_corr'].mean()
    best = {}
    for g, sub in df.groupby('group', sort=False):
        res = optimize.minimize_scalar(
            lambda q: ((sub['JT_meas/(K/MPa)'] + q * sub['dQ_per_W/(K/MPa/W)'] + sub['dk/(K/MPa)']
                        - sub['JT_eos/(K/MPa)'])**2).sum(),
            bounds=(-0.2, 0.2), method='bounded')
        best[g] = res.x * 1e3
    summary['Q_best/mW'] = pd.Series(best)
    summary.to_csv(os.path.join(JT_DIR, 'corrections_summary.csv'))
    pd.set_option('display.width', 160)
    print(summary.round(3))


if __name__ == '__main__':
    main()
