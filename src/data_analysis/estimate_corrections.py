"""
Heat-loss and kinetic-energy corrections to the measured path coefficient,

    mu_JT = (dT/dp)_path - Q / (m c_p dp) + de_k / (c_p dp),

with Q the heat gained by the fluid and dp = p_out - p_in (eq 7 of the paper).

evaluated for every point used in the mu_JT evaluation.

The mass flow rate m is not measured; it is interpolated at the measured outlet
pressure from the capillary design calculations of each isenthalp
(data/derived_data/capillary_mass_flow, D = 0.365 mm, L = 6 m). The heat
exchange is a constant in-leak Q = 27 mW, and
de_k = (u2^2 - u1^2) / 2 with u = m / (rho A). c_p is the mean of the inlet and
outlet values.

Outputs (data/derived_data/jt_coeffs):
  * corrections_per_point.csv  per-point corrections in K/MPa and % of |mu_JT|,
                               and the offset of the outlet temperature from the
                               EOS isenthalp through the inlet state
  * corrections_summary.csv    mean |corrections| per fluid group, the values of
                               CORRECTIONS_PERC in calculate_jt_coefficients.py
  * heat_exchange_fit.csv      per isenthalp, a fit of the outlet-temperature
                               offsets (diagnostic, not used in U, see below)
and prints the deviation from the EOS with and without the corrections, and the
heat exchange that best fits each fluid group.

Heat-exchange diagnostic. For every point of an isenthalp, the measured outlet
temperature minus the EOS temperature on the isenthalp through the measured
inlet state is fitted as
    delta * c_p,in / c_p,out + Q / (m c_p,out),
an offset delta of the inlet thermometer TT101 and a constant heat Q gained by
the fluid. TT102 sits in a copper block downstream of the capillary, where the
gas has slowed down, so it reads close to the stagnation temperature and the
kinetic-energy term of the capillary exit is not part of the fit. The fit
reports delta, Q, the RMS of its residuals, and the AAD of mu_JT after removing
the slope of the fitted terms.

Run after calculate_jt_coefficients.py.
"""
import os
import numpy as np
import pandas as pd
from scipy import optimize
from FluidProps import FluidProps
from repo_root import REPO_ROOT

Q_IN = 27e-3                        # W, heat gained by the fluid
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
            dp = (p2 - p1) * 1e6   # p_out - p_in, Pa
            u1 = m / fp.get_rhomass(p1, t1) / AREA
            u2 = m / fp.get_rhomass(p2, t2) / AREA
            # offset of the measured outlet temperature from the EOS isenthalp
            # through the measured inlet state of the point
            t2_eos = fp.get_T(p=p2, h=fp.get_hmass(p=p1, T=t1))
            rows.append({
                'isenthalp': stem, 'group': group_name(fluid, x_he),
                'p/MPa': p2, 'm/(g/s)': m * 1e3,
                'JT_meas/(K/MPa)': r['JT_meas/(K/MPa)'], 'JT_eos/(K/MPa)': r['JT_eos/(K/MPa)'],
                # correction terms of eq 7, added to the path coefficient:
                # heat term per watt of Q, and kinetic-energy term
                'dQ_per_W/(K/MPa/W)': -1e6 / (m * cp * dp),
                'dk/(K/MPa)': (u2**2 - u1**2) / 2 / (cp * dp) * 1e6,
                'T_offset/K': t2 - t2_eos,
            })
    df = pd.DataFrame(rows)
    df['dQ/(K/MPa)'] = Q_IN * df['dQ_per_W/(K/MPa/W)']
    df['dQ/%'] = 100 * df['dQ/(K/MPa)'].abs() / df['JT_meas/(K/MPa)'].abs()
    df['dk/%'] = 100 * df['dk/(K/MPa)'].abs() / df['JT_meas/(K/MPa)'].abs()
    return df


def heat_exchange_fit():
    """ Per-isenthalp fit of the outlet-temperature offsets (module docstring). """
    summary = pd.read_csv(os.path.join(JT_DIR, 'jt_coefficients_summary.csv')).set_index('filename')
    rows = []
    for stem, grid_file in GRIDS.items():
        fluid = stem.split('_')[0]
        pt = pd.read_csv(os.path.join(PT_DIR, f'{stem}.csv')).sort_values('PT102/MPa').reset_index(drop=True)
        jt = pd.read_csv(os.path.join(JT_DIR, f'{stem}_JT.csv')).sort_values('p/MPa').reset_index(drop=True)
        grid = pd.read_csv(os.path.join(MF_DIR, grid_file)).sort_values('p2/MPa')
        fp = FluidProps(FLUIDS[fluid])
        if fluid in ('Helium-Neon', 'Nitrogen-Helium'):
            fp.set_composition_from_1st_fraction(summary.loc[f'{stem}.csv', 'x1'])
        p1, t1 = pt['PT101/MPa'].values, pt['TT101/K'].values
        p2, t2 = pt['PT102/MPa'].values, pt['TT102/K'].values
        m = np.interp(p2, grid['p2/MPa'], grid['m/g_per_s']) / 1e3
        cp_in = np.array([fp.get_cp(a, b) for a, b in zip(p1, t1)])
        cp_out = np.array([fp.get_cp(a, b) for a, b in zip(p2, t2)])
        t2_eos = np.array([fp.get_T(p=b, h=fp.get_hmass(p=a, T=c)) for a, b, c in zip(p1, p2, t1)])
        y = t2 - t2_eos
        A = np.column_stack([cp_in / cp_out, 1.0 / (m * cp_out)])
        (delta, q), *_ = np.linalg.lstsq(A, y, rcond=None)
        fitted = A @ np.array([delta, q])
        # slope of the fitted offset terms along the isenthalp (K/MPa)
        bias = np.gradient(fitted, p2)
        used = jt['used'].values
        mu, mu_eos = jt['JT_meas/(K/MPa)'].values, jt['JT_eos/(K/MPa)'].values
        mu_fit = mu - bias
        mixture = '-' in fluid
        rows.append({
            'isenthalp': stem, 'n_points': len(p2),
            'delta/K': delta, 'Q/mW': 1e3 * q,
            'rms/mK': 1e3 * np.sqrt(np.mean((y - fitted)**2)),
            'AAD_perc': 100 * np.mean(np.abs((mu - mu_eos) / (mu if mixture else mu_eos))[used]),
            'AAD_fit_removed_perc': 100 * np.mean(np.abs((mu_fit - mu_eos) / (mu_fit if mixture else mu_eos))[used]),
        })
    return pd.DataFrame(rows)


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
    df['dev_corr'] = deviation(df, Q_IN).abs()
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

    fit = heat_exchange_fit()
    fit.to_csv(os.path.join(JT_DIR, 'heat_exchange_fit.csv'), index=False)
    print(fit.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
