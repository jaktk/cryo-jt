"""
Heat-exchange and kinetic-energy corrections to mu_JT, evaluated for every
point used in the mu_JT evaluation. Neither is applied to the data; both enter
U(mu_JT) as uncorrected biases (calculate_jt_coefficients.py).

The measured path coefficient relates to mu_JT as

    mu_JT = (dT/dp)_path - Q / (m c_p dp) + de_k / (c_p dp),

with Q the heat gained by the fluid and dp = p_out - p_in (eq 7 of the paper).
The mass flow rate m is not measured; it is interpolated at the measured outlet
pressure from the capillary design calculation of each isenthalp
(heat_exchange.GRIDS, D = 0.365 mm, L = 6 m).

  * Heat exchange: Q is modelled per point from the measured shield and outlet
    temperatures (heat_exchange.py). Because mu_JT is the local slope of the
    fitted isenthalp, its bias is the slope of the outlet warming Q/(m c_p)
    along the isenthalp, refitted with the degree and weights of the data
    (heat_exchange.heat_bias), not the secant term of eq 7.
  * Kinetic energy: de_k = (u2^2 - u1^2) / 2 with u = m / (rho A), and c_p the
    mean of the inlet and outlet values (secant term of eq 7).

Outputs (data/derived_data/jt_coeffs):
  * corrections_per_point.csv  per point: modelled heat gain, heat-exchange bias
                               and kinetic-energy correction in K/MPa and % of
                               |mu_JT|, and the offset of the outlet temperature
                               from the EOS isenthalp through the inlet state
  * corrections_summary.csv    per fluid group: range of the heat gain and mean
                               |corrections|; dk/% is KINETIC_PERC in
                               calculate_jt_coefficients.py
  * heat_exchange_fit.csv      per isenthalp, a fit of the outlet-temperature
                               offsets against the EOS (diagnostic, not used in U)

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
from FluidProps import FluidProps
from chebyshev_fit import chebyshev_fit_weights
from heat_exchange import GRIDS, heat_bias, mass_flow
from measurement_uncertainty import K_TARGET, point_uncertainties
from repo_root import REPO_ROOT

AREA = np.pi * (0.365e-3)**2 / 4    # m^2, equivalent internal diameter

ROOT = REPO_ROOT
JT_DIR = os.path.join(ROOT, 'data', 'derived_data', 'jt_coeffs')
PT_DIR = os.path.join(ROOT, 'data', 'derived_data', 'p_T_pairs')

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
    summary = pd.read_csv(os.path.join(JT_DIR, 'jt_coefficients_summary.csv')).set_index('filename')
    rows = []
    for stem in GRIDS:
        fluid = stem.split('_')[0]
        meta = summary.loc[f'{stem}.csv']
        jt = pd.read_csv(os.path.join(JT_DIR, f'{stem}_JT.csv'))
        pt = pd.read_csv(os.path.join(PT_DIR, f'{stem}.csv'))
        fp = FluidProps(FLUIDS[fluid])
        x_he = None
        if fluid in ('Helium-Neon', 'Nitrogen-Helium'):
            fp.set_composition_from_1st_fraction(meta['x1'])
            x_he = meta['x1'] if fluid == 'Helium-Neon' else 1.0 - meta['x1']
        p_unc, t_unc = point_uncertainties(pt)
        weights = chebyshev_fit_weights(pt['PT102/MPa'].values, pt['TT102/K'].values,
                                        p_unc / K_TARGET, t_unc / K_TARGET)
        bias_q, q_heat = heat_bias(stem, pt, fp, int(meta['polynomial_degree']), weights)
        for i, q in pt.iterrows():
            r = jt.loc[(jt['p/MPa'] - q['PT102/MPa']).abs().idxmin()]
            if not r['used']:
                continue
            p1, t1, p2, t2 = q['PT101/MPa'], q['TT101/K'], q['PT102/MPa'], q['TT102/K']
            m = mass_flow(stem, p2)
            cp = np.mean([fp.get_cp(p1, t1), fp.get_cp(p2, t2)])
            dp = (p2 - p1) * 1e6   # p_out - p_in, Pa
            u1 = m / fp.get_rhomass(p1, t1) / AREA
            u2 = m / fp.get_rhomass(p2, t2) / AREA
            # offset of the measured outlet temperature from the EOS isenthalp
            # through the measured inlet state of the point
            t2_eos = fp.get_T(p=p2, h=fp.get_hmass(p=p1, T=t1))
            rows.append({
                'isenthalp': stem, 'group': group_name(fluid, x_he),
                'p/MPa': p2, 'm/(g/s)': m * 1e3, 'Q/mW': 1e3 * q_heat[i],
                'JT_meas/(K/MPa)': r['JT_meas/(K/MPa)'], 'JT_eos/(K/MPa)': r['JT_eos/(K/MPa)'],
                'dQ/(K/MPa)': bias_q[i],
                'dk/(K/MPa)': (u2**2 - u1**2) / 2 / (cp * dp) * 1e6,
                'T_offset/K': t2 - t2_eos,
            })
    df = pd.DataFrame(rows)
    df['dQ/%'] = 100 * df['dQ/(K/MPa)'].abs() / df['JT_meas/(K/MPa)'].abs()
    df['dk/%'] = 100 * df['dk/(K/MPa)'].abs() / df['JT_meas/(K/MPa)'].abs()
    return df


def heat_exchange_fit():
    """ Per-isenthalp fit of the outlet-temperature offsets (module docstring). """
    summary = pd.read_csv(os.path.join(JT_DIR, 'jt_coefficients_summary.csv')).set_index('filename')
    rows = []
    for stem in GRIDS:
        fluid = stem.split('_')[0]
        pt = pd.read_csv(os.path.join(PT_DIR, f'{stem}.csv')).sort_values('PT102/MPa').reset_index(drop=True)
        jt = pd.read_csv(os.path.join(JT_DIR, f'{stem}_JT.csv')).sort_values('p/MPa').reset_index(drop=True)
        fp = FluidProps(FLUIDS[fluid])
        if fluid in ('Helium-Neon', 'Nitrogen-Helium'):
            fp.set_composition_from_1st_fraction(summary.loc[f'{stem}.csv', 'x1'])
        p1, t1 = pt['PT101/MPa'].values, pt['TT101/K'].values
        p2, t2 = pt['PT102/MPa'].values, pt['TT102/K'].values
        m = mass_flow(stem, p2)
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


def main():
    df = per_point_corrections()
    df.to_csv(os.path.join(JT_DIR, 'corrections_per_point.csv'), index=False)

    summary = df.groupby('group', sort=False).agg(**{
        'Q_min/mW': ('Q/mW', 'min'),
        'Q_max/mW': ('Q/mW', 'max'),
        'dQ/(K/MPa)': ('dQ/(K/MPa)', lambda s: s.abs().mean()),
        'dQ/%': ('dQ/%', 'mean'),
        'dk/(K/MPa)': ('dk/(K/MPa)', lambda s: s.abs().mean()),
        'dk/%': ('dk/%', 'mean'),
    })
    summary.to_csv(os.path.join(JT_DIR, 'corrections_summary.csv'))
    pd.set_option('display.width', 160)
    print(summary.round(3))

    fit = heat_exchange_fit()
    fit.to_csv(os.path.join(JT_DIR, 'heat_exchange_fit.csv'), index=False)
    print(fit.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
