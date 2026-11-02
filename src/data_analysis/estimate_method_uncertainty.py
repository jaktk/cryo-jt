"""
Cross-check of the method uncertainty u_m against the pure-fluid validation.

u_m (U_METHOD_REL in calculate_jt_coefficients.py) is set by the synthetic
validation in validate_degree_criterion.py, without the measured deviations.
This script checks it against the N2 and Ar points, whose deviations from their
reference equations of state (Span et al. 2000; Tegeler et al. 1999) are
compared with the modelled uncertainty terms. The excess not explained by the
Monte Carlo, bias (heat exchange, kinetic energy, and purity), and EOS terms,

    u_excess^2 = RMS(dev)^2 - mean(u_MC^2 + u_sys^2 + u_EOS^2),

with all terms relative and at k = 1, should not exceed u_m. Neither EOS
publication states an uncertainty for mu_JT; u_EOS is taken from their isobaric
heat capacity tolerance diagrams (Span et al. Fig. 41, Tegeler et al. Fig. 44),
since mu_JT = (T (dv/dT)_p - v) / c_p.

The script also counts the N2 and Ar points within the expanded uncertainty
and bootstraps u_excess over the isenthalps.

Run after calculate_jt_coefficients.py with the final constants.
"""
import glob
import os
import numpy as np
import pandas as pd
from repo_root import REPO_ROOT
from calculate_jt_coefficients import K_TARGET, U_METHOD_REL

N_BOOTSTRAP = 4000
SEED = 20210101


def eos_rel_std_uncertainty(fluid, p_mpa):
    """ Relative standard uncertainty (k = 1) of c_p from the EOS tolerance diagrams (expanded values / 2). """
    if fluid == 'Nitrogen':
        return 0.8 / 100 / 2   # region C, (100 to 160) K below 30 MPa
    if fluid == 'Argon':
        return (2.0 if p_mpa >= 5.5 else 0.3) / 100 / 2   # region C above about 5.5 MPa at (155 to 175) K, else A
    raise ValueError(fluid)


def excess(d):
    """ Signed sqrt of RMS(dev)^2 minus the mean modelled variance (negative if over-explained). """
    v = np.mean(d['dev']**2) - np.mean(d['u_mc']**2 + d['u_sys']**2 + d['u_eos']**2)
    return np.sign(v) * np.sqrt(abs(v))


def main():
    root = REPO_ROOT
    rows = []
    for fluid in ('Nitrogen', 'Argon'):
        for f in sorted(glob.glob(os.path.join(root, 'data', 'derived_data', 'jt_coeffs', f'{fluid}_*_JT.csv'))):
            df = pd.read_csv(f)
            df = df[df['used'] == True]
            mu_eos = df['JT_eos/(K/MPa)'].values
            mu = df['JT_meas/(K/MPa)'].values
            rows.append(pd.DataFrame({
                'isenthalp': os.path.basename(f)[:-len('_JT.csv')],
                'dev': (mu - mu_eos) / mu_eos,
                'u_mc': df['JT_UNC_MC/(K/MPa)'].values / K_TARGET / np.abs(mu),
                'u_sys': np.hypot(df['JT_UNC_Q/(K/MPa)'].values, df['JT_UNC_PUR/(K/MPa)'].values) / K_TARGET / np.abs(mu),
                'u_eos': [eos_rel_std_uncertainty(fluid, p) for p in df['p/MPa'].values],
                'U_rel': df['JT_UNC/(K/MPa)'].values / np.abs(mu),
            }))
    d = pd.concat(rows, ignore_index=True)
    rms = np.sqrt(np.mean(d['dev']**2))
    explained = np.mean(d['u_mc']**2 + d['u_sys']**2 + d['u_eos']**2)
    print(f"N2 + Ar points: {len(d)}")
    print(f"mean relative deviation: {100 * d['dev'].mean():+.2f} %")
    print(f"RMS relative deviation:  {100 * rms:.2f} %")
    print(f"RMS of modelled terms:   {100 * np.sqrt(explained):.2f} % "
          f"(MC {100 * np.sqrt(np.mean(d['u_mc']**2)):.2f}, bias {100 * np.sqrt(np.mean(d['u_sys']**2)):.2f}, "
          f"EOS {100 * np.sqrt(np.mean(d['u_eos']**2)):.2f})")
    print(f"unexplained excess (k = 1): {100 * excess(d):+.2f} %  (negative: over-explained)")

    rng = np.random.default_rng(SEED)
    isenthalps = d['isenthalp'].unique()
    boot = np.array([excess(pd.concat([d[d['isenthalp'] == i]
                                       for i in rng.choice(isenthalps, len(isenthalps))]))
                     for _ in range(N_BOOTSTRAP)])
    lo, hi = np.percentile(boot, [5, 95])
    print(f"bootstrap over isenthalps: 90 % interval ({100 * lo:+.2f} to {100 * hi:+.2f}) %, "
          f"excess <= 0 in {100 * np.mean(boot <= 0):.0f} % of resamples")

    within = np.abs(d['dev']) <= d['U_rel']
    print(f"pooled u_m (U_METHOD_REL, k = 1): {100 * U_METHOD_REL:.2f} %")
    print(f"points within U (k = 2): {within.sum()} of {len(d)}; "
          f"mean U_r {100 * d['U_rel'].mean():.2f} %")


if __name__ == '__main__':
    main()
