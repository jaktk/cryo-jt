"""
Hysteresis check: the nitrogen isenthalp at 150 K/6 MPa was measured with
decreasing outlet pressure (isenthalp 6 of extraction_windows.json, the
published one) and then with increasing outlet pressure (the
"hysteresis_check" stamps, sharing the lowest point).

Both runs are averaged over their steady-state windows as in get_pT_pairs.py,
fitted with the pipeline rule (chebyshev_fit.select_degree), and compared:
the outlet temperatures at the same set points, and mu_JT of the two fits at
the interior points of the downward run.

Output: data/derived_data/jt_coeffs/hysteresis_check.csv
No REFPROP needed.
"""
import json
import os
import numpy as np
import pandas as pd
from repo_root import REPO_ROOT
from chebyshev_fit import select_degree
from measurement_uncertainty import K_TARGET, point_uncertainties

RAW_DIR = os.path.join(REPO_ROOT, 'data', 'raw_data')
JT_DIR = os.path.join(REPO_ROOT, 'data', 'derived_data', 'jt_coeffs')
COLUMNS = ['PT102/MPa', 'PT101/MPa', 'TT101/K', 'TT102/K']


def average_windows(stamps):
    """ Mean, standard deviation and expanded uncertainty of the mean (k = 2)
    over each window, as in get_pT_pairs.py. """
    rows = []
    for stamp in stamps:
        df = pd.read_csv(os.path.join(RAW_DIR, f"{stamp['date/yyyy-mm-dd']}.csv"))
        mp = df.loc[(df['time/hh:mm:ss'] >= stamp['time1/hh:mm:ss'])
                    & (df['time/hh:mm:ss'] <= stamp['time2/hh:mm:ss'])]
        row = {}
        for col in COLUMNS:
            row[col] = mp[col].mean()
            row[f'{col}_STD'] = mp[col].std()
            row[f'{col}_EXP_UNC'] = 2.0 * row[f'{col}_STD'] / np.sqrt(len(mp))
        rows.append(row)
    return pd.DataFrame(rows).sort_values('PT102/MPa').reset_index(drop=True)


def fit(df):
    p_unc, t_unc = point_uncertainties(df)
    return select_degree(df['PT102/MPa'].values, df['TT102/K'].values,
                         p_unc / K_TARGET, t_unc / K_TARGET)


def main():
    windows = json.load(open(os.path.join(RAW_DIR, 'extraction_windows.json')))
    check = windows['hysteresis_check']
    down_stamps = next(i['stamps'] for i in windows['isenthalps'] if i['id'] == check['downward_isenthalp_id'])
    down = average_windows(down_stamps)
    up = average_windows(check['upward_stamps'])
    fit_down, fit_up = fit(down), fit(up)

    p = down['PT102/MPa'].values
    mu_down = fit_down.polynomial.deriv()(p)
    mu_up = fit_up.polynomial.deriv()(p)
    interior = np.zeros(len(p), dtype=bool)
    interior[1:-1] = True
    out = pd.DataFrame({
        'p_down/MPa': p,
        'p_up/MPa': up['PT102/MPa'].values,
        'T_down/K': down['TT102/K'].values,
        'T_up/K': up['TT102/K'].values,
        'T_in_down/K': down['TT101/K'].values,
        'T_in_up/K': up['TT101/K'].values,
        'JT_down/(K/MPa)': mu_down,
        'JT_up/(K/MPa)': mu_up,
        'used': interior,
    })
    out['dT/mK'] = 1e3 * (out['T_up/K'] - out['T_down/K'])
    out['dJT_rel_perc'] = 100 * (out['JT_up/(K/MPa)'] - out['JT_down/(K/MPa)']) / out['JT_down/(K/MPa)']
    out.to_csv(os.path.join(JT_DIR, 'hysteresis_check.csv'), index=False)

    print(f"degree down {fit_down.degree} (chi2_red {fit_down.chi2_red:.2f}), "
          f"up {fit_up.degree} (chi2_red {fit_up.chi2_red:.2f})")
    print(out.round(4).to_string(index=False))
    used = out[out['used']]
    print(f"\nmax |T_up - T_down|: {used['dT/mK'].abs().max():.1f} mK; "
          f"max |relative difference of mu_JT| at the interior points: {used['dJT_rel_perc'].abs().max():.2f} %")


if __name__ == '__main__':
    main()
