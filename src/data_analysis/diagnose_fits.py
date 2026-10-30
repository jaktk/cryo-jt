"""
Diagnostics of the Chebyshev fits quoted in SI Section S4. Neither enters the
results.

1. Sensitivity of the chi-square degree rule to the noise model. The rule
   treats the whole sensor-chain uncertainty as independent per point. Part
   of it, the Cernox calibration and the calibration polynomial T(R), is a
   smooth error shared by the points of an isenthalp. The rule is repeated
   with only the window statistics and the acquisition-module (CABTR) term in
   u(T), and only the window statistics in u(p).
2. EOS-only residual: the EOS isenthalp through the mean inlet state at the
   measured pressures, fitted with the selected degree and the same weights.
   Its residual RMS is the truncation error of the polynomial alone.

Output: data/derived_data/jt_coeffs/fit_diagnostics.csv
Run after calculate_jt_coefficients.py.
"""
import os
import numpy as np
import pandas as pd
from FluidProps import FluidProps
from repo_root import REPO_ROOT
from chebyshev_fit import effective_t_uncertainty, fit_chebyshev, select_degree
from measurement_uncertainty import (K_TARGET, STUDENT_SCALE, cabtr_uncertainty,
                                     point_uncertainties)

JT_DIR = os.path.join(REPO_ROOT, 'data', 'derived_data', 'jt_coeffs')
PT_DIR = os.path.join(REPO_ROOT, 'data', 'derived_data', 'p_T_pairs')


def main():
    summary = pd.read_csv(os.path.join(JT_DIR, 'jt_coefficients_summary.csv'))
    rows = []
    for _, meta in summary.iterrows():
        df = pd.read_csv(os.path.join(PT_DIR, meta['filename']))
        p = df['PT102/MPa'].to_numpy(dtype=float)
        T = df['TT102/K'].to_numpy(dtype=float)
        p_unc, t_unc = point_uncertainties(df)
        u_p, u_T = p_unc / K_TARGET, t_unc / K_TARGET
        fit = select_degree(p, T, u_p, u_T)

        # uncorrelated part only (standard uncertainties)
        u_p_rand = STUDENT_SCALE * df['PT102/MPa_EXP_UNC'].to_numpy(dtype=float) / 2
        u_T_rand = np.sqrt((STUDENT_SCALE * df['TT102/K_EXP_UNC'].to_numpy(dtype=float) / 2)**2
                           + (cabtr_uncertainty(T) * 1e-3 / 2)**2)
        fit_rand = select_degree(p, T, np.maximum(u_p_rand, 1e-9), u_T_rand)

        # EOS-only residual at the selected degree
        fluid = meta['fluid']
        fp = FluidProps(tuple(fluid.split('-')) if '-' in fluid else fluid)
        if '-' in fluid:
            fp.set_composition_from_1st_fraction(meta['x1'])
        h_in = fp.get_hmass(p=df['PT101/MPa'].mean(), T=df['TT101/K'].mean())
        T_eos = np.array([fp.get_T(p=pi, h=h_in) for pi in p])
        u_eff = effective_t_uncertainty(p, T_eos, u_p, u_T)
        poly = fit_chebyshev(p, T_eos, fit.degree, 1.0 / u_eff)
        s_t_eos = np.sqrt(np.mean((T_eos - poly(p))**2))

        rows.append({
            'isenthalp': meta['filename'][:-4], 'n_points': len(p),
            'degree': fit.degree, 'chi2_red': fit.chi2_red, 'chi2_pass': fit.passed,
            'degree_uncorr': fit_rand.degree, 'chi2_red_uncorr': fit_rand.chi2_red,
            'chi2_pass_uncorr': fit_rand.passed,
            's_T/mK': 1e3 * fit.s_T, 's_T_eos/mK': 1e3 * s_t_eos,
        })
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(JT_DIR, 'fit_diagnostics.csv'), index=False)
    pd.set_option('display.width', 200)
    print(out.round(2).to_string(index=False))


if __name__ == '__main__':
    main()
