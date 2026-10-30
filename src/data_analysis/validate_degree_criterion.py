"""
Synthetic validation of the Chebyshev degree rule and estimate of the method
uncertainty u_m.

For every published isenthalp, the measured pressures are placed on the EOS
isenthalp through the mean inlet state. Gaussian noise at the standard (k = 1)
measurement uncertainty of each point is added to p and T, the degree is
selected, and mu_JT at the interior points is compared with the exact EOS
slope. This is repeated over N_REAL realisations for three degree rules:

  * chi2   the rule of chebyshev_fit.select_degree used by the pipeline;
  * r2     the earlier rule, lowest degree (capped at N_s - 1) with R^2 > 0.999,
           with the earlier normalised 1/u^2 weights;
  * ftest  a nested F-test, raising n while the next degree lowers the weighted
           residual sum of squares at the 5 % level (n <= N_s - 2).

U_MC(n) is the k = 2 Monte Carlo uncertainty at fixed degree n, computed as in
the pipeline. The method uncertainty u_m is the smallest relative standard
uncertainty for which U = 2 sqrt(u_MC^2 + u_m^2) covers the exact slope at
COVERAGE_GOAL of all interior points under the chi2 rule. It covers the fit and
differentiation error that the fixed-degree Monte Carlo misses. The measured
deviations from the EOS are not used, so the N2 and Ar data remain an
independent test of U.

A second part applies the three rules to the measured data and reports the
selected degree and the interior deviation from the EOS, for the sensitivity
statement in the SI. That part is not used to choose anything.

Outputs (data/derived_data/jt_coeffs):
  * degree_validation.csv   per isenthalp and rule: degrees selected, bias,
                            scatter, RMS of the relative error, coverage
  * degree_sensitivity.csv  per isenthalp and rule: degree and AAD on the
                            measured data
  * degree_pooled.csv       pooled synthetic results per fluid group and rule
  * si_degree_tables.tex    the two tables of SI Section S4
and prints the pooled results and u_m (the constant U_METHOD_REL in
calculate_jt_coefficients.py).

Run after calculate_jt_coefficients.py (reads jt_coefficients_summary.csv for
the list of isenthalps). Takes about 10 minutes.
"""
import os
import numpy as np
import pandas as pd
from numpy.polynomial.chebyshev import Chebyshev
from scipy import stats
from sklearn.metrics import r2_score

from FluidProps import FluidProps
from repo_root import REPO_ROOT
from chebyshev_fit import (MAX_DEGREE, chebyshev_fit_weights, effective_t_uncertainty,
                           fit_chebyshev, select_degree)
from calculate_jt_coefficients import K_TARGET, TempUncertainty

JT_DIR = os.path.join(REPO_ROOT, 'data', 'derived_data', 'jt_coeffs')
PT_DIR = os.path.join(REPO_ROOT, 'data', 'derived_data', 'p_T_pairs')
N_REAL = 1000        # noise realisations per isenthalp
N_MC = 2000          # Monte Carlo samples for U_MC(n), as in the pipeline
SEED = 20210101
COVERAGE_GOAL = 0.95
RULES = ['chi2', 'r2', 'ftest']


def degree_r2(p, T, u_p, u_T):
    """ Earlier rule: lowest degree (1 to 6, capped at N_s - 1) with R^2 > 0.999. """
    w = 1.0 / effective_t_uncertainty(p, T, u_p, u_T)**2
    w = w / w.sum()
    best = (None, -np.inf)
    for n in range(1, min(MAX_DEGREE, len(p) - 1) + 1):
        poly = fit_chebyshev(p, T, n, w)
        r2 = r2_score(T, poly(p))
        if r2 > best[1]:
            best = (poly, r2)
        if r2 > 0.999:
            break
    return best[0]


def degree_ftest(p, T, u_p, u_T, alpha=0.05):
    """ Nested F-test on the weighted residual sum of squares, n <= N_s - 2. """
    u_eff = effective_t_uncertainty(p, T, u_p, u_T)
    w = 1.0 / u_eff
    cap = max(1, min(MAX_DEGREE, len(p) - 2))
    poly = fit_chebyshev(p, T, 1, w)
    rss = np.sum(((T - poly(p)) / u_eff)**2)
    for n in range(2, cap + 1):
        nxt = fit_chebyshev(p, T, n, w)
        rss_nxt = np.sum(((T - nxt(p)) / u_eff)**2)
        dof = len(p) - n - 1
        if (rss - rss_nxt) / (rss_nxt / dof) <= stats.f.ppf(1 - alpha, 1, dof):
            break
        poly, rss = nxt, rss_nxt
    return poly


def fit_rule(rule, p, T, u_p, u_T):
    if rule == 'chi2':
        return select_degree(p, T, u_p, u_T).polynomial
    if rule == 'r2':
        return degree_r2(p, T, u_p, u_T)
    if rule == 'ftest':
        return degree_ftest(p, T, u_p, u_T)
    raise ValueError(rule)


def load_isenthalp(filename, fluid):
    """ Measured points and standard (k = 1) uncertainties, as in calculate_jt_coefficients.py. """
    df = pd.read_csv(os.path.join(PT_DIR, filename))
    p = df['PT102/MPa'].values.astype(float)
    T = df['TT102/K'].values.astype(float)
    u_p = np.array([max(u, K_TARGET * 0.01 / 100 * 13.7) for u in df['PT102/MPa_EXP_UNC']]) / K_TARGET
    t_chain = TempUncertainty('X93303')
    u_T = np.array([max(u, t_chain(Ti)) for u, Ti in zip(df['TT102/K_EXP_UNC'], T)]) / K_TARGET
    fp = FluidProps(tuple(fluid.split('-')) if '-' in fluid else fluid)
    if '-' in fluid:
        x1 = df['x1'].values.astype(float)
        x2 = df['x2'].values.astype(float)
        valid = np.abs(x1 + x2 - 1.0) < 0.02
        if not valid.any():
            valid[:] = True
        fp.set_composition_from_1st_fraction(x1[valid].mean())
    interior = np.ones(len(p), dtype=bool)
    order = np.argsort(p)
    interior[order[0]] = interior[order[-1]] = False
    return p, T, u_p, u_T, df['PT101/MPa'].mean(), df['TT101/K'].mean(), fp, interior


def mc_uncertainty(p, T_eos, u_p, u_T, degree, rng):
    """ k = 2 Monte Carlo uncertainty of mu_JT at a fixed degree. """
    w = chebyshev_fit_weights(p, T_eos, u_p, u_T) if degree < len(p) - 1 else None
    samples = np.empty((N_MC, len(p)))
    for k in range(N_MC):
        poly = fit_chebyshev(p + rng.normal(0, u_p), T_eos + rng.normal(0, u_T), degree, w)
        samples[k] = poly.deriv()(p)
    return K_TARGET * samples.std(axis=0, ddof=1)


def required_u_m(err, u_mc, goal=COVERAGE_GOAL):
    """ Smallest u_m (step 0.0005) with |err| <= 2 sqrt(u_mc^2 + u_m^2) for a fraction goal. """
    for u_m in np.arange(0.0, 0.1, 0.0005):
        if np.mean(np.abs(err) <= K_TARGET * np.sqrt(u_mc**2 + u_m**2)) >= goal:
            return u_m
    return np.nan


SI_GROUPS = [('N2/Ar', r'nitrogen, argon'), ('mixtures', 'mixtures'),
             ('He', 'helium-4'), ('all', 'all')]
SI_FLUIDS = {'Nitrogen': 'nitrogen', 'Argon': 'argon', 'Helium': 'helium-4',
             'Helium-Neon': r'${\rm ^4He-Ne}$', 'Nitrogen-Helium': r'${\rm ^4He-N_2}$'}
SI_ORDER = ['Nitrogen', 'Argon', 'Helium', 'Helium-Neon', 'Nitrogen-Helium']


def write_si_tables(pooled, sens, u_m):
    """ LaTeX tables of SI Section S4: the pooled synthetic results and the
    degree and AAD of the measured isenthalps under each rule. """
    tex = [r"\begin{table}[h]", r"  \centering", r"  \footnotesize",
           r"  \caption{\label{tab:si_synth}Synthetic test of the degree rules: RMS of the relative "
           r"error of $\mu_{\rm JT}$ against the exact EOS slope at the interior points, and the "
           r"fraction of points for which the error lies within $U = 2\sqrt{u_{\rm MC}^2 + "
           rf"(u_{{\rm m}}\,\mu_{{\rm JT}})^2}}$ with $u_{{\rm m}} = {u_m:.3f}$ "
           rf"({N_REAL} realizations per isenthalp).}}",
           r"  \begin{tabular}{lcccccc}", r"    \toprule",
           r"    & \multicolumn{2}{c}{chi-square} & \multicolumn{2}{c}{$R^2 > 0.999$} & "
           r"\multicolumn{2}{c}{F-test} \\",
           r"    \cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
           r"    isenthalps & RMS/\% & within $U$ & RMS/\% & within $U$ & RMS/\% & within $U$ \\",
           r"    \midrule"]
    for key, label in SI_GROUPS:
        cells = []
        for rule in RULES:
            r = pooled.loc[(key, rule)]
            cells += [f"{r['rms_perc']:.1f}", f"{r['coverage_mc_um']:.2f}"]
        tex.append(f"    {label} & " + " & ".join(cells) + r" \\")
    tex += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}", ""]

    tex += [r"\begin{table}[h]", r"  \centering", r"  \footnotesize",
            r"  \caption{\label{tab:si_sens}Selected degree $n$ and AAD of the interior points "
            r"(as defined in the main text) of each measured isenthalp under the three degree rules. "
            r"The chi-square rule is the one used for all results.}",
            r"  \begin{tabular}{lccccrcrcr}", r"    \toprule",
            r"    & & & & \multicolumn{2}{c}{chi-square} & \multicolumn{2}{c}{$R^2 > 0.999$} & "
            r"\multicolumn{2}{c}{F-test} \\",
            r"    \cmidrule(lr){5-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}",
            r"    fluid & $T_{\rm in}$/K & $p_{\rm in}$/MPa & $N_{\rm s}$ & $n$ & AAD/\% & "
            r"$n$ & AAD/\% & $n$ & AAD/\% \\",
            r"    \midrule"]
    sens = sens.assign(order=sens['fluid'].map(SI_ORDER.index))
    for iso, g in sens.sort_values(['order', 'isenthalp']).groupby('isenthalp', sort=False):
        g = g.set_index('rule')
        _, t_in, p_in = iso.rsplit('_', 2)
        cells = []
        for rule in RULES:
            cells += [f"{int(g.loc[rule, 'degree'])}", f"{g.loc[rule, 'AAD_perc']:.1f}"]
        tex.append(f"    {SI_FLUIDS[g['fluid'].iloc[0]]} & {t_in[:-1]} & {p_in[:-3]} & "
                   f"{int(g['n_points'].iloc[0])} & " + " & ".join(cells) + r" \\")
    tex += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}"]
    with open(os.path.join(JT_DIR, 'si_degree_tables.tex'), 'w') as f:
        f.write("\n".join(tex) + "\n")


def main():
    summary = pd.read_csv(os.path.join(JT_DIR, 'jt_coefficients_summary.csv'))
    rng = np.random.default_rng(SEED)
    points, sens_rows = [], []
    for _, meta in summary.iterrows():
        fname, fluid = meta['filename'], meta['fluid']
        p, T, u_p, u_T, p_in, T_in, fp, interior = load_isenthalp(fname, fluid)
        h_in = fp.get_hmass(p=p_in, T=T_in)
        T_eos = np.array([fp.get_T(p=pi, h=h_in) for pi in p])
        mu_true = np.array([fp.get_JT_coefficient(p=pi, T=Ti) for pi, Ti in zip(p, T_eos)])
        u_mc = {n: mc_uncertainty(p, T_eos, u_p, u_T, n, rng) / K_TARGET / np.abs(mu_true)
                for n in range(1, min(MAX_DEGREE, len(p) - 1) + 1)}

        # Measured data: degree and interior deviation from the EOS for each
        # rule, relative to the EOS for pure fluids and to the measurement for
        # mixtures, as in the paper (eq AAD_res)
        mu_eos = np.array([fp.get_JT_coefficient(p=pi, T=Ti) for pi, Ti in zip(p, T)])
        for rule in RULES:
            poly = fit_rule(rule, p, T, u_p, u_T)
            mu = poly.deriv()(p)
            dev = (mu - mu_eos) / (mu if '-' in fluid else mu_eos)
            sens_rows.append({'isenthalp': fname[:-4], 'fluid': fluid, 'rule': rule,
                              'n_points': len(p), 'degree': poly.degree(),
                              'AAD_perc': 100 * np.mean(np.abs(dev[interior]))})

        # Synthetic data on the EOS isenthalp
        for _ in range(N_REAL):
            p_s = p + rng.normal(0, u_p)
            T_s = T_eos + rng.normal(0, u_T)
            for rule in RULES:
                poly = fit_rule(rule, p_s, T_s, u_p, u_T)
                n = poly.degree()
                err = (poly.deriv()(p) - mu_true) / np.abs(mu_true)
                for i in np.where(interior)[0]:
                    points.append((fname[:-4], fluid, rule, n, err[i], u_mc[n][i]))
        print(f'{fname[:-4]} done', flush=True)

    pts = pd.DataFrame(points, columns=['isenthalp', 'fluid', 'rule', 'degree', 'err', 'u_mc'])
    chi2 = pts[pts['rule'] == 'chi2']
    u_m = required_u_m(chi2['err'].values, chi2['u_mc'].values)

    def summarise(g):
        U = K_TARGET * np.sqrt(g['u_mc']**2 + u_m**2)
        counts = g.groupby('degree').size() / len(g)
        return pd.Series({
            'degree_mode': int(counts.idxmax()),
            'degree_mode_frac': counts.max(),
            'degree_min': int(g['degree'].min()),
            'degree_max': int(g['degree'].max()),
            'bias_perc': 100 * g['err'].mean(),
            'sd_perc': 100 * g['err'].std(),
            'rms_perc': 100 * np.sqrt(np.mean(g['err']**2)),
            'coverage_mc': np.mean(np.abs(g['err']) <= K_TARGET * g['u_mc']),
            'coverage_mc_um': np.mean(np.abs(g['err']) <= U),
        })

    degree_cols = ['degree_mode', 'degree_min', 'degree_max']
    per_iso = pts.groupby(['isenthalp', 'rule']).apply(summarise, include_groups=False).reset_index()
    per_iso[degree_cols] = per_iso[degree_cols].astype(int)
    per_iso.insert(2, 'n_points', per_iso['isenthalp'].map(
        {f[:-4]: n for f, n in zip(summary['filename'], summary['n_points'])}))
    per_iso.to_csv(os.path.join(JT_DIR, 'degree_validation.csv'), index=False)
    pd.DataFrame(sens_rows).to_csv(os.path.join(JT_DIR, 'degree_sensitivity.csv'), index=False)

    group = {'Nitrogen': 'N2/Ar', 'Argon': 'N2/Ar', 'Helium': 'He'}
    pts['group'] = pts['fluid'].map(lambda f: group.get(f, 'mixtures'))
    pooled = pts.groupby(['group', 'rule']).apply(summarise, include_groups=False)
    pooled_all = pts.groupby('rule').apply(summarise, include_groups=False)
    pooled_all.index = pd.MultiIndex.from_product([['all'], pooled_all.index], names=['group', 'rule'])
    pooled = pd.concat([pooled, pooled_all])
    pooled[degree_cols] = pooled[degree_cols].astype(int)
    pooled.to_csv(os.path.join(JT_DIR, 'degree_pooled.csv'))
    write_si_tables(pooled, pd.DataFrame(sens_rows), u_m)
    cols = ['rms_perc', 'bias_perc', 'coverage_mc', 'coverage_mc_um']
    print(f'\nSynthetic validation, {N_REAL} realisations per isenthalp, interior points')
    print(pooled[cols].round(3).to_string())
    print(f'\nu_m (k = 1) for {COVERAGE_GOAL:.0%} coverage under the chi2 rule: {u_m:.4f}')


if __name__ == '__main__':
    main()
