"""
Extract the per-isenthalp Chebyshev polynomial fit and report its coefficients,
the reduced chi-square, and the RMS of the fit residuals, so that the JCED
Supporting Information can tabulate the correlating-equation parameters
(mu_JT = dT/dp).

The coefficients refer to the scaled pressure
    x = (2 p - p_max - p_min) / (p_max - p_min),
which maps the measured pressure range of the isenthalp (rounded to 1 kPa)
onto [-1, 1], so that
    T/K = sum_k a_k C_k(x),   mu_JT = 2 / (p_max - p_min) * sum_k a_k C_k'(x)
with the Chebyshev polynomials C_k of the first kind (SI Section S4). The
scaled basis keeps the coefficients well conditioned at degree 5.

The fit is the one of calculate_jt_coefficients.py (chebyshev_fit.select_degree,
same uncertainties from measurement_uncertainty.py, same weights). No REFPROP
needed: only the measured p-T pairs are required.

Verification: the derivative of the polynomial rebuilt from the printed
coefficients is compared with the published JT_meas/(K/MPa) column of
data/derived_data/jt_coeffs/*_JT.csv.
"""
import os
import numpy as np
import pandas as pd
from numpy.polynomial.chebyshev import Chebyshev
from repo_root import REPO_ROOT
from chebyshev_fit import select_degree
from measurement_uncertainty import K_TARGET, point_uncertainties

ROOT = REPO_ROOT
PT_DIR = os.path.join(ROOT, "data", "derived_data", "p_T_pairs")
JT_DIR = os.path.join(ROOT, "data", "derived_data", "jt_coeffs")

DECIMALS = 5      # printed decimals of a_k in K
MAX_ORDER = 5     # highest coefficient column in the table


def scaled_series(coef, p_min, p_max):
    """ Chebyshev series in p/MPa from the coefficients in the scaled pressure x. """
    return Chebyshev(coef, domain=[p_min, p_max])


def main():
    # Process only the published isenthalps (those with a *_JT.csv produced by
    # the pipeline).
    published = sorted(f.replace("_JT.csv", ".csv") for f in os.listdir(JT_DIR)
                       if f.endswith("_JT.csv"))

    rows = []
    for fname in published:
        df = pd.read_csv(os.path.join(PT_DIR, fname))
        p = df["PT102/MPa"].values.astype(float)
        T = df["TT102/K"].values.astype(float)
        if len(p) < 3:
            continue
        p_unc, t_unc = point_uncertainties(df)

        fit = select_degree(p, T, p_unc / K_TARGET, t_unc / K_TARGET)
        # The scaling uses the printed (rounded) pressure range, so that the
        # table reproduces the fit exactly
        p_min, p_max = round(p.min(), 3), round(p.max(), 3)
        a = fit.polynomial.convert(domain=[p_min, p_max]).coef
        a_printed = np.round(a, DECIMALS)

        # Verification against the published JT_meas, at full precision and
        # with the printed coefficients
        jt = pd.read_csv(os.path.join(JT_DIR, fname.replace(".csv", "_JT.csv")))
        jt_pub = jt["JT_meas/(K/MPa)"].values
        pp = jt["p/MPa"].values
        max_dev = np.max(np.abs(scaled_series(a, p_min, p_max).deriv()(pp) - jt_pub))
        max_dev_printed = np.max(np.abs(scaled_series(a_printed, p_min, p_max).deriv()(pp) - jt_pub))

        rows.append({"file": fname, "N": len(p), "n": fit.degree, "p_min": p_min, "p_max": p_max,
                     "chi2_red": fit.chi2_red, "passed": fit.passed, "rms_mK": 1e3 * fit.s_T,
                     "a": a_printed, "maxdev_JT": max_dev, "maxdev_printed": max_dev_printed})

    # Report
    print(f"{'isenthalp':28s} N  n  chi2_red  s_T/mK  max|dJT|  max|dJT|printed  a_k (scaled p)")
    print("-" * 120)
    for r in rows:
        astr = "  ".join(f"{ai:+.{DECIMALS}f}" for ai in r["a"])
        print(f"{r['file'][:-4]:28s} {r['N']:2d} {r['n']:2d}  {r['chi2_red']:8.2f}  {r['rms_mK']:6.1f}  "
              f"{r['maxdev_JT']:.1e}   {r['maxdev_printed']:.1e}      {astr}")
    print("\nMax verification dev (full coeffs):    "
          f"{np.nanmax([r['maxdev_JT'] for r in rows]):.2e} K/MPa")
    print(f"Max verification dev (printed table, {DECIMALS} decimals): "
          f"{np.nanmax([r['maxdev_printed'] for r in rows]):.2e} K/MPa")

    # ---- emit LaTeX table for the Supporting Information ----
    label = {"Argon": "argon", "Helium": "helium-4", "Nitrogen": "nitrogen",
             "Helium-Neon": r"${\rm ^4He-Ne}$", "Nitrogen-Helium": r"${\rm ^4He-N_2}$"}

    def parse(fname):
        parts = fname[:-4].split("_")
        fluid = "_".join(parts[:-2]).replace("_", "-")
        return label.get(fluid, fluid), parts[-2].replace("K", ""), parts[-1].replace("MPa", "")

    coef_cols = " & ".join(rf"\multicolumn{{1}}{{c}}{{$a_{k}$}}" for k in range(MAX_ORDER + 1))
    tex = [
        r"\begin{sidewaystable}",
        r"  \centering",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{4pt}",
        r"  \caption{\label{tab:si_cheb}Chebyshev fit $T/{\rm K} = \sum_{k=0}^{n} a_k\, C_k(x)$ "
        r"(eq~\ref{eq:si_cheb}) of each measured isenthalp in the scaled pressure "
        r"$x = (2p - p_{\max} - p_{\min})/(p_{\max} - p_{\min})$: number of points $N_{\rm s}$, "
        r"selected degree $n$, pressure range, reduced chi-square $\chi^2_\nu$ with "
        r"$\nu = N_{\rm s} - n - 1$, root mean square $s_T$ of the temperature residuals, and "
        r"coefficients $a_k$. $\chi^2_\nu$ values above the 0.95 quantile of $\chi^2_\nu/\nu$ are "
        r"marked with an asterisk; for these isenthalps no degree up to the cap $n = N_{\rm s} - 2$ "
        r"passes the test. The three-point isenthalps are interpolated ($n = 2$, $\nu = 0$). The "
        r"Joule-Thomson coefficient is the analytic derivative (eq~\ref{eq:si_chebderiv}). Empty "
        r"cells denote coefficients absent at the selected degree.}",
        r"  \begin{tabular}{lccccccrr" + "D{.}{.}{3.5}" + "D{.}{.}{2.5}" * MAX_ORDER + "}",
        r"    \toprule",
        r"    fluid & $T_{\rm in}$/K & $p_{\rm in}$/MPa & $N_{\rm s}$ & $n$ & $p_{\min}$/MPa & "
        r"$p_{\max}$/MPa & \multicolumn{1}{c}{$\chi^2_\nu$} & \multicolumn{1}{c}{$s_T$/mK} & "
        + coef_cols + r" \\",
        r"    \midrule",
    ]
    order = ["Nitrogen_", "Argon_", "Helium_", "Helium-Neon_", "Nitrogen-Helium_"]
    rows_sorted = sorted(rows, key=lambda r: (next((i for i, o in enumerate(order)
                         if r["file"].startswith(o)), 99), r["file"]))
    for r in rows_sorted:
        fluid, t_in, p_in = parse(r["file"])
        a = list(r["a"]) + [None] * (MAX_ORDER + 1 - len(r["a"]))
        cells = ["" if ai is None else f"{ai:.{DECIMALS}f}" for ai in a]
        if r["passed"] is None:
            chi2 = "--"
        else:
            chi2 = f"{r['chi2_red']:.1f}" + ("" if r["passed"] else r"$^{*}$")
        tex.append(f"    {fluid} & {t_in} & {p_in} & {r['N']} & {r['n']} & {r['p_min']:.3f} & "
                   f"{r['p_max']:.3f} & {chi2} & {r['rms_mK']:.0f} & " + " & ".join(cells) + r" \\")
    tex += [r"    \bottomrule", r"  \end{tabular}", r"\end{sidewaystable}"]

    out_path = os.path.join(JT_DIR, "si_cheb_table.tex")
    with open(out_path, "w") as f:
        f.write("\n".join(tex) + "\n")
    print(f"\nLaTeX table written to {out_path}")


if __name__ == "__main__":
    main()
