"""
Extract the per-isenthalp Chebyshev polynomial fit and report its coefficients
T/K = sum_k a_k C_k(p/MPa) and the RMS of the fit residuals, so that the JCED
Supporting Information can tabulate the correlating-equation parameters
(mu_JT = dT/dp).

The fit reproduces calculate_jt_coefficients._fit_polynomial_weighted exactly
(same weights, same degree-selection rule). No REFPROP needed: only the measured
p-T pairs and the closed-form sensor-chain temperature uncertainty are required.

Verification: the analytic derivative of the recovered polynomial is compared with
the published JT_meas/(K/MPa) column of data/derived_data/jt_coeffs/*_JT.csv.
"""
import os
import numpy as np
import pandas as pd
from numpy.polynomial.chebyshev import Chebyshev
from numpy.polynomial import Polynomial
from sklearn.metrics import r2_score
from repo_root import REPO_ROOT

ROOT = REPO_ROOT
PT_DIR = os.path.join(ROOT, "data", "derived_data", "p_T_pairs")
JT_DIR = os.path.join(ROOT, "data", "derived_data", "jt_coeffs")

# Coverage factors: source specifications are expanded at K_SOURCE, the analysis
# reports everything at K_TARGET (approx. 95 % level of confidence).
K_TARGET = 2.0
K_SOURCE = 2.0

# Sensor-chain temperature uncertainty, X93303 (closed form, == TempUncertainty.__call__)
_POLY_BANDS = [
    {"Tmin": 1.4, "Tmax": 14.1, "N": 31, "n": 9, "DTrms": 0.93},
    {"Tmin": 14.1, "Tmax": 80.0, "N": 31, "n": 7, "DTrms": 1.60},
    {"Tmin": 80.0, "Tmax": 325.0, "N": 32, "n": 8, "DTrms": 5.84},
]


def t_unc_chain(T):
    cernox = -5e-05 * T**2 + 0.1281 * T + 6.4603
    cabtr = 2e-05 * T**2 + 0.0695 * T - 0.1001
    poly = 0.0
    for b in _POLY_BANDS:
        if b["Tmin"] <= T < b["Tmax"]:
            poly = K_SOURCE * np.sqrt(b["N"] / (b["N"] - b["n"]) * b["DTrms"] ** 2)
            break
    # Contributions are independent and expanded at K_SOURCE; combine them in
    # quadrature and rescale the chain to K_TARGET used throughout the analysis.
    return np.sqrt(cernox**2 + cabtr**2 + poly**2) * 1e-3 * (K_TARGET / K_SOURCE)  # K (k=K_TARGET)


def fit_isenthalp(p, T, p_unc, t_unc, max_degree=6):
    """Reproduce _fit_polynomial_weighted."""
    max_degree = min(max_degree, len(p) - 1)
    order = np.argsort(p)
    slope = np.empty_like(T, dtype=float)
    slope[order] = np.gradient(T[order], p[order])
    w = 1.0 / (t_unc**2 + (p_unc * slope) ** 2)
    w = w / np.sum(w)

    best = (None, 1, -np.inf)
    for deg in range(1, max_degree + 1):
        fit = Chebyshev.fit(p, T, deg, w=w)
        poly = Chebyshev(fit.convert().coef)  # raw-p Chebyshev (matches eq S5)
        r2 = r2_score(T, poly(p))
        if r2 > best[2]:
            best = (poly, deg, r2)
        if r2 > 0.999:
            break
    return best


# Process only the published isenthalps (those with a *_JT.csv produced by the
# pipeline).
published = sorted(f.replace("_JT.csv", ".csv") for f in os.listdir(JT_DIR)
                   if f.endswith("_JT.csv"))

rows = []
for fname in published:
    df = pd.read_csv(os.path.join(PT_DIR, fname))
    p = df["PT102/MPa"].values.astype(float)
    T = df["TT102/K"].values.astype(float)
    if len(p) < 3:
        continue
    p_unc = np.array([max(u, K_TARGET * 0.01 / 100 * 13.7) for u in df["PT102/MPa_EXP_UNC"].values])
    t_unc = np.array([max(us, t_unc_chain(Ti)) for us, Ti in zip(df["TT102/K_EXP_UNC"].values, T)])

    poly, deg, r2 = fit_isenthalp(p, T, p_unc, t_unc)

    # a_k : Chebyshev coefficients in raw p, i.e. T/K = sum_k a_k C_k(p/MPa) with
    # C_0=1, C_1=p, C_{k+1}=2 p C_k - C_{k-1}  (eqs S5-S6 of the SI, verbatim).
    a = poly.coef

    # Verification against published JT_meas, at full and at 8-sig-fig rounded coeffs
    jt_csv = os.path.join(JT_DIR, fname.replace(".csv", "_JT.csv"))
    max_dev = max_dev_round = np.nan

    def round_sig(x, s=8):
        out = np.zeros_like(x)
        for i, v in enumerate(x):
            out[i] = 0.0 if v == 0 else round(v, -int(np.floor(np.log10(abs(v)))) + (s - 1))
        return out

    if os.path.exists(jt_csv):
        jt = pd.read_csv(jt_csv)
        jt_pub = jt["JT_meas/(K/MPa)"].values
        pp = jt["p/MPa"].values
        jt_here = poly.deriv()(pp)
        jt_round = Chebyshev(round_sig(a)).deriv()(pp)
        max_dev = np.max(np.abs(jt_here - jt_pub))
        max_dev_round = np.max(np.abs(jt_round - jt_pub))

    # RMS of the temperature residuals over all N_s points; zero when the
    # polynomial interpolates (n = N_s - 1)
    rms_mK = 1e3 * np.sqrt(np.mean((T - poly(p)) ** 2))

    rows.append({"file": fname, "N": len(p), "n": deg, "R2": r2, "rms_mK": rms_mK,
                 "a": a, "maxdev_JT": max_dev, "maxdev_round": max_dev_round})

# Report
print(f"{'isenthalp':28s} N  n   R^2       max|dJT|  max|dJT|8sf  a_k (Chebyshev, raw p/MPa)")
print("-" * 120)
for r in rows:
    astr = "  ".join(f"{ai:+.7e}" for ai in r["a"])
    print(f"{r['file'][:-4]:28s} {r['N']:2d} {r['n']:2d}  {r['R2']:.6f}  {r['maxdev_JT']:.1e}  {r['maxdev_round']:.1e}    {astr}")

print("\nMax verification dev (full coeffs):     "
      f"{np.nanmax([r['maxdev_JT'] for r in rows]):.2e} K/MPa")
print("Max verification dev (8 sig-fig coeffs): "
      f"{np.nanmax([r['maxdev_round'] for r in rows]):.2e} K/MPa")


# ---- emit LaTeX table for the Supporting Information ----
LABEL = {"Argon": "argon", "Helium": "helium-4", "Nitrogen": "nitrogen",
         "Helium-Neon": r"${\rm ^4He-Ne}$", "Nitrogen-Helium": r"${\rm ^4He-N_2}$"}


def parse(fname):
    stem = fname[:-4]
    parts = stem.split("_")
    fluid = "_".join(parts[:-2]).replace("_", "-")
    Tin = parts[-2].replace("K", "")
    pin = parts[-1].replace("MPa", "")
    return LABEL.get(fluid, fluid), Tin, pin


# Worst-case mu_JT reproduction at the printed 7-decimal precision
worst7 = 0.0
for r in rows:
    a7 = np.array([float(f"{ai:.7f}") for ai in r["a"]])
    jt = pd.read_csv(os.path.join(JT_DIR, r["file"].replace(".csv", "_JT.csv")))
    worst7 = max(worst7, float(np.max(np.abs(
        Chebyshev(a7).deriv()(jt["p/MPa"].values) - jt["JT_meas/(K/MPa)"].values))))
print(f"Max verification dev (7-decimal printed table): {worst7:.2e} K/MPa")

tex = []
tex.append(r"\begin{table}[h]")
tex.append(r"  \centering")
tex.append(r"  \footnotesize")
tex.append(r"  \setlength{\tabcolsep}{4pt}")
tex.append(r"  \caption{\label{tab:si_cheb}Coefficients $a_k$ of the Chebyshev fit "
           r"$T/{\rm K} = \sum_{k=0}^{n} a_k\, C_k(p/{\rm MPa})$ (eq~\ref{eq:si_cheb}) "
           r"for each measured isenthalp, with the number of points $N_{\rm s}$ and the "
           r"fitted degree $n$ (all fits reach $R^2 > 0.999$), and the root mean square "
           r"$s_T$ of the temperature residuals (zero where the polynomial interpolates, "
           r"$n = N_{\rm s} - 1$). The Joule-Thomson "
           r"coefficient is the analytic derivative (eq~\ref{eq:si_chebderiv}). "
           r"Empty cells denote coefficients absent at the fitted degree.}")
tex.append(r"  \begin{tabular}{lcccrD{.}{.}{3.7}D{.}{.}{2.7}D{.}{.}{1.7}D{.}{.}{1.7}}")
tex.append(r"    \toprule")
tex.append(r"    fluid & $T_{\rm in}$/K & $p_{\rm in}$/MPa & ($N_{\rm s}$, $n$) & "
           r"$s_T$/mK & \multicolumn{1}{c}{$a_0$} & \multicolumn{1}{c}{$a_1$} & "
           r"\multicolumn{1}{c}{$a_2$} & \multicolumn{1}{c}{$a_3$} \\")
tex.append(r"    \midrule")
order = ["Nitrogen_", "Argon_", "Helium_", "Helium-Neon_", "Nitrogen-Helium_"]
rows_sorted = sorted(rows, key=lambda r: (next((i for i, o in enumerate(order)
                     if r["file"].startswith(o)), 99), r["file"]))
for r in rows_sorted:
    fluid, Tin, pin = parse(r["file"])
    a = list(r["a"]) + [None] * (4 - len(r["a"]))
    cells = ["" if ai is None else f"{ai:.7f}" for ai in a]
    rms = "0" if r["n"] == r["N"] - 1 else f"{r['rms_mK']:.0f}"
    tex.append(f"    {fluid} & {Tin} & {pin} & ({r['N']}, {r['n']}) & {rms} & "
               + " & ".join(cells) + r" \\")
tex.append(r"    \bottomrule")
tex.append(r"  \end{tabular}")
tex.append(r"\end{table}")

out_path = os.path.join(JT_DIR, "si_cheb_table.tex")
with open(out_path, "w") as f:
    f.write("\n".join(tex) + "\n")
print(f"\nLaTeX table written to {out_path}")
