"""
Generate the LaTeX bodies of the result tables in the manuscript from the
derived JT-coefficient data, so the tables are fully reproducible from the
committed CSVs and stay consistent with the analysis pipeline.

Outputs (printed to stdout, and written to data/derived_data/jt_coeffs/paper_tables.tex):
  * summary table of the pure-fluid isenthalps          (tab:res_pure)
  * summary table of the mixture isenthalps             (tab:res_mix)
  * per-point longtable for the pure fluids             (tab:res_pure_data)
  * per-point longtable for the He-Ne isenthalps        (tab:res_HeNe_data)
  * per-point longtable for the He-N2 isenthalps        (tab:res_HeN2_data)

Run after calculate_jt_coefficients.py.
"""
import os
import numpy as np
import pandas as pd
from get_git_root import get_git_root

JT_DIR = os.path.join(get_git_root(os.getcwd()), "data", "derived_data", "jt_coeffs")
SUMMARY = os.path.join(JT_DIR, "jt_coefficients_summary.csv")

# Display order of the isenthalps in each table, by derived-data filename stem.
PURE_DATA_ORDER = [
    "Nitrogen_150K_6MPa", "Nitrogen_160K_12MPa", "Nitrogen_160K_5MPa",
    "Nitrogen_160K_6MPa", "Nitrogen_160K_7MPa", "Nitrogen_160K_9MPa",
    "Argon_180K_12MPa", "Argon_180K_5MPa",
    "Helium_141K_10MPa", "Helium_65K_5MPa", "Helium_65K_7MPa",
]
PURE_SUMMARY_ORDER = [
    "Nitrogen_150K_6MPa", "Nitrogen_160K_5MPa", "Nitrogen_160K_6MPa",
    "Nitrogen_160K_7MPa", "Nitrogen_160K_9MPa", "Nitrogen_160K_12MPa",
    "Argon_180K_5MPa", "Argon_180K_12MPa",
    "Helium_65K_5MPa", "Helium_65K_7MPa", "Helium_141K_10MPa",
]
HENE_DATA_ORDER = [
    "Helium-Neon_65K_10MPa", "Helium-Neon_65K_5MPa", "Helium-Neon_65K_7MPa",
    "Helium-Neon_65K_8MPa", "Helium-Neon_66K_10MPa", "Helium-Neon_80K_10MPa",
    "Helium-Neon_80K_7MPa", "Helium-Neon_80K_8MPa",
]
HEN2_DATA_ORDER = [
    "Nitrogen-Helium_140K_8MPa", "Nitrogen-Helium_160K_5MPa",
    "Nitrogen-Helium_160K_8MPa",
]
MIX_SUMMARY_ORDER = [
    "Helium-Neon_65K_5MPa", "Helium-Neon_65K_7MPa", "Helium-Neon_65K_8MPa",
    "Helium-Neon_65K_10MPa", "Helium-Neon_66K_10MPa", "Helium-Neon_80K_7MPa",
    "Helium-Neon_80K_8MPa", "Helium-Neon_80K_10MPa",
    "Nitrogen-Helium_140K_8MPa", "Nitrogen-Helium_160K_5MPa",
    "Nitrogen-Helium_160K_8MPa",
]

# Coverage factor of the expanded uncertainties stored in the derived data
# (== calculate_jt_coefficients.K_TARGET). The per-point tables report the
# standard (k = 1) uncertainties u(p) and u(T) and the expanded U(mu_JT).
K_TARGET = 2.0

PURE_DISPLAY = {"Nitrogen": "nitrogen", "Argon": "argon", "Helium": "helium-4"}

# The two dilute He-N2 isenthalps come from the same nominal 14 % helium fill;
# their in-situ means (0.1449 and 0.1454) straddle the two-decimal rounding
# boundary, so the summary table reports the shared nominal composition for a
# consistent label. The precise per-point compositions remain in tab:res_HeN2_data.
NOMINAL_XHE = {"Nitrogen-Helium_160K_8MPa": 0.14}


def sig2(x, sig=2):
    """Format a value to `sig` significant figures (for uncertainties)."""
    if x == 0:
        return "0"
    d = sig - 1 - int(np.floor(np.log10(abs(x))))
    return f"{round(x, d):.{max(d, 0)}f}"


def parse_name(stem):
    """`Nitrogen_150K_6MPa` -> (fluid, Tin, pin)."""
    fluid, t, p = stem.split("_")
    return fluid, int(t[:-1]), int(p[:-3])


def load(stem):
    return pd.read_csv(os.path.join(JT_DIR, f"{stem}_JT.csv"))


def used_stats(stem):
    """Number of points used in the AAD (extremities and inlet excluded) and
    their mean relative expanded uncertainty U(mu_JT) in %."""
    df = load(stem)
    df = df[df["used"] == True]
    u_rel = 100 * df["JT_UNC/(K/MPa)"] / df["JT_meas/(K/MPa)"].abs()
    return len(df), u_rel.mean()


def fmt_rows_pure(df):
    df = df.sort_values("p/MPa")
    out = []
    for _, r in df.iterrows():
        out.append(
            f"  {r['p/MPa']:7.3f} & {sig2(r['p_UNC/MPa']/K_TARGET*1e3):>5} & {r['T/K']:8.3f} & "
            f"{sig2(r['T_UNC/K']/K_TARGET*1e3):>4} & {r['JT_meas/(K/MPa)']:7.3f} & "
            f"{r['JT_UNC/(K/MPa)']:6.3f} & {r['JT_eos/(K/MPa)']:7.3f} \\\\"
        )
    return "\n".join(out)


def fmt_rows_mix(df, he_is_first):
    df = df.sort_values("p/MPa")
    out = []
    for _, r in df.iterrows():
        x_he = r["x1"] if he_is_first else 1.0 - r["x1"]
        out.append(
            f"  {r['p/MPa']:7.3f} & {sig2(r['p_UNC/MPa']/K_TARGET*1e3):>5} & {r['T/K']:8.3f} & "
            f"{sig2(r['T_UNC/K']/K_TARGET*1e3):>4} & {x_he:6.4f} & {r['JT_meas/(K/MPa)']:7.3f} & "
            f"{r['JT_UNC/(K/MPa)']:6.3f} & {r['JT_eos/(K/MPa)']:7.3f} \\\\"
        )
    return "\n".join(out)


def pure_data_table():
    blocks = []
    for stem in PURE_DATA_ORDER:
        fluid, t, p = parse_name(stem)
        hdr = (f"  \\multicolumn{{7}}{{c}}{{\\textit{{{PURE_DISPLAY[fluid]}: "
               f"$T_{{\\rm in}} = {t}$~K, $p_{{\\rm in}} = {p}$~MPa}}}} \\\\")
        blocks.append(hdr + "\n" + fmt_rows_pure(load(stem)))
    return "\n".join(blocks)


def mix_data_table(order, he_is_first):
    summ = pd.read_csv(SUMMARY).set_index("filename")
    blocks = []
    for stem in order:
        fluid, t, p = parse_name(stem)
        x1 = summ.loc[f"{stem}.csv", "x1"]
        xbar_he = x1 if he_is_first else 1.0 - x1
        hdr = (f"  \\multicolumn{{8}}{{c}}{{\\textit{{$T_{{\\rm in}} = {t}$~K, "
               f"$p_{{\\rm in}} = {p}$~MPa, $\\bar{{x}}_{{\\rm He}} \\approx "
               f"{xbar_he:.3f}$}}}} \\\\")
        blocks.append(hdr + "\n" + fmt_rows_mix(load(stem), he_is_first))
    return "\n".join(blocks)


def pure_summary_table():
    summ = pd.read_csv(SUMMARY).set_index("filename")
    rows = []
    for stem in PURE_SUMMARY_ORDER:
        fluid, t, p = parse_name(stem)
        r = summ.loc[f"{stem}.csv"]
        n_used, u_mean = used_stats(stem)
        rows.append(f"    {PURE_DISPLAY[fluid]} & {t} & {p} & {int(r['n_points'])} & "
                    f"{int(r['polynomial_degree'])} & {n_used} & {u_mean:.1f} & "
                    f"{r['mean_rel_abs_err_perc']:.1f} \\\\")
    return "\n".join(rows)


def mix_summary_table():
    summ = pd.read_csv(SUMMARY).set_index("filename")
    rows = []
    for stem in MIX_SUMMARY_ORDER:
        fluid, t, p = parse_name(stem)
        r = summ.loc[f"{stem}.csv"]
        n_used, u_mean = used_stats(stem)
        he_is_first = fluid == "Helium-Neon"
        x1 = r["x1"]
        x_he = x1 if he_is_first else 1.0 - x1
        x_he = NOMINAL_XHE.get(stem, x_he)  # use nominal label where defined
        if he_is_first:
            mix = "${\\rm ^4He-Ne}$"
        else:
            mix = "${\\rm ^4He-N_2}$"
        rows.append(f"    {mix} & {x_he:.2f} & {1-x_he:.2f} & {t} & {p} & "
                    f"{int(r['n_points'])} & {int(r['polynomial_degree'])} & "
                    f"{n_used} & {u_mean:.1f} & {r['mean_rel_abs_err_perc']:.1f} \\\\")
    return "\n".join(rows)


def main():
    parts = {
        "% ==== tab:res_pure (summary, pure fluids) ====": pure_summary_table(),
        "% ==== tab:res_mix (summary, mixtures) ====": mix_summary_table(),
        "% ==== tab:res_pure_data (per-point, pure fluids) ====": pure_data_table(),
        "% ==== tab:res_HeNe_data (per-point, He-Ne) ====": mix_data_table(HENE_DATA_ORDER, True),
        "% ==== tab:res_HeN2_data (per-point, He-N2) ====": mix_data_table(HEN2_DATA_ORDER, False),
    }
    text = "\n\n".join(f"{k}\n{v}" for k, v in parts.items())
    out_path = os.path.join(JT_DIR, "paper_tables.tex")
    with open(out_path, "w") as f:
        f.write(text + "\n")
    print(text)
    print(f"\n% written to {out_path}")


if __name__ == "__main__":
    main()
