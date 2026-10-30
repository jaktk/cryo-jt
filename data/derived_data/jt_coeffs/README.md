# `data/derived_data/jt_coeffs/`

Derived Joule-Thomson coefficients with uncertainties, one isenthalp per
file, plus a summary table over all isenthalps and the same content as JSON
for programmatic use. The directory also holds the Monte Carlo bound on the
achievable uncertainty, the heat-loss and kinetic-energy correction
estimates, and the LaTeX tables generated for the paper.

The full column dictionary is in
[`../../metadata/codebook.csv`](../../metadata/codebook.csv).

## Files

| File | Written by | Content |
| --- | --- | --- |
| `Fluid_TK_PMPa_JT.csv` | `calculate_jt_coefficients.py` | One isenthalp: per-point pressure, temperature, and composition with their expanded uncertainties, the derived $\mu_{\rm JT}$ with its Monte Carlo part and its combined expanded uncertainty (k = 2), the EOS reference at each point, and the `used` flag. |
| `jt_coefficients_summary.csv` | `calculate_jt_coefficients.py` | One row per isenthalp: fluid, composition, fit metadata, mean relative expanded uncertainty, and mean deviation from the EOS. |
| `jt_coefficients_detailed.json` | `calculate_jt_coefficients.py` | The same results in nested JSON form including the per-point arrays. |
| `Fluid_TK_PMPa_MC_uncertainty.csv` | `monte_carlo_uncertainty.py` | One isenthalp: the Monte Carlo bound on the achievable $\mu_{\rm JT}$ uncertainty given the sensor floor (and, for mixtures, the composition uncertainty), evaluated on the EOS isenthalp through the mean inlet condition. |
| `mc_uncertainty_summary.csv` | `monte_carlo_uncertainty.py` | One row per isenthalp: mean and maximum Monte Carlo relative expanded uncertainty (and, for mixtures, the composition contribution). Written only by a full run, not with `--isenthalp`. |
| `corrections_per_point.csv` | `estimate_corrections.py` | Heat-loss (constant 27 mW) and kinetic-energy corrections for every used point, in K/MPa and in % of $\lvert\mu_{\rm JT}\rvert$. |
| `corrections_summary.csv` | `estimate_corrections.py` | Mean corrections per fluid group (He-N2 split by composition), deviation from the EOS with and without the corrections, and the best-fit heat exchange. Source of `CORRECTIONS_PERC` in `calculate_jt_coefficients.py`. |
| `paper_tables.tex` | `generate_paper_tables.py` | LaTeX bodies of the result tables of the paper. |
| `si_cheb_table.tex` | `extract_fit_coefficients.py` | LaTeX table of the Chebyshev fit coefficients and the RMS of the fit residuals $s_T$ for the Supporting Information. |

All scripts are in [`../../../src/data_analysis/`](../../../src/data_analysis).

## Filename convention

`<Fluid>_<T>K_<P>MPa_JT.csv` and `<Fluid>_<T>K_<P>MPa_MC_uncertainty.csv`,
where:

| Token | Meaning |
| --- | --- |
| `Fluid` | `Nitrogen`, `Argon`, `Helium`, `Helium-Neon`, `Nitrogen-Helium` (for mixtures the order matches the `x1`/`x2` columns) |
| `T` | mean inlet temperature, rounded to K |
| `P` | mean inlet pressure, rounded to MPa |

The same stem is used for the matching p-T input in
[`../p_T_pairs/`](../p_T_pairs).

## Uncertainty conventions

All uncertainties in the `*_JT.csv` files are expanded with k = 2 (about
95 % level of confidence).

- `p_UNC/MPa` and `T_UNC/K` are the larger of the statistical uncertainty of
  the steady-state window and the sensor-chain uncertainty.
- `JT_UNC_MC/(K/MPa)` is the Monte Carlo part: each point is perturbed with a
  Gaussian error at the standard (k = 1) measurement uncertainty, the data are
  refitted with the same weighted Chebyshev polynomial, and the spread of the
  derivative over 2000 realisations is expanded with k = 2.
- `JT_UNC/(K/MPa)` combines in quadrature the Monte Carlo standard
  uncertainty, the uncorrected heat-loss and kinetic-energy biases (treated as
  rectangular, u = δ/√3, with δ from `corrections_summary.csv`), and the
  relative method uncertainty u_m = 2.36 % (k = 1, from
  `estimate_method_uncertainty.py`), and expands the result with k = 2.
- `mean_unc_perc` is the mean of 100 `JT_UNC`/|`JT_meas`| over the used
  points of the isenthalp; it is repeated on every row of the `*_JT.csv` file
  and also given in the summary.

The isenthalp end points (lowest and highest pressure, `used` = False) are
fitted but excluded from `mean_unc_perc`, the deviation statistics, and the
paper tables.

The `*_MC_uncertainty.csv` files are a separate, idealised bound: the measured
pressures are placed on the EOS isenthalp through the mean inlet condition and
perturbed within the sensor floor only (2000 realisations, k = 2).
