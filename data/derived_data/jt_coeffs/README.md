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
| `si_cheb_table.tex` | `extract_fit_coefficients.py` | LaTeX table of the Chebyshev fit coefficients (in the scaled pressure), the pressure range, $\chi^2_\nu$, and the RMS of the fit residuals $s_T$ for the Supporting Information. |
| `heat_exchange_fit.csv` | `estimate_corrections.py` | Per isenthalp, a fit of the outlet-temperature offsets from the EOS isenthalp as an inlet-thermometer offset plus a constant heat exchange (diagnostic, not used in U). |
| `fit_diagnostics.csv` | `diagnose_fits.py` | Degree rule with only the uncorrelated uncertainty, and the residual RMS of the EOS isenthalp at the selected degree (diagnostics for the SI). |
| `hysteresis_check.csv` | `hysteresis_check.py` | Nitrogen 150 K/6 MPa measured with decreasing and increasing outlet pressure: outlet temperatures and μ_JT of both fits. |
| `hene_eos_refit.csv` | `refit_hene_eos.py` | Refit of the He-Ne EOS parameters to the measured $\mu_{\rm JT}$: parameters, AAD, points outside U, and the change of the bubble points and densities. |
| `degree_validation.csv` | `validate_degree_criterion.py` | Synthetic test, one row per isenthalp and degree rule (`chi2`, `r2`, `ftest`): degrees selected, bias, scatter, and RMS of the relative error of $\mu_{\rm JT}$ against the exact EOS slope, and the coverage of the uncertainty. |
| `degree_pooled.csv` | `validate_degree_criterion.py` | The synthetic test pooled per fluid group and rule. |
| `degree_sensitivity.csv` | `validate_degree_criterion.py` | The three degree rules applied to the measured data: degree and AAD per isenthalp. |
| `si_degree_tables.tex` | `validate_degree_criterion.py` | LaTeX tables of the synthetic test and the sensitivity for the Supporting Information. |

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

- `p_UNC/MPa` and `T_UNC/K` combine in quadrature the statistical
  uncertainty of the steady-state window (with a Student factor for the five
  samples), the sensor-chain uncertainty, and, for points that had not settled, a
  drift term (the drift over 60 s before the window × 60 s / √3)
  (`measurement_uncertainty.py`).
- `JT_UNC_MC/(K/MPa)` is the Monte Carlo part: each point is perturbed with a
  Gaussian error at the standard (k = 1) measurement uncertainty, the data are
  refitted with the same weighted Chebyshev polynomial (same degree), and the spread of the
  derivative over 2000 realisations is expanded with k = 2. Where the fit fails
  the χ² test, it is scaled by the Birge ratio √χ²_ν (`birge_ratio` in the
  summary).
- `JT_UNC_X/(K/MPa)` is the composition part for mixtures (k = 2):
  2 |∂μ_JT/∂x1| u(x1), with the sensitivity from the EOS at the measured state
  and u(x1) the larger of 0.0005 and the standard deviation of the point
  compositions of the isenthalp (`u_x1` in the summary); zero for pure fluids.
- `JT_UNC/(K/MPa)` combines in quadrature the Monte Carlo standard
  uncertainty, the composition term, the uncorrected heat-loss and kinetic-energy biases (treated as
  rectangular, u = δ/√3, with δ from `corrections_summary.csv`), and the
  long-term stability of the Cernox sensor (6 × 10⁻⁴ of μ_JT / √3), the
  relative method uncertainty u_m = 0.65 % (k = 1, the fit error from the
  synthetic test of `validate_degree_criterion.py`), and expands the result
  with k = 2.
- `mean_unc_perc` is the mean of 100 `JT_UNC`/|`JT_meas`| over the used
  points of the isenthalp; it is repeated on every row of the `*_JT.csv` file
  and also given in the summary.

The isenthalp end points (lowest and highest pressure, `used` = False) are
fitted but excluded from `mean_unc_perc`, the deviation statistics, and the
paper tables.

The `*_MC_uncertainty.csv` files are a separate, idealised bound: the measured
pressures are placed on the EOS isenthalp through the mean inlet condition and
perturbed within the sensor floor only (and, for mixtures, the composition
uncertainty; 2000 realisations, k = 2).
