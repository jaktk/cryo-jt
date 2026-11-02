# Data Analysis

This directory contains the Python scripts that turn the raw measurements
into Joule-Thomson coefficients with uncertainties, and that produce the
tables and result figures of the paper.

## Overview

The pipeline extracts steady-state pressure-temperature pairs from the raw
time series, fits each isenthalp with a Chebyshev polynomial, and takes
μ_JT as its derivative. The methodology is based on PhD thesis Chapter 4:
"Isenthalpic Joule-Thomson coefficient measurements."

The scripts locate the repository root from their own file location
([`repo_root.py`](repo_root.py)), so they run from a git clone or from the
unpacked Zenodo archive; git is not needed. Run them from any directory as
`python src/data_analysis/<script>.py`.

## Run order

| Step | Script | REFPROP | Reads | Writes |
| --- | --- | --- | --- | --- |
| 1 | `get_pT_pairs.py` | no | `data/raw_data/*.csv`, `extraction_windows.json` | `p_T_pairs/*.csv` |
| 1b | `analyze_purity.py` | no | `data/metadata/campaigns.csv`, `data/raw_data/*.csv` | `purity/purity_per_day.csv`, `purity/purity_per_gas.csv`, `purity/si_purity_tables.tex` |
| 2 | `calculate_jt_coefficients.py` | yes | `p_T_pairs/*.csv`, `p_T_pairs/index.json`, Cernox calibration, `capillary_mass_flow/*.csv`, `purity/purity_per_gas.csv`, `jt_coeffs/method_uncertainty.csv` | `jt_coeffs/*_JT.csv`, `jt_coefficients_summary.csv`, `jt_coefficients_detailed.json` |
| 3 | `monte_carlo_uncertainty.py` | yes | `p_T_pairs/*.csv`, `jt_coefficients_summary.csv` | `jt_coeffs/*_MC_uncertainty.csv`, `mc_uncertainty_summary.csv` |
| 4 | `estimate_corrections.py` | yes | `*_JT.csv`, `jt_coefficients_summary.csv`, `p_T_pairs/*.csv`, `capillary_mass_flow/*.csv` | `corrections_per_point.csv`, `corrections_summary.csv`, `heat_exchange_fit.csv` |
| 5 | `validate_degree_criterion.py` | yes | `p_T_pairs/*.csv`, `jt_coefficients_summary.csv` | `degree_validation.csv`, `degree_sensitivity.csv`, `degree_pooled.csv`, `si_degree_tables.tex`, `method_uncertainty.csv`; prints u_m |
| 6 | `calculate_jt_coefficients.py` again, only if the constants below or `method_uncertainty.csv` changed | yes | as step 2 | as step 2 |
| 7 | `estimate_method_uncertainty.py` | no | N2 and Ar `*_JT.csv` | prints the cross-check of u_m |
| 7b | `diagnose_fits.py`, `hysteresis_check.py` | yes / no | `p_T_pairs/*.csv`, `jt_coefficients_summary.csv`, `data/raw_data/2021-03-04.csv` | `fit_diagnostics.csv`, `hysteresis_check.csv` |
| 8 | `refit_hene_eos.py` | yes | He-Ne `*_JT.csv`, `jt_coefficients_summary.csv` | `hene_eos_refit.csv` |
| 10 | `generate_paper_tables.py`, `extract_fit_coefficients.py`, figure scripts | no | `jt_coeffs/`, `p_T_pairs/` | `paper_tables.tex`, `si_cheb_table.tex`, `img/*.pdf` |

Paths in the table are relative to `data/derived_data/` unless given in full.
`p_T_pairs/index.json` is maintained by hand: `get_pT_pairs.py` does not
write it, so add or remove its entries when the set of isenthalps changes.

Two constants in `calculate_jt_coefficients.py` are hard-coded from the
outputs of later steps:

- `KINETIC_PERC`: the mean relative kinetic-energy correction δ_k in % per
  fluid group, copied from the `dk/%` column of `corrections_summary.csv`
  (step 4).
- `U_METHOD_REL`: the pooled relative method uncertainty u_m printed by
  `validate_degree_criterion.py` (step 5), from the synthetic test only. The
  per-isenthalp values of `method_uncertainty.csv` (step 5) are read directly.

The two constants and `method_uncertainty.csv` only affect `JT_UNC/(K/MPa)`,
`JT_UNC_Q/(K/MPa)`, and `mean_unc_perc`.
`JT_meas`, `JT_eos`, and `JT_UNC_MC`, which steps 3 to 5 read, do not depend
on them, so step 6 does not require rerunning steps 3 to 5. Step 7 reads
`JT_UNC`, so it runs after step 6.

## Scripts

### `get_pT_pairs.py`

Cuts each steady-state window listed in
`data/raw_data/extraction_windows.json` out of the raw time series and
averages it (5 samples at 1 s). For every column `X` it writes the mean, the
standard deviation `X_STD`, and the expanded uncertainty of the mean
`X_EXP_UNC` = 2 `X_STD`/√N (k = 2), and for the outlet pressure and
temperature the drift `X_DRIFT` (slope per second of a straight line fitted
over the 60 s before the window). The
output name is `<fluid>_<T_in>K_<p_in>MPa.csv` from the rounded mean inlet
conditions; the script raises an error if two isenthalps map to the same name.

### `calculate_jt_coefficients.py`

Main script for the Joule-Thomson coefficients.

1. **Measurement uncertainties** (expanded, k = 2, from
   [`measurement_uncertainty.py`](measurement_uncertainty.py)): `p_UNC` and
   `T_UNC` combine in quadrature the statistical uncertainty of the window
   mean (`*_EXP_UNC`, scaled by t₀.₉₇₅,₄/2 = 1.39 for the 5-sample window),
   a drift term (`*_DRIFT` × 60 s / √3, for points that had not settled),
   and the sensor chain: 0.01 % of 13.7 MPa (k = 2) for pressure; for TT102
   the Cernox calibration (uncertainty table of the X93303 calibration report,
   CX-1050 column), the CABTR module, and the calibration polynomial, in
   quadrature.
2. **Composition** (mixtures): gas-analyser readings with |x1 + x2 - 1| ≥ 0.02
   are rejected and replaced by the isenthalp mean.
3. **Fit**: weighted least-squares Chebyshev fit of T(p) at the outlet that
   minimises χ² = Σ (T - T_fit)²/u_eff², with u_eff² = u_T² + (u_p·dT/dp)² at
   k = 1. The degree is the lowest (1 to 6, at most N - 2) whose reduced
   χ²_ν lies below the 0.95 quantile of χ²_ν/ν (ν = N - n - 1); if none
   passes, the degree with the lowest χ²_ν is used and `chi2_pass` is False. Three-point
   isenthalps are interpolated with n = 2. The rule is in
   [`chebyshev_fit.py`](chebyshev_fit.py), shared by all scripts that fit.
4. **μ_JT**: analytic derivative of the fitted polynomial at each point.
5. **Monte Carlo** (`JT_UNC_MC`): each point is perturbed with a Gaussian error
   at the standard (k = 1) measurement uncertainty, refitted with the same
   weighted Chebyshev degree, and re-differentiated over 2000 seeded
   realisations; `JT_UNC_MC` is 2 × the sample standard deviation, scaled by
   the Birge ratio √χ²_ν (`birge_ratio` in the summary) where the fit fails
   the χ² test.
6. **Combined uncertainty** (`JT_UNC`): the Monte Carlo standard uncertainty,
   for mixtures the composition term `JT_UNC_X` (2 |∂μ_JT/∂x1| u(x1) at k = 2, with the
   sensitivity from the EOS at the measured state; u(x1), `u_x1` in the
   summary, is the larger of 0.0005 and the standard deviation of the point
   compositions of the isenthalp, combined with the slope bias if the
   composition varied along the isenthalp as read, |bias|/√3), the bias term
   `JT_UNC_Q` (the local slope bias from the heat gain modelled in
   [`heat_exchange.py`](heat_exchange.py) and the kinetic-energy correction
   from `KINETIC_PERC`, He-N2 split at x_He = 0.3, each rectangular, u = δ/√3),
   for pure fluids the purity term `JT_UNC_PUR` (change of the EOS μ_JT for
   the in-situ impurity of `purity_per_gas.csv`, contaminant O₂ for N₂ and N₂
   for Ar and He, rectangular), the method uncertainty u_m (k = 1, `u_m` in the
   summary, from `method_uncertainty.csv`), and the long-term stability of the
   Cernox sensor (6 × 10⁻⁴ of μ_JT / √3) are added in quadrature and expanded
   with k = 2. The summary flags `steady` (every point drifts by at most
   10 mK/min and 10 kPa/min before its window) and `well_conditioned`
   (N ≥ 4, χ² test passed, and steady).
7. **Comparison with the EOS**: μ_JT from REFPROP at each measured (p, T) and
   the relative absolute deviation, referred to the EOS value for pure fluids
   and to the measured value for mixtures.

The lowest- and highest-pressure points of each isenthalp are fitted but get
`used` = False and are excluded from `mean_unc_perc`, the deviation
statistics, the figures, and the paper tables.

### `monte_carlo_uncertainty.py`

Monte Carlo bound on the achievable μ_JT uncertainty for each measured
isenthalp. It places the measured pressures on the EOS isenthalp through the
mean inlet condition, perturbs each point within its sensor floor (and, for
mixtures, the composition alone within 0.1 mol-%, k = 2), refits with the
degree from `jt_coefficients_summary.csv`, and reports the spread of the
derived μ_JT (2000 realisations, k = 2).

```bash
python src/data_analysis/monte_carlo_uncertainty.py   # all isenthalps
# one isenthalp:
python src/data_analysis/monte_carlo_uncertainty.py --isenthalp Helium-Neon_65K_5MPa -v
```

It writes `data/derived_data/jt_coeffs/<stem>_MC_uncertainty.csv` for each
isenthalp it runs. `mc_uncertainty_summary.csv` is written only by a full run;
with `--isenthalp` the summary is left unchanged.

### `estimate_corrections.py`

Heat-exchange and kinetic-energy corrections for every used point, with the
mass flow interpolated from `data/derived_data/capillary_mass_flow/`. The heat
gain is modelled per point in `heat_exchange.py`, and its effect on μ_JT is
the local slope bias of the refitted isenthalp. Writes
`corrections_per_point.csv` and `corrections_summary.csv`. The corrections
are not applied to the data; they enter `JT_UNC` as uncorrected biases. It also fits the offsets of the outlet temperatures from the EOS
isenthalp through each point's inlet state as an inlet-thermometer offset plus
a constant heat exchange (`heat_exchange_fit.csv`); this is a diagnostic and
is not used in U.

### `validate_degree_criterion.py`

Synthetic test of the degree rule. For every isenthalp it places the measured
pressures on the EOS isenthalp through the mean inlet state, adds Gaussian
noise at the standard measurement uncertainties (1000 seeded realisations),
applies three degree rules (the χ² rule of `chebyshev_fit.py`, the earlier
R² > 0.999 rule, and a nested F-test), and compares μ_JT at the interior
points with the exact EOS slope. Isenthalps whose EOS curve enters the
two-phase region (nitrogen 150 K/6 MPa) are left out of u_m and of the pooled
results. The method uncertainty u_m is the smallest
value for which 2 √(u_MC² + u_m²) covers the exact slope for 95 % of the
interior points under the χ² rule; copy it to `U_METHOD_REL`. Where the pooled
value covers less than 95 % of the points of one isenthalp, that isenthalp
gets the value that reaches 95 % on its own points (`method_uncertainty.csv`). The measured
deviations from the EOS are not used. The script also applies the three rules
to the measured data (degree and AAD per isenthalp, for the sensitivity
statement in the SI). Writes `degree_validation.csv`,
`degree_sensitivity.csv`, `degree_pooled.csv`, `si_degree_tables.tex`, and
`method_uncertainty.csv`.
Takes a few minutes.

### `diagnose_fits.py`

Diagnostics of the fits quoted in SI Section S4 (`fit_diagnostics.csv`): the
χ² degree rule repeated with only the uncorrelated part of the uncertainty
(window statistics and CABTR), and the residual RMS of the EOS isenthalp
fitted at the selected degree.

### `hysteresis_check.py`

Compares the published nitrogen isenthalp at 150 K/6 MPa (decreasing outlet
pressure) with the return run at increasing pressure recorded the same day
(`hysteresis_check` stamps in `extraction_windows.json`): outlet
temperatures and μ_JT of the two fits (`hysteresis_check.csv`).

### `refit_hene_eos.py`

Refits the binary parameters of the He-Ne model (β_T, γ_T, β_v, γ_v, and the
departure-function weight F) to the measured He-Ne μ_JT, weighted by U/2,
free and with the phase envelope of the published model held (bubble points
at (30 to 42) K as pseudo-data), each to all He-Ne points and to the
well-conditioned isenthalps. The parameters are changed in memory through
REFPROP's `SETKTVdll`; no REFPROP file is modified. Writes
`hene_eos_refit.csv` with the AAD, the points outside U, and the change of the
bubble points and densities for each fit.

### `analyze_purity.py`

In-situ purity of the pure fluids from the gas analyser: drops null frames
(impurity 0) and frozen frames (runs of 300 or more identical readings), and
writes the per-day and per-gas purities and the SI tables to
`data/derived_data/purity/`.

### `estimate_method_uncertainty.py`

Cross-check of u_m against the N2 and Ar points: the excess
u_ex² = RMS(dev)² - mean(u_MC² + u_sys² + u_EOS²), all relative and at k = 1,
with u_EOS taken from the c_p tolerance diagrams of the reference equations of
state, a bootstrap of u_ex over the isenthalps, and the number of points
within `JT_UNC`. It only prints. It makes no REFPROP calls but imports
`calculate_jt_coefficients.py`, so `ctREFPROP` must be installed.

### `generate_paper_tables.py`

Writes the LaTeX bodies of the result tables to
`data/derived_data/jt_coeffs/paper_tables.tex` and prints them. The per-point
tables give the standard (k = 1) u(p) and u(T) and the expanded U(μ_JT); μ_JT
and U are left empty for the isenthalp end points (`used` = False).

### `extract_fit_coefficients.py`

Reproduces the Chebyshev fit of every isenthalp without REFPROP, checks that
the derivative of the printed coefficients matches `JT_meas`, and writes the
Supporting Information table to
`data/derived_data/jt_coeffs/si_cheb_table.tex`: the coefficients a_k in the
scaled pressure x = (2p - p_max - p_min)/(p_max - p_min), the pressure range,
the number of points and degree, χ²_ν, and the RMS of the fit residuals s_T
(mK).

### `theoretical_jt_uncertainty.py`

Closed-form (GUM) and Monte Carlo propagation of the sensor-limited
uncertainty of μ_JT, parametrised by isenthalp slope, curvature, and number
of points. Writes `img/theoretical_uncertainty_combined.pdf`, the lower bound
on the achievable uncertainty shown in the paper.

### `compare_theory_vs_experiment.py`

Overlays the measured μ_JT deviations from the EOS (used points only) on the
sensor-floor band from `theoretical_jt_uncertainty.py`. Writes
`img/theory_vs_experiment_uncertainty.pdf`.

### `plot_paper_results.py` and `plot_toc_graphic.py`

`plot_paper_results.py` writes the per-fluid result figures
(`measurements_pures_errorbars.pdf`, `measurements_HeNe_errorbars.pdf`,
`measurements_HeN2_errorbars.pdf`) with `JT_UNC` error bars for the used
points. `plot_toc_graphic.py` writes the table-of-contents graphic
`img/toc.pdf` at 3.25 in × 1.75 in.

### Supporting files

- `FluidProps.py`: REFPROP wrapper for thermodynamic property calculations.
  It finds REFPROP through the `RPPREFIX` environment variable or a default
  per-platform location.
- `measurement_uncertainty.py`: the per-point pressure and temperature
  uncertainties (sensor chain and window statistics). No REFPROP needed.
- `chebyshev_fit.py`: the weighted Chebyshev fit, its weights, and the χ²
  degree rule, shared by the scripts above. No REFPROP needed.
- `heat_exchange.py`: the heat gained by the fluid at each point, radiation
  from the shield (TT009) onto the capillary, treated as bare stainless steel
  (ε = 0.1, no credit for the multi-layer insulation), plus conduction down
  the outlet tube from TT010 (4 mm × 3 mm, 1 m, as in the design
  calculation), and the resulting local slope bias of μ_JT. Used by
  `calculate_jt_coefficients.py` and `estimate_corrections.py`.
- `repo_root.py`: defines `REPO_ROOT`, the repository root computed from the
  location of this file.

## Data structure

Full column definitions for every file are in
[`data/metadata/codebook.csv`](../../data/metadata/codebook.csv).

### Input: `p_T_pairs` CSV files

```
PT102/MPa,PT102/MPa_STD,PT102/MPa_EXP_UNC,PT101/MPa,PT101/MPa_STD,PT101/MPa_EXP_UNC,
TT009/K,TT009/K_STD,TT009/K_EXP_UNC,TT010/K,TT010/K_STD,TT010/K_EXP_UNC,
TT008/K,TT008/K_STD,TT008/K_EXP_UNC,TT101/K,TT101/K_STD,TT101/K_EXP_UNC,
RT101/ohm,RT101/ohm_STD,RT101/ohm_EXP_UNC,TT102/K,TT102/K_STD,TT102/K_EXP_UNC,
RT102/ohm,RT102/ohm_STD,RT102/ohm_EXP_UNC,Q008/W,Q008/W_STD,Q008/W_EXP_UNC,
x1,x1_STD,x1_EXP_UNC,x2,x2_STD,x2_EXP_UNC,impurity,impurity_STD,impurity_EXP_UNC,
err(x),err(x)_STD,err(x)_EXP_UNC,PT102/MPa_DRIFT,TT102/K_DRIFT
```

The `Q008/W` columns are absent where the raw data have no heater power.

**Key variables:**
- `PT102/MPa`, `TT102/K`: downstream pressure and temperature (measurement point)
- `PT101/MPa`, `TT101/K`: upstream pressure and temperature (inlet condition)
- `*_EXP_UNC`: expanded (k = 2) uncertainty of the window mean
- `PT102/MPa_DRIFT`, `TT102/K_DRIFT`: drift (per second) over the 60 s before the window
- `x1`, `x2`: mixture composition (binary mixtures)

### Output: `jt_coefficients_summary.csv`

`filename`, `fluid`, `is_mixture`, `x1`, `x2`, `TT101_mean/K`, `n_points`,
`polynomial_degree`, `chi2_red`, `chi2_crit`, `chi2_pass`, `s_T/mK`, `steady`, `well_conditioned`, `birge_ratio`, `u_x1`,
`mean_rel_abs_err_perc`,
`std_rel_abs_err_perc`, `rms_rel_abs_err_perc`, `mean_unc_perc`. The
deviation statistics and `mean_unc_perc` are taken over the used points.

### Output: `<stem>_JT.csv`

| Column | Content |
| --- | --- |
| `p/MPa`, `T/K` | outlet pressure and temperature of the point |
| `p_UNC/MPa`, `T_UNC/K` | expanded (k = 2) measurement uncertainties |
| `x1` | first-component mole fraction at the point (mixtures only) |
| `JT_meas/(K/MPa)` | measured μ_JT |
| `JT_UNC_MC/(K/MPa)` | Monte Carlo part of the uncertainty (k = 2, with the Birge ratio) |
| `JT_UNC_X/(K/MPa)` | composition part of the uncertainty (k = 2, mixtures) |
| `JT_UNC/(K/MPa)` | combined expanded uncertainty (k = 2) |
| `JT_eos/(K/MPa)` | μ_JT from the EOS at the same point |
| `mean_unc_perc` | mean of 100 `JT_UNC`/\|`JT_meas`\| over the used points (same on every row) |
| `rel_abs_err_perc` | relative absolute deviation from the EOS |
| `used` | False for the isenthalp end points |

### Output: `<stem>_MC_uncertainty.csv`

`p/MPa`, `T_eos/K`, `jt_eos/(K/MPa)`, `jt_mc_mean/(K/MPa)`,
`jt_mc_std/(K/MPa)` (k = 1), `jt_mc_unc/(K/MPa)` (k = 2),
`jt_mc_rel_unc_perc`, and for mixtures `jt_mc_cmp_rel_unc_perc` (composition
contribution).

## Experimental conditions

**Fluids measured:**
- **Pure fluids:** nitrogen, argon, helium-4
- **Binary mixtures:** helium-neon, helium-nitrogen

**Inlet temperature range:** 65 K to 180 K
**Pressure range:** 0.1 MPa to 12 MPa

**Measurement uncertainty** (k = 2, computed by `theoretical_jt_uncertainty.combined_temperature_uncertainty`, using the outlet thermometer TT102 / sensor X93303):
- **Temperature (Cernox sensor only, calibration report of X93303):** ±12 mK at 50 K to ±29 mK at 180 K
- **Temperature (full chain, sensor + CABTR + calibration polynomial, combined in quadrature):** ±15 mK at 65 K to ±34 mK at 180 K
- **Pressure:** ±0.01 % of full scale (13.7 MPa), k = 2
- **Composition:** ±0.1 mol-% (binary mixtures)

## Theoretical background

The Joule-Thomson coefficient is defined as:
```
μ_JT = (∂T/∂p)_h,x̄
```

Where:
- T = temperature
- p = pressure
- h = specific enthalpy (constant)
- x̄ = composition vector (constant)

**Physical significance:**
- Measures temperature change during isenthalpic expansion
- Positive μ_JT: cooling upon expansion
- Negative μ_JT: heating upon expansion (He, Ne at ambient conditions)
- Zero μ_JT: inversion temperature

**Measurement principle:**
1. Fluid undergoes isenthalpic expansion through capillary
2. Pressure and temperature measured upstream and downstream
3. Multiple points measured along single isenthalpic line
4. Polynomial fitted to p-T data
5. μ_JT calculated as derivative dT/dp

## Requirements

**Python packages:** `numpy`, `pandas`, `scipy`, `scikit-learn` (R² of the
earlier degree rule, in `validate_degree_criterion.py`),
`matplotlib`, and `ctREFPROP` (bindings to REFPROP). See
[`requirements.txt`](../../requirements.txt).

**External dependencies:**
- **REFPROP** for the scripts that evaluate the equation of state (steps 2,
  3, 4, 5, and the `theory_figures/` scripts except the schematic).
- The custom mixture models for He-Ne, He-Ar, Ne-Ar, and He-N2 in
  [`REFPROP/FLUIDS/HMX.BNC`](../../REFPROP/FLUIDS/HMX.BNC), merged into the
  local REFPROP installation (see [`REFPROP/README.md`](../../REFPROP/README.md)).
- **LaTeX** (`latex` + `dvipng`, with the packages `amsmath`, `helvet`, and
  `newtxsf`) for the figure scripts: `matplotlib`'s `text.usetex` renders
  text in Helvetica with sans-serif math, as set in the shared style sheet.
  Set `text.usetex: False` in `jced.mplstyle` to skip LaTeX (the fonts will
  then differ from the published figures).

## Figure regeneration

All paper figures are generated by the scripts below and written to the
repository [`img/`](../../img) directory. The scripts share a single
matplotlib style sheet, [`jced.mplstyle`](../jced.mplstyle), loaded with
`plt.style.use(...)`; each script then sets its own figure-specific layout.
Figures are drawn at their printed size.

| Script | Output(s) |
| --- | --- |
| [`plot_paper_results.py`](plot_paper_results.py) | `img/measurements_pures_errorbars.pdf`, `img/measurements_HeNe_errorbars.pdf`, `img/measurements_HeN2_errorbars.pdf` |
| [`theoretical_jt_uncertainty.py`](theoretical_jt_uncertainty.py) | `img/theoretical_uncertainty_combined.pdf` |
| [`compare_theory_vs_experiment.py`](compare_theory_vs_experiment.py) | `img/theory_vs_experiment_uncertainty.pdf` |
| [`plot_toc_graphic.py`](plot_toc_graphic.py) | `img/toc.pdf` |
| [`generate_paper_tables.py`](generate_paper_tables.py) | `data/derived_data/jt_coeffs/paper_tables.tex` |
| [`extract_fit_coefficients.py`](extract_fit_coefficients.py) | `data/derived_data/jt_coeffs/si_cheb_table.tex` |

These scripts read the committed derived data and do not require REFPROP.
The theory figures are listed in
[`../theory_figures/README.md`](../theory_figures/README.md) and in
[`img/README.md`](../../img/README.md).

```bash
python src/data_analysis/plot_paper_results.py
python src/data_analysis/compare_theory_vs_experiment.py
python src/data_analysis/theoretical_jt_uncertainty.py
python src/data_analysis/plot_toc_graphic.py
python src/data_analysis/generate_paper_tables.py
python src/data_analysis/extract_fit_coefficients.py
```

## Error handling

- Files listed in `index.json` but missing on disk are skipped with a warning.
- Isenthalps with fewer than 3 points are skipped.
- Failed Monte Carlo refits are dropped (and counted in
  `monte_carlo_uncertainty.py`); in the composition Monte Carlo, realisations
  whose EOS μ_JT is non-finite or more than 5 times off the nominal value are
  also dropped.
- EOS failures at a point give NaN for `JT_eos`.
- `get_pT_pairs.py` stops with an error if two isenthalps map to the same
  output file name.

## Data quality

**Fitting:**
- The selected degrees are 2 to 4. Eleven fits pass the χ² test, 9 do not
  (they take the degree with the lowest χ²_ν), and the two three-point
  isenthalps are interpolated. Six isenthalps are well-conditioned.
- The weights account for the measurement uncertainties.

**Mean relative absolute deviation from the EOS (used points):**
- Nitrogen and argon: 0.3 % to 2.2 %
- Helium-4: 16 % to 53 %
- Helium-neon: 3.1 % to 19 %
- Helium-nitrogen: 0.7 % to 11 %

The mean relative expanded uncertainty (`mean_unc_perc`) is 1.6 % to 2.5 %
for nitrogen and argon, 2.0 % for dilute He-N2 (5.0 % for the He-N2
isenthalp with about 50 % helium), 6.5 % to 9.6 % for He-Ne, and 32 % to 68 %
for helium-4; the larger values belong to fits that fail the χ² test, whose
Monte Carlo part is scaled by the Birge ratio. See
`jt_coefficients_summary.csv` for the values per isenthalp.
