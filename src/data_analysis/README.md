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
| 2 | `calculate_jt_coefficients.py` | yes | `p_T_pairs/*.csv`, `p_T_pairs/index.json`, Cernox calibration | `jt_coeffs/*_JT.csv`, `jt_coefficients_summary.csv`, `jt_coefficients_detailed.json` |
| 3 | `monte_carlo_uncertainty.py` | yes | `p_T_pairs/*.csv`, `jt_coefficients_summary.csv` | `jt_coeffs/*_MC_uncertainty.csv`, `mc_uncertainty_summary.csv` |
| 4 | `estimate_corrections.py` | yes | `*_JT.csv`, `p_T_pairs/*.csv`, `capillary_mass_flow/*.csv` | `corrections_per_point.csv`, `corrections_summary.csv` |
| 5 | `estimate_method_uncertainty.py` | no | N2 and Ar `*_JT.csv`, `CORRECTIONS_PERC` | prints u_m |
| 6 | `calculate_jt_coefficients.py` again, only if the constants below changed | yes | as step 2 | as step 2 |
| 7 | `generate_paper_tables.py`, `extract_fit_coefficients.py`, figure scripts | no | `jt_coeffs/`, `p_T_pairs/` | `paper_tables.tex`, `si_cheb_table.tex`, `img/*.pdf` |

Paths in the table are relative to `data/derived_data/` unless given in full.
`p_T_pairs/index.json` is maintained by hand: `get_pT_pairs.py` does not
write it, so add or remove its entries when the set of isenthalps changes.

Two constants in `calculate_jt_coefficients.py` are hard-coded from the
outputs of later steps:

- `CORRECTIONS_PERC`: the mean relative heat-loss and kinetic-energy
  corrections (δ_Q, δ_k in %) per fluid group, copied from the `dQ/%` and
  `dk/%` columns of `corrections_summary.csv` (step 4).
- `U_METHOD_REL` (0.0236): the relative method uncertainty u_m printed by
  `estimate_method_uncertainty.py` (step 5). That script uses
  `CORRECTIONS_PERC`, so update `CORRECTIONS_PERC` first.

The two constants only affect `JT_UNC/(K/MPa)` and `mean_unc_perc`.
`JT_meas`, `JT_eos`, and `JT_UNC_MC`, which steps 3 to 5 read, do not depend
on them, so step 6 does not require rerunning steps 3 to 5.

## Scripts

### `get_pT_pairs.py`

Cuts each steady-state window listed in
`data/raw_data/extraction_windows.json` out of the raw time series and
averages it (5 samples at 1 s). For every column `X` it writes the mean, the
standard deviation `X_STD`, and the expanded uncertainty of the mean
`X_EXP_UNC` = 2 `X_STD`/√N (k = 2). Isenthalps with an `"excluded"` field
(13, 15, 23, 24, 25: not used in the published analysis) are skipped. The
output name is `<fluid>_<T_in>K_<p_in>MPa.csv` from the rounded mean inlet
conditions; the script raises an error if two isenthalps map to the same name.

### `calculate_jt_coefficients.py`

Main script for the Joule-Thomson coefficients.

1. **Measurement uncertainties** (expanded, k = 2): `p_UNC` is the larger of
   `PT102/MPa_EXP_UNC` and the sensor floor 2 × 0.01 % of 13.7 MPa; `T_UNC` is
   the larger of `TT102/K_EXP_UNC` and the sensor-chain uncertainty of TT102
   (Cernox sensor, CABTR module, calibration polynomial, in quadrature).
2. **Composition** (mixtures): gas-analyser readings with |x1 + x2 - 1| ≥ 0.02
   are rejected and replaced by the isenthalp mean.
3. **Fit**: weighted least-squares Chebyshev fit of T(p) at the outlet
   (weights `1/(u_T² + (u_p·dT/dp)²)`), choosing the lowest degree (1 to 6,
   capped at N - 1) that reaches R² > 0.999.
4. **μ_JT**: analytic derivative of the fitted polynomial at each point.
5. **Monte Carlo** (`JT_UNC_MC`): each point is perturbed with a Gaussian error
   at the standard (k = 1) measurement uncertainty, refitted with the same
   weighted Chebyshev degree, and re-differentiated over 2000 seeded
   realisations; `JT_UNC_MC` is 2 × the sample standard deviation.
6. **Combined uncertainty** (`JT_UNC`): the Monte Carlo standard uncertainty,
   the uncorrected heat-loss and kinetic-energy biases (rectangular,
   u = δ/√3, from `CORRECTIONS_PERC`; He-N2 split at x_He = 0.3), and the
   method uncertainty u_m = 2.36 % (k = 1) are added in quadrature and
   expanded with k = 2.
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

Heat-loss and kinetic-energy corrections for every used point, for a constant
heat exchange of 27 mW, with the mass flow interpolated from
`data/derived_data/capillary_mass_flow/`. Writes `corrections_per_point.csv`
and `corrections_summary.csv`, and prints the deviation from the EOS with and
without the corrections and the best-fit heat exchange per fluid group. The
corrections are not applied to the data; they enter `JT_UNC` as uncorrected
biases.

### `estimate_method_uncertainty.py`

Estimates the method uncertainty from the N2 and Ar points as
u_m² = RMS(dev)² - mean(u_MC² + u_sys² + u_EOS²), all relative and at k = 1,
with u_EOS taken from the c_p tolerance diagrams of the reference equations of
state. It only prints the result; copy it to `U_METHOD_REL`. It makes no
REFPROP calls but imports `calculate_jt_coefficients.py`, so `ctREFPROP`
must be installed.

### `generate_paper_tables.py`

Writes the LaTeX bodies of the result tables to
`data/derived_data/jt_coeffs/paper_tables.tex` and prints them. The per-point
tables give the standard (k = 1) u(p) and u(T) and the expanded U(μ_JT); μ_JT
and U are left empty for the isenthalp end points (`used` = False).

### `extract_fit_coefficients.py`

Reproduces the Chebyshev fit of every isenthalp without REFPROP, checks that
its derivative matches `JT_meas`, and writes the Supporting Information table
of the coefficients a_k, the number of points and degree, and the RMS of the
fit residuals s_T (mK) to `data/derived_data/jt_coeffs/si_cheb_table.tex`.

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
err(x),err(x)_STD,err(x)_EXP_UNC
```

The `Q008/W` columns are absent where the raw data have no heater power.

**Key variables:**
- `PT102/MPa`, `TT102/K`: downstream pressure and temperature (measurement point)
- `PT101/MPa`, `TT101/K`: upstream pressure and temperature (inlet condition)
- `*_EXP_UNC`: expanded (k = 2) uncertainty of the window mean
- `x1`, `x2`: mixture composition (binary mixtures)

### Output: `jt_coefficients_summary.csv`

`filename`, `fluid`, `is_mixture`, `x1`, `x2`, `TT101_mean/K`, `n_points`,
`polynomial_degree`, `r2_score`, `mean_rel_abs_err_perc`,
`std_rel_abs_err_perc`, `rms_rel_abs_err_perc`, `mean_unc_perc`. The
deviation statistics and `mean_unc_perc` are taken over the used points.

### Output: `<stem>_JT.csv`

| Column | Content |
| --- | --- |
| `p/MPa`, `T/K` | outlet pressure and temperature of the point |
| `p_UNC/MPa`, `T_UNC/K` | expanded (k = 2) measurement uncertainties |
| `x1` | first-component mole fraction at the point (mixtures only) |
| `JT_meas/(K/MPa)` | measured μ_JT |
| `JT_UNC_MC/(K/MPa)` | Monte Carlo part of the uncertainty (k = 2) |
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
- **Temperature (Cernox sensor only):** ±13 mK at 50 K to ±28 mK at 180 K
- **Temperature (full chain, sensor + CABTR + calibration polynomial, combined in quadrature):** ±16 mK at 65 K to ±34 mK at 180 K
- **Pressure:** ±0.01 % of full scale (13.7 MPa)
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

**Python packages:** `numpy`, `pandas`, `scipy`, `scikit-learn` (R²),
`matplotlib`, and `ctREFPROP` (bindings to REFPROP). See
[`requirements.txt`](../../requirements.txt).

**External dependencies:**
- **REFPROP** for the scripts that evaluate the equation of state (steps 2,
  3, 4, and the `theory_figures/` scripts except the schematic).
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
- Every fit reaches R² > 0.999; the selected degrees are 1 to 3.
- The weights account for the measurement uncertainties.

**Mean relative absolute deviation from the EOS (used points):**
- Nitrogen and argon: 0.6 % to 3.7 %
- Helium-4: 16 % to 22 %
- Helium-neon: 1.2 % to 33 %
- Helium-nitrogen: 0.7 % to 11 %

The mean relative expanded uncertainty (`mean_unc_perc`) is about 5 % for
nitrogen, argon, and dilute He-N2 (6.7 % for the He-N2 isenthalp with about
50 % helium), 8 % to 10 % for He-Ne, and 25 % to 26 % for helium-4. See
`jt_coefficients_summary.csv` for the values per isenthalp.
