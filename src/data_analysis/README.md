# Data Analysis

This directory contains Python scripts for analyzing Joule-Thomson coefficient measurements from cryogenic fluid experiments.

## Overview

The data analysis pipeline processes experimental pressure-temperature pairs from isenthalpic expansions to calculate Joule-Thomson coefficients with full uncertainty propagation. The methodology is based on PhD thesis Chapter 4: "Isenthalpic Joule-Thomson coefficient measurements."

## Key Scripts

### `calculate_jt_coefficients.py`

Main script for calculating Joule-Thomson coefficients from processed experimental data.

**Features:**
- **Full uncertainty propagation** through polynomial fitting and differentiation
- **Weighted least squares fitting** using measurement uncertainties
- **Monte Carlo uncertainty** of the derived coefficient (Gaussian perturbations, expanded to k = 1.96)
- **Automatic polynomial degree selection**: lowest degree reaching R² > 0.999
- **Comprehensive output** with individual and summary results

**Usage:**
```bash
python src/data_analysis/calculate_jt_coefficients.py
```

**Input Data:**
- `data/derived_data/p_T_pairs/*.csv` - Processed pressure-temperature pairs with uncertainties
- `data/derived_data/p_T_pairs/index.json` - Manifest of per-isenthalp CSVs (filename, fluid, mean inlet temperature)

**Output:**
- `data/derived_data/jt_coeffs/jt_coefficients_summary.csv` - Summary results table
- `data/derived_data/jt_coeffs/jt_coefficients_detailed.json` - Complete results with arrays
- `data/derived_data/jt_coeffs/*_JT.csv` - Individual measurement files

**Methodology:**

1. **Polynomial Fitting:**
   - Fits first-kind Chebyshev polynomials to isenthalpic p-T data
   - Uses weighted least squares with measurement uncertainties
     (weights `1/(u_T² + (u_p·dT/dp)²)`; weights do not affect an exact interpolation)
   - Selects the lowest degree (1-6, capped at N_points-1) reaching R² > 0.999

2. **JT Coefficient Calculation:**
   - Computes μ_JT = dT/dp as the analytic derivative of the fitted polynomial

3. **Uncertainty Analysis (Monte Carlo):**
   - Perturbs each measured point with a **Gaussian** error whose standard
     deviation is the **standard (k = 1)** uncertainty (expanded value / 1.96)
   - Refits with the same weighted Chebyshev fit and re-differentiates, over
     **2000** realizations (seeded for reproducibility)
   - Reports the **expanded** uncertainty: k = 1.96 × the sample standard deviation
   - The closed-form GUM propagation and sensor-floor bound are computed
     separately by `theoretical_jt_uncertainty.py`

4. **Validation:**
   - Compares measured vs. theoretical JT coefficients using REFPROP
   - Calculates relative errors and goodness-of-fit metrics
   - Provides detailed quality assessment

### `theoretical_jt_uncertainty.py`

Closed-form (GUM) and Monte Carlo propagation of the sensor-limited uncertainty of $\mu_{\rm JT}$, parametrised by isenthalp slope, curvature, and number of points. Writes `img/theoretical_uncertainty_combined.{pdf,svg}` and is used in the paper as the lower bound on achievable uncertainty.

### `compare_theory_vs_experiment.py`

Overlays the measured $\mu_{\rm JT}$ deviations from the EOS on the sensor-floor band from `theoretical_jt_uncertainty.py`. Writes `img/theory_vs_experiment_uncertainty.pdf`.

### `monte_carlo_uncertainty.py`

Monte Carlo bound on the achievable $\mu_{\rm JT}$ uncertainty *for each actually-measured isenthalp*. Reads `data/derived_data/p_T_pairs/*.csv` and `data/derived_data/jt_coeffs/jt_coefficients_summary.csv` (for the polynomial degree), places the measured pressures on the EOS-derived isenthalp through the mean inlet condition, perturbs each point within its sensor-floor (and, for mixtures, composition) uncertainty, refits Chebyshev, and reports the per-point spread of the derived $\mu_{\rm JT}$. No hard-coded $p$, $T$, $x$ values.

```bash
python src/data_analysis/monte_carlo_uncertainty.py             # all isenthalps
python src/data_analysis/monte_carlo_uncertainty.py \           # one isenthalp
    --isenthalp Helium-Neon_65K_5MPa --n-samples 2000 -v
```

Writes `data/derived_data/jt_coeffs/<stem>_MC_uncertainty.csv` per isenthalp and a combined `mc_uncertainty_summary.csv`. Requires REFPROP for the EOS lookups.

### Supporting Files

- `FluidProps.py` - REFPROP wrapper for thermodynamic property calculations
- `get_git_root.py` - Utility for finding the repository root directory

## Data Structure

### Input Data Format

**p_T_pairs CSV files:**
```
PT102/MPa,PT102/MPa_STD,PT102/MPa_EXP_UNC,PT101/MPa,PT101/MPa_STD,PT101/MPa_EXP_UNC,
TT009/K,TT009/K_STD,TT009/K_EXP_UNC,TT010/K,TT010/K_STD,TT010/K_EXP_UNC,
TT008/K,TT008/K_STD,TT008/K_EXP_UNC,TT101/K,TT101/K_STD,TT101/K_EXP_UNC,
RT101/ohm,RT101/ohm_STD,RT101/ohm_EXP_UNC,TT102/K,TT102/K_STD,TT102/K_EXP_UNC,
RT102/ohm,RT102/ohm_STD,RT102/ohm_EXP_UNC,Q008/W,Q008/W_STD,Q008/W_EXP_UNC,
x1,x1_STD,x1_EXP_UNC,x2,x2_STD,x2_EXP_UNC,impurity,impurity_STD,impurity_EXP_UNC,
err(x),err(x)_STD,err(x)_EXP_UNC
```

**Key variables:**
- `PT102/MPa` - Downstream pressure (measurement point)
- `TT102/K` - Downstream temperature (measurement point)  
- `PT101/MPa` - Upstream pressure (reference condition)
- `TT101/K` - Upstream temperature (reference condition)
- `*_EXP_UNC` - Expanded measurement uncertainties
- `x1`, `x2` - Mixture composition (for binary mixtures)

### Output Data Format

**Summary CSV columns:**
- `filename` - Original data file name
- `fluid` - Fluid type (e.g., "Nitrogen", "Helium-Neon")
- `is_mixture` - Boolean indicating mixture vs. pure fluid
- `composition` - Mole fraction of first component (mixtures only)
- `n_points` - Number of measurement points
- `polynomial_degree` - Optimal polynomial degree used
- `r2_score` - Goodness of fit (R² value)
- `mean_relative_error_percent` - Mean deviation from EOS
- `rms_relative_error_percent` - RMS deviation from EOS
- `mean_uncertainty_percent` - Mean measurement uncertainty

**Individual measurement files:**
- `pressure_MPa` - Measurement pressures
- `temperature_K` - Measurement temperatures
- `pressure_uncertainty_MPa` - Pressure measurement uncertainties
- `temperature_uncertainty_K` - Temperature measurement uncertainties
- `jt_measured_K_per_MPa` - Measured JT coefficients
- `jt_theoretical_K_per_MPa` - Theoretical JT coefficients from EOS
- `jt_uncertainty_analytical_K_per_MPa` - Analytical uncertainty propagation
- `jt_uncertainty_mc_K_per_MPa` - Monte Carlo uncertainty estimation
- `jt_uncertainty_K_per_MPa` - Combined uncertainty (maximum of both methods)
- `relative_error_percent` - Relative error vs. theoretical values

## Experimental Conditions

**Fluids measured:**
- **Pure fluids:** Nitrogen, Argon, Helium-4
- **Binary mixtures:** Helium-Neon, Helium-Nitrogen

**Temperature range:** 65K to 180K
**Pressure range:** 0.1 MPa to 12 MPa

**Measurement uncertainty** (k=1.96, computed by `theoretical_jt_uncertainty.combined_temperature_uncertainty`, using the outlet thermometer TT102 / sensor X93303):
- **Temperature (Cernox sensor only):** ±13 mK at 50 K to ±28 mK at 180 K
- **Temperature (full chain, sensor + CABTR + calibration polynomial):** ±23 mK at 65 K to ±54 mK at 180 K
- **Pressure:** ±0.01% of full scale (13.7 MPa)
- **Composition:** ±0.1 mol% (binary mixtures)

## Theoretical Background

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

**Python packages:**
- `numpy` - Numerical computations
- `pandas` - Data manipulation
- `scipy` - Scientific computing
- `sklearn` - R² calculation
- `matplotlib` - Plotting (for legacy scripts)

**External dependencies:**
- **REFPROP** - Thermodynamic property calculations
- Custom equations of state for He-Ne, He-Ar, Ne-Ar mixtures (see `REFPROP/` directory)
- **LaTeX** (`latex` + `dvipng`) is required by the figure scripts for
  `matplotlib`'s `text.usetex` (the Computer Modern rendering used by the
  shared style sheet). Set `text.usetex: False` in `jced.mplstyle` to skip.

A pinned dependency list is provided at the repository root in
[`requirements.txt`](../../requirements.txt).

## Figure Regeneration

All paper figures are generated by the scripts below and written to the
repository [`img/`](../../img) directory. The scripts share a single
declarative matplotlib style sheet, [`jced.mplstyle`](../jced.mplstyle), which is
loaded explicitly with `plt.style.use(...)`; each script then sets its own
figure-specific layout. There is no shared style function.

| Script | Output(s) in `img/` |
| --- | --- |
| [`plot_paper_results.py`](plot_paper_results.py) | `measurements_pures_errorbars.pdf`, `measurements_HeNe_errorbars.pdf`, `measurements_HeN2_errorbars.pdf` |
| [`theoretical_jt_uncertainty.py`](theoretical_jt_uncertainty.py) | `theoretical_uncertainty_combined.pdf` |
| [`compare_theory_vs_experiment.py`](compare_theory_vs_experiment.py) | `theory_vs_experiment_uncertainty.pdf` |
| [`plot_toc_graphic.py`](plot_toc_graphic.py) | `toc.pdf` |
| [`generate_paper_tables.py`](generate_paper_tables.py) | `data/derived_data/jt_coeffs/paper_tables.tex` (LaTeX bodies of the result tables) |

`plot_paper_results.py`, `compare_theory_vs_experiment.py`, and
`plot_toc_graphic.py` read the derived JT CSVs in
[`data/derived_data/jt_coeffs/`](../../data/derived_data/jt_coeffs) and do not
require REFPROP. `theoretical_jt_uncertainty.py` propagates the equipment
uncertainties analytically and by Monte Carlo and also does not require
REFPROP. The full data-reduction pipeline (`calculate_jt_coefficients.py`)
does require REFPROP.

Regenerate the full figure set with:

```bash
python src/data_analysis/plot_paper_results.py
python src/data_analysis/compare_theory_vs_experiment.py
python src/data_analysis/theoretical_jt_uncertainty.py
python src/data_analysis/plot_toc_graphic.py
```

## Error Handling

The analysis includes robust error handling for:
- **Missing data files** - Graceful skipping with warnings
- **Insufficient data points** - Minimum 3 points required for fitting
- **Numerical instabilities** - Regularization and fallback methods
- **EOS calculation failures** - NaN handling and error reporting
- **Uncertainty propagation errors** - Fallback to simplified estimates

## Data Quality Metrics

**Fitting quality:**
- R² > 0.999 reached by every retained fit
- Polynomial degree: lowest degree (1-6) reaching R² > 0.999
- Weighted fitting accounts for measurement uncertainties

**Measurement accuracy:**
- Typical relative errors: 1-5% vs. EOS predictions
- Best measurements: <1% relative error
- Uncertainty analysis validates measurement quality

## Publication Notes

This analysis forms the foundation for:
- PhD thesis Chapter 4: "Isenthalpic Joule-Thomson coefficient measurements"
- Planned journal publication on He-Ne mixture measurements
- Validation of equations of state for cryogenic applications
- Data repository for reproducible research

The methodology represents a significant advancement in JT coefficient measurements through:
- Modern cryogenic instrumentation
- Rigorous uncertainty analysis
- Comprehensive data validation
- Open-source analysis tools