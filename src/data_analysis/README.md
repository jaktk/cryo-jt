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
- **Monte Carlo validation** of analytical uncertainty calculations
- **Automatic polynomial degree selection** based on R² optimization
- **Comprehensive output** with individual and summary results

**Usage:**
```bash
python src/data_analysis/calculate_jt_coefficients.py
```

**Input Data:**
- `data/derived_data/p_T_pairs/*.csv` - Processed pressure-temperature pairs with uncertainties
- `data/metadata/p_T_pairs.json` - Metadata with fluid types and reference temperatures

**Output:**
- `data/derived_data/JT_coeffs/jt_coefficients_summary.csv` - Summary results table
- `data/derived_data/JT_coeffs/jt_coefficients_detailed.json` - Complete results with arrays
- `data/derived_data/JT_coeffs/*_JT.csv` - Individual measurement files

**Methodology:**

1. **Polynomial Fitting:**
   - Fits Chebyshev polynomials to isenthalpic p-T data
   - Uses weighted least squares with measurement uncertainties
   - Automatically selects optimal polynomial degree (1-6)
   - Calculates covariance matrix for coefficient uncertainties

2. **JT Coefficient Calculation:**
   - Computes μ_JT = dT/dp by taking polynomial derivative
   - Propagates uncertainties through differentiation using covariance matrix
   - Implements proper mathematical uncertainty propagation

3. **Uncertainty Analysis:**
   - **Analytical propagation:** Uses covariance matrix of polynomial coefficients
   - **Monte Carlo validation:** 1000 samples with perturbed measurements
   - **Conservative approach:** Uses maximum of both uncertainty estimates

4. **Validation:**
   - Compares measured vs. theoretical JT coefficients using REFPROP
   - Calculates relative errors and goodness-of-fit metrics
   - Provides detailed quality assessment

### `analyze_and_plot_JT.py`

Original analysis script from thesis work containing plotting and validation functions.

**Key Classes:**
- `CernoxCal` - Cernox temperature sensor calibration
- `StdTempUncertainty` - Temperature measurement uncertainty calculations
- `Isenthalpic` - Isenthalpic transformation calculations

**Functions:**
- `plot_isenthalps_and_JT_coefs()` - Comprehensive plotting of results
- `validate_nitrogen()`, `validate_argon()`, `validate_helium()` - Pure fluid validation
- `JT_helium_nitrogen()`, `JT_helium_neon()` - Mixture measurements

### `monte_carlo_uncertainty_analysis.py`

Monte Carlo uncertainty analysis for JT coefficient measurements.

**Key Functions:**
- `monte_carlo_pures()` - Pure fluid uncertainty analysis
- `monte_carlo_mixtures()` - Mixture uncertainty with composition effects
- `monte_carlo_*_general()` - General uncertainty estimation routines

### Supporting Files

- `FluidProps.py` - REFPROP wrapper for thermodynamic property calculations
- `get_git_root.py` - Utility for finding repository root directory
- `extract_file_info.py` - Utility for extracting metadata from data files

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

**Measurement uncertainty:**
- **Temperature:** ±12-17 mK (Cernox sensors)
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

## Error Handling

The analysis includes robust error handling for:
- **Missing data files** - Graceful skipping with warnings
- **Insufficient data points** - Minimum 3 points required for fitting
- **Numerical instabilities** - Regularization and fallback methods
- **EOS calculation failures** - NaN handling and error reporting
- **Uncertainty propagation errors** - Fallback to simplified estimates

## Data Quality Metrics

**Fitting quality:**
- R² score > 0.95 typically achieved
- Polynomial degree selection based on error minimization
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