# Derived data

Data products computed from the raw measurements in [`../raw_data`](../raw_data)
by the analysis pipeline in [`../../src/data_analysis`](../../src/data_analysis).
Column-level details are in
[`../metadata/codebook.csv`](../metadata/codebook.csv); the pipeline and its
run order are described in
[`../../src/data_analysis/README.md`](../../src/data_analysis/README.md).

## `p_T_pairs/`

The steady-state pressure-temperature(-composition) points extracted from the
raw time series, one CSV per isenthalpic line. File names follow
`Fluid_TK_PMPa.csv`, where `T` is the mean inlet temperature and `P` the
mean inlet pressure, both rounded (e.g. `Nitrogen_160K_6MPa.csv`). These are
the inputs to the Joule-Thomson coefficient extraction. See
[`p_T_pairs/README.md`](p_T_pairs/README.md).

## `jt_coeffs/`

The Joule-Thomson coefficients derived from each isenthalp by Chebyshev
fitting and differentiation, with uncertainties, plus the Monte Carlo bound,
the correction estimates, and the generated LaTeX tables. See
[`jt_coeffs/README.md`](jt_coeffs/README.md).

- `Fluid_TK_PMPa_JT.csv`: one file per isenthalp, with the columns:

  | Column | Unit | Description |
  | --- | --- | --- |
  | `p/MPa` | MPa | pressure of the measured point |
  | `p_UNC/MPa` | MPa | expanded (k = 2) uncertainty in pressure |
  | `T/K` | K | temperature of the measured point |
  | `T_UNC/K` | K | expanded (k = 2) uncertainty in temperature |
  | `x1` | mole fraction | first-component mole fraction at the point (mixtures only) |
  | `JT_meas/(K/MPa)` | K MPa⁻¹ | measured μ_JT (derivative of the fit) |
  | `JT_UNC_MC/(K/MPa)` | K MPa⁻¹ | Monte Carlo part of the expanded (k = 2) uncertainty |
  | `JT_UNC/(K/MPa)` | K MPa⁻¹ | expanded (k = 2) uncertainty of μ_JT: Monte Carlo, bias, and method terms combined |
  | `JT_eos/(K/MPa)` | K MPa⁻¹ | μ_JT from the reference EOS at the same point |
  | `mean_unc_perc` | % | mean of 100 `JT_UNC`/\|`JT_meas`\| over the used points of the isenthalp |
  | `rel_abs_err_perc` | % | relative absolute deviation from the EOS |
  | `used` | - | False for the isenthalp end points, which are fitted but not used |

- `jt_coefficients_summary.csv`: one row per isenthalp: fluid, mixture flag,
  composition (`x1`, `x2`), mean inlet temperature, number of points,
  polynomial degree, fit `r2_score`, the mean/std/rms relative absolute
  deviation from the EOS, and `mean_unc_perc`.
- `jt_coefficients_detailed.json`: the same results in nested form, including
  per-point arrays, for programmatic use.

## `capillary_mass_flow/`

Capillary pressure-drop calculations (not measurements) of mass flow rate versus
outlet pressure for each measured isenthalp (D = 0.365 mm, L = 6 m), used by
`estimate_corrections.py` to estimate the mass flow of each point. File names
encode the inlet pressure, temperature, and, for mixtures, the helium fraction.
No script in this repository generates them.

## `calculated_mass_flow_rate/`

Design calculations (not measurements) of the downstream pressure as a function
of mass flow rate through the expansion capillary, used to size the capillary
and to plan the measurement conditions. File names encode the capillary inner
diameter, length, composition, inlet pressure, and temperature, e.g.
`Helium-Neon_D0.365_L6_x0.167_p10_T100.csv` (D in mm, L in m, p in MPa, T in K).
