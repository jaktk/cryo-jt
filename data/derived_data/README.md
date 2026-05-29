# Derived data

Data products computed from the raw measurements in [`../raw_data`](../raw_data)
by the analysis pipeline in [`../../src/data_analysis`](../../src/data_analysis).
Column-level details of the input/output formats are documented in
[`../../src/data_analysis/README.md`](../../src/data_analysis/README.md).

## `p_T_pairs/`

The steady-state pressure–temperature(–composition) points extracted from the
raw time series, one CSV per isenthalpic line. File names follow
`Fluid_TK_PMPa.csv`, where `T` is the nominal inlet temperature and `P` the
nominal inlet pressure (e.g. `Nitrogen_160K_6MPa.csv`). These are the inputs to
the Joule-Thomson coefficient extraction.

## `jt_coeffs/`

The Joule-Thomson coefficients derived from each isenthalp by Chebyshev
fitting and differentiation, with uncertainties.

- `Fluid_TK_PMPa_JT.csv` — one file per isenthalp, with the columns:

  | Column | Unit | Description |
  | --- | --- | --- |
  | `p/MPa` | MPa | pressure of the measured point |
  | `p_UNC/MPa` | MPa | standard uncertainty in pressure |
  | `T/K` | K | temperature of the measured point |
  | `T_UNC/K` | K | standard uncertainty in temperature |
  | `JT_meas/(K/MPa)` | K MPa⁻¹ | measured μ_JT (derivative of the fit) |
  | `JT_UNC/(K/MPa)` | K MPa⁻¹ | uncertainty in the measured μ_JT |
  | `JT_eos/(K/MPa)` | K MPa⁻¹ | μ_JT from the reference EOS at the same point |
  | `mean_unc_perc` | % | conventional relative expanded uncertainty |
  | `rel_abs_err_perc` | % | relative absolute deviation from the EOS |

- `jt_coefficients_summary.csv` — one row per isenthalp: fluid, mixture flag,
  composition (`x1`, `x2`), mean inlet temperature, number of points,
  polynomial degree, fit `r2_score`, and the mean/std/rms relative absolute
  deviation from the EOS together with the mean conventional uncertainty.
- `jt_coefficients_detailed.json` — the same results in nested form, including
  per-point arrays, for programmatic use.

## `calculated_mass_flow_rate/`

Design calculations (not measurements) of the downstream pressure as a function
of mass flow rate through the expansion capillary, used to size the capillary
and to plan the measurement conditions. File names encode the capillary inner
diameter, length, composition, inlet pressure, and temperature, e.g.
`Helium-Neon_D0.365_L6_x0.167_p10_T100.csv` (D in mm, L in m, p in MPa, T in K).
